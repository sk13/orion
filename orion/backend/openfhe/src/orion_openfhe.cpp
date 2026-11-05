// OpenFHE implementation of the Orion backend C ABI (orion/backend/c_api/orion_backend.h).
//
// Mapping notes
//  * CKKS with FIXEDMANUAL scaling: Orion decides when to rescale, like with Lattigo.
//  * OpenFHE tracks scales inexactly in this mode (it assumes 2^p), so this backend keeps
//    the exact scale of every ciphertext/plaintext in OpenFHE's scaling-factor field and
//    corrects values on decode. Plaintexts are encoded at arbitrary scales by pre-scaling.
//  * Levels use Lattigo's convention (level = #RNS limbs - 1).
//  * ConjugateInvariant rings do not exist in OpenFHE: a conjugate-invariant ring of
//    degree N (N real slots) is emulated by a standard ring of degree 2N (N complex slots,
//    real values in the real lane), which gives the same slot count and rotation semantics.
//  * Rotation-key ids exposed to Python are plain rotation amounts.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "openfhe.h"
#include "scheme/ckksrns/ckksrns-fhe.h"
#include "orion_backend.h"

using namespace lbcrypto;
using CT = Ciphertext<DCRTPoly>;

namespace {

[[noreturn]] void fatal(const char* fn, const char* what) {
    std::fprintf(stderr, "[orion-openfhe] %s failed: %s\n", fn, what);
    std::fflush(stderr);
    std::abort();
}

struct ProfEntry { double ms = 0; long calls = 0; };
std::map<std::string, ProfEntry>& prof_table() { static std::map<std::string, ProfEntry> t; return t; }
bool prof_on() { static const bool on = std::getenv("ORION_OPENFHE_PROFILE") != nullptr; return on; }
struct ProfScope {
    const char* name; std::chrono::steady_clock::time_point t0;
    explicit ProfScope(const char* n) : name(n), t0(std::chrono::steady_clock::now()) {}
    ~ProfScope() {
        if (!prof_on()) return;
        auto& e = prof_table()[name];
        e.ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        e.calls++;
    }
};
#define PROF(name) ProfScope _prof_scope(name)
#define API_BEGIN try {
#define API_END(name) } catch (const std::exception& e) { fatal(name, e.what()); } catch (...) { fatal(name, "unknown error"); }
[[noreturn]] void unsupported(const char* fn) { fatal(fn, "not implemented in the OpenFHE backend yet"); }

// Handle allocator: reuses the lowest freed id, like the Go backend.
template <class T>
struct Table {
    std::unordered_map<int, T> m;
    std::priority_queue<int, std::vector<int>, std::greater<int>> freed;
    int next = 0;
    int add(T v) {
        int id;
        if (!freed.empty()) { id = freed.top(); freed.pop(); } else { id = next++; }
        m[id] = std::move(v);
        return id;
    }
    T& get(int id) {
        auto it = m.find(id);
        if (it == m.end()) throw std::runtime_error("heap object not found: " + std::to_string(id));
        return it->second;
    }
    void del(int id) { if (m.erase(id)) freed.push(id); }
    void reset() { m.clear(); freed = {}; next = 0; }
    std::vector<int> keys() const { std::vector<int> k; for (auto& kv : m) k.push_back(kv.first); return k; }
};

struct PtEntry {
    Plaintext pt;
    std::vector<std::complex<double>> vals;  // original values (for re-encoding at a lower level)
    bool real = true;
    double scale = 1.0;     // exact scale of the encoded values
    double decodeCorr = 1.0; // multiplier applied to OpenFHE's decoded values
};

struct LinTransform {
    int level = 0;
    int n1 = 1;
    std::vector<int> babyIdx;     // distinct baby-step rotations (excluding 0)
    std::vector<int> giants;      // ascending giant offsets G (multiples of n1)
    std::set<int> rotKeys;        // every rotation amount the evaluation needs
    // For each giant offset: (baby offset, plaintext in the extended basis Q_l*P)
    std::vector<std::vector<std::pair<int, Plaintext>>> terms;
    double scale = 1.0;
};

struct PolyEntry { std::vector<double> coeffs; bool chebyshev = false; };

struct Scheme {
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> kp;
    uint32_t depth = 0;      // max level (Lattigo convention)
    uint32_t slots = 0;
    uint32_t scaleBits = 0;  // p, plaintext scaling exponent
    std::vector<uint64_t> Q, P;
    std::set<int> rotKeys;   // rotation amounts with generated keys
    bool relinGenerated = false;

    void reset() {
        cc = nullptr;
        kp = KeyPair<DCRTPoly>();
        depth = slots = scaleBits = 0;
        Q.clear(); P.clear(); rotKeys.clear();
        relinGenerated = false;
    }
} S;

Table<PtEntry> pts;
Table<CT> cts;
Table<LinTransform> lts;
Table<PolyEntry> polys;

int normRot(int r) {
    int n = static_cast<int>(S.slots);
    return ((r % n) + n) % n;
}

int level_of(const CT& ct) { return static_cast<int>(ct->GetElements()[0].GetNumOfElements()) - 1; }
int level_of(const Plaintext& pt) { return static_cast<int>(pt->GetElement<DCRTPoly>().GetNumOfElements()) - 1; }

Plaintext encode(const std::vector<std::complex<double>>& v, bool real, int level, double scale) {
    if (level < 0 || level > static_cast<int>(S.depth)) throw std::runtime_error("encode: invalid level");
    uint32_t ofheLevel = S.depth - level;
    double p2 = std::pow(2.0, S.scaleBits);
    double f = scale / p2;
    auto make = [&](bool zeros) {
        if (real) {
            std::vector<double> x(v.size());
            for (size_t i = 0; i < v.size(); i++) x[i] = zeros ? 0.0 : v[i].real() * f;
            return S.cc->MakeCKKSPackedPlaintext(x, 1, ofheLevel, nullptr, S.slots);
        }
        std::vector<std::complex<double>> x(v.size());
        for (size_t i = 0; i < v.size(); i++) x[i] = zeros ? 0.0 : v[i] * f;
        return S.cc->MakeCKKSPackedPlaintext(x, 1, ofheLevel, nullptr, S.slots);
    };
    try {
        return make(false);
    } catch (const std::exception& e) {
        // Values too small for the scale all round to zero; OpenFHE refuses to encode those.
        if (std::strstr(e.what(), "too small")) return make(true);
        throw;
    }
}

// Like Lattigo, operations between a ciphertext and a plaintext encoded at a higher
// level use the lower level: re-encode the plaintext at the ciphertext's level.
Plaintext pt_at_level(PtEntry& e, int level) {
    if (level_of(e.pt) == level) return e.pt;
    if (level > level_of(e.pt) || e.vals.empty())
        throw std::runtime_error("plaintext op: cannot adapt plaintext level " + std::to_string(level_of(e.pt)) +
                                 " to ciphertext level " + std::to_string(level));
    return encode(e.vals, e.real, level, e.scale);
}

CT ct_at_level(const CT& ct, int level) {
    int diff = level_of(ct) - level;
    if (diff == 0) return ct;
    if (diff < 0) throw std::runtime_error("cannot raise ciphertext level");
    return S.cc->LevelReduce(ct, nullptr, diff);
}

int push_ct(CT ct, double scale) {
    ct->SetScalingFactor(scale);
    return cts.add(std::move(ct));
}

double scale_of(const CT& ct) { return ct->GetScalingFactor(); }

// Replace the ciphertext behind `id` (in-place semantics of the Go backend).
int replace_ct(int id, CT ct, double scale) {
    ct->SetScalingFactor(scale);
    cts.m[id] = std::move(ct);
    return id;
}

template <class F>
int unary(int id, bool inPlace, F f) {  // f: (const CT&) -> pair<CT,double scale>
    auto r = f(cts.get(id));
    return inPlace ? replace_ct(id, r.first, r.second) : push_ct(r.first, r.second);
}

// In OpenFHE ops with a Python-side in-place flag
using Res = std::pair<CT, double>;

ArrayResultFloat float_array(const std::vector<double>& v) {
    ArrayResultFloat r{nullptr, 0};
    if (v.empty()) return r;
    r.Data = static_cast<float*>(std::malloc(v.size() * sizeof(float)));
    for (size_t i = 0; i < v.size(); i++) r.Data[i] = static_cast<float>(v[i]);
    r.Length = v.size();
    return r;
}
ArrayResultInt int_array(const std::vector<int>& v) {
    ArrayResultInt r{nullptr, 0};
    if (v.empty()) return r;
    r.Data = static_cast<int*>(std::malloc(v.size() * sizeof(int)));
    std::copy(v.begin(), v.end(), r.Data);
    r.Length = v.size();
    return r;
}
ArrayResultUInt64 u64_array(const std::vector<uint64_t>& v) {
    ArrayResultUInt64 r{nullptr, 0};
    if (v.empty()) return r;
    r.Data = static_cast<unsigned long long*>(std::malloc(v.size() * sizeof(unsigned long long)));
    std::copy(v.begin(), v.end(), r.Data);
    r.Length = v.size();
    return r;
}

void ensure_rot_key(int rot) {
    rot = normRot(rot);
    if (rot == 0 || S.rotKeys.count(rot)) return;
    S.cc->EvalAtIndexKeyGen(S.kp.secretKey, {rot});
    S.rotKeys.insert(rot);
}

// Rotation by +r moves slot r+i to slot i (same convention as Lattigo).
CT rotate(const CT& ct, int amount) {
    int r = normRot(amount);
    if (r == 0) return ct->Clone();
    ensure_rot_key(r);
    return S.cc->EvalRotate(ct, r);
}

}  // namespace

extern "C" {

// ---------------------------------------------------------------- scheme --
void NewScheme(int logN, int* logQ, int lenQ, int* logP, int lenP, int logScale, int h,
               char* ringType, char* keysPath, char* ioMode) {
    (void)h; (void)keysPath; (void)ioMode;
    API_BEGIN
    if (lenQ < 1) throw std::runtime_error("LogQ must not be empty");
    for (int i = 1; i < lenQ; i++)
        if (logQ[i] != logScale)
            throw std::runtime_error("OpenFHE backend (FIXEDMANUAL) requires every LogQ entry after the first "
                                     "to equal LogScale");
    bool ci = std::string(ringType) != "standard";
    uint32_t ringDim = 1u << (logN + (ci ? 1 : 0));

    CCParams<CryptoContextCKKSRNS> p;
    p.SetSecurityLevel(HEStd_NotSet);
    p.SetRingDim(ringDim);
    p.SetBatchSize(ringDim / 2);
    p.SetMultiplicativeDepth(lenQ - 1);
    p.SetScalingModSize(logScale);
    p.SetFirstModSize(logQ[0]);
    p.SetScalingTechnique(FIXEDMANUAL);
    p.SetKeySwitchTechnique(HYBRID);
    p.SetSecretKeyDist(UNIFORM_TERNARY);
    if (lenP > 0) p.SetNumLargeDigits((lenQ + lenP - 1) / lenP);

    S.reset();
    S.cc = GenCryptoContext(p);
    S.cc->Enable(PKE);
    S.cc->Enable(KEYSWITCH);
    S.cc->Enable(LEVELEDSHE);
    S.depth = lenQ - 1;
    S.slots = ringDim / 2;
    S.scaleBits = logScale;

    auto cp = S.cc->GetCryptoParameters();
    for (auto& t : cp->GetElementParams()->GetParams()) S.Q.push_back(t->GetModulus().ConvertToInt());
    if (auto pp = std::dynamic_pointer_cast<CryptoParametersRNS>(cp)->GetParamsP())
        for (auto& t : pp->GetParams()) S.P.push_back(t->GetModulus().ConvertToInt());
    API_END("NewScheme")
}

void DeleteScheme(void) {
    if (prof_on()) {
        for (auto& kv : prof_table())
            std::fprintf(stderr, "[PROF] %-26s %8.1f ms  %6ld calls\n", kv.first.c_str(), kv.second.ms, kv.second.calls);
        prof_table().clear();
    }
    pts.reset(); cts.reset(); lts.reset(); polys.reset();
    S.reset();
    CryptoContextFactory<DCRTPoly>::ReleaseAllContexts();
}

void FreeCArray(void* ptr) { std::free(ptr); }

// --------------------------------------------------------------- tensors --
void DeletePlaintext(int id) { pts.del(id); }
void DeleteCiphertext(int id) { cts.del(id); }
unsigned long GetPlaintextScale(int id) { return static_cast<unsigned long>(std::llround(pts.get(id).scale)); }
unsigned long GetCiphertextScale(int id) { return static_cast<unsigned long>(std::llround(scale_of(cts.get(id)))); }
void SetPlaintextScale(int id, unsigned long s) { pts.get(id).scale = static_cast<double>(s); }
void SetCiphertextScale(int id, unsigned long s) { cts.get(id)->SetScalingFactor(static_cast<double>(s)); }
int GetPlaintextLevel(int id) { return level_of(pts.get(id).pt); }
int GetCiphertextLevel(int id) { return level_of(cts.get(id)); }
int GetPlaintextSlots(int) { return static_cast<int>(S.slots); }
int GetCiphertextSlots(int) { return static_cast<int>(S.slots); }
int GetCiphertextDegree(int id) { return static_cast<int>(cts.get(id)->NumberCiphertextElements()) - 1; }
ArrayResultUInt64 GetModuliChain(void) { return u64_array(S.Q); }
ArrayResultUInt64 GetAuxModuliChain(void) { return u64_array(S.P); }
ArrayResultInt GetLivePlaintexts(void) { return int_array(pts.keys()); }
ArrayResultInt GetLiveCiphertexts(void) { return int_array(cts.keys()); }

// ----------------------------------------------------------- key generation --
void NewKeyGenerator(void) {}
void GenerateSecretKey(void) {
    API_BEGIN
    S.kp = S.cc->KeyGen();  // OpenFHE generates the secret/public pair together
    API_END("GenerateSecretKey")
}
void GeneratePublicKey(void) {}
void GenerateRelinearizationKey(void) {
    API_BEGIN
    S.cc->EvalMultKeyGen(S.kp.secretKey);
    S.relinGenerated = true;
    API_END("GenerateRelinearizationKey")
}
void GenerateEvaluationKeys(void) {}
ArrayResultByte SerializeSecretKey(void) { unsupported("SerializeSecretKey"); }
void LoadSecretKey(unsigned char*, unsigned long) { unsupported("LoadSecretKey"); }

// ------------------------------------------------------------------ encoder --
void NewEncoder(void) {}

int Encode(float* values, int len, int level, unsigned long scale) {
    API_BEGIN PROF("Encode");
    std::vector<std::complex<double>> v(len);
    for (int i = 0; i < len; i++) v[i] = values[i];
    PtEntry e;
    e.pt = encode(v, true, level, static_cast<double>(scale));
    e.vals = v; e.real = true;
    e.scale = static_cast<double>(scale);
    return pts.add(std::move(e));
    API_END("Encode")
}

ArrayResultFloat Decode(int id) {
    API_BEGIN
    auto& e = pts.get(id);
    auto vals = e.pt->GetRealPackedValue();
    for (auto& x : vals) x *= e.decodeCorr;
    vals.resize(S.slots, 0.0);
    return float_array(vals);
    API_END("Decode")
}

int EncodeComplex(float* interleaved, int len, int level, unsigned long scale) {
    API_BEGIN
    std::vector<std::complex<double>> v(len / 2);
    for (int i = 0; i < len / 2; i++) v[i] = {interleaved[2 * i], interleaved[2 * i + 1]};
    PtEntry e;
    e.pt = encode(v, false, level, static_cast<double>(scale));
    e.vals = v; e.real = false;
    e.scale = static_cast<double>(scale);
    return pts.add(std::move(e));
    API_END("EncodeComplex")
}

ArrayResultFloat DecodeComplex(int id) {
    API_BEGIN
    auto& e = pts.get(id);
    auto vals = e.pt->GetCKKSPackedValue();
    std::vector<double> out(2 * S.slots, 0.0);
    for (size_t i = 0; i < vals.size() && i < S.slots; i++) {
        out[2 * i] = vals[i].real() * e.decodeCorr;
        out[2 * i + 1] = vals[i].imag() * e.decodeCorr;
    }
    return float_array(out);
    API_END("DecodeComplex")
}

// ---------------------------------------------------------------- encryptor --
void NewEncryptor(void) {}
void NewDecryptor(void) {}

int Encrypt(int id) {
    API_BEGIN PROF("Encrypt");
    auto& e = pts.get(id);
    CT ct = S.cc->Encrypt(S.kp.publicKey, e.pt);
    return push_ct(ct, e.scale);
    API_END("Encrypt")
}

int Decrypt(int id) {
    API_BEGIN PROF("Decrypt");
    CT ct = cts.get(id);
    if (ct->GetNoiseScaleDeg() != 1)
        throw std::runtime_error("decrypt: ciphertext must be rescaled first (noise scale degree " +
                                 std::to_string(ct->GetNoiseScaleDeg()) + ")");
    Plaintext pt;
    S.cc->Decrypt(S.kp.secretKey, ct, &pt);
    PtEntry e;
    e.pt = pt;
    e.scale = scale_of(ct);
    // OpenFHE decodes assuming scale 2^p; correct to the true scale.
    e.decodeCorr = std::pow(2.0, S.scaleBits) / e.scale;
    return pts.add(std::move(e));
    API_END("Decrypt")
}

// ---------------------------------------------------------------- evaluator --
void NewEvaluator(void) {
    API_BEGIN
    for (uint32_t i = 1; i < S.slots; i *= 2) ensure_rot_key(static_cast<int>(i));
    API_END("NewEvaluator")
}
void AddRotationKey(int rotation) {
    API_BEGIN
    ensure_rot_key(rotation);
    API_END("AddRotationKey")
}

int Negate(int id) {
    API_BEGIN
    auto& ct = cts.get(id);
    return push_ct(S.cc->EvalNegate(ct), scale_of(ct));
    API_END("Negate")
}

static int do_rotate(int id, int amount, bool inPlace) {
    return unary(id, inPlace, [&](const CT& ct) { return Res{rotate(ct, amount), scale_of(ct)}; });
}
int Rotate(int id, int amount) { API_BEGIN PROF("Rotate"); return do_rotate(id, amount, true); API_END("Rotate") }
int RotateNew(int id, int amount) { API_BEGIN PROF("RotateNew"); return do_rotate(id, amount, false); API_END("RotateNew") }

static int do_rescale(int id, bool inPlace) {
    return unary(id, inPlace, [&](const CT& ct) {
        int lvl = level_of(ct);
        if (lvl < 1) throw std::runtime_error("rescale: no level left");
        double s = scale_of(ct) / static_cast<double>(S.Q[lvl]);
        return Res{S.cc->ModReduce(ct), s};
    });
}
int Rescale(int id) { API_BEGIN PROF("Rescale"); return do_rescale(id, true); API_END("Rescale") }
int RescaleNew(int id) { API_BEGIN PROF("RescaleNew"); return do_rescale(id, false); API_END("RescaleNew") }

static int scalar_add(int id, double s, bool inPlace, bool sub) {
    return unary(id, inPlace, [&](const CT& ct) {
        return Res{sub ? S.cc->EvalSub(ct, s) : S.cc->EvalAdd(ct, s), scale_of(ct)};
    });
}
int AddScalar(int id, float s) { API_BEGIN return scalar_add(id, s, true, false); API_END("AddScalar") }
int AddScalarNew(int id, float s) { API_BEGIN return scalar_add(id, s, false, false); API_END("AddScalarNew") }
int SubScalar(int id, float s) { API_BEGIN return scalar_add(id, s, true, true); API_END("SubScalar") }
int SubScalarNew(int id, float s) { API_BEGIN return scalar_add(id, s, false, true); API_END("SubScalarNew") }

// Integer scalars multiply the RNS limbs directly: no scale or level change.
static int scalar_mul_int(int id, int s, bool inPlace) {
    return unary(id, inPlace, [&](const CT& ct) {
        CT out = ct->Clone();
        for (auto& el : out->GetElements()) {
            el = el.Times(NativeInteger(static_cast<uint64_t>(std::llabs(static_cast<long long>(s)))));
            if (s < 0) el = el.Negate();
        }
        return Res{out, scale_of(ct)};
    });
}
int MulScalarInt(int id, int s) { API_BEGIN return scalar_mul_int(id, s, true); API_END("MulScalarInt") }
int MulScalarIntNew(int id, int s) { API_BEGIN return scalar_mul_int(id, s, false); API_END("MulScalarIntNew") }

// Float scalars are encoded at scale 2^p; the caller rescales afterwards.
static int scalar_mul_float(int id, double s, bool inPlace) {
    return unary(id, inPlace, [&](const CT& ct) {
        return Res{S.cc->EvalMult(ct, s), scale_of(ct) * std::pow(2.0, S.scaleBits)};
    });
}
int MulScalarFloat(int id, float s) { API_BEGIN return scalar_mul_float(id, s, true); API_END("MulScalarFloat") }
int MulScalarFloatNew(int id, float s) { API_BEGIN return scalar_mul_float(id, s, false); API_END("MulScalarFloatNew") }

static int pt_op(int cid, int pid, bool inPlace, int kind) {  // 0 add, 1 sub, 2 mul
    return unary(cid, inPlace, [&](const CT& ct) {
        auto& e = pts.get(pid);
        Plaintext pt = pt_at_level(e, level_of(ct));
        if (kind == 0) return Res{S.cc->EvalAdd(ct, pt), scale_of(ct)};
        if (kind == 1) return Res{S.cc->EvalSub(ct, pt), scale_of(ct)};
        return Res{S.cc->EvalMult(ct, pt), scale_of(ct) * e.scale};
    });
}
int AddPlaintext(int c, int p) { API_BEGIN PROF("AddPlaintext"); return pt_op(c, p, true, 0); API_END("AddPlaintext") }
int AddPlaintextNew(int c, int p) { API_BEGIN PROF("AddPlaintextNew"); return pt_op(c, p, false, 0); API_END("AddPlaintextNew") }
int SubPlaintext(int c, int p) { API_BEGIN return pt_op(c, p, true, 1); API_END("SubPlaintext") }
int SubPlaintextNew(int c, int p) { API_BEGIN return pt_op(c, p, false, 1); API_END("SubPlaintextNew") }
int MulPlaintext(int c, int p) { API_BEGIN PROF("MulPlaintext"); return pt_op(c, p, true, 2); API_END("MulPlaintext") }
int MulPlaintextNew(int c, int p) { API_BEGIN return pt_op(c, p, false, 2); API_END("MulPlaintextNew") }

static int ct_op(int a, int b, bool inPlace, int kind) {  // 0 add, 1 sub, 2 mul+relin
    return unary(a, inPlace, [&](const CT& x) {
        int lvl = std::min(level_of(x), level_of(cts.get(b)));
        CT xa = ct_at_level(x, lvl), ya = ct_at_level(cts.get(b), lvl);
        if (kind == 0) return Res{S.cc->EvalAdd(xa, ya), scale_of(x)};
        if (kind == 1) return Res{S.cc->EvalSub(xa, ya), scale_of(x)};
        return Res{S.cc->EvalMult(xa, ya), scale_of(x) * scale_of(ya)};
    });
}
int AddCiphertext(int a, int b) { API_BEGIN PROF("AddCiphertext"); return ct_op(a, b, true, 0); API_END("AddCiphertext") }
int AddCiphertextNew(int a, int b) { API_BEGIN return ct_op(a, b, false, 0); API_END("AddCiphertextNew") }
int SubCiphertext(int a, int b) { API_BEGIN return ct_op(a, b, true, 1); API_END("SubCiphertext") }
int SubCiphertextNew(int a, int b) { API_BEGIN return ct_op(a, b, false, 1); API_END("SubCiphertextNew") }
int MulRelinCiphertext(int a, int b) { API_BEGIN PROF("MulRelinCiphertext"); return ct_op(a, b, true, 2); API_END("MulRelinCiphertext") }
int MulRelinCiphertextNew(int a, int b) { API_BEGIN PROF("MulRelinCiphertextNew"); return ct_op(a, b, false, 2); API_END("MulRelinCiphertextNew") }

// ------------------------------------------------------ polynomial evaluator --
void NewPolynomialEvaluator(void) {}
int GenerateMonomial(float* coeffs, int len) {
    API_BEGIN
    PolyEntry p; p.coeffs.assign(coeffs, coeffs + len);
    return polys.add(std::move(p));
    API_END("GenerateMonomial")
}
int GenerateChebyshev(float* coeffs, int len) {
    API_BEGIN
    PolyEntry p; p.coeffs.assign(coeffs, coeffs + len); p.chebyshev = true;
    return polys.add(std::move(p));
    API_END("GenerateChebyshev")
}
int EvaluatePolynomial(int, int, unsigned long) { unsupported("EvaluatePolynomial"); }
ArrayResultDouble GenerateMinimaxSignCoeffs(int*, int, int, int, int, int) { unsupported("GenerateMinimaxSignCoeffs"); }

// ---------------------------------------------------------- linear transforms --
void NewLinearTransformEvaluator(void) {}

// Baby-step/giant-step with double hoisting (as in OpenFHE's own bootstrapping linear
// transforms): diagonal d = G + i is applied as rot_G( sum_i rot_{-G}(diag_d) * rot_i(ct) ).
// Baby rotations share one digit decomposition and stay in the extended basis Q_l*P, the
// giant steps are accumulated Horner-style in the extended basis, and a single mod-down
// is done at the end.
int GenerateLinearTransform(int* diagIdxs, int lenIdxs, float* diagData, int lenData, int level,
                            float bsgsRatio, char* ioMode) {
    API_BEGIN
    if (std::string(ioMode) != "none") unsupported("GenerateLinearTransform (io_mode save/load)");
    (void)bsgsRatio;
    const int n = static_cast<int>(S.slots);
    if (lenData != lenIdxs * n) throw std::runtime_error("linear transform: unexpected diagonal data length");
    if (level < 1 || level > static_cast<int>(S.depth)) throw std::runtime_error("linear transform: invalid level");

    std::vector<int> idx(lenIdxs);
    for (int i = 0; i < lenIdxs; i++) idx[i] = normRot(diagIdxs[i]);

    // Choose the baby-step size: hoisted baby rotations cost about half of a giant step.
    const double babyCost = 0.5, giantCost = 1.0;
    int bestN1 = 1; double bestCost = 1e300;
    for (int n1 = 1; n1 <= n; n1++) {
        std::set<int> b, g;
        for (int d : idx) { if (d % n1) b.insert(d % n1); g.insert((d / n1) * n1); }
        double cost = babyCost * b.size() + giantCost * (g.size() - 1 + (*g.begin() > 0 ? 1 : 0));
        if (cost < bestCost) { bestCost = cost; bestN1 = n1; }
    }

    LinTransform lt;
    lt.level = level;
    lt.n1 = bestN1;
    lt.scale = static_cast<double>(S.Q[level]);

    std::set<int> babies, giantSet;
    for (int d : idx) { if (d % bestN1) babies.insert(d % bestN1); giantSet.insert((d / bestN1) * bestN1); }
    lt.babyIdx.assign(babies.begin(), babies.end());
    lt.giants.assign(giantSet.begin(), giantSet.end());
    lt.terms.resize(lt.giants.size());
    lt.rotKeys = babies;
    for (size_t k = 0; k + 1 < lt.giants.size(); k++) lt.rotKeys.insert(lt.giants[k + 1] - lt.giants[k]);
    if (lt.giants.front() > 0) lt.rotKeys.insert(lt.giants.front());

    // Element parameters of the extended basis Q_l * P.
    auto cp = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(S.cc->GetCryptoParameters());
    auto elementParams = *(cp->GetElementParams());
    for (uint32_t i = 0; i < S.depth - level; i++) elementParams.PopLastParam();
    const auto& paramsQ = elementParams.GetParams();
    const auto& paramsP = cp->GetParamsP()->GetParams();
    std::vector<NativeInteger> moduli, roots;
    for (auto& t : paramsQ) { moduli.push_back(t->GetModulus()); roots.push_back(t->GetRootOfUnity()); }
    for (auto& t : paramsP) { moduli.push_back(t->GetModulus()); roots.push_back(t->GetRootOfUnity()); }
    auto paramsQlP = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(S.cc->GetCyclotomicOrder(), moduli, roots);

    struct Job { size_t g; int i; int k; };
    std::vector<Job> jobs;
    for (int k = 0; k < lenIdxs; k++) {
        int d = idx[k], G = (d / bestN1) * bestN1;
        size_t gi = std::lower_bound(lt.giants.begin(), lt.giants.end(), G) - lt.giants.begin();
        jobs.push_back({gi, d % bestN1, k});
        lt.terms[gi].push_back({d % bestN1, nullptr});
    }
    // fill plaintexts (each term slot matches its job in insertion order per giant)
    std::vector<size_t> fill(lt.giants.size(), 0);
    std::vector<std::pair<size_t, size_t>> where(jobs.size());
    for (size_t j = 0; j < jobs.size(); j++) where[j] = {jobs[j].g, fill[jobs[j].g]++};

    const double f = lt.scale / std::pow(2.0, S.scaleBits);
    const uint32_t ofheLevel = S.depth - level;
    const int njobs = static_cast<int>(jobs.size());
    std::string err;
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(njobs))
    for (int j = 0; j < njobs; j++) {
        try {
            int G = lt.giants[jobs[j].g];
            const float* diag = diagData + static_cast<size_t>(jobs[j].k) * n;
            // plaintext = diag rotated by -G so that the giant rotation by G realigns it
            std::vector<std::complex<double>> v(n);
            for (int t = 0; t < n; t++) v[t] = static_cast<double>(diag[((t - G) % n + n) % n]) * f;
            lt.terms[where[j].first][where[j].second].second =
                    FHECKKSRNS::MakeAuxPlaintext(*S.cc, paramsQlP, v, 1, ofheLevel, n);
        } catch (const std::exception& e) {
#pragma omp critical
            err = e.what();
        }
    }
    if (!err.empty()) throw std::runtime_error(err);
    return lts.add(std::move(lt));
    API_END("GenerateLinearTransform")
}

int EvaluateLinearTransform(int tid, int cid) {
    API_BEGIN PROF("EvaluateLinearTransform");
    auto& lt = lts.get(tid);
    const CT& ct = cts.get(cid);
    if (level_of(ct) != lt.level)
        throw std::runtime_error("linear transform: input level " + std::to_string(level_of(ct)) +
                                 " != transform level " + std::to_string(lt.level));
    using clk = std::chrono::steady_clock;
    static const bool prof = std::getenv("ORION_OPENFHE_PROFILE") != nullptr;
    auto ms = [](auto a, auto b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    auto t0 = clk::now();

    // Hoisted baby steps in the extended basis (one digit decomposition for all of them).
    for (int i : lt.babyIdx) ensure_rot_key(i);
    for (int k : lt.rotKeys) ensure_rot_key(k);
    std::vector<CT> fast(lt.n1);
    fast[0] = S.cc->KeySwitchExt(ct, true);
    if (!lt.babyIdx.empty()) {
        auto digits = S.cc->EvalFastRotationPrecompute(ct);
        const int nb = static_cast<int>(lt.babyIdx.size());
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(nb))
        for (int b = 0; b < nb; b++)
            fast[lt.babyIdx[b]] = S.cc->EvalFastRotationExt(ct, lt.babyIdx[b], digits, true);
    }
    auto t1 = clk::now();

    // Inner sums per giant offset (extended basis).
    const int ng = static_cast<int>(lt.giants.size());
    std::vector<CT> inner(ng);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(ng))
    for (int g = 0; g < ng; g++) {
        auto& terms = lt.terms[g];
        // acc = sum_t fast[b_t] * pt_t, without copying the plaintexts and reusing one scratch poly
        auto firstEl = [&](const CT& a, const Plaintext& pt) {
            CT r = a->Clone();
            const DCRTPoly& pe = pt->GetElement<DCRTPoly>();
            for (auto& c : r->GetElements()) c *= pe;
            return r;
        };
        CT acc = firstEl(fast[terms[0].first], terms[0].second);
        DCRTPoly tmp;
        for (size_t t = 1; t < terms.size(); t++) {
            const DCRTPoly& pe = terms[t].second->GetElement<DCRTPoly>();
            const auto& src = fast[terms[t].first]->GetElements();
            auto& dst = acc->GetElements();
            for (size_t e = 0; e < dst.size(); e++) {
                tmp = src[e];
                tmp *= pe;
                dst[e] += tmp;
            }
        }
        acc->SetNoiseScaleDeg(fast[0]->GetNoiseScaleDeg() + 1);
        inner[g] = std::move(acc);
    }
    auto t2 = clk::now();

    // Horner accumulation of the giant steps: outer = rot(outer) + inner_g, descending.
    auto giantRotate = [&](const CT& outer, int amount) {
        uint32_t autoIndex = 0;
        std::vector<uint32_t> map;
        auto key = FHECKKSRNS::GetGiantStepRotation(ct, amount, autoIndex, map);
        return FHECKKSRNS::EvalHornerGiantRotate(outer, autoIndex, map, key);
    };
    CT outer = inner[ng - 1];
    for (int g = ng - 2; g >= 0; g--) {
        outer = giantRotate(outer, lt.giants[g + 1] - lt.giants[g]);
        FHECKKSRNS::EvalAddExtInPlace(outer, inner[g]);
    }
    if (lt.giants.front() > 0) outer = giantRotate(outer, lt.giants.front());
    CT out = S.cc->KeySwitchDown(outer);
    if (prof) {
        auto t3 = clk::now();
        std::fprintf(stderr, "[LT] n1=%d babies=%zu giants=%d | hoisted baby rot %.1f ms, mul+add %.1f ms, giant+down %.1f ms\n",
                     lt.n1, lt.babyIdx.size(), ng, ms(t0, t1), ms(t1, t2), ms(t2, t3));
    }
    return push_ct(out, scale_of(ct) * lt.scale);
    API_END("EvaluateLinearTransform")
}

void DeleteLinearTransform(int id) { lts.del(id); }

ArrayResultInt GetLinearTransformRotationKeys(int tid) {
    API_BEGIN
    auto& lt = lts.get(tid);
    return int_array(std::vector<int>(lt.rotKeys.begin(), lt.rotKeys.end()));
    API_END("GetLinearTransformRotationKeys")
}

void GenerateLinearTransformRotationKey(int keyId) {
    API_BEGIN
    ensure_rot_key(keyId);
    API_END("GenerateLinearTransformRotationKey")
}
ArrayResultByte GenerateAndSerializeRotationKey(int) { unsupported("GenerateAndSerializeRotationKey"); }
void LoadRotationKey(unsigned char*, unsigned long, unsigned long) { unsupported("LoadRotationKey"); }
ArrayResultByte SerializeDiagonal(int, int) { unsupported("SerializeDiagonal"); }
void LoadPlaintextDiagonal(unsigned char*, unsigned long, int, unsigned long) { unsupported("LoadPlaintextDiagonal"); }
void RemovePlaintextDiagonals(int) { unsupported("RemovePlaintextDiagonals"); }
void RemoveRotationKeys(void) { unsupported("RemoveRotationKeys"); }

// ------------------------------------------------------------- bootstrapping --
void NewBootstrapper(int*, int, int) { unsupported("NewBootstrapper"); }
int Bootstrap(int, int) { unsupported("Bootstrap"); }
void DeleteBootstrappers(void) {}

}  // extern "C"
