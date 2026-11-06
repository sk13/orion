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
    // I/O modes (save/load): where each diagonal lives in `terms`, the element parameters of Q_l*P and a
    // plaintext carrying all metadata but no polynomial, used to rebuild a diagonal from its serialized form.
    std::unordered_map<int, std::pair<size_t, size_t>> diagLoc;
    std::shared_ptr<ILDCRTParams<DCRTPoly::Integer>> params;
    Plaintext shell;
};

struct PolyEntry { std::vector<double> coeffs; bool chebyshev = false; };

struct Scheme {
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> kp;
    uint32_t depth = 0;      // total depth of the modulus chain (Lattigo level convention: level = towers-1)
    uint32_t userDepth = 0;  // max level visible to Orion (= len(LogQ)-1); the rest is bootstrapping headroom
    bool bootEnabled = false;
    int off = 0;             // Orion level l <-> internal level l+off (one hidden bottom limb when bootstrapping)
    std::vector<uint32_t> levelBudget{3, 3};
    std::set<int> bootSlots;
    uint32_t slots = 0;
    uint32_t scaleBits = 0;  // p, plaintext scaling exponent (internal scale 2^p)
    double pyDefaultScale = 0;  // Orion's default scale 2^LogScale, mapped to 2^p when they differ
    std::vector<uint64_t> Q, P;
    std::set<int> rotKeys;   // rotation amounts with generated keys
    std::set<int> transientKeys;     // rotation keys loaded from disk, dropped again by RemoveRotationKeys
    std::set<uint32_t> permAuto;     // automorphism indices that must never be dropped (bootstrapping keys)
    bool relinGenerated = false;

    void reset() {
        cc = nullptr;
        kp = KeyPair<DCRTPoly>();
        depth = userDepth = slots = scaleBits = 0;
        pyDefaultScale = 0;
        off = 0;
        bootSlots.clear();
        Q.clear(); P.clear(); rotKeys.clear(); transientKeys.clear(); permAuto.clear();
        relinGenerated = false;
    }
} S;

struct BootConfig { bool enabled = false; uint32_t enc = 3, dec = 3; } g_boot;

Table<PtEntry> pts;
Table<CT> cts;
Table<LinTransform> lts;
Table<PolyEntry> polys;

// Orion passes its default scale 2^LogScale to mean "the default scale"; map it to the internal one.
double map_scale(unsigned long scale) {
    double s = static_cast<double>(scale);
    return (S.pyDefaultScale != 0 && s == S.pyDefaultScale) ? std::pow(2.0, S.scaleBits) : s;
}

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

// ---- raw (de)serialization of polynomials, used by the save/load I/O modes ----------------
// A polynomial is stored as: u32 towers, u32 ring dimension, u8 format, u64 fingerprint of the
// moduli, then the coefficients of every tower. The fingerprint catches data that was written
// for a different parameter set (OpenFHE picks the primes deterministically from the parameters).
using ElemParams = ILDCRTParams<DCRTPoly::Integer>;
static_assert(sizeof(NativeInteger) == sizeof(uint64_t), "NativeInteger must be a plain 64-bit word");

uint64_t params_fingerprint(const ElemParams& ps) {
    uint64_t h = 1469598103934665603ULL;  // FNV-1a
    auto mix = [&](uint64_t v) { for (int i = 0; i < 8; i++) { h ^= (v >> (8 * i)) & 0xff; h *= 1099511628211ULL; } };
    mix(ps.GetCyclotomicOrder());
    for (const auto& t : ps.GetParams()) mix(t->GetModulus().ConvertToInt());
    return h;
}

void put_raw(std::string& out, const void* data, size_t n) { out.append(static_cast<const char*>(data), n); }
template <typename T> void put(std::string& out, T v) { put_raw(out, &v, sizeof(T)); }

struct Reader {
    const unsigned char* cur;
    const unsigned char* end;
    void need(size_t n) const {
        if (static_cast<size_t>(end - cur) < n) throw std::runtime_error("serialized data is truncated");
    }
    template <typename T> T get() { need(sizeof(T)); T v; std::memcpy(&v, cur, sizeof(T)); cur += sizeof(T); return v; }
};

void write_poly(std::string& out, const DCRTPoly& p) {
    const auto& towers = p.GetAllElements();
    if (towers.empty()) throw std::runtime_error("cannot serialize an empty polynomial");
    const uint32_t n = static_cast<uint32_t>(towers[0].GetLength());
    put<uint32_t>(out, static_cast<uint32_t>(towers.size()));
    put<uint32_t>(out, n);
    put<uint8_t>(out, static_cast<uint8_t>(p.GetFormat()));
    put<uint64_t>(out, params_fingerprint(*p.GetParams()));
    out.reserve(out.size() + towers.size() * n * sizeof(uint64_t));
    for (const auto& t : towers) put_raw(out, &t.GetValues()[0], n * sizeof(uint64_t));
}

// Reads a polynomial written by write_poly; `candidates` are the element parameters it may live over.
DCRTPoly read_poly(Reader& r, const std::vector<std::shared_ptr<ElemParams>>& candidates) {
    const uint32_t towers = r.get<uint32_t>(), n = r.get<uint32_t>();
    const auto fmt = static_cast<Format>(r.get<uint8_t>());
    const uint64_t fp = r.get<uint64_t>();
    std::shared_ptr<ElemParams> ps;
    for (const auto& c : candidates)
        if (c && towers == c->GetParams().size() && n == c->GetRingDimension() && fp == params_fingerprint(*c)) {
            ps = c;
            break;
        }
    if (!ps)
        throw std::runtime_error("serialized data was written for a different parameter set (" +
                                 std::to_string(towers) + " towers, ring dimension " + std::to_string(n) + ")");
    r.need(static_cast<size_t>(towers) * n * sizeof(uint64_t));
    DCRTPoly poly(ps, fmt, true);
    for (uint32_t i = 0; i < towers; i++) {
        NativeVector v(n, ps->GetParams()[i]->GetModulus());
        std::memcpy(&v[0], r.cur, static_cast<size_t>(n) * sizeof(uint64_t));
        r.cur += static_cast<size_t>(n) * sizeof(uint64_t);
        poly.GetAllElements()[i].SetValues(std::move(v), fmt);
    }
    return poly;
}

ArrayResultByte byte_array(const std::string& s) {
    ArrayResultByte r{nullptr, 0};
    if (s.empty()) return r;
    r.Data = static_cast<char*>(std::malloc(s.size()));
    std::memcpy(r.Data, s.data(), s.size());
    r.Length = s.size();
    return r;
}

// Element parameters a key polynomial can live over: Q*P, the public-key basis or Q alone.
std::vector<std::shared_ptr<ElemParams>> key_params() {
    auto cp = std::dynamic_pointer_cast<CryptoParametersRNS>(S.cc->GetCryptoParameters());
    return {cp->GetParamsQP(), cp->GetParamsPK(), cp->GetElementParams()};
}

// ---- rotation keys -------------------------------------------------------------------------
using AutoKeyMap = std::map<uint32_t, EvalKey<DCRTPoly>>;

uint32_t auto_index(int rot) { return S.cc->FindAutomorphismIndex(static_cast<uint32_t>(rot)); }

// Automorphism keys of the current secret key (nullptr while there are none).
AutoKeyMap* auto_keys() {
    auto& all = CryptoContextImpl<DCRTPoly>::GetAllEvalAutomorphismKeys();
    auto it = all.find(S.kp.secretKey->GetKeyTag());
    return it == all.end() ? nullptr : it->second.get();
}

bool key_present(int rot) {
    auto* m = auto_keys();
    return m && m->count(auto_index(rot));
}

void ensure_rot_key(int rot) {
    rot = normRot(rot);
    if (rot == 0 || S.rotKeys.count(rot)) return;
    S.cc->EvalAtIndexKeyGen(S.kp.secretKey, {rot});
    S.rotKeys.insert(rot);
    if (prof_on() && S.rotKeys.size() % 50 == 0) std::fprintf(stderr, "[KEYS] %zu rotation keys\n", S.rotKeys.size());
}

LinTransform& lt_with_diag(int tid, int diagIdx, std::pair<size_t, size_t>& loc) {
    auto& lt = lts.get(tid);
    auto it = lt.diagLoc.find(normRot(diagIdx));
    if (it == lt.diagLoc.end()) throw std::runtime_error("linear transform has no diagonal " + std::to_string(diagIdx));
    loc = it->second;
    return lt;
}

// Rotation by +r moves slot r+i to slot i (same convention as Lattigo).
CT rotate(const CT& ct, int amount) {
    int r = normRot(amount);
    if (r == 0) return ct->Clone();
    ensure_rot_key(r);
    return S.cc->EvalRotate(ct, r);
}

// ---- ciphertext helpers used by polynomial evaluation and bootstrapping ----------------

CT rescale_ct(const CT& ct) {
    int lvl = level_of(ct);
    if (lvl < 1) throw std::runtime_error("rescale: no level left");
    double s = scale_of(ct) / static_cast<double>(S.Q[lvl]);
    CT out = S.cc->ModReduce(ct);
    out->SetScalingFactor(s);
    return out;
}

CT align_level(const CT& ct, int level) {
    int diff = level_of(ct) - level;
    if (diff < 0) throw std::runtime_error("cannot raise ciphertext level");
    return diff == 0 ? ct : S.cc->LevelReduce(ct, nullptr, diff);
}

// ct * c where the constant is encoded at scale targetScale/scale(ct), so that the product has
// exactly the scale targetScale (the multiplication is not rescaled).
CT mul_const_scale(const CT& ct, double c, double targetScale) {
    double constScale = targetScale / scale_of(ct);
    CT out = S.cc->EvalMult(ct, c * constScale / std::pow(2.0, S.scaleBits));
    out->SetScalingFactor(targetScale);
    return out;
}

// ct + c, with c encoded at the exact scale of ct.
CT add_const_exact(const CT& ct, double c) {
    CT out = ct->Clone();
    auto& cv = out->GetElements();
    const uint32_t N = cv[0].GetLength();
    long double v = std::fabs(static_cast<long double>(c)) * static_cast<long double>(scale_of(ct));
    unsigned __int128 mag = static_cast<unsigned __int128>(v + 0.5L);
    DCRTPoly elems(cv[0].GetParams(), Format::COEFFICIENT, true);
    const uint32_t towers = elems.GetNumOfElements();
    for (uint32_t i = 0; i < towers; i++) {
        auto element = cv[0].GetElementAtIndex(i);
        auto modulus = element.GetModulus();
        NativeVector vec(N, modulus);
        uint64_t q = modulus.ConvertToInt();
        uint64_t r = static_cast<uint64_t>(mag % static_cast<unsigned __int128>(q));
        if (c < 0 && r != 0) r = q - r;
        vec[0] = NativeInteger(r);
        element.SetValues(std::move(vec), Format::COEFFICIENT);
        elems.SetElementAtIndex(i, std::move(element));
    }
    elems.SetFormat(Format::EVALUATION);
    cv[0] += elems;
    return out;
}

int ceil_log2(int d) { int m = 0; while ((1 << m) < d) m++; return m; }

// Depth-optimal polynomial evaluation (Paterson-Stockmeyer style, as in Lattigo): intermediate
// results stay "unrescaled" (scale ~ 2^(2p)); a node of the recursion is evaluated as
//     rescale(q) * T_k + r
// at the level where T_k lives, so that the whole polynomial of degree d consumes exactly
// bits(d) = ceil(log2(d+1)) levels, with every constant encoded at ~2^p precision.
CT evaluate_polynomial(const CT& x, const PolyEntry& poly, double outScale) {
    std::vector<double> coeffs = poly.coeffs;
    while (coeffs.size() > 1 && coeffs.back() == 0.0) coeffs.pop_back();
    const int d = static_cast<int>(coeffs.size()) - 1;
    const int L = level_of(x);
    const int D = d == 0 ? 1 : 64 - __builtin_clzll(static_cast<unsigned long long>(d));
    if (D > L) throw std::runtime_error("polynomial of degree " + std::to_string(d) + " needs " + std::to_string(D) +
                                        " levels but the input has level " + std::to_string(L));
    const bool cheb = poly.chebyshev;
    constexpr int LEAF_MAX = 16;

    std::map<int, CT> T;
    T[1] = x;
    std::function<const CT&(int)> power = [&](int i) -> const CT& {
        auto it = T.find(i);
        if (it != T.end()) return it->second;
        int a = (i + 1) / 2, b = i / 2;
        CT Ta = power(a), Tb = power(b);
        int lvl = std::min(level_of(Ta), level_of(Tb));
        CT A = align_level(Ta, lvl), B = align_level(Tb, lvl);
        CT prod = S.cc->EvalMult(A, B);
        double sp = scale_of(A) * scale_of(B);
        prod->SetScalingFactor(sp);
        if (cheb) {  // T_{a+b} = 2 T_a T_b - T_{|a-b|}
            prod = S.cc->EvalAdd(prod, prod);
            prod->SetScalingFactor(sp);
            if (a == b) {
                prod = add_const_exact(prod, -1.0);
            } else {
                CT sub = mul_const_scale(align_level(power(a - b), lvl), 1.0, sp);
                prod = S.cc->EvalSub(prod, sub);
                prod->SetScalingFactor(sp);
            }
        }
        return T[i] = rescale_ct(prod);
    };

    auto leaf = [&](const std::vector<double>& cf, int lvl, double sc) -> CT {
        CT res;
        for (size_t i = 1; i < cf.size(); i++) {
            if (cf[i] == 0.0) continue;
            CT t = mul_const_scale(align_level(power(static_cast<int>(i)), lvl), cf[i], sc);
            if (res) { res = S.cc->EvalAdd(res, t); res->SetScalingFactor(sc); } else res = std::move(t);
        }
        if (!res) res = mul_const_scale(align_level(power(1), lvl), 0.0, sc);
        if (cf[0] != 0.0) res = add_const_exact(res, cf[0]);
        return res;
    };

    std::function<CT(std::vector<double>, int, double)> node = [&](std::vector<double> cf, int lvl, double sc) -> CT {
        while (cf.size() > 1 && cf.back() == 0.0) cf.pop_back();
        const int dd = static_cast<int>(cf.size()) - 1;
        const int m = dd >= 2 ? ceil_log2(dd) : 0;
        const bool leafFeasible = lvl <= L - m;
        if (dd <= 1 || (leafFeasible && dd <= LEAF_MAX)) return leaf(cf, lvl, sc);
        if (lvl + 1 > L) throw std::runtime_error("polynomial evaluation: not enough levels");
        int k = 1;
        while (2 * k < dd) k *= 2;  // largest power of two below the degree
        std::vector<double> q(dd - k + 1, 0.0), r;
        std::vector<double> work = cf;
        if (cheb) {
            for (int i = dd; i > k; i--) {
                q[i - k] += 2.0 * work[i];
                work[std::abs(2 * k - i)] -= work[i];
                work[i] = 0.0;
            }
            q[0] += work[k];
        } else {
            for (int i = dd; i >= k; i--) q[i - k] = work[i];
        }
        r.assign(work.begin(), work.begin() + k);
        const CT& X = power(k);
        double sx = scale_of(X);
        double sb = sc * static_cast<double>(S.Q[lvl + 1]) / sx;
        CT b = rescale_ct(node(q, lvl + 1, sb));
        CT prod = S.cc->EvalMult(b, align_level(X, lvl));
        prod->SetScalingFactor(scale_of(b) * sx);
        CT a = node(r, lvl, sc);
        CT out = S.cc->EvalAdd(prod, a);
        out->SetScalingFactor(sc);
        return out;
    };

    const int rootLevel = L - D + 1;
    CT root = node(coeffs, rootLevel, outScale * static_cast<double>(S.Q[rootLevel]));
    return rescale_ct(root);
}

}  // namespace

extern "C" {

// ---------------------------------------------------------------- scheme --
void NewScheme(int logN, int* logQ, int lenQ, int* logP, int lenP, int logScale, int h,
               char* ringType, char* keysPath, char* ioMode) {
    (void)h; (void)keysPath; (void)ioMode;
    API_BEGIN
    if (lenQ < 1) throw std::runtime_error("LogQ must not be empty");
    if (!g_boot.enabled)
        for (int i = 1; i < lenQ; i++)
            if (logQ[i] != logScale)
                throw std::runtime_error("OpenFHE backend (FIXEDMANUAL) requires every LogQ entry after the first "
                                         "to equal LogScale");
    bool ci = std::string(ringType) != "standard";
    uint32_t ringDim = 1u << (logN + (ci ? 1 : 0));

    uint32_t bootDepth = 0;
    std::vector<uint32_t> levelBudget{g_boot.enc, g_boot.dec};
    if (g_boot.enabled) bootDepth = FHECKKSRNS::GetBootstrapDepth(levelBudget, SPARSE_TERNARY);

    CCParams<CryptoContextCKKSRNS> p;
    p.SetSecurityLevel(HEStd_NotSet);
    p.SetRingDim(ringDim);
    p.SetBatchSize(ringDim / 2);
    const int off = g_boot.enabled ? 1 : 0;
    p.SetMultiplicativeDepth(lenQ - 1 + off + bootDepth);
    // OpenFHE's bootstrapping is only accurate with a scale of about 2^52 or more and a 60-bit q0
    // (it fails outright at 2^40). When bootstrapping is enabled the backend therefore uses its own
    // internal scale and maps Orion's default scale onto it. 2^59 gives ~17 bits of precision after
    // a bootstrap, 2^53 only ~11 bits, which is too coarse for Orion's composite-sign ReLU.
    int scaleBits = logScale;
    if (g_boot.enabled) {
        const char* e = std::getenv("ORION_OPENFHE_BOOT_SCALE");
        scaleBits = std::max(logScale, e ? std::atoi(e) : 59);
    }
    p.SetScalingModSize(scaleBits);
    // OpenFHE's fixed-scaling bootstrapping requires q0 within a few bits of the scale.
    uint32_t firstMod = logQ[0];
    if (g_boot.enabled) {
        const char* e = std::getenv("ORION_OPENFHE_BOOT_Q0");
        firstMod = std::max<uint32_t>(firstMod, e ? static_cast<uint32_t>(std::atoi(e)) : 60);
    }
    p.SetFirstModSize(firstMod);
    p.SetScalingTechnique(FIXEDMANUAL);
    p.SetKeySwitchTechnique(HYBRID);
    p.SetSecretKeyDist(g_boot.enabled ? SPARSE_TERNARY : UNIFORM_TERNARY);
    if (lenP > 0 && !g_boot.enabled) p.SetNumLargeDigits((lenQ + lenP - 1) / lenP);

    S.reset();
    S.cc = GenCryptoContext(p);
    S.cc->Enable(PKE);
    S.cc->Enable(KEYSWITCH);
    S.cc->Enable(LEVELEDSHE);
    S.cc->Enable(ADVANCEDSHE);
    if (g_boot.enabled) S.cc->Enable(FHE);
    S.bootEnabled = g_boot.enabled;
    S.levelBudget = levelBudget;
    S.off = off;
    S.userDepth = lenQ - 1 + off;  // max internal level visible to Orion
    S.depth = lenQ - 1 + off + bootDepth;
    S.slots = ringDim / 2;
    S.scaleBits = scaleBits;
    S.pyDefaultScale = std::pow(2.0, logScale);

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
int GetPlaintextLevel(int id) { return level_of(pts.get(id).pt) - S.off; }
int GetCiphertextLevel(int id) { return level_of(cts.get(id)) - S.off; }
int GetPlaintextSlots(int) { return static_cast<int>(S.slots); }
int GetCiphertextSlots(int) { return static_cast<int>(S.slots); }
int GetCiphertextDegree(int id) { return static_cast<int>(cts.get(id)->NumberCiphertextElements()) - 1; }
ArrayResultUInt64 GetModuliChain(void) {
    // Orion's Q[l] is the prime dropped when rescaling at its level l (internal level l+off).
    return u64_array(std::vector<uint64_t>(S.Q.begin() + S.off, S.Q.end()));
}
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
ArrayResultByte SerializeSecretKey(void) {
    API_BEGIN
    std::string buf;
    write_poly(buf, S.kp.secretKey->GetPrivateElement());
    return byte_array(buf);
    API_END("SerializeSecretKey")
}
// Restores the secret key written by SerializeSecretKey. No public key is derived from it (OpenFHE
// has no API for that): encryption then uses the secret key, which is equivalent for Orion.
void LoadSecretKey(unsigned char* data, unsigned long len) {
    API_BEGIN
    Reader r{data, data + len};
    DCRTPoly poly = read_poly(r, key_params());
    auto sk = std::make_shared<PrivateKeyImpl<DCRTPoly>>(S.cc);
    sk->SetPrivateElement(std::move(poly));
    S.kp.secretKey = sk;
    S.kp.publicKey = nullptr;
    API_END("LoadSecretKey")
}

// ------------------------------------------------------------------ encoder --
void NewEncoder(void) {}

int Encode(float* values, int len, int level, unsigned long scale) {
    API_BEGIN PROF("Encode");
    std::vector<std::complex<double>> v(len);
    for (int i = 0; i < len; i++) v[i] = values[i];
    PtEntry e;
    e.scale = map_scale(scale);
    e.pt = encode(v, true, level + S.off, e.scale);
    e.vals = v; e.real = true;
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
    e.scale = map_scale(scale);
    e.pt = encode(v, false, level + S.off, e.scale);
    e.vals = v; e.real = false;
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
    CT ct = S.kp.publicKey ? S.cc->Encrypt(S.kp.publicKey, e.pt) : S.cc->Encrypt(S.kp.secretKey, e.pt);
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
    // Eagerly generate the power-of-two rotation keys for small rings; for big rings they
    // are generated on first use (keys span the whole modulus chain and are large).
    if (S.slots <= 8192)
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
    return unary(cid, inPlace, [&](const CT& ct0) {
        auto& e = pts.get(pid);
        // Like Lattigo, operate at the lower of the two levels.
        int lvl = std::min(level_of(ct0), level_of(e.pt));
        CT ct = align_level(ct0, lvl);
        Plaintext pt = pt_at_level(e, lvl);
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
int EvaluatePolynomial(int ctId, int polyId, unsigned long outScale) {
    API_BEGIN PROF("EvaluatePolynomial");
    double sc = outScale ? map_scale(outScale) : std::pow(2.0, S.scaleBits);
    CT out = evaluate_polynomial(cts.get(ctId), polys.get(polyId), sc);
    return push_ct(out, scale_of(out));
    API_END("EvaluatePolynomial")
}
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
    const std::string mode(ioMode);
    if (mode != "none" && mode != "save" && mode != "load")
        throw std::runtime_error("linear transform: unknown io_mode '" + mode + "'");
    (void)bsgsRatio;
    level += S.off;  // Orion level -> internal level
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
    for (size_t j = 0; j < jobs.size(); j++) {
        where[j] = {jobs[j].g, fill[jobs[j].g]++};
        lt.diagLoc[idx[jobs[j].k]] = where[j];
    }
    lt.params = paramsQlP;
    if (mode != "none") {
        // metadata-only plaintext to rebuild diagonals from their serialized polynomial later
        lt.shell = FHECKKSRNS::MakeAuxPlaintext(*S.cc, paramsQlP, std::vector<std::complex<double>>(n, 0.0), 1,
                                                S.depth - level, n);
        lt.shell->GetElement<DCRTPoly>() = DCRTPoly();
    }
    if (mode == "load") return lts.add(std::move(lt));  // diagonals are loaded from disk before each evaluation

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
    if (level_of(cts.get(cid)) < lt.level)
        throw std::runtime_error("linear transform: input level " + std::to_string(level_of(cts.get(cid))) +
                                 " < transform level " + std::to_string(lt.level));
    for (const auto& giant : lt.terms)
        for (const auto& term : giant)
            if (!term.second) throw std::runtime_error("linear transform: diagonals are not loaded (io_mode load/save)");
    const CT ct = align_level(cts.get(cid), lt.level);  // a higher-level input is dropped to the transform's level
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
// Generates the key, hands it to Python as bytes and forgets it again unless it was there before.
ArrayResultByte GenerateAndSerializeRotationKey(int keyId) {
    API_BEGIN
    const int rot = normRot(keyId);
    if (rot == 0) throw std::runtime_error("no rotation key for the identity rotation");
    const uint32_t ai = auto_index(rot);
    const bool had = key_present(rot);
    if (!had) S.cc->EvalAtIndexKeyGen(S.kp.secretKey, {rot});
    auto* m = auto_keys();
    auto it = m ? m->find(ai) : AutoKeyMap::iterator();
    if (!m || it == m->end()) throw std::runtime_error("rotation key was not generated");
    const EvalKey<DCRTPoly>& key = it->second;
    std::string buf;
    const auto& A = key->GetAVector();
    const auto& B = key->GetBVector();
    put<uint32_t>(buf, static_cast<uint32_t>(A.size()));
    for (const auto& a : A) write_poly(buf, a);
    for (const auto& b : B) write_poly(buf, b);
    if (!had) m->erase(ai);
    return byte_array(buf);
    API_END("GenerateAndSerializeRotationKey")
}

// Installs a serialized rotation key until RemoveRotationKeys is called.
void LoadRotationKey(unsigned char* data, unsigned long len, unsigned long keyId) {
    API_BEGIN
    const int rot = normRot(static_cast<int>(keyId));
    if (rot == 0 || key_present(rot)) return;  // already available (e.g. generated on demand)
    Reader r{data, data + len};
    const auto ps = key_params();
    const uint32_t digits = r.get<uint32_t>();
    std::vector<DCRTPoly> A, B;
    for (uint32_t i = 0; i < digits; i++) A.push_back(read_poly(r, ps));
    for (uint32_t i = 0; i < digits; i++) B.push_back(read_poly(r, ps));
    auto key = std::make_shared<EvalKeyRelinImpl<DCRTPoly>>(S.cc);
    key->SetAVector(std::move(A));
    key->SetBVector(std::move(B));
    key->SetKeyTag(S.kp.secretKey->GetKeyTag());
    auto one = std::make_shared<AutoKeyMap>();
    (*one)[auto_index(rot)] = key;
    CryptoContextImpl<DCRTPoly>::InsertEvalAutomorphismKey(one, S.kp.secretKey->GetKeyTag());
    S.rotKeys.insert(rot);
    S.transientKeys.insert(rot);
    API_END("LoadRotationKey")
}

// Drops every rotation key installed by LoadRotationKey (bootstrapping keys stay).
void RemoveRotationKeys(void) {
    API_BEGIN
    auto* m = auto_keys();
    for (int rot : S.transientKeys) {
        const uint32_t ai = auto_index(rot);
        if (S.permAuto.count(ai)) continue;
        if (m) m->erase(ai);
        S.rotKeys.erase(rot);
    }
    S.transientKeys.clear();
    API_END("RemoveRotationKeys")
}

// Returns the encoded diagonal and frees it: it is reloaded right before the transform is evaluated.
ArrayResultByte SerializeDiagonal(int tid, int diagIdx) {
    API_BEGIN
    std::pair<size_t, size_t> loc;
    auto& lt = lt_with_diag(tid, diagIdx, loc);
    Plaintext& pt = lt.terms[loc.first][loc.second].second;
    if (!pt) throw std::runtime_error("diagonal " + std::to_string(diagIdx) + " is not in memory");
    std::string buf;
    write_poly(buf, pt->GetElement<DCRTPoly>());
    pt = nullptr;
    return byte_array(buf);
    API_END("SerializeDiagonal")
}

void LoadPlaintextDiagonal(unsigned char* data, unsigned long len, int tid, unsigned long diagIdx) {
    API_BEGIN
    std::pair<size_t, size_t> loc;
    auto& lt = lt_with_diag(tid, static_cast<int>(diagIdx), loc);
    if (!lt.shell) throw std::runtime_error("linear transform was not created in an I/O mode");
    Reader r{data, data + len};
    auto pt = std::make_shared<CKKSPackedEncoding>(static_cast<const CKKSPackedEncoding&>(*lt.shell));
    pt->GetElement<DCRTPoly>() = read_poly(r, {lt.params});
    lt.terms[loc.first][loc.second].second = pt;
    API_END("LoadPlaintextDiagonal")
}

void RemovePlaintextDiagonals(int tid) {
    API_BEGIN
    for (auto& giant : lts.get(tid).terms)
        for (auto& term : giant) term.second = nullptr;
    API_END("RemovePlaintextDiagonals")
}

// ------------------------------------------------------------- bootstrapping --
// Backend-specific extension (not exported by the Lattigo library): must be called before
// NewScheme to reserve bootstrapping levels at the top of the modulus chain.
void OrionOpenFHEConfigureBootstrapping(int enable, int levelBudgetEnc, int levelBudgetDec) {
    g_boot.enabled = enable != 0;
    g_boot.enc = levelBudgetEnc > 0 ? levelBudgetEnc : 3;
    g_boot.dec = levelBudgetDec > 0 ? levelBudgetDec : 3;
}

void NewBootstrapper(int*, int, int numSlots) {
    API_BEGIN PROF("NewBootstrapper");
    if (!S.bootEnabled) throw std::runtime_error("bootstrapping was not enabled before NewScheme");
    if (S.bootSlots.count(numSlots)) return;
    // The first modulus is larger than the scale; OpenFHE's correction factor must cover the difference.
    S.cc->EvalBootstrapSetup(S.levelBudget, {0, 0}, numSlots);
    S.cc->EvalBootstrapKeyGen(S.kp.secretKey, numSlots);
    if (auto* m = auto_keys())  // bootstrapping needs these for good: RemoveRotationKeys must not drop them
        for (const auto& kv : *m) S.permAuto.insert(kv.first);
    S.bootSlots.insert(numSlots);
    API_END("NewBootstrapper")
}

// Same semantics as the Lattigo backend: the input holds `numSlots` values in its first slots
// (zeros elsewhere); the output has the values replicated with period numSlots over all slots.
int Bootstrap(int id, int numSlots) {
    API_BEGIN PROF("Bootstrap");
    if (!S.bootSlots.count(numSlots)) throw std::runtime_error("no bootstrapper for this slot count");
    CT x = cts.get(id)->Clone();
    const double sIn = scale_of(x);
    if (x->GetNoiseScaleDeg() != 1) throw std::runtime_error("bootstrap: ciphertext must be rescaled first");
    int lvl = level_of(x);
    if (lvl < 1) throw std::runtime_error("bootstrap: input needs at least two limbs");
    if (lvl > 1) x = S.cc->LevelReduce(x, nullptr, lvl - 1);  // OpenFHE consumes one level before raising the modulus
    for (int k = numSlots; k < static_cast<int>(S.slots); k *= 2) {  // fold the zero-padded vector into its periodic form
        CT r = rotate(x, k);
        x = S.cc->EvalAdd(x, r);
    }
    x->SetSlots(numSlots);
    if (prof_on()) std::fprintf(stderr, "[BTP] input towers=%d deg=%zu slots=%u\n", level_of(x) + 1, x->GetNoiseScaleDeg(), x->GetSlots());
    CT out = S.cc->EvalBootstrap(x);
    out->SetSlots(S.slots);
    if (prof_on()) std::fprintf(stderr, "[BTP] output towers=%d deg=%zu\n", level_of(out) + 1, out->GetNoiseScaleDeg());
    while (out->GetNoiseScaleDeg() > 1) out = S.cc->ModReduce(out);
    int ol = level_of(out);
    if (ol > static_cast<int>(S.userDepth)) out = S.cc->LevelReduce(out, nullptr, ol - static_cast<int>(S.userDepth));
    return push_ct(out, sIn);
    API_END("Bootstrap")
}
void DeleteBootstrappers(void) {}

}  // extern "C"
