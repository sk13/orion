/*
 * Orion backend C ABI.
 *
 * This header is the contract between Orion's Python layer (ctypes, see
 * orion/backend/lattigo/bindings.py) and an FHE library backend. The Lattigo
 * backend (Go, cgo) exports exactly these symbols; the OpenFHE backend
 * (C++, orion/backend/openfhe) implements the same functions, so the Python
 * side is shared by both.
 *
 * Conventions
 *  - Plaintexts, ciphertexts, linear transforms and polynomials live in the
 *    backend and are referred to by integer handles. "...New" functions
 *    return a fresh handle; the others modify the object in place and return
 *    the same handle.
 *  - Levels use Lattigo's convention: a ciphertext with k RNS limbs has
 *    level k-1 and a freshly encoded top-level plaintext has the maximum
 *    level (len(LogQ)-1).
 *  - Scales are positive integers (the double scale rounded to uint64).
 *  - Rotation-key ids returned by GetLinearTransformRotationKeys and taken
 *    by GenerateLinearTransformRotationKey / GenerateAndSerializeRotationKey
 *    / LoadRotationKey are opaque to Python. (Lattigo: Galois elements,
 *    OpenFHE: rotation amounts.)
 *  - Arrays returned by value ("ArrayResult*") are malloc'ed by the backend
 *    and must be released with FreeCArray.
 */
#ifndef ORION_BACKEND_H
#define ORION_BACKEND_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { int*                Data; unsigned long Length; } ArrayResultInt;
typedef struct { float*              Data; unsigned long Length; } ArrayResultFloat;
typedef struct { double*             Data; unsigned long Length; } ArrayResultDouble;
typedef struct { unsigned long long* Data; unsigned long long Length; } ArrayResultUInt64;
typedef struct { char*               Data; unsigned long Length; } ArrayResultByte;

/* ---- scheme ------------------------------------------------------------ */
void NewScheme(int logN, int* logQ, int lenQ, int* logP, int lenP, int logScale,
               int h, char* ringType, char* keysPath, char* ioMode);
void DeleteScheme(void);
void FreeCArray(void* ptr);

/* ---- tensors ----------------------------------------------------------- */
void DeletePlaintext(int id);
void DeleteCiphertext(int id);
unsigned long GetPlaintextScale(int id);
unsigned long GetCiphertextScale(int id);
void SetPlaintextScale(int id, unsigned long scale);
void SetCiphertextScale(int id, unsigned long scale);
int GetPlaintextLevel(int id);
int GetCiphertextLevel(int id);
int GetPlaintextSlots(int id);
int GetCiphertextSlots(int id);
int GetCiphertextDegree(int id);
ArrayResultUInt64 GetModuliChain(void);
ArrayResultUInt64 GetAuxModuliChain(void);
ArrayResultInt GetLivePlaintexts(void);
ArrayResultInt GetLiveCiphertexts(void);

/* ---- key generation ---------------------------------------------------- */
void NewKeyGenerator(void);
void GenerateSecretKey(void);
void GeneratePublicKey(void);
void GenerateRelinearizationKey(void);
void GenerateEvaluationKeys(void);
ArrayResultByte SerializeSecretKey(void);
void LoadSecretKey(unsigned char* data, unsigned long len);

/* ---- encoder ----------------------------------------------------------- */
void NewEncoder(void);
int Encode(float* values, int len, int level, unsigned long scale);
ArrayResultFloat Decode(int plaintextId);
int EncodeComplex(float* interleaved, int len, int level, unsigned long scale);
ArrayResultFloat DecodeComplex(int plaintextId);

/* ---- encryptor --------------------------------------------------------- */
void NewEncryptor(void);
void NewDecryptor(void);
int Encrypt(int plaintextId);
int Decrypt(int ciphertextId);

/* ---- evaluator --------------------------------------------------------- */
void NewEvaluator(void);
void AddRotationKey(int rotation);
int Negate(int ct);
int Rotate(int ct, int amount);
int RotateNew(int ct, int amount);
int Rescale(int ct);
int RescaleNew(int ct);
int AddScalar(int ct, float s);
int AddScalarNew(int ct, float s);
int SubScalar(int ct, float s);
int SubScalarNew(int ct, float s);
int MulScalarInt(int ct, int s);
int MulScalarIntNew(int ct, int s);
int MulScalarFloat(int ct, float s);
int MulScalarFloatNew(int ct, float s);
int AddPlaintext(int ct, int pt);
int AddPlaintextNew(int ct, int pt);
int SubPlaintext(int ct, int pt);
int SubPlaintextNew(int ct, int pt);
int MulPlaintext(int ct, int pt);
int MulPlaintextNew(int ct, int pt);
int AddCiphertext(int ct0, int ct1);
int AddCiphertextNew(int ct0, int ct1);
int SubCiphertext(int ct0, int ct1);
int SubCiphertextNew(int ct0, int ct1);
int MulRelinCiphertext(int ct0, int ct1);
int MulRelinCiphertextNew(int ct0, int ct1);

/* ---- polynomial evaluator ---------------------------------------------- */
void NewPolynomialEvaluator(void);
int GenerateMonomial(float* coeffs, int len);
int GenerateChebyshev(float* coeffs, int len);
int EvaluatePolynomial(int ct, int poly, unsigned long outScale);
ArrayResultDouble GenerateMinimaxSignCoeffs(int* degrees, int lenDegrees,
                                            int prec, int logalpha, int logerr, int debug);

/* ---- linear transforms ------------------------------------------------- */
void NewLinearTransformEvaluator(void);
int GenerateLinearTransform(int* diagIdxs, int lenIdxs, float* diagData, int lenData,
                            int level, float bsgsRatio, char* ioMode);
int EvaluateLinearTransform(int transformId, int ctId);
void DeleteLinearTransform(int id);
ArrayResultInt GetLinearTransformRotationKeys(int transformId);
void GenerateLinearTransformRotationKey(int keyId);
ArrayResultByte GenerateAndSerializeRotationKey(int keyId);
void LoadRotationKey(unsigned char* data, unsigned long len, unsigned long keyId);
ArrayResultByte SerializeDiagonal(int transformId, int diagIdx);
void LoadPlaintextDiagonal(unsigned char* data, unsigned long len, int transformId,
                           unsigned long diagIdx);
void RemovePlaintextDiagonals(int transformId);
void RemoveRotationKeys(void);

/* ---- bootstrapping ----------------------------------------------------- */
void NewBootstrapper(int* logPs, int lenLogPs, int numSlots);
int Bootstrap(int ct, int numSlots);
void DeleteBootstrappers(void);

#ifdef __cplusplus
}
#endif
#endif /* ORION_BACKEND_H */
