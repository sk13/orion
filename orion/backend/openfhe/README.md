# OpenFHE backend

A second FHE runtime for Orion, selected with `backend: openfhe` in the config
(see `configs/mlp_openfhe.yml`). It implements the C ABI declared in
`orion/backend/c_api/orion_backend.h`, which is the same ABI the Lattigo (Go)
library exports, so the Python layers are shared by both backends.

## Build

Needs an OpenFHE >= 1.3 installation (developed against v1.6.1) and CMake.

```
cmake -S orion/backend/openfhe -B build -DCMAKE_PREFIX_PATH=/path/to/openfhe
cmake --build build
cp build/openfhe-linux.so orion/backend/openfhe/
```

`ORION_OPENFHE_LIB=/path/to/openfhe-linux.so` loads the library from elsewhere.

## Status

Working: encode/decode (real and complex), encrypt/decrypt, add/sub/mul with
scalars, plaintexts and ciphertexts, rotations, rescale, and BSGS linear
transforms. The MLP example runs end to end (`examples/run_mlp.py` with the
OpenFHE config).

Not implemented yet (the functions abort with a clear message): polynomial
evaluation (activations other than `Quad`), bootstrapping, minimax sign
coefficients, and the `save`/`load` I/O modes (key and diagonal serialization).

## Differences from Lattigo

* CKKS with `FIXEDMANUAL` scaling, so Orion keeps deciding when to rescale.
  OpenFHE tracks scales approximately in this mode; the backend keeps the exact
  scale of each ciphertext itself and corrects values on decode.
* Every `LogQ` entry after the first must equal `LogScale` (OpenFHE uses one
  scaling-modulus size).
* OpenFHE has no conjugate-invariant ring. `RingType: ConjugateInvariant` with
  `LogN` is emulated with a standard ring of degree `2^(LogN+1)`, which gives the
  same number of slots (`2^LogN`) and the same rotation semantics.
* The secret-key Hamming weight (`H`) is ignored (uniform ternary key).
* Rotation-key ids exposed to Python are plain rotation amounts, not Galois
  elements.
