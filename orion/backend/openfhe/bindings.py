import os
import ctypes

from orion.backend.lattigo.bindings import (
    LattigoLibrary, LattigoFunction, ArrayResultDouble)


class OpenFHELibrary(LattigoLibrary):
    """OpenFHE backend.

    The shared library (built from orion/backend/openfhe/src with CMake)
    exports the same C ABI as the Lattigo library (see
    orion/backend/c_api/orion_backend.h), so all ctypes declarations and the
    Python layers above are shared with the Lattigo backend. Only the library
    that is loaded differs. Set ORION_OPENFHE_LIB to load it from elsewhere.

    Bootstrapping needs extra levels at the top of OpenFHE's modulus chain, so
    it must be enabled before the scheme is created. It is enabled when the
    configuration has a `boot_params` section (LogP different from the
    regular LogP); ORION_OPENFHE_BOOTSTRAP=0/1 overrides this.
    """
    LIB_NAME = "openfhe-linux.so"

    # In the save/load I/O modes the shim keeps the raw float32 diagonals (the ones `save_transforms`
    # writes and `load_transforms` reads) and encodes them while evaluating, so the evaluator has no
    # per-plaintext files to write or read. (Encoded diagonals are ~50x larger than raw ones.)
    plaintexts_on_disk = False

    def _load_library(self):
        path = os.environ.get("ORION_OPENFHE_LIB") or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), self.LIB_NAME)
        try:
            return ctypes.CDLL(path)
        except OSError as e:
            raise RuntimeError(
                f"Failed to load the OpenFHE backend library '{path}': {e}. "
                f"Build it with `cmake -S orion/backend/openfhe -B build && "
                f"cmake --build build` (needs an OpenFHE installation) and "
                f"copy {self.LIB_NAME} next to orion/backend/openfhe/bindings.py.")

    def setup_scheme(self, orion_params):
        enable = orion_params.get_boot_logp() != orion_params.get_logp()
        override = os.environ.get("ORION_OPENFHE_BOOTSTRAP")
        if override is not None:
            enable = override not in ("0", "", "false", "False")

        configure = self.lib.OrionOpenFHEConfigureBootstrapping
        configure.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        configure.restype = None
        configure(int(enable), 3, 3)

        super().setup_scheme(orion_params)

    def setup_poly_evaluator(self):
        super().setup_poly_evaluator()

        # Minimax sign coefficients are pure math (a multi-precision Remez
        # algorithm) that only the Lattigo library implements. Reuse it for
        # this one function; it does not touch any scheme state.
        self._lattigo_for_minimax = None
        backend = self

        class _Minimax:
            def __call__(self, *args):
                if backend._lattigo_for_minimax is None:
                    try:
                        lgo = LattigoLibrary()
                    except RuntimeError as e:
                        raise RuntimeError(
                            "Generating minimax sign polynomials (ReLU, sign) on the "
                            "OpenFHE backend currently needs the Lattigo library "
                            "(orion/backend/lattigo/lattigo-linux.so, built with Go) "
                            f"for the coefficient computation: {e}")
                    fn = LattigoFunction(
                        lgo.lib.GenerateMinimaxSignCoeffs,
                        argtypes=[
                            ctypes.POINTER(ctypes.c_int), ctypes.c_int,
                            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
                        ],
                        restype=ArrayResultDouble)
                    backend._lattigo_for_minimax = (lgo, fn)
                return backend._lattigo_for_minimax[1](*args)

        self.GenerateMinimaxSignCoeffs = _Minimax()
