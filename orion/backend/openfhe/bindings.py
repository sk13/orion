import os
import ctypes

from orion.backend.lattigo.bindings import LattigoLibrary


class OpenFHELibrary(LattigoLibrary):
    """OpenFHE backend.

    The shared library (built from orion/backend/openfhe/src with CMake)
    exports the same C ABI as the Lattigo library (see
    orion/backend/c_api/orion_backend.h), so all ctypes declarations and the
    Python layers above are shared with the Lattigo backend. Only the library
    that is loaded differs. Set ORION_OPENFHE_LIB to load it from elsewhere.
    """
    LIB_NAME = "openfhe-linux.so"

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
