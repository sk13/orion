"""Smoke tests for the OpenFHE backend (skipped when the library is not built)."""
import os

import pytest
import torch

import orion

LIB = os.environ.get("ORION_OPENFHE_LIB") or os.path.join(
    os.path.dirname(orion.__file__), "backend", "openfhe", "openfhe-linux.so")
pytestmark = pytest.mark.skipif(not os.path.exists(LIB), reason="OpenFHE backend not built")

CONFIG = {
    "ckks_params": {"LogN": 13, "LogQ": [29, 26, 26, 26, 26, 26], "LogP": [29, 29],
                    "LogScale": 26, "H": 8192, "RingType": "ConjugateInvariant"},
    "orion": {"margin": 2, "embedding_method": "hybrid", "backend": "openfhe",
              "fuse_modules": True, "debug": False, "diags_path": "", "keys_path": "",
              "io_mode": "none"},
}


@pytest.fixture(scope="module")
def scheme():
    s = orion.init_scheme(CONFIG)
    yield s
    s.delete_scheme()


def _roundtrip(scheme, x, level=None):
    ct = scheme.encrypt(scheme.encode(x, level))
    return ct, ct.decrypt().decode()


def test_encode_encrypt_decrypt(scheme):
    x = torch.rand(8192, dtype=torch.float64) - 0.5
    _, y = _roundtrip(scheme, x)
    assert torch.allclose(x, y, atol=3e-2)


def test_rotation_matches_lattigo_convention(scheme):
    x = torch.rand(8192, dtype=torch.float64)
    ct, _ = _roundtrip(scheme, x)
    r = scheme.evaluator.rotate(ct.ids[0], 5, in_place=False)
    out = orion.backend.python.tensors.CipherTensor(scheme, r, x.shape).decrypt().decode()
    assert torch.allclose(out, torch.roll(x, -5), atol=3e-2)


def test_mult_rescale_and_add(scheme):
    x = torch.rand(8192, dtype=torch.float64)
    ct, _ = _roundtrip(scheme, x)
    sq = (ct * ct).decrypt().decode()
    assert torch.allclose(sq, x * x, atol=3e-2)
    s = (ct + ct).decrypt().decode()
    assert torch.allclose(s, 2 * x, atol=3e-2)


def _poly_check(scheme, coeffs, chebyshev):
    import numpy as np
    pe = scheme.poly_evaluator
    x = torch.rand(8192, dtype=torch.float64) * 2 - 1
    ct = scheme.encrypt(scheme.encode(x, level=5))
    if chebyshev:
        pid = pe.generate_chebyshev(coeffs)
        ref = np.polynomial.chebyshev.chebval(x.numpy(), coeffs)
    else:
        pid = pe.generate_monomial(coeffs[::-1])
        ref = np.polynomial.polynomial.polyval(x.numpy(), coeffs)
    out = pe.evaluate_polynomial(ct, pid)
    # Orion plans exactly bits(degree) = ceil(log2(degree+1)) levels for a polynomial.
    assert scheme.backend.GetCiphertextLevel(out.ids[0]) == 5 - (len(coeffs) - 1).bit_length()
    return np.abs(out.decrypt().decode().numpy() - ref).max()


@pytest.mark.parametrize("degree", [1, 2, 3, 7, 8, 15])
@pytest.mark.parametrize("chebyshev", [True, False])
def test_polynomial_evaluation_depth_and_values(scheme, degree, chebyshev):
    import numpy as np
    coeffs = (np.random.default_rng(degree).standard_normal(degree + 1) / (np.arange(degree + 1) + 1)).tolist()
    assert _poly_check(scheme, coeffs, chebyshev) < 0.3


@pytest.mark.skipif(not os.environ.get("ORION_OPENFHE_TEST_BOOTSTRAP"),
                    reason="slow (N=2^16); set ORION_OPENFHE_TEST_BOOTSTRAP=1")
def test_sparse_bootstrap_matches_lattigo_semantics():
    cfg = {
        "ckks_params": {"LogN": 16, "LogQ": [55] + [40] * 10, "LogP": [61, 61, 61],
                        "LogScale": 40, "H": 192, "RingType": "standard"},
        "boot_params": {"LogP": [61] * 8},
        "orion": {"margin": 2, "embedding_method": "hybrid", "backend": "openfhe",
                  "fuse_modules": True, "debug": False, "diags_path": "", "keys_path": "",
                  "io_mode": "none"},
    }
    from orion.backend.python.tensors import CipherTensor
    s = orion.init_scheme(cfg)
    try:
        slots = 8
        s.bootstrapper.generate_bootstrapper(slots)
        x = torch.zeros(s.params.get_slots(), dtype=torch.float64)
        x[:slots] = torch.rand(slots, dtype=torch.float64) - 0.5
        ct = s.encrypt(s.encode(x, level=1))
        out = CipherTensor(s, s.bootstrapper.bootstrap(ct.ids[0], slots), x.shape)
        assert out.level() == s.params.get_max_level()
        y = out.decrypt().decode()
        # values come back replicated with period `slots`, like with Lattigo
        for k in range(0, 4 * slots, slots):
            assert torch.allclose(y[k:k + slots], x[:slots], atol=1e-3)
    finally:
        s.delete_scheme()
