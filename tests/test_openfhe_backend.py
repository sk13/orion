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
    assert torch.allclose(x, y, atol=1e-2)


def test_rotation_matches_lattigo_convention(scheme):
    x = torch.rand(8192, dtype=torch.float64)
    ct, _ = _roundtrip(scheme, x)
    r = scheme.evaluator.rotate(ct.ids[0], 5, in_place=False)
    out = orion.backend.python.tensors.CipherTensor(scheme, r, x.shape).decrypt().decode()
    assert torch.allclose(out, torch.roll(x, -5), atol=1e-2)


def test_mult_rescale_and_add(scheme):
    x = torch.rand(8192, dtype=torch.float64)
    ct, _ = _roundtrip(scheme, x)
    sq = (ct * ct).decrypt().decode()
    assert torch.allclose(sq, x * x, atol=1e-2)
    s = (ct + ct).decrypt().decode()
    assert torch.allclose(s, 2 * x, atol=1e-2)
