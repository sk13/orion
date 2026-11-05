# Running Orion in an Occlum (Intel SGX) enclave

Everything (Python env, Orion code, Lattigo `.so`, configs, MNIST data, HELRM model)
is packed into the enclave's protected file system.

```
docker build -t orion-occlum tools/occlum
mkdir -p /tmp/occ_w
# first run builds the enclave image (~15 min); later runs reuse it (REBUILD=1 forces a rebuild)
docker run --rm --entrypoint bash -e REBUILD=1 \
  --device /dev/sgx_enclave --device /dev/sgx_provision \
  -v $PWD:/src:ro -v $PWD/tools/occlum:/encl:ro -v /tmp/occ_w:/tmp/w \
  orion-occlum /encl/build_and_run.sh examples/run_mlp.py /work/examples
```

Arguments: `<script relative to repo> <working dir inside enclave>`
(`run_lola.py` and `run_helrm.py` expect `/work`, `run_mlp.py` and `run_resnet.py` expect `/work/examples`).
Needs `data/` (MNIST/CIFAR) in the repo dir. `USS` sets the enclave user space (default 40GB;
HELRM needs more than 12GB). `launch.py` overrides `platform.system()` because Occlum's
`uname` does not report "Linux".
