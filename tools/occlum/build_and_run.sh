#!/bin/bash
# Everything (python env, Orion code, Lattigo .so, data, model) is baked into the enclave's protected FS.
set -e
export PATH=/usr/local/go/bin:$PATH GOCACHE=/tmp/gocache GOPATH=/tmp/gopath HOME=/root
SCRIPT=${1:-examples/run_mlp.py}; CWD=${2:-/work/examples}
W=/tmp/w; I=$W/occlum_instance; mkdir -p $W
if [ ! -d $I/build ] || [ "$REBUILD" = 1 ]; then
  rm -rf $W/stage $I; mkdir -p $W/stage/model
  rsync -a --exclude .venv --exclude .git /src/ $W/orion-src/
  (cd $W/orion-src && /opt/python-orion/bin/python tools/build_lattigo.py | tail -1)
  for d in orion examples configs data; do cp -r $W/orion-src/$d $W/stage/$d; done
  cp $W/orion-src/model_state/helrm.pth $W/stage/model/; cp /encl/launch.py $W/stage/
  cd $W && occlum new occlum_instance >/dev/null && cd $I
  copy_bom -f /encl/orion.yaml --root image --include-dir /opt/occlum/etc/template
  jq '.resource_limits.user_space_size = "'${USS:-40GB}'" |
    .resource_limits.kernel_space_heap_size = "512MB" |
    .resource_limits.max_num_of_threads = 256 |
    .env.default += ["PYTHONHOME=/opt/python-orion","LD_LIBRARY_PATH=/opt/python-orion/lib","HOME=/root","USER=root","LOGNAME=root","OMP_NUM_THREADS=8","MKL_NUM_THREADS=8","TORCHINDUCTOR_CACHE_DIR=/tmp/ti"]' Occlum.json > o.json && mv o.json Occlum.json
  occlum build ${OCCLUM_BUILD_ARGS}
fi
cd $I
time occlum run /bin/python3 /work/launch.py /work/$SCRIPT $CWD
