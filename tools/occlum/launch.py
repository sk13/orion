import os, sys, runpy, platform
platform.system = lambda: "Linux"   # Occlum's uname reports a different sysname
script, cwd = sys.argv[1], sys.argv[2]
sys.path.insert(0, "/work"); os.chdir(cwd)
sys.argv = [script]
runpy.run_path(script, run_name="__main__")
