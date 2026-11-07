import os, sys, runpy
script, cwd = sys.argv[1], sys.argv[2]
sys.path.insert(0, "/work"); os.chdir(cwd)
sys.argv = [script]
runpy.run_path(script, run_name="__main__")
