# Test ase.utils.eos without ASE (eoswoase)
# This test will run if numpy is in the default sys.path
# because PYTHONPATH is overwritten under os.system call!

import os

import ase

dir = os.path.abspath(os.path.dirname(ase.__file__))
utilsdir = os.path.join(dir, 'utils')

sjeos = os.path.join(utilsdir, 'sjeos.py')
assert os.system("PYTHONPATH=" + utilsdir + " python " + sjeos) == 0
