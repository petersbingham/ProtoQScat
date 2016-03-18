from DoubleChannel import *
import argparse

parentArgs = argparse.ArgumentParser(add_help=False)
parentArgs.add_argument("r0_", help="Width of Well", type=float)
parentArgs.add_argument("v1_", help="Potential 1", type=float)
parentArgs.add_argument("v2_", help="Potential 2", type=float)
parentArgs.add_argument("t1_", help="Threshold 1", type=float)
parentArgs.add_argument("t2_", help="Threshold 2", type=float)
parentArgs.add_argument("lam_", help="Coupling Factor", type=float)