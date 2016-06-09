import argparse

parentArgs = argparse.ArgumentParser(description="Numercal Data Fit Routine", add_help=False)
parentArgs.add_argument("subN_", help="Sub Set Number", type=int, nargs='?', default=-1)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=-1)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=-1)