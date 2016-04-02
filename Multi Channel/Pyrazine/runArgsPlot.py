from runBase import *

args_range = argparse.ArgumentParser(description="Pyrazine Fit Plot", parents=[parentArgs], add_help=False)
args_range.add_argument("plotStart_", help="Plot Start Energy", type=float)
args_range.add_argument("plotEnd_", help="Plot End Energy", type=float)
args_range.add_argument("plotComplex_", help="Plot Complex Energy Offset", type=float)
args_range.add_argument("steps_", help="Number of steps", type=int)
args = args_range.parse_args()