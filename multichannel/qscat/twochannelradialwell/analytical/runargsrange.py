from runargscommon import *

tcp_range = argparse.ArgumentParser(parents=[parentArgs], add_help=False)
tcp_range.add_argument("eneStart_", help="Energy Start",type=QSfloat)
tcp_range.add_argument("eneEnd_", help="Energy End",type=QSfloat)
tcp_range.add_argument("eneComplex_", help="Energy Complex Offset",type=QSfloat)
tcp_range.add_argument("eneSteps_", help="Energy Steps",type=int)
