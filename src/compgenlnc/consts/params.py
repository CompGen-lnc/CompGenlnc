import sys
from multiprocessing import cpu_count

# molecule type
LNCRNA = "lnc"
MIRNA = "mir"
PRECURSOR = "pre"
MATURE = "mat"

# record mode
SEQ_MODE = "seq"
SS_MODE = "ss"

# workers
WORKERS = cpu_count() - 2 if "pytest" not in sys.modules else 0
