import sys
from multiprocessing import cpu_count

# molecule type
LNCRNA_TYPE = "lnc"
MIRNA_TYPE = "mir"

# record mode
SEQ_MODE = "seq"
SS_MODE = "ss"

# workers
WORKERS = cpu_count() - 2 if "pytest" not in sys.modules else 0
