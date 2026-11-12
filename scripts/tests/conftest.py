"""The scripts of the set import the library as `coreset` (scripts/coreset/): put scripts/ on the path."""

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SCRIPTS))
