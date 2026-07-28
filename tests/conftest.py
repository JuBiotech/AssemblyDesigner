import os
import sys
import warnings

import matplotlib  # must be imported before non-import statements

# Headless plotting for matplotlib (do this immediately after import)
matplotlib.use("Agg")

# Make the repo root importable so `assembly_designer` can be found
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Silence the pandas usecols FutureWarning #anoying
warnings.filterwarnings(
    "ignore",
    message="Defining usecols with out of bounds indices is deprecated",
    category=FutureWarning,
)
