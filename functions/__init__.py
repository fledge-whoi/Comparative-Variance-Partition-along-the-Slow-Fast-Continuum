"""Functions called by main.py and plotresults_articlescript.py.

One module per function of the original MATLAB project (same names), so that
functions/<name>.py corresponds to functions/<name>.m.
"""
from .combinations import combinations
from .popHeteromatfact import popHeteromatfact
from .det_outputfactorial import det_outputfactorial
from .lifetimeoutput import lifetimeoutput
from .reprooutputfactorial import reprooutputfactorial
from .varcalc import varcalc
from .vardecomposition import vardecomposition
from .minusparamfinder_ConstantG import minusparamfinder_ConstantG
from .minusparamfinder_ConstantLambda import minusparamfinder_ConstantLambda

__all__ = [
    "combinations", "popHeteromatfact", "det_outputfactorial", "lifetimeoutput",
    "reprooutputfactorial", "varcalc", "vardecomposition",
    "minusparamfinder_ConstantG", "minusparamfinder_ConstantLambda",
]
