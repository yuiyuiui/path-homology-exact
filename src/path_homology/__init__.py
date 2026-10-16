"""Exact non-regular first path homology of finite loopless digraphs."""

from .core import Digraph, Homology, Limits, ResourceLimitError, betti1, compute_h1
from .verify import verify_certificate

__version__ = "0.1.0"
__all__ = [
    "Digraph",
    "Homology",
    "Limits",
    "ResourceLimitError",
    "betti1",
    "compute_h1",
    "verify_certificate",
]
