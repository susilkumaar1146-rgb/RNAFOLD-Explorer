"""
Shared pytest setup.

* Adds the project folder to sys.path so tests can `import modules...`
  when pytest is run from the project root.
* Provides a brute-force structure enumerator used to PROVE that both
  dynamic programming algorithms find the true optimum on short sequences.
"""

import os
import sys
from functools import lru_cache

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from modules.pairing import MIN_HAIRPIN, can_pair  # noqa: E402


def enumerate_structures(sequence):
    """
    Return EVERY valid nested structure of `sequence` as a list of pair
    tuples. Exponential - only for sequences of about 16 nt or fewer.

    Recurrence: in the segment i..j, base i is either unpaired, or it pairs
    with some k and splits the segment into (i+1..k-1) and (k+1..j).
    """
    n = len(sequence)

    @lru_cache(maxsize=None)
    def structures(i, j):
        if i > j:
            return ((),)
        result = list(structures(i + 1, j))                 # i unpaired
        for k in range(i + MIN_HAIRPIN + 1, j + 1):
            if can_pair(sequence[i], sequence[k]):
                for inside in structures(i + 1, k - 1):
                    for outside in structures(k + 1, j):
                        result.append(((i, k),) + inside + outside)
        return tuple(result)

    return [sorted(s) for s in structures(0, n - 1)]
