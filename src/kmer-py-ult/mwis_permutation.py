"""
Maximum Weight Increasing Subsequence — permutation variant.

Given:
  x  : strictly increasing array of n distinct values
  y  : a permutation of x
  w  : keep-weight for each index

Goal: choose a subset of indices S (kept in both arrays) that maximises
sum(w[i] for i in S), subject to y[i] being strictly increasing over S.
(x[i] is automatically increasing since x is already sorted.)

Key insight: since y is a permutation of sorted x, rank(y[i]) is just the
position of y[i] in x — no separate coordinate compression needed.
The problem is a 1D Maximum Weight Increasing Subsequence on y.

Algorithm (O(n log n)):
  Process indices left to right.
  dp[i] = w[i] + max(dp[j]  for j < i  where y[j] < y[i])
         = w[i] + BIT.query(rank(y[i]) - 1)
  Then update BIT at rank(y[i]) with dp[i].

Contrast with the unweighted case:
  Without weights, patience sort (bisect) suffices — O(n log n), simpler.
  With weights, we need a max-BIT because we optimise total weight, not count.
"""

import bisect
import sys
from typing import List, Tuple


# ---------------------------------------------------------------------------
# Shared: Max Fenwick tree
# ---------------------------------------------------------------------------

class _MaxBIT:
    def __init__(self, size: int):
        self.n = size
        self.val = [0] * (size + 1)
        self.src = [-1] * (size + 1)

    def query(self, pos: int) -> Tuple[int, int]:
        """Max (val, src) in prefix [1, pos]."""
        best_val, best_src = 0, -1
        while pos > 0:
            if self.val[pos] > best_val:
                best_val = self.val[pos]
                best_src = self.src[pos]
            pos -= pos & (-pos)
        return best_val, best_src

    def update(self, pos: int, val: int, src: int) -> None:
        while pos <= self.n:
            if val > self.val[pos]:
                self.val[pos] = val
                self.src[pos] = src
            pos += pos & (-pos)


# ---------------------------------------------------------------------------
# Unweighted — O(n log n) patience sort (baseline, no weights)
# ---------------------------------------------------------------------------

def min_removals_unweighted(x: List[int], y: List[int]) -> Tuple[int, List[int], List[int]]:
    """
    Patience sort / binary search.
    tails[k] = smallest tail value among all increasing subsequences of length k+1.
    Works only because all y values are distinct (permutation of sorted x).
    """
    n = len(x)
    tails: List[int] = []
    tail_src: List[int] = []   # which index produced each tail
    parent = [-1] * n
    pile = [-1] * n

    for i, v in enumerate(y):
        p = bisect.bisect_left(tails, v)
        if p == len(tails):
            tails.append(v)
            tail_src.append(i)
        else:
            tails[p] = v
            tail_src[p] = i
        pile[i] = p
        parent[i] = tail_src[p - 1] if p > 0 else -1

    kept = []
    cur = tail_src[-1]
    while cur != -1:
        kept.append(cur)
        cur = parent[cur]
    kept.reverse()

    removed = [i for i in range(n) if i not in set(kept)]
    return n - len(kept), kept, removed


# ---------------------------------------------------------------------------
# Weighted — O(n log n) max-BIT
# ---------------------------------------------------------------------------

def max_weight_kept(
    x: List[int], y: List[int], w: List[int]
) -> Tuple[int, List[int], List[int]]:
    """
    Maximum Weight Increasing Subsequence on y using a max-Fenwick tree.

    rank(y[i]) = bisect_left(x, y[i]) + 1   (1-indexed position in sorted x)
    Since y is a permutation of x, this is always an exact hit — O(log n) lookup.

    Returns (max_kept_weight, kept_indices, removed_indices).
    """
    n = len(x)
    bit = _MaxBIT(n)
    dp = [0] * n
    parent = [-1] * n

    for i in range(n):
        r = bisect.bisect_left(x, y[i]) + 1          # rank in x, 1-indexed
        prev_val, prev_src = bit.query(r - 1)         # best predecessor with y[j] < y[i]
        dp[i] = w[i] + prev_val
        parent[i] = prev_src
        bit.update(r, dp[i], i)

    best_end = max(range(n), key=lambda i: dp[i])

    kept = []
    cur = best_end
    while cur != -1:
        kept.append(cur)
        cur = parent[cur]
    kept.reverse()

    removed = [i for i in range(n) if i not in set(kept)]
    return dp[best_end], kept, removed


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    examples = [
        {
            "desc": "Two local conflicts in y",
            "x": [1, 2, 3, 4, 5, 6, 7],
            "y": [1, 2, 4, 3, 5, 7, 6],
        },
        {
            "desc": "Descending run [5,4,3] in y",
            "x": [1, 2, 3, 4, 5, 6, 7],
            "y": [1, 2, 5, 4, 3, 6, 7],
        },
    ]

    weight_cases = [
        ("uniform",           [1, 1, 1, 1, 1, 1, 1]),
        ("prefer high index", [1, 1, 1, 10, 1, 10, 1]),
        ("prefer low index",  [1, 1, 10, 1, 1, 1, 10]),
    ]

    for ex in examples:
        x, y = ex["x"], ex["y"]
        print(f"=== {ex['desc']} ===")
        print(f"x = {x}")
        print(f"y = {y}")

        # Unweighted baseline
        rm, k0, r0 = min_removals_unweighted(x, y)
        print(f"  [unweighted]  min_removals={rm}  kept={k0}  removed={r0}")

        # Weighted
        for label, w in weight_cases:
            best, kept, removed = max_weight_kept(x, y, w)
            print(f"  [{label}]  w={w}")
            print(f"    max_kept_weight={best}  kept={kept}  removed={removed}")
        print()
