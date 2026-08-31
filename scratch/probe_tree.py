import numpy as np, time
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform
from lrg_eegfc.utils.metrics.tree import pair_merge_level, pair_tree_octave, merge_height_pair_counts

rng = np.random.default_rng(0)
n = 118
X = rng.random((n, 5))
D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
np.fill_diagonal(D, 0)
Z = linkage(squareform(D, checks=False), method="average")
lev = pair_merge_level(Z)
oct_ = pair_tree_octave(Z)
cnt = merge_height_pair_counts(Z)
C = cophenet(Z)
print("M", lev.size, n * (n - 1) // 2, "levels", lev.min(), lev.max())
print("octave counts", np.bincount(oct_))
assert np.allclose(C, Z[lev - 1, 2]), "merge level inconsistent with cophenet"
assert np.array_equal(cnt, np.bincount(lev, minlength=n)[1:]), "counts mismatch"
print("OK: pair_merge_level agrees with cophenet; counts agree; sum", cnt.sum())
t = time.time()
for _ in range(50):
    pair_tree_octave(Z)
print("pair_tree_octave %.2f ms" % ((time.time() - t) / 50 * 1000))
