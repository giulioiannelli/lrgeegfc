"""Heat-kernel multiscale cophenetic analysis of a weighted graph.

The LRG diffusion primitive, factored so a single eigendecomposition drives an
entire diffusion-time (``tau``) sweep. On a **connected** graph (e.g. a
percolation backbone) the combinatorial Laplacian heat kernel ``K=e^{-tau L}``
integrates all path lengths; sweeping the dimensionless scale ``s = tau*lambda_max``
from ``1`` (fastest mode) upward resolves structure at successively coarser
scales that a single adjacency threshold cannot.

Conventions
-----------
- Combinatorial Laplacian ``L = D - W`` (never the random-walk normalisation;
  see ``feedback_combinatorial_laplacian_only``).
- ``K = V exp(-tau Lambda) V^T >= 0`` exactly (``-L`` is Metzler); tiny negative
  round-off / underflow is floored at ``+RHO_FLOOR`` before inverting.
- Communication distance ``D_ij = (1 - delta_ij) / rho_ij``; ``rho = K / Tr K``.
- Cophenetic distances from average-linkage (UPGMA) UPGMA of ``D``.
- Normalised von Neumann entropy ``Shat = -Tr[rho ln rho] / ln N in [0, 1]``;
  specific heat ``C = -dShat / d log10 tau``.
- Cross-phase trace ``rho_sym`` = symmetric split-half Spearman of cophenetic
  reorganisations (removes the arbitrary half assignment of bare ``rho_split``).

Nothing here is dataset-specific: it takes matrices / eigendecompositions and
returns arrays. FC loading and cohort constants stay script-side.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.signal import find_peaks
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

__all__ = [
    "RHO_FLOOR",
    "laplacian_eig",
    "entropy_specific_heat",
    "specific_heat_peaks",
    "scale_grid",
    "cophenetic_at_tau",
    "cophenetic_at_scale",
    "linkage_at_scale",
    "rho_sym",
    "rho_sym_over_scales",
    "effective_cluster_count",
    "communication_neighbourhood_size",
    "scale_resolution",
    "CROSS_PHASE_ROLES",
    "cross_phase_functionals",
    "cross_phase_functionals_over_scales",
]

RHO_FLOOR = 1e-30          # heat-kernel underflow floor (K>=0 exactly; far pairs)


# --------------------------------------------------------------------------- #
# spectrum
# --------------------------------------------------------------------------- #
def laplacian_eig(W: NDArray) -> tuple[NDArray, NDArray]:
    """Eigendecomposition of the combinatorial Laplacian ``L = D - W``.

    Returns ``(eigenvalues, eigenvectors)`` with eigenvalues clipped at 0
    (numerical noise on the trivial ``lambda_1 = 0`` mode).
    """
    W = np.asarray(W, float)
    deg = W.sum(1)
    ev, V = np.linalg.eigh(np.diag(deg) - W)
    return np.maximum(ev, 0.0), V


def _entropy_nats(ev: NDArray, tau: float) -> float:
    """von Neumann entropy ``-Tr[rho ln rho]`` in nats via the Boltzmann form."""
    log_boltz = -tau * ev
    log_Z = np.logaddexp.reduce(log_boltz)
    log_p = log_boltz - log_Z
    p = np.exp(log_p)
    return float(-np.sum(p * np.where(p > 1e-300, log_p, 0.0)))


# --------------------------------------------------------------------------- #
# entropy / specific heat with the "Shat from 1 to 0" range rule
# --------------------------------------------------------------------------- #
def entropy_specific_heat(ev: NDArray, n: int = 300,
                          s_hi: float = 0.999, s_lo: float = 0.001) -> dict:
    """Normalised entropy ``Shat(tau)`` and specific heat ``C(tau)``.

    The ``tau`` grid is auto-ranged so ``Shat`` spans ~``[s_lo, s_hi]`` -- i.e.
    the full 0->1 sigmoid -- regardless of the spectrum (robust when
    ``lambda_2 -> 0`` makes ``1/lambda_2`` explode). ``Shat -> 1`` as ``tau->0``
    (maximally mixed) and ``Shat -> 0`` as ``tau->inf`` (collapse to the
    connected ground state).

    Returns dict: ``tau, s (=tau*lambda_max), log10_tau, S (=Shat), C, tau_min
    (=1/lambda_max), tau_max (=1/lambda_2), lambda_max, lambda_2, N``.
    """
    ev = np.maximum(np.asarray(ev, float), 0.0)
    N = ev.size
    logN = np.log(N)
    lam_max = ev[-1]
    lam2 = ev[1] if N > 1 else np.nan
    tau0 = 1.0 / lam_max

    t = tau0                                   # bracket small-tau (high entropy) end
    for _ in range(400):
        if _entropy_nats(ev, t) / logN >= s_hi:
            break
        t *= 0.5
    tau_lo = t
    t = tau0                                   # bracket large-tau (low entropy) end
    for _ in range(400):
        if _entropy_nats(ev, t) / logN <= s_lo:
            break
        t *= 1.5
    tau_hi = t

    tau = np.logspace(np.log10(tau_lo), np.log10(tau_hi), n)
    S = np.array([_entropy_nats(ev, x) / logN for x in tau])
    log10_tau = np.log10(tau)
    C = -np.gradient(S, log10_tau)
    return dict(tau=tau, s=tau * lam_max, log10_tau=log10_tau, S=S, C=C,
                tau_min=tau0, tau_max=(1.0 / lam2 if lam2 > 1e-30 else np.nan),
                lambda_max=float(lam_max), lambda_2=float(lam2), N=int(N))


def specific_heat_peaks(res: dict, rel_prominence: float = 0.05) -> dict:
    """Interior peaks of ``C(tau)`` -> characteristic scale(s).

    1 peak = collapsed single scale (Wigner-semicircle regime); >1 = multiscale
    ladder (Villegas 2025 App. A/F). Takes the dict from
    :func:`entropy_specific_heat`. Returns ``n_peaks``, peak ``s`` values, peak
    ``tau`` values, and dominant-peak ``s``.
    """
    C = res["C"]
    cmax = float(C.max()) if C.size else 0.0
    prom = max(rel_prominence * cmax, 1e-9)
    idx, _ = find_peaks(C, prominence=prom)
    if idx.size == 0 and C.size:
        idx = np.array([int(np.argmax(C))])
    s_peaks = res["s"][idx] if idx.size else np.array([])
    tau_peaks = res["tau"][idx] if idx.size else np.array([])
    dom_s = float(res["s"][idx[np.argmax(C[idx])]]) if idx.size else np.nan
    return dict(n_peaks=int(idx.size), s_peaks=s_peaks, tau_peaks=tau_peaks,
                dominant_s=dom_s)


def scale_grid(ev: NDArray, n: int = 16, s_floor: float = 0.02,
               s_min: float = 1.0, s_cap: float = 200.0) -> NDArray:
    """Dimensionless scale grid ``s = tau*lambda_max`` for the trace arc.

    From ``s_min=1`` (``tau=1/lambda_max``) up to the scale where ``Shat`` drops
    to ``s_floor`` (graph collapsed to one cluster), capped at ``s_cap``. This is
    the finite proxy for ``tau_max`` (``1/lambda_min = 1/0`` is undefined).
    Log-spaced, ``n`` points. Range is set by the reference spectrum ``ev``
    (typically phase A) so a whole cell shares one comparable ``s`` axis.
    """
    ev = np.maximum(np.asarray(ev, float), 0.0)
    N = ev.size
    logN = np.log(N)
    lam_max = ev[-1]
    tau = 1.0 / lam_max
    s = 1.0
    for _ in range(400):
        if _entropy_nats(ev, tau) / logN <= s_floor or s >= s_cap:
            break
        tau *= 1.3
        s = tau * lam_max
    s_max = min(max(s, s_min * 2), s_cap)
    return np.logspace(np.log10(s_min), np.log10(s_max), n)


# --------------------------------------------------------------------------- #
# cophenetic distances at a diffusion time
# --------------------------------------------------------------------------- #
def _linkage_at_tau(ev: NDArray, V: NDArray, tau: float,
                    rho_floor: float = RHO_FLOOR) -> NDArray:
    """UPGMA linkage ``Z`` of the diffusion communication tree at time ``tau``.

    ``rho = e^{-tau L}/Tr``, floored at ``rho_floor`` (underflow only), ``D=1/rho``,
    average linkage. Shared body of :func:`cophenetic_at_tau` and
    :func:`linkage_at_scale` (single source for the tree geometry).
    """
    rho = (V * np.exp(-tau * ev)) @ V.T
    rho /= np.trace(rho)
    rho = np.where(rho > rho_floor, rho, rho_floor)
    T = 1.0 / rho
    np.fill_diagonal(T, 0.0)
    T = np.maximum(T, T.T)
    return linkage(squareform(T, checks=False), method="average")


def cophenetic_at_tau(ev: NDArray, V: NDArray, tau: float,
                      rho_floor: float = RHO_FLOOR) -> NDArray:
    """Condensed cophenetic distances of the UPGMA tree at diffusion time ``tau``.

    ``rho = e^{-tau L}/Tr``, floored at ``rho_floor`` (underflow only), ``D=1/rho``,
    average linkage, cophenet. Identical to the established ``ultra``/``diff_coph``
    at ``tau=1/lambda_max`` (the floor never triggers on a full graph).
    """
    return cophenet(_linkage_at_tau(ev, V, tau, rho_floor))


def cophenetic_at_scale(ev: NDArray, V: NDArray, s: float,
                        rho_floor: float = RHO_FLOOR) -> NDArray:
    """:func:`cophenetic_at_tau` addressed by the dimensionless scale ``s=tau*lambda_max``."""
    return cophenetic_at_tau(ev, V, s / ev[-1], rho_floor)


def linkage_at_scale(ev: NDArray, V: NDArray, s: float,
                     rho_floor: float = RHO_FLOOR) -> NDArray:
    """UPGMA linkage ``Z`` of the diffusion tree at dimensionless scale ``s=tau*lambda_max``.

    The tree behind :func:`cophenetic_at_scale`, returned as a SciPy ``linkage``
    matrix so callers can draw the dendrogram / cut clades at scale ``s``
    (tanglegrams, chord layouts, ribbon trees, the tau-morph). ``tau = s / lambda_max``.
    """
    return _linkage_at_tau(ev, V, s / ev[-1], rho_floor)


# --------------------------------------------------------------------------- #
# how much structure a scale actually resolves (the unit of the s axis)
# --------------------------------------------------------------------------- #
def effective_cluster_count(ev: NDArray, s: float) -> float:
    """Effective number of resolved components ``N_eff(s) = exp(S_nats(tau))``.

    The exponential of the (un-normalised) von Neumann entropy of
    ``rho = e^{-tau L}/Tr``, at ``tau = s/lambda_max``. This is the standard
    entropy-based effective rank: ``N_eff = N`` when ``rho`` is maximally mixed
    (``tau -> 0``, every node its own component) and ``N_eff -> 1`` when the
    diffusion has collapsed the graph to its connected ground state
    (``tau -> inf``). Continuous, threshold-free, and defined by the spectrum
    alone.

    Its purpose is to give the dimensionless scale axis ``s`` an interpretable
    unit. A scale is never described as fine / meso / coarse without a number;
    ``N_eff(s)`` is that number (see ``feedback_no_unquantified_scale_labels``).
    """
    ev = np.maximum(np.asarray(ev, float), 0.0)
    return float(np.exp(_entropy_nats(ev, s / ev[-1])))


def communication_neighbourhood_size(ev: NDArray, V: NDArray, s: float) -> float:
    """Effective number of nodes one node communicates with at scale ``s``.

    Row ``i`` of the heat kernel, normalised to a probability vector
    ``p_ij = K_ij / sum_j K_ij``, is where the heat released at ``i`` sits at
    time ``tau = s/lambda_max``. Its participation ratio
    ``m_i = 1 / sum_j p_ij^2`` is the effective size of that neighbourhood:
    ``m_i -> 1`` when the heat has not left ``i`` and ``m_i -> N`` at
    equilibrium. The returned value is the mean of ``m_i`` over nodes.

    Threshold-free companion to :func:`effective_cluster_count`, read from the
    kernel rather than the spectrum: the two together say "at ``s`` the process
    resolves ``N_eff`` groups and each node mixes with ``m`` others". Use them
    to attach a number to any scale that is described in words.
    """
    ev = np.maximum(np.asarray(ev, float), 0.0)
    K = (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
    K = np.maximum(K, 0.0)
    row = K.sum(1, keepdims=True)
    P = np.divide(K, row, out=np.zeros_like(K), where=row > 0)
    denom = np.sum(P * P, axis=1)
    m = np.divide(1.0, denom, out=np.full(denom.shape, np.nan), where=denom > 0)
    return float(np.nanmean(m))


def scale_resolution(ev: NDArray, V: NDArray, s_grid: NDArray) -> dict:
    """Both resolution readouts over a scale grid.

    Returns ``dict(s, n_eff, m_comm, N)`` with ``n_eff`` from
    :func:`effective_cluster_count` and ``m_comm`` from
    :func:`communication_neighbourhood_size`. Attach to every reported scale so
    the ``s`` axis carries an interpretable unit rather than a bare number.

    A UPGMA-tree cut was tried as a third readout and **discarded**: on
    ``D = 1/rho`` the merge heights span many decades and every fixed or
    relative cut (at ``D = N``, at ``N/e``, at the largest log-height gap) is
    either saturated at ``N`` over the whole usable range or jumps
    non-monotonically at the coarse end. The spectral and kernel readouts above
    are monotone by construction and are the ones to report.
    """
    s_grid = np.asarray(s_grid, float)
    return dict(
        s=s_grid,
        n_eff=np.array([effective_cluster_count(ev, s) for s in s_grid]),
        m_comm=np.array([communication_neighbourhood_size(ev, V, s) for s in s_grid]),
        N=int(ev.size),
    )


# --------------------------------------------------------------------------- #
# cross-phase trace
# --------------------------------------------------------------------------- #
def rho_sym(DA: NDArray, DB: NDArray, Dt: NDArray, Dp: NDArray) -> tuple[float, float]:
    """Symmetric split-half cross-phase trace ``(rho_sym, rho_split)``.

    ``rho_sym = 1/2[ Spearman(Dt-DA, Dp-DB) + Spearman(Dt-DB, Dp-DA) ]`` on
    condensed cophenetic distances of the four phases (A,B split-half rest_pre;
    ``t`` = task_test; ``p`` = rest_post). Removes the arbitrary half assignment
    of bare ``rho_split`` (audit_149). Byte-equivalent to audit_150 ``rho_sym_split``.
    """
    r_ab, _ = spearmanr(Dt - DA, Dp - DB)
    r_ba, _ = spearmanr(Dt - DB, Dp - DA)
    return float(0.5 * (r_ab + r_ba)), float(r_ab)


def rho_sym_over_scales(eig_by_phase: dict, s_grid: NDArray,
                        rho_floor: float = RHO_FLOOR) -> NDArray:
    """``rho_sym(s)`` swept over the scale grid, one eigendecomposition per phase.

    ``eig_by_phase`` maps each phase in ``{A, B, task_test, rest_post}`` to its
    ``(eigenvalues, eigenvectors)``. Returns ``rho_sym`` for each ``s`` in
    ``s_grid`` (the cross-phase trace as a function of diffusion scale -- the
    "role of tau" curve). ``NaN`` where a cophenetic tree is degenerate.
    """
    phases = ("A", "B", "task_test", "rest_post")
    out = np.full(len(s_grid), np.nan)
    for i, s in enumerate(s_grid):
        try:
            D = {ph: cophenetic_at_scale(*eig_by_phase[ph], s, rho_floor)
                 for ph in phases}
            out[i], _ = rho_sym(D["A"], D["B"], D["task_test"], D["rest_post"])
        except Exception:
            out[i] = np.nan
    return out


# --------------------------------------------------------------------------- #
# cross-phase functionals over an arbitrary phase set
# --------------------------------------------------------------------------- #
#: Semantic role -> default phase key. The functionals below are defined on
#: ROLES, not on dataset-specific phase names, so a study with different labels
#: only has to remap this dict. ``encode`` is optional: drop it from the phase
#: set and only the ``probe``-based functional is returned, so a four-phase and
#: a five-phase caller share one code path.
CROSS_PHASE_ROLES = {
    "baseline_a": "A",
    "baseline_b": "B",
    "encode": "task_learn",
    "probe": "task_test",
    "follow": "rest_post",
}


def _average_ranks(x: NDArray) -> NDArray:
    """Average ranks (ties shared). Pearson on these equals Spearman."""
    from scipy.stats import rankdata
    return rankdata(x)


def _partial_from_corr(r_xy: float, r_xz: float, r_yz: float) -> float:
    """First-order partial correlation ``r(x, y | z)`` from the three pairwise r."""
    den = np.sqrt(max(0.0, (1.0 - r_xz ** 2) * (1.0 - r_yz ** 2)))
    return float((r_xy - r_xz * r_yz) / den) if den > 0 else float("nan")


def cross_phase_functionals(D: dict, roles: dict | None = None) -> dict:
    """Symmetric split-half cross-phase functionals from cophenetic distances.

    ``D`` maps a phase name to its condensed cophenetic distance vector (all on
    the same pair index set). ``roles`` maps the five semantic roles --
    ``baseline_a``, ``baseline_b``, ``encode``, ``probe``, ``follow`` -- to the
    phase names present; it defaults to :data:`CROSS_PHASE_ROLES`.

    With a **two-stage task** (an ``encode`` phase that establishes a structure
    and a ``probe`` phase that applies it) the reorganisation vectors are

    ``e  = D_encode - D_A``      (encoding, arm A)      ``e2 = D_encode - D_B``
    ``f  = D_probe  - D_encode`` (probe-specific; arm-invariant by construction)
    ``g  = D_probe  - D_A``      (total task change)    ``g2 = D_probe  - D_B``
    ``p  = D_follow - D_B``      (persistence)          ``p2 = D_follow - D_A``

    and the returned functionals, each symmetrised over the arbitrary A/B arm
    assignment and each **cross-baseline** (an A-referenced change is always
    paired with a B-referenced persistence, so no arm shares a baseline with
    itself -- shared-baseline pairing inflates the correlation):

    ``T_probe``        ``1/2[rho(g, p) + rho(g2, p2)]``  -- the standard trace.
    ``T_encode``       ``1/2[rho(e, p) + rho(e2, p2)]``  -- does the *encoding*
                       reorganisation persist?
    ``T_probespec``    ``1/2[rho(f, p) + rho(f, p2)]``   -- does the change the
                       probe adds *on top of* encoding persist?
    ``T_probespec_pe`` ``1/2[pr(f, p | e) + pr(f, p2 | e2)]`` -- the same with
                       the encoding component partialled out.

    ``rho`` is Spearman throughout, obtained from a single rank correlation
    matrix per call rather than by repeated pairwise ``spearmanr``. If no
    ``encode`` phase is present only ``T_probe`` is returned.

    Note: ``T_probe`` is numerically identical to :func:`rho_sym` on the same
    four phases.
    """
    roles = dict(CROSS_PHASE_ROLES if roles is None else roles)
    A, B = D[roles["baseline_a"]], D[roles["baseline_b"]]
    P, F = D[roles["probe"]], D[roles["follow"]]
    enc_key = roles.get("encode")
    has_enc = enc_key is not None and enc_key in D

    vecs = {"g": P - A, "g2": P - B, "p": F - B, "p2": F - A}
    if has_enc:
        E = D[enc_key]
        vecs.update({"e": E - A, "e2": E - B, "f": P - E})

    names = list(vecs)
    Rmat = np.corrcoef(np.vstack([_average_ranks(vecs[k]) for k in names]))
    ix = {k: i for i, k in enumerate(names)}

    def r(a, b):
        return float(Rmat[ix[a], ix[b]])

    out = {"T_probe": 0.5 * (r("g", "p") + r("g2", "p2"))}
    if has_enc:
        out["T_encode"] = 0.5 * (r("e", "p") + r("e2", "p2"))
        out["T_probespec"] = 0.5 * (r("f", "p") + r("f", "p2"))
        out["T_probespec_pe"] = 0.5 * (
            _partial_from_corr(r("f", "p"), r("f", "e"), r("p", "e"))
            + _partial_from_corr(r("f", "p2"), r("f", "e2"), r("p2", "e2"))
        )
    return out


def cross_phase_functionals_over_scales(eig_by_phase: dict, s_grid: NDArray,
                                        roles: dict | None = None,
                                        rho_floor: float = RHO_FLOOR) -> dict:
    """:func:`cross_phase_functionals` swept over a diffusion-scale grid.

    ``eig_by_phase`` maps each phase present to its ``(eigenvalues,
    eigenvectors)``. The phase set is taken as **data**: whichever phases are in
    the dict are used, so adding or removing the ``encode`` phase needs no
    change here or in any caller. Returns ``{functional: array(len(s_grid))}``,
    ``NaN`` wherever a cophenetic tree is degenerate.
    """
    roles = dict(CROSS_PHASE_ROLES if roles is None else roles)
    need = [roles[k] for k in ("baseline_a", "baseline_b", "probe", "follow")]
    missing = [ph for ph in need if ph not in eig_by_phase]
    if missing:
        raise KeyError(f"cross-phase functionals need phases {missing}")
    has_enc = roles.get("encode") in eig_by_phase
    phases = need + ([roles["encode"]] if has_enc else [])
    keys = ["T_probe"] + (["T_encode", "T_probespec", "T_probespec_pe"] if has_enc else [])
    out = {k: np.full(len(s_grid), np.nan) for k in keys}
    for i, s in enumerate(s_grid):
        try:
            D = {ph: cophenetic_at_scale(*eig_by_phase[ph], s, rho_floor) for ph in phases}
            vals = cross_phase_functionals(D, roles)
        except Exception:
            continue
        for k in keys:
            out[k][i] = vals.get(k, np.nan)
    return out
