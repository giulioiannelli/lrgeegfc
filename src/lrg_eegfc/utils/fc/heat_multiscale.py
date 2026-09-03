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
    "partial_spearman",
    "CROSS_PHASE_FUNCTIONALS",
    "CROSS_PHASE_FUNCTIONALS_DRIFT",
    "CROSS_PHASE_VECTOR_NAMES",
    "cross_phase_vectors",
    "cross_phase_rank_corr",
    "cross_phase_rank_corr_over_scales",
    "partial_from_corr_matrix",
    "cross_phase_functionals_from_corr",
    "cross_phase_functionals_from_corr_stack",
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


def partial_spearman(a: NDArray, b: NDArray, c: NDArray) -> float:
    """Partial Spearman correlation of ``a`` and ``b`` controlling for ``c``.

    ``(r_ab - r_ac r_bc) / sqrt((1 - r_ac^2)(1 - r_bc^2))`` on the rank
    correlations. ``NaN`` when the denominator vanishes.

    .. warning::
       A partial correlation is **not** automatically zero-centred under a
       no-signal input: the conditioning is entangled with the estimator, and a
       construction that manufactures shared structure between ``a``, ``b`` and
       ``c`` can return a large positive partial from data containing no effect.
       This is not hypothetical -- a windowed sham arc built entirely from
       pre-task ``rest_pre`` returned ``+0.243`` for beta (real value ``+0.091``,
       sham-above-zero ``p = 0.007``); see
       ``.agents/preprint/supplementary/S2_drift_controls.md``. Any null used with
       :data:`CROSS_PHASE_FUNCTIONALS` ``T_infspec_pe`` must therefore be
       CALIBRATED on a no-signal input first, and its p-values must not be
       reported until it is shown to sit near zero there.
    """
    r_ab, _ = spearmanr(a, b)
    r_ac, _ = spearmanr(a, c)
    r_bc, _ = spearmanr(b, c)
    den = np.sqrt(max(0.0, (1.0 - r_ac ** 2) * (1.0 - r_bc ** 2)))
    return float((r_ab - r_ac * r_bc) / den) if den > 0 else float("nan")


#: The five-phase cross-phase functionals, keyed by name. Each maps a dict of
#: condensed cophenetic distance vectors -- keys ``A``, ``B``, ``task_learn``,
#: ``task_test``, ``rest_post`` -- to a scalar. Phase labels are passed in by the
#: caller, so the same functionals apply to a real arc, a pseudo-phase arc from a
#: block permutation, or a sham arc built inside a single recording.
CROSS_PHASE_FUNCTIONALS = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")

#: Drift-controlled counterparts: the same four functionals with the task-free
#: drift direction ``d`` partialled out of both arguments (see
#: :func:`cross_phase_functionals_from_corr`).
CROSS_PHASE_FUNCTIONALS_DRIFT = ("T_test_d", "T_learn_d", "T_infspec_d",
                                 "T_infspec_ped")

#: Reorganisation vectors of a five-phase arc, in the row order used by
#: :func:`cross_phase_rank_corr`. Everything the cross-phase functionals measure
#: is a (partial) correlation between two of these, so storing their rank
#: correlation matrix stores every functional -- present and future -- exactly.
#:
#: ``d = D_baseline_b - D_baseline_a`` is the **task-free drift** vector. The two
#: baselines are contiguous halves of ONE pre-task recording, so ``d`` is the
#: reorganisation the system performs over a comparable timespan with no task in
#: it: same subject, same session, same estimator. It is the natural covariate
#: for any claim that a cross-phase change is task-driven rather than elapsed
#: time, and it costs nothing -- it is already implied by every arc.
CROSS_PHASE_VECTOR_NAMES = ("e", "e2", "f", "g", "g2", "p", "p2", "d")

#: Two naming schemes for the same four functionals arrived from two independent
#: Wave-0 implementations, verified to agree to full float precision before being
#: unified here. The role-neutral names are canonical in the library; the
#: task-specific names are canonical in the manuscript and in every script written
#: against ``05_enc_inf_arc``. Both are emitted so neither breaks.
_FUNCTIONAL_ALIASES = {
    "T_probe": "T_test",
    "T_encode": "T_learn",
    "T_probespec": "T_infspec",
    "T_probespec_pe": "T_infspec_pe",
}






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
    for _k, _v in list(out.items()):
        out.setdefault(_FUNCTIONAL_ALIASES[_k], _v)
    return out


def cross_phase_functionals_over_scales(eig_by_phase: dict, s_grid: NDArray,
                                        roles: dict | None = None,
                                        rho_floor: float = RHO_FLOOR,
                                        phases: tuple | None = None) -> dict:
    """:func:`cross_phase_functionals` swept over a diffusion-scale grid.

    ``eig_by_phase`` maps each phase present to its ``(eigenvalues,
    eigenvectors)``. The phase set is taken as **data**: whichever phases are in
    the dict are used, so adding or removing the ``encode`` phase needs no
    change here or in any caller. Returns ``{functional: array(len(s_grid))}``,
    ``NaN`` wherever a cophenetic tree is degenerate.
    """
    if phases is not None and roles is None:
        roles = dict(zip(("baseline_a", "baseline_b", "encode", "probe", "follow"), phases))
    roles = dict(CROSS_PHASE_ROLES if roles is None else roles)
    need = [roles[k] for k in ("baseline_a", "baseline_b", "probe", "follow")]
    missing = [ph for ph in need if ph not in eig_by_phase]
    if missing:
        raise KeyError(f"cross-phase functionals need phases {missing}")
    has_enc = roles.get("encode") in eig_by_phase
    phases = need + ([roles["encode"]] if has_enc else [])
    keys = ["T_probe"] + (["T_encode", "T_probespec", "T_probespec_pe"] if has_enc else [])
    keys = keys + [_FUNCTIONAL_ALIASES[k] for k in keys]
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


# --------------------------------------------------------------------------- #
# rank-correlation tensor: every cross-phase functional as a post-hoc derivation
# --------------------------------------------------------------------------- #
def cross_phase_vectors(D: dict, roles: dict | None = None) -> dict:
    """Reorganisation vectors of a five-phase arc, keyed by
    :data:`CROSS_PHASE_VECTOR_NAMES`.

    ``D`` maps phase name -> condensed cophenetic distance vector, all on the
    same pair index set; ``roles`` maps the semantic roles to phase names and
    defaults to :data:`CROSS_PHASE_ROLES`. See :func:`cross_phase_functionals`
    for what each vector means; ``d`` is the task-free baseline drift.
    """
    roles = dict(CROSS_PHASE_ROLES if roles is None else roles)
    A, B = D[roles["baseline_a"]], D[roles["baseline_b"]]
    P, F, E = D[roles["probe"]], D[roles["follow"]], D[roles["encode"]]
    return {"e": E - A, "e2": E - B, "f": P - E, "g": P - A, "g2": P - B,
            "p": F - B, "p2": F - A, "d": B - A}


def cross_phase_rank_corr(D: dict, roles: dict | None = None) -> NDArray:
    """Spearman correlation matrix of :func:`cross_phase_vectors`, ``(8, 8)``.

    This is the sufficient statistic for the whole cross-phase family. Every
    functional in :data:`CROSS_PHASE_FUNCTIONALS`, its drift-controlled
    counterpart, the encode/probe role swap and any contrast between them is an
    exact function of this matrix, so a pipeline that stores it can be
    re-interrogated with a new estimator without recomputing a single
    eigendecomposition. Rows are in :data:`CROSS_PHASE_VECTOR_NAMES` order.
    """
    V = cross_phase_vectors(D, roles)
    return np.corrcoef(np.vstack([_average_ranks(V[k])
                                  for k in CROSS_PHASE_VECTOR_NAMES]))


def partial_from_corr_matrix(R: NDArray, i: int, j: int, controls=()) -> float:
    """Partial correlation ``r(x_i, x_j | x_controls)`` from a correlation matrix.

    Any order, via the precision matrix of the relevant submatrix:
    ``r_ij.rest = -P_ij / sqrt(P_ii P_jj)``. With no controls this returns
    ``R[i, j]``; with one it agrees with the closed-form first-order expression
    to floating point. ``NaN`` if the submatrix is singular (a control that is a
    linear combination of the others).
    """
    ctrl = [int(c) for c in controls if int(c) not in (int(i), int(j))]
    if not ctrl:
        return float(R[i, j])
    idx = [int(i), int(j)] + ctrl
    S = np.asarray(R, float)[np.ix_(idx, idx)]
    if not np.all(np.isfinite(S)):
        return float("nan")
    try:
        P = np.linalg.inv(S)
    except np.linalg.LinAlgError:
        return float("nan")
    den = np.sqrt(P[0, 0] * P[1, 1])
    return float(-P[0, 1] / den) if np.isfinite(den) and den > 0 else float("nan")


#: Encode<->probe role swap as a SIGNED PERMUTATION of the vectors, not a
#: recomputation. Exchanging which phase is "encode" and which is "probe" maps
#: ``e -> g``, ``e2 -> g2``, ``g -> e``, ``g2 -> e2``, ``f -> -f`` and leaves
#: ``p``, ``p2``, ``d`` untouched -- so the swapped arc's rank correlations are
#: obtained from the observed ones by relabelling and a sign, exactly. This is
#: what makes the swap a construction-preserving null: it changes only the
#: semantic assignment of two phases and nothing about the graphs.
_SWAP_ENCODE_PROBE = {"e": ("g", 1), "e2": ("g2", 1), "f": ("f", -1),
                      "g": ("e", 1), "g2": ("e2", 1), "p": ("p", 1),
                      "p2": ("p2", 1), "d": ("d", 1)}


def _swapped_corr(R: NDArray) -> NDArray:
    ix = {k: i for i, k in enumerate(CROSS_PHASE_VECTOR_NAMES)}
    perm = [ix[_SWAP_ENCODE_PROBE[k][0]] for k in CROSS_PHASE_VECTOR_NAMES]
    sgn = np.array([_SWAP_ENCODE_PROBE[k][1] for k in CROSS_PHASE_VECTOR_NAMES],
                   float)
    return np.asarray(R, float)[np.ix_(perm, perm)] * np.outer(sgn, sgn)


def cross_phase_functionals_from_corr(R: NDArray, *, swap: bool = False,
                                      drift: bool = True) -> dict:
    """Every cross-phase functional derived from a stored rank-correlation matrix.

    ``R`` is the ``(8, 8)`` matrix from :func:`cross_phase_rank_corr`. Returns
    :data:`CROSS_PHASE_FUNCTIONALS` (both naming schemes) and, when ``drift``,
    :data:`CROSS_PHASE_FUNCTIONALS_DRIFT`.

    The drift-controlled variants partial the task-free drift direction ``d`` out
    of both arguments of every correlation. Motivation: the two baselines are the
    contiguous halves of one resting recording, so a session that drifts
    monotonically makes *every* later phase resemble *every* other later phase,
    and the arc inherits that similarity with no task involved. The symmetric
    A/B arm average cancels shared-*baseline* bias but NOT drift -- both arms
    pick up the drift component with the same sign -- which is why an arc carved
    out of pure rest, in temporal order, returns positive functionals.
    Conditioning on ``d`` removes the component of each change vector that lies
    along the measured drift direction.

    Two caveats that belong with any number this produces. ``d`` is estimated
    from half-length recordings, so it is noisier than the full-phase vectors,
    and partialling out a noisy regressor UNDER-corrects (regression dilution):
    a surviving drift-controlled value is an upper bound on the drift-free
    effect, not an unbiased estimate of it. And ``d`` spans one rest recording
    while the arc spans much longer, so if drift is non-linear in time ``d`` has
    the right direction but the wrong magnitude. Neither is a reason not to
    condition; both are reasons to validate the corrected statistic on a
    no-signal arc rather than trust it.

    ``swap`` evaluates the encode<->probe role exchange via
    :data:`_SWAP_ENCODE_PROBE`, an exact signed relabelling of the same
    vectors -- no second pass over the data.
    """
    R = _swapped_corr(R) if swap else np.asarray(R, float)
    ix = {k: i for i, k in enumerate(CROSS_PHASE_VECTOR_NAMES)}

    def pr(a, b, ctrl=()):
        return partial_from_corr_matrix(R, ix[a], ix[b], [ix[c] for c in ctrl])

    out = {}
    for suffix, base in ((("", ()),) + ((("_d", ("d",)),) if drift else ())):
        out["T_probe" + suffix] = 0.5 * (pr("g", "p", base) + pr("g2", "p2", base))
        out["T_encode" + suffix] = 0.5 * (pr("e", "p", base) + pr("e2", "p2", base))
        out["T_probespec" + suffix] = 0.5 * (pr("f", "p", base) + pr("f", "p2", base))
        key = "T_probespec_pe" + ("d" if suffix else "")
        out[key] = 0.5 * (pr("f", "p", ("e",) + base) + pr("f", "p2", ("e2",) + base))
    for k, v in list(out.items()):
        base_k = k[:-2] if k.endswith("_d") else (k[:-1] if k.endswith("_ped") else k)
        alias = _FUNCTIONAL_ALIASES.get(base_k)
        if alias:
            out.setdefault(alias + k[len(base_k):], v)
    return out


def cross_phase_rank_corr_over_scales(eig_by_phase: dict, s_grid: NDArray,
                                      roles: dict | None = None,
                                      rho_floor: float = RHO_FLOOR) -> NDArray:
    """:func:`cross_phase_rank_corr` swept over a diffusion-scale grid.

    ``eig_by_phase`` maps each of the five phases to its ``(eigenvalues,
    eigenvectors)``. Returns ``(len(s_grid), 8, 8)``, ``NaN`` at any scale whose
    cophenetic tree is degenerate (non-finite, or constant so the rank
    correlation is undefined).

    This is the storage form a knob-integrated pipeline should write instead of
    collapsed functionals: one pass over the eigendecompositions yields an object
    from which every cross-phase estimator -- including ones invented after the
    compute finished -- is recovered exactly by
    :func:`cross_phase_functionals_from_corr`.
    """
    roles = dict(CROSS_PHASE_ROLES if roles is None else roles)
    phases = [roles[k] for k in ("baseline_a", "baseline_b", "encode", "probe", "follow")]
    n = len(CROSS_PHASE_VECTOR_NAMES)
    out = np.full((len(s_grid), n, n), np.nan)
    for j, s in enumerate(s_grid):
        try:
            D = {ph: cophenetic_at_scale(*eig_by_phase[ph], s, rho_floor)
                 for ph in phases}
            if any((not np.all(np.isfinite(v))) or np.std(v) == 0 for v in D.values()):
                continue
            out[j] = cross_phase_rank_corr(D, roles)
        except Exception:                                          # noqa: BLE001
            continue
    return out


def _partial1(r_xy, r_xz, r_yz):
    """Vectorised first-order partial ``r(x, y | z)`` on broadcastable arrays."""
    den = np.sqrt(np.clip(1.0 - r_xz ** 2, 0.0, None)
                  * np.clip(1.0 - r_yz ** 2, 0.0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (r_xy - r_xz * r_yz) / den, np.nan)


def cross_phase_functionals_from_corr_stack(R: NDArray, *, swap: bool = False,
                                            drift: bool = True) -> dict:
    """:func:`cross_phase_functionals_from_corr` over a whole stack at once.

    ``R`` is ``(..., 8, 8)``; returns ``{functional: array(...)}``. Identical
    values to the per-matrix function (verified elementwise), obtained from the
    recursive first-order partial identity rather than a matrix inverse, so a
    surrogate ensemble of millions of matrices is evaluated in one pass. Both
    naming schemes are emitted, as elsewhere in this module.
    """
    R = _swapped_corr_stack(R) if swap else np.asarray(R, float)
    ix = {k: i for i, k in enumerate(CROSS_PHASE_VECTOR_NAMES)}

    def r(a, b):
        return R[..., ix[a], ix[b]]

    out = {
        "T_probe": 0.5 * (r("g", "p") + r("g2", "p2")),
        "T_encode": 0.5 * (r("e", "p") + r("e2", "p2")),
        "T_probespec": 0.5 * (r("f", "p") + r("f", "p2")),
        "T_probespec_pe": 0.5 * (
            _partial1(r("f", "p"), r("f", "e"), r("p", "e"))
            + _partial1(r("f", "p2"), r("f", "e2"), r("p2", "e2"))),
    }
    if drift:
        def pd_(a, b):
            return _partial1(r(a, b), r(a, "d"), r(b, "d"))
        out["T_probe_d"] = 0.5 * (pd_("g", "p") + pd_("g2", "p2"))
        out["T_encode_d"] = 0.5 * (pd_("e", "p") + pd_("e2", "p2"))
        out["T_probespec_d"] = 0.5 * (pd_("f", "p") + pd_("f", "p2"))
        # second order: condition on the encoding vector, then on drift
        halves = []
        for pk, ek in (("p", "e"), ("p2", "e2")):
            a = _partial1(r("f", pk), r("f", ek), r(pk, ek))       # r(f,p|e)
            b = _partial1(r("f", "d"), r("f", ek), r("d", ek))     # r(f,d|e)
            c = _partial1(r(pk, "d"), r(pk, ek), r("d", ek))       # r(p,d|e)
            halves.append(_partial1(a, b, c))
        out["T_probespec_ped"] = 0.5 * (halves[0] + halves[1])
    for k, v in list(out.items()):
        base_k = k[:-2] if k.endswith("_d") else (k[:-1] if k.endswith("_ped") else k)
        alias = _FUNCTIONAL_ALIASES.get(base_k)
        if alias:
            out.setdefault(alias + k[len(base_k):], v)
    return out


def _swapped_corr_stack(R: NDArray) -> NDArray:
    ix = {k: i for i, k in enumerate(CROSS_PHASE_VECTOR_NAMES)}
    perm = [ix[_SWAP_ENCODE_PROBE[k][0]] for k in CROSS_PHASE_VECTOR_NAMES]
    sgn = np.array([_SWAP_ENCODE_PROBE[k][1] for k in CROSS_PHASE_VECTOR_NAMES],
                   float)
    return np.asarray(R, float)[..., perm, :][..., :, perm] * np.outer(sgn, sgn)
