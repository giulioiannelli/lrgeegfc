#!/usr/bin/env python3
"""Audit 83 — white-matter-node-stratified cross-phase trace (cophenetic + raw).

White-matter analogue of audit_77. Re-runs the split-baseline cross-phase
trace across six WM views — full graph, gray↔gray pairs, WM-excluded subgraph,
gray↔WM interface pairs, WM↔WM pairs, WM-only subgraph — for **two substrates**
(LRG cophenetic ρ_split and raw |ImCoh| ρ_split^raw), all six bands, the full
n=10 cohort, into ONE tidy cache feeding the WM-exclusion report.

Both substrates are read off the SAME matched-strength eigendecomposition
ensemble per (config, patient, band, phase): the cophenetic vector via
``cophenetic_condensed_from_eigs`` and the raw adjacency via
``adjacency_from_laplacian_eigs`` (= the off-diagonal of −L). So adding the raw
substrate costs no extra surrogate generation.

Critical preamble (per CLAUDE.md rule, before any code)
=======================================================
(1) **Claim.** The established split-baseline trace (positive = TRACE;
    ρ_split = Spearman(D_task − D_preA, D_post − D_preB)) is (a) robust to
    removing white-matter contacts from the montage (exclude_wm — the literal
    "drop WM from the timeseries" view), (b) present or absent within the
    white-matter subgraph (wm_only), and (c) carried disproportionately by a
    particular WM-stratified pair class — gray↔gray, gray↔WM interface (cross),
    or WM↔WM.

(2) **Null.** Matched-strength 4-cycle ±δ surrogate, R=200. Subgraph configs
    regenerate the ensemble on the submatrix W[keep,keep] (own cache dir, never
    clobbering the canonical full-graph cache). Pair-class configs use the
    FULL-graph ensemble with ρ_split recomputed on the restricted pair set
    (holds the global per-node strength sequence fixed, asks whether the
    class-restricted co-movement survives).

(3) **Strongest plausible alternatives.** (i) Class/montage size: gray↔gray and
    the full graph have the most pairs → a more stable Spearman, independent of
    any localization; conversely exclude_wm/wm_only lose 30–57 % of nodes →
    *weaker* Spearman from sample size alone, which could masquerade as the
    trace "needing" WM. (ii) The global strength sequence rather than tissue
    topology. (iii) WM contacts are physiologically suspect (volume-conducted,
    low-amplitude); a trace they carry could be a recording-geometry artifact
    rather than neural reorganization — exclude_wm is exactly the test of that.

(4) **Null's mechanical reach.** Matched-strength fixes per-node strength
    exactly, so a surviving statistic is NOT explained by (ii). It does NOT
    control for (i): node/pair-count differences are reported (n_nodes, n_pairs)
    and read alongside the verdict; a *weakening* under exclude_wm is therefore
    ambiguous between "trace was WM-carried" and "fewer nodes → noisier
    Spearman" — disambiguated only by whether wm_only and gray↔gray point the
    same way. The pair-class null is global-rewiring, not within-class.

(5) **Falsification + limitations.** "Trace lives in gray matter / survives WM
    removal" requires exclude_wm = persist (verdict separated, Wilcoxon
    p<0.05). "Trace is WM-carried" requires exclude_wm = weaken AND
    (wm_only separated OR wm_wm/cross carrying it). "Trace is interface" requires
    cross separated AND gray_gray weakening. Limitations: WM is a large fraction
    so exclude_wm halves some graphs (power loss is real, not just statistical
    convenience); wm_only is the weaker arm but — unlike epi_only — has 36–67
    nodes so its matched-strength null is NOT degenerate; Pat_15 has the most WM
    (67, 57 %) and is the known β-LRG anti-aligned patient, so a
    leave-Pat_15-out cohort p is reported.

Outputs
-------
``data/audit/wm_stratified/``
    cophenetic_raw_per_patient.csv   tidy long (one row per pat×band×config×substrate)
    cophenetic_raw_cohort.csv        per (substrate, band, config) cohort verdict
    README.md

Surrogate caches (per config, no clobber)
-----------------------------------------
- full / pair-class → ``data/cache/matched_strength_surrogate_lrg/`` (seed 20260511, reused)
- exclude_wm        → ``data/cache/matched_strength_surrogate_wm_excluded_lrg/`` (seed 20260608)
- wm_only           → ``data/cache/matched_strength_surrogate_wm_only_lrg/`` (seed 20260608)
"""
from __future__ import annotations

import argparse
import gc
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.surrogate import (
    adjacency_from_laplacian_eigs,
    cophenetic_condensed_from_eigs,
    load_or_compute_eigs_at_path,
)

# Reuse audit_63 cophenetic + half-FC pipeline and the shared WM vocabulary.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs,
    load_phase_fc,
    lrg_ultrametric_condensed,
)
import _wm_stratify as ws  # type: ignore


SUBSTRATES = ("cophenetic", "raw")
OUT = ROOT / "data" / "audit" / "wm_stratified"
OUT.mkdir(parents=True, exist_ok=True)

LO_PATIENT = "Pat_15"   # max WM (57%) + known β-LRG anti-aligned → robustness line


# ---------------------------------------------------------------------------
# Condensed vectors
# ---------------------------------------------------------------------------
def raw_condensed(W: np.ndarray) -> np.ndarray:
    """Off-diagonal upper-triangle of the raw |ImCoh| adjacency (audit_74)."""
    W = W.copy()
    np.fill_diagonal(W, 0.0)
    W = np.maximum(W, W.T)
    return squareform(W, checks=False)


def _surr_stacks(W: np.ndarray, cache_path: Path, rng: np.random.Generator,
                 n_surr: int, swap_factor: int, verbose: bool
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(coph_stack, raw_stack)`` each ``(R, n_pairs)`` for the
    matched-strength ensemble of ``W`` — cophenetic condensed and raw-adjacency
    condensed, both reconstructed from the SAME cached eigendecomposition.
    Strength-violating surrogate rows stay NaN for downstream masking."""
    evals, evecs = load_or_compute_eigs_at_path(
        cache_path, W, n_surr, swap_factor, rng, verbose=verbose)
    R = evals.shape[0]
    N = W.shape[0]
    n_pairs = N * (N - 1) // 2
    C = np.full((R, n_pairs), np.nan, dtype=np.float64)
    Rw = np.full((R, n_pairs), np.nan, dtype=np.float64)
    for r in range(R):
        if not np.isfinite(evals[r]).all():
            continue
        C[r] = cophenetic_condensed_from_eigs(evals[r], evecs[r])
        A = adjacency_from_laplacian_eigs(evals[r], evecs[r])
        np.fill_diagonal(A, 0.0)
        A = np.maximum(A, A.T)
        Rw[r] = squareform(A, checks=False)
    return C, Rw


def _rho(a: np.ndarray, b: np.ndarray) -> float:
    rho, _ = spearmanr(a, b)
    return float(rho)


# ---------------------------------------------------------------------------
# Row builders (substrate-parameterized)
# ---------------------------------------------------------------------------
def _row(pat: str, band: str, config: str, substrate: str, obs: float,
         surr: np.ndarray | None, n_nodes: int, n_pairs: int,
         surrogate_source: str, n_swaps: int, seed: int, defined: bool) -> dict:
    if surr is not None and surr.size > 0:
        smean = float(np.mean(surr)); sstd = float(np.std(surr, ddof=1))
        qs = np.quantile(surr, [0.05, 0.25, 0.50, 0.75, 0.95])
        z = (obs - smean) / sstd if sstd > 0 else float("nan")
        p_up = float(np.mean(surr >= obs))
        n_s = int(surr.size)
    else:
        smean = sstd = float("nan")
        qs = [float("nan")] * 5
        z = p_up = float("nan"); n_s = 0
    return {
        "patient": pat, "band": band, "config": config,
        "config_class": ws.CONFIG_CLASS[config], "substrate": substrate,
        "n_nodes": int(n_nodes), "n_pairs": int(n_pairs),
        "obs_stat": float(obs) if obs == obs else float("nan"),
        "surr_mean": smean, "surr_std": sstd,
        "surr_p5": float(qs[0]), "surr_p25": float(qs[1]),
        "surr_p50": float(qs[2]), "surr_p75": float(qs[3]),
        "surr_p95": float(qs[4]),
        "obs_z": z, "obs_p_one_sided": p_up,
        "n_surrogates": n_s, "n_swaps_per_surrogate": int(n_swaps),
        "surrogate_source": surrogate_source, "seed": int(seed),
        "defined": bool(defined),
    }


def _stats_row(pat: str, band: str, config: str, substrate: str,
               dT_obs: np.ndarray, dR_obs: np.ndarray, dT_surr: np.ndarray,
               dR_surr: np.ndarray, pair_keep: np.ndarray | None, n_nodes: int,
               surrogate_source: str, n_swaps: int, seed: int) -> dict:
    sel = slice(None) if pair_keep is None else pair_keep
    n_pairs = int(dT_obs.size) if pair_keep is None else int(pair_keep.sum())
    obs = _rho(dT_obs[sel], dR_obs[sel])
    R = dT_surr.shape[0]
    srho = np.full(R, np.nan)
    for r in range(R):
        a, b = dT_surr[r][sel], dR_surr[r][sel]
        if np.isfinite(a).all() and np.isfinite(b).all():
            srho[r] = _rho(a, b)
    surr = srho[np.isfinite(srho)]
    return _row(pat, band, config, substrate, obs, surr, n_nodes, n_pairs,
                surrogate_source, n_swaps, seed, defined=True)


def _undefined_rows(pat: str, band: str, config: str, n_nodes: int) -> list[dict]:
    return [_row(pat, band, config, s, float("nan"), None, n_nodes, 0,
                 "undefined", 0, ws.config_seed(config), defined=False)
            for s in SUBSTRATES]


# ---------------------------------------------------------------------------
# Per (patient, band): all requested configs, both substrates, shared eigs
# ---------------------------------------------------------------------------
def _diffs(D: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    return (D["task_test"] - D["rest_pre_A"],
            D["rest_post"] - D["rest_pre_B"])


def per_patient_band(pat: str, band: str, configs: list[str], n_surr: int,
                     swap_factor: int, verbose: bool) -> list[dict]:
    try:
        Ws = {ph: load_phase_fc(pat, ph, band) for ph in ws.PHASES_4}
    except Exception as e:
        if verbose:
            print(f"[audit_83] SKIP {pat}/{band}: load failed: {e}")
        return []
    N = Ws["rest_pre_A"].shape[0]
    if any(W.shape != (N, N) for W in Ws.values()):
        if verbose:
            print(f"[audit_83] SKIP {pat}/{band}: phase shape mismatch")
        return []
    wm = ws.wm_mask_for(pat, n_expected=N)
    rows: list[dict] = []

    # ---- Full-graph: serves 'full' + pair-class configs (both substrates) ----
    pc = [c for c in configs if c in ws.PAIRCLASS_CONFIGS]
    if "full" in configs or pc:
        Dc = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in ws.PHASES_4}
        Dr = {ph: raw_condensed(Ws[ph]) for ph in ws.PHASES_4}
        dTc, dRc = _diffs(Dc)
        dTr, dRr = _diffs(Dr)
        sc, sr = {}, {}
        for ph in ws.PHASES_4:
            sc[ph], sr[ph] = _surr_stacks(
                Ws[ph], ws.surr_eig_path("full", pat, band, ph, n_surr, swap_factor),
                ws.cell_rng(pat, band, ph, "full"), n_surr, swap_factor, verbose)
        dTcs, dRcs = sc["task_test"] - sc["rest_pre_A"], sc["rest_post"] - sc["rest_pre_B"]
        dTrs, dRrs = sr["task_test"] - sr["rest_pre_A"], sr["rest_post"] - sr["rest_pre_B"]
        n_swaps_full = swap_factor * (N * (N - 1)) // 2
        seed_full = ws.config_seed("full")

        def _emit_full(config, pk, source):
            rows.append(_stats_row(pat, band, config, "cophenetic", dTc, dRc,
                                   dTcs, dRcs, pk, N, source, n_swaps_full, seed_full))
            rows.append(_stats_row(pat, band, config, "raw", dTr, dRr,
                                   dTrs, dRrs, pk, N, source, n_swaps_full, seed_full))

        if "full" in configs:
            _emit_full("full", None, "fullgraph")
        for c in pc:
            pk = ws.pair_keep_for_config(wm, c)
            if pk is None:
                rows.extend(_undefined_rows(pat, band, c, N))
            else:
                _emit_full(c, pk, "fullgraph_restricted")
        del sc, sr, dTcs, dRcs, dTrs, dRrs
        gc.collect()

    # ---- Subgraph configs: exclude_wm, wm_only ----
    for c in [c for c in configs if c in ("exclude_wm", "wm_only")]:
        keep = ws.node_mask_for_config(wm, c)
        if keep is None:
            rows.extend(_undefined_rows(pat, band, c, 0))
            continue
        idx = np.ix_(keep, keep)
        Wk = {ph: np.ascontiguousarray(Ws[ph][idx]) for ph in ws.PHASES_4}
        Nk = int(keep.sum())
        Dc = {ph: lrg_ultrametric_condensed(Wk[ph]) for ph in ws.PHASES_4}
        Dr = {ph: raw_condensed(Wk[ph]) for ph in ws.PHASES_4}
        dTc, dRc = _diffs(Dc)
        dTr, dRr = _diffs(Dr)
        sc, sr = {}, {}
        for ph in ws.PHASES_4:
            sc[ph], sr[ph] = _surr_stacks(
                Wk[ph], ws.surr_eig_path(c, pat, band, ph, n_surr, swap_factor),
                ws.cell_rng(pat, band, ph, c), n_surr, swap_factor, verbose)
        dTcs, dRcs = sc["task_test"] - sc["rest_pre_A"], sc["rest_post"] - sc["rest_pre_B"]
        dTrs, dRrs = sr["task_test"] - sr["rest_pre_A"], sr["rest_post"] - sr["rest_pre_B"]
        n_swaps = swap_factor * (Nk * (Nk - 1)) // 2
        seed = ws.config_seed(c)
        rows.append(_stats_row(pat, band, c, "cophenetic", dTc, dRc, dTcs, dRcs,
                               None, Nk, "submatrix", n_swaps, seed))
        rows.append(_stats_row(pat, band, c, "raw", dTr, dRr, dTrs, dRrs,
                               None, Nk, "submatrix", n_swaps, seed))
        del sc, sr, dTcs, dRcs, dTrs, dRrs
        gc.collect()

    return rows


# ---------------------------------------------------------------------------
# Cohort summary (per substrate, band, config)
# ---------------------------------------------------------------------------
def _band_order(df: pd.DataFrame) -> list[str]:
    present = set(df.band.unique())
    return [b for b in ws.ALL_BANDS if b in present]


def cohort_summary(per_pat: pd.DataFrame) -> pd.DataFrame:
    out = []
    bands = _band_order(per_pat)
    for substrate in SUBSTRATES:
        sp = per_pat[per_pat.substrate == substrate]
        full_p: dict[str, float] = {}
        for band in bands:
            d = sp[(sp.band == band) & (sp.config == "full") & sp.defined]
            n_above = int((d.obs_p_one_sided < 0.05).sum())
            full_p[band] = ws.cohort_verdict(
                d.obs_stat.values, d.surr_p50.values, n_above)["wilcoxon_p"]
        for band in bands:
            present = set(sp[sp.band == band].config)
            for config in [c for c in ws.ALL_CONFIGS if c in present]:
                sub = sp[(sp.band == band) & (sp.config == config)]
                d = sub[sub.defined]
                n_total = len(sub)
                n_above = int((d.obs_p_one_sided < 0.05).sum())
                v = ws.cohort_verdict(d.obs_stat.values, d.surr_p50.values, n_above)
                d2 = d[d.patient != LO_PATIENT]
                n_above2 = int((d2.obs_p_one_sided < 0.05).sum())
                v2 = ws.cohort_verdict(d2.obs_stat.values, d2.surr_p50.values, n_above2)
                flag = ("baseline" if config == "full"
                        else ws.sensitivity_flag(
                            full_p.get(band, np.nan), v["wilcoxon_p"],
                            cfg_defined=(v["verdict"] != "undefined")))
                out.append({
                    "substrate": substrate, "band": band, "config": config,
                    "config_class": ws.CONFIG_CLASS[config],
                    "n_patients": n_total, "n_defined": v["n_defined"],
                    "obs_median": v["med_obs"],
                    "surr_median_median": v["med_surr"],
                    "median_n_nodes": float(np.median(d.n_nodes)) if len(d) else np.nan,
                    "median_n_pairs": float(np.median(d.n_pairs)) if len(d) else np.nan,
                    "n_above_own_surrogate": f"{n_above}/{v['n_defined']}",
                    "paired_wilcoxon_z": v["wilcoxon_z"],
                    "paired_wilcoxon_p": v["wilcoxon_p"],
                    "lopat15_wilcoxon_p": v2["wilcoxon_p"],
                    "verdict": v["verdict"], "sensitivity_flag": flag,
                })
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# README
# ---------------------------------------------------------------------------
def write_readme(per_pat: pd.DataFrame, cohort: pd.DataFrame,
                 n_surr: int, swap_factor: int, runtime_s: float) -> None:
    L: list[str] = []
    a = L.append
    a("---")
    a("name: wm_stratified_cophenetic_raw")
    a("scope: white_matter_node_stratified_rho_split_trace_all_bands_2substrates")
    a("era: COHORT_N10 / IMCOH_ABS")
    a(f"date: {time.strftime('%Y-%m-%d')}")
    a("status: first_pass")
    a("build_script: scripts/01_compute/audit/audit_83_wm_stratified_cophenetic.py")
    a("---")
    a("")
    a("# White-matter-node-stratified ρ_split trace (cophenetic + raw, n=10)")
    a("")
    a("**Head.** The split-baseline trace (positive = TRACE) re-run across six "
      "white-matter views — full graph, gray↔gray pairs, WM-excluded subgraph "
      "(gray-only montage), gray↔WM interface pairs, WM↔WM pairs, WM-only "
      "subgraph — for two substrates (LRG cophenetic, raw |ImCoh|), all six "
      f"bands, under R={n_surr} matched-strength surrogacy. Subgraph configs "
      "regenerate the null on the submatrix; pair-class configs restrict the "
      "full-graph null. WM is 30–57 % of every montage, so `exclude_wm` is a "
      "structural change — the literal 'drop WM channels before computing FC' "
      "test (exact: ImCoh is pairwise, so the gray-only FC = the gray-only "
      "submatrix; only the node-coupled LRG step changes).")
    a("")
    for substrate in SUBSTRATES:
        a(f"## {substrate} — cohort verdict (band × config)")
        a("")
        a("| band | config | n_def | obs med | surr med | n>surr | Wilcoxon p | "
          "LO-P15 p | med nodes | verdict | vs full |")
        a("|---|---|---|---|---|---|---|---|---|---|---|")
        cs = cohort[cohort.substrate == substrate]
        for band in _band_order(per_pat):
            for _, r in cs[cs.band == band].iterrows():
                mn = "—" if not np.isfinite(r["median_n_nodes"]) else f"{int(r['median_n_nodes'])}"
                a(f"| {BRAIN_BAND_TEX_DICT.get(band, band)} | {r['config']} "
                  f"| {r['n_defined']} | {r['obs_median']:+.3f} "
                  f"| {r['surr_median_median']:+.3f} | {r['n_above_own_surrogate']} "
                  f"| {r['paired_wilcoxon_p']:.4f} | {r['lopat15_wilcoxon_p']:.4f} "
                  f"| {mn} | {r['verdict']} | **{r['sensitivity_flag']}** |")
        a("")
    a("## Sensitivity labels (config vs full, same substrate)")
    a("")
    a("- **persist**: full p<0.05 AND config p<0.05 — trace present in this view.")
    a("- **weaken**: full p<0.05 AND config p≥0.05 — trace needs the removed "
      "pairs/nodes (ambiguous vs power loss; read med nodes / n_pairs).")
    a("- **emerge**: full p≥0.05 AND config p<0.05 — view reveals a masked trace.")
    a("- **absent**: neither significant. **undefined**: <3 defined patients.")
    a("")
    a("## Caveats")
    a("- **WM is a large fraction** (30–57 %, cohort ≈40 %). `exclude_wm` halves "
      "some graphs, so a `weaken` is ambiguous between 'trace was WM-carried' and "
      "'fewer nodes → noisier Spearman'. Disambiguate with wm_only + gray_gray.")
    a("- **wm_only is the weaker arm but NOT degenerate**: 36–67 nodes (vs "
      "epi_only's 6–30), so the matched-strength null is well-posed here.")
    a(f"- **{LO_PATIENT}** has the most WM (57 %) and is the known β-LRG "
      "anti-aligned patient → the LO-P15 column is the leave-one-out robustness "
      "line.")
    a("- **Pair-class null is global-rewiring**, not within-class; class sizes "
      "are unequal — read med n_pairs.")
    a("")
    a("## Median nodes / pairs per config (size context)")
    a("")
    a("| substrate | band | " + " | ".join(ws.ALL_CONFIGS) + " | (nodes) |")
    a("|---|---|" + "|".join(["---"] * (len(ws.ALL_CONFIGS) + 1)) + "|")
    for substrate in SUBSTRATES:
        cs = cohort[cohort.substrate == substrate]
        for band in _band_order(per_pat):
            cells = []
            for c in ws.ALL_CONFIGS:
                row = cs[(cs.band == band) & (cs.config == c)]
                cells.append("—" if row.empty or not np.isfinite(row.iloc[0]["median_n_pairs"])
                             else f"{int(row.iloc[0]['median_n_pairs'])}")
            sub = cs[(cs.band == band) & (cs.config == "exclude_wm")]
            nodes = "—" if sub.empty or not np.isfinite(sub.iloc[0]["median_n_nodes"]) \
                else f"{int(sub.iloc[0]['median_n_nodes'])}"
            a(f"| {substrate} | {band} | " + " | ".join(cells) + f" | {nodes} |")
    a("")
    a("## Provenance")
    a(f"- N_surrogates = {n_surr}, SWAP_FACTOR = {swap_factor}")
    a("- Cohort: " + ", ".join(ws.COHORT))
    a("- Phases: " + ", ".join(ws.PHASES_4))
    a("- FC method: imcoh_abs; τ = 1/λ_max; average-linkage ultrametric.")
    a("- WM = dominant Desikan-Killiany tissue == 'Wm' (atlas argmax).")
    a("- Seeds: full/pair-class 20260511 (reused canonical), "
      "exclude_wm/wm_only 20260608.")
    a(f"- Wall-clock: {runtime_s:.1f} s")
    a("")
    a("## Files")
    a("- `cophenetic_raw_per_patient.csv` — tidy long per pat×band×config×substrate")
    a("- `cophenetic_raw_cohort.csv` — per (substrate, band, config) verdict")
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=list(ws.ALL_BANDS))
    ap.add_argument("--patients", nargs="+", default=list(ws.COHORT))
    ap.add_argument("--configs", nargs="+", default=list(ws.ALL_CONFIGS))
    ap.add_argument("--n-surrogates", type=int, default=ws.N_SURROGATES)
    ap.add_argument("--swap-factor", type=int, default=ws.SWAP_FACTOR)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    print("[audit_83] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, args.bands)

    t0 = time.time()
    all_rows: list[dict] = []
    for band in args.bands:
        for pat in args.patients:
            tc = time.time()
            rows = per_patient_band(pat, band, args.configs,
                                    args.n_surrogates, args.swap_factor,
                                    args.verbose)
            all_rows.extend(rows)
            bits = " ".join(
                f"{r['config'][:6]}/{r['substrate'][:4]}:{r['obs_stat']:+.2f}"
                f"(p{r['obs_p_one_sided']:.2f})"
                for r in rows if r["defined"])
            print(f"[audit_83] {pat}/{band}: {bits} ({time.time()-tc:.1f}s)")

    runtime = time.time() - t0
    per_pat = pd.DataFrame(all_rows)
    per_pat.to_csv(OUT / "cophenetic_raw_per_patient.csv", index=False)
    cohort = cohort_summary(per_pat)
    cohort.to_csv(OUT / "cophenetic_raw_cohort.csv", index=False)
    write_readme(per_pat, cohort, args.n_surrogates, args.swap_factor, runtime)
    print(f"[audit_83] {len(per_pat)} rows, {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
