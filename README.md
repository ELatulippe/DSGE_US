# Replication package — St. Louis Fed DSGE model (Faria-e-Castro, 2026)

Independent replication in Dynare on US data, 1959Q1–2026Q1: re-estimation, historical
decompositions, FEVD, impulse responses, natural-rate (r\*) comparisons, unconditional and
conditional forecasts, and a companion note (with full derivations). This README lists every file,
the run order, and — importantly — the **naming rules** that must be respected for the pipeline to work.

---

## 1. Prerequisites

- **Dynare 7.1** on **MATLAB** (tested on Windows). A C compiler for `use_dll` (MATLAB's `mex -setup`).
- MATLAB **Statistics** and (for parallel MCMC) **Parallel Computing** toolboxes.
- A **LaTeX** distribution (pdflatex) to build the note.
- Internet access once, to download the HLW natural-rate file for Figure 6 (see §6).

---

## 2. What's in the package

```
frbstl-replication/
├── README.md
├── model/      Dynare .mod files, parallel config, saved posterior mode
├── data/       observables (.mat + .csv) and the data-build script
├── scripts/    MATLAB post-processing + diagnostics; external r* data
└── note/       LaTeX companion note + derivations appendix + figures/
```

### model/
| File | Purpose |
|---|---|
| `frbstl_us_est_fastlr_lean_ESTIM.mod` | **Estimation** run (STEP B): random-walk Metropolis–Hastings, 12 chains × 250k draws, loads the mode, writes the posterior draws. |
| `frbstl_us_est_fastlr_lean_POSTEST.mod` | **Post-estimation** run: adds the flexible-price / r\* block, reuses the draws (`load_mh_file`), and produces the smoother, historical decompositions, forecasts (`forecast=40`), and the conditional scenarios. |
| `fevd_run.mod` | **Standalone fast run** (no MCMC): loads the mode and produces the FEVD arrays and the 40-quarter impulse responses in ~1–2 min. |
| `cluster.ini` | Local parallel-pool config for the MH chains. |
| `frbstl_us_est_fastlr_mode.mat` | Saved posterior **mode** (input to all three `.mod` files via `mode_file=`). |

> The two `_ESTIM` / `_POSTEST` suffixes are **storage names only**. At run time both must be renamed to
> the same base name — see §4. This is the single most important thing in this README.

### data/
| File | Purpose |
|---|---|
| `usdata_1959_2026.mat` | The 15 quarterly observables read by every `.mod` (`datafile=`). |
| `usdata_1959_2026.csv` | Human-readable copy of the same data. |
| `build_usdata.py` | Script that constructs the `.mat` from source series (edit here to refresh the data). |

### scripts/  (run in MATLAB, after the Dynare runs)
| File | Produces | Reads |
|---|---|---|
| `make_paper_figures.m` | Figs 4–8 + Table 5 (`table5_forecast.csv`, `table5_growth_decomp.csv`) | POSTEST results |
| `make_detailed_decomp.m` | Figs 9–10 (12-shock decompositions) | POSTEST results |
| `make_fevd.m` | `fevd_*.png` + FEVD table | `fevd_run` (or POSTEST) results |
| `make_irf_panel.m` | `fig_irf_page1/2.png` (12 shocks × 4 vars) | `fevd_run` results |
| `make_counterfactual.m` | `cf_mp_hike_100bp.png` (+ higher-for-longer if present) | `fevd_run` results |
| `diag_gap.m` | gap sanity diagnostic | POSTEST results |
| `diag_forecast.m` | forecast/r\*/inflation-target diagnostic | POSTEST results |
| `add_hlw_to_external.m` | fills the HLW column of `rstar_external.csv` | official HLW `.xlsx` (you download) |
| `rstar_external.csv` | external r\* series for Fig 6 (Lubik–Matthes pre-filled; HLW added by the script above) | — |

### note/
| File | Purpose |
|---|---|
| `frbstl_replication_note.tex` | The companion note (methodology, estimation, results, replication assumptions, references). |
| `app_derivations.tex` | The full model derivations, `\input` by the note as the appendix. |
| `figures/` | Where the note looks for its figures (drop the generated PNGs here — see §5). |

---

## 3. End-to-end run order

```
(0) data       python build_usdata.py                  ->  data/usdata_1959_2026.mat   [optional; only to refresh data]
(1) estimate   dynare frbstl_us_est_fastlr_lean        ->  posterior draws              [~34 h, 12 chains]
(2) postest    dynare frbstl_us_est_fastlr_lean        ->  smoother/forecast/decomp     [~2 h, reuses draws]
(3) fast run   dynare fevd_run                         ->  FEVD + IRFs                  [~2 min]
(4) figures    (MATLAB) make_paper_figures; make_detailed_decomp; make_fevd;
                        make_irf_panel; make_counterfactual; diag_gap; diag_forecast
(5) external   download HLW xlsx; (MATLAB) add_hlw_to_external
(6) note       copy figures -> note/figures/;  pdflatex x2
```

Steps (1) and (2) run **the same base name** `frbstl_us_est_fastlr_lean` (see §4). Step (3) is
independent and cheap — it's the quickest way to regenerate FEVD/IRFs without touching the 2 h smoother.

---

## 4. Naming rules — read this before renaming anything

The pipeline hard-codes a few names. Change one and you must change its matches, or the run breaks.

**(a) The estimation and post-estimation files must share the base name `frbstl_us_est_fastlr_lean`.**
Dynare writes the MH draws to a folder named after the running `.mod` (`frbstl_us_est_fastlr_lean/metropolis/`),
and the post-estimation step finds them via `load_mh_file`, which looks in a folder with the **same base
name**. The two files are stored here as `..._ESTIM.mod` and `..._POSTEST.mod` only to keep them apart on
disk. To run:

```
# 1) estimation
copy  frbstl_us_est_fastlr_lean_ESTIM.mod    frbstl_us_est_fastlr_lean.mod
dynare frbstl_us_est_fastlr_lean            # writes frbstl_us_est_fastlr_lean/metropolis/

# 2) post-estimation (same folder, same base name; draws are reused)
copy  frbstl_us_est_fastlr_lean_POSTEST.mod  frbstl_us_est_fastlr_lean.mod   # overwrite
dynare frbstl_us_est_fastlr_lean            # load_mh_file reads the draws from step 1
```

Equivalent alternative: keep them in two folders, each holding a file literally named
`frbstl_us_est_fastlr_lean.mod`, and copy the generated `frbstl_us_est_fastlr_lean/` results folder from
the estimation folder into the post-estimation folder before step 2.

If you want a different base name, rename **consistently** everywhere: the two `.mod` files, and there is
nothing else to change inside them (they refer to their own base name implicitly) — but the results folder
and any hard path you use in the MATLAB scripts must follow.

**(b) `mode_file=frbstl_us_est_fastlr_mode`.** All three `.mod` files load the mode by this name. The
`.mat` must sit in the run folder and be named `frbstl_us_est_fastlr_mode.mat`. Rename the mode → update
the `mode_file=` option in `*_ESTIM.mod`, `*_POSTEST.mod`, and `fevd_run.mod`.

**(c) `datafile='usdata_1959_2026.mat'`.** Same idea — appears in the `estimation(...)` call of all three
`.mod`. Rename the data → update `datafile=` in each.

**(d) Results path read by the MATLAB scripts.** The scripts auto-detect
`frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat` (and a couple of fallbacks). If
you changed the base name in (a), pass the new results path explicitly, e.g.
`make_paper_figures('MYNAME/Output/MYNAME_results.mat')`. Every script accepts the results path as its
first argument.

**(e) Figure filenames ↔ the note.** The scripts write PNGs into `figures_paper/` with fixed names, and
the note includes them by those exact names via `\incfig{...}`. If you rename a PNG, update the matching
`\incfig{}` line in `frbstl_replication_note.tex`. The full mapping is in §5.

**(f) Conditional-scenario files.** The POSTEST run saves `cond_constrate.mat`, `cond_infltarget.mat`, and
`cond_tighter.mat` in the run folder; `make_paper_figures` (Fig 8) and `make_counterfactual` look for them
by those names in the current MATLAB folder. Keep them alongside the scripts, or `cd` into the run folder.

---

## 5. Figures the note expects (drop these into `note/figures/`)

| Note figure | Filename | Produced by |
|---|---|---|
| Fig 4a / 4b | `fig4a_core_pce_yoy.png`, `fig4b_headline_pce_yoy.png` | `make_paper_figures` |
| Fig 5a / 5b | `fig5a_output_gap.png`, `fig5b_rstar_spot.png` | `make_paper_figures` |
| Fig 6 (r\*) | `fig6_rstar.png` | `make_paper_figures` (+ `add_hlw_to_external`) |
| Fig 7 / 8 | `fig7_uncond_forecast.png`, `fig8_cond_scenarios.png` | `make_paper_figures` |
| Figs 9–10 | `fig9a_core_detailed.png`, `fig9b_headline_detailed.png`, `fig10a_gap_detailed.png`, `fig10b_rstar_detailed.png` | `make_detailed_decomp` |
| FEVD | `fevd_gdp.png`, `fevd_obs_pi.png`, `fevd_obs_ffr.png`, `fevd_rstar.png`, `fevd_ogap.png`, `fevd_iv.png` | `make_fevd` |
| IRFs | `fig_irf_page1.png`, `fig_irf_page2.png` | `make_irf_panel` |
| Counterfactual | `cf_mp_hike_100bp.png` | `make_counterfactual` |

The scripts save into `figures_paper/`. Copy that folder's PNGs into `note/figures/` (same names), then
build the note:

```
cd note
pdflatex frbstl_replication_note.tex
pdflatex frbstl_replication_note.tex     # second pass for the ToC and cross-references
```

The note uses `\graphicspath{{figures/}}`, so only the folder name `figures/` matters, not the location of
`figures_paper/`.

---

## 6. External r\* data for Figure 6

`rstar_external.csv` ships with the **Lubik–Matthes** median already filled (1985Q1–2026Q1). For the
**Holston–Laubach–Williams** line, download the official file once and run the loader:

1. Get `Holston_Laubach_Williams_current_estimates.xlsx` from the New York Fed's r\* page and put it in the
   scripts/run folder.
2. In MATLAB: `add_hlw_to_external` (fills the `HLW_rstar` column of `rstar_external.csv`).
3. Re-run `make_paper_figures` — Fig 6 now shows both external series.

---

## 7. Reproducibility notes

- **Software:** Dynare 7.1; results were produced with `use_dll` and a 12-worker local pool. Different
  Dynare/BLAS versions can shift low digits.
- **Runtime:** estimation ≈ 34 h (12 × 250k), post-estimation ≈ 2 h, `fevd_run` ≈ 2 min, figures seconds.
- **What's an assumption, not a model output:** the ≈0.4pp population add-on in Table 5 / the growth
  decomposition (`POP_GROWTH` in `make_paper_figures.m`) is a flat placeholder for the Census
  population-growth trend; the flexible-price/r\* block and the gap-centering convention are described in
  the note's "Assumptions" section. Read that section before interpreting r\*, the gap, or the forecast
  table.
- **Three published-equation tensions** (utilization tax wedge, the firm's nominal SDF, the steady-state
  wage markup) are flagged in the note's appendix; they matter for any recalibration.

---

## 8. Quick start (if the estimation is already done)

If you already have the posterior draws (`frbstl_us_est_fastlr_lean/metropolis/`) and the mode:

```
# post-estimation objects
copy frbstl_us_est_fastlr_lean_POSTEST.mod frbstl_us_est_fastlr_lean.mod
dynare frbstl_us_est_fastlr_lean
dynare fevd_run
# MATLAB
make_paper_figures; make_detailed_decomp; make_fevd; make_irf_panel; make_counterfactual
diag_gap; diag_forecast
# note
copy figures_paper\*.png note\figures\
cd note && pdflatex frbstl_replication_note.tex && pdflatex frbstl_replication_note.tex
```
