# Replication package — St. Louis Fed DSGE model (Faria-e-Castro, 2026)

Independent replication in Dynare on US data, 1959Q1–2026Q1: re-estimation of the 63 free
parameters, historical decompositions, forecast-error variance decompositions, impulse responses,
a model-consistent natural rate (r\*), the unconditional and conditional forecasts including the
paper's Table 5, and a companion note with full derivations.

The package adds a flexible-price / flexible-wage block that generates r\* and the output gap. It is
block-recursive and likelihood-neutral: it changes nothing in the estimation, so the posterior mode
and the MH draws are the same objects with or without it.

This README lists every file, the run order with the runtimes actually measured, and the naming
rules the pipeline depends on.

---

## 1. Prerequisites

- **Dynare 7.1** on **MATLAB** (tested on Windows). A C compiler for `use_dll` (`mex -setup`).
- MATLAB **Statistics** and, for the parallel MCMC, **Parallel Computing** toolboxes.
- **On Windows only:** Dynare's parallel engine launches its workers through Sysinternals
  **PsExec**. Install PSTools, put `psexec.exe` on the `PATH`, and run it once by hand to accept the
  EULA. Without that, `dynare ... parallel conffile=cluster.ini` reports
  `'psexec' is not recognized`, prints acceptance ratios of zero, and produces no draws.
- A **LaTeX** distribution (pdflatex) to build the note.
- For the data build only: Python with `fredapi pandas numpy statsmodels scipy openpyxl` and a free
  FRED API key.

---

## 2. Layout

```
DSGE_US/
├── README.md
├── model/     .mod files, parallel config, mode files, data, run inputs   <- THE RUN FOLDER
├── scripts/   MATLAB post-processing and diagnostics
├── data/      data-build scripts and a human-readable copy of the observables
├── note/      LaTeX companion note, derivations appendix, figures/
└── results/   estimation outputs, so the figures regenerate without Dynare
```

### Working folder

**`model/` is the run folder.** Dynare needs the `.mod`, the data file and the mode file in one
directory and writes its output there, and the MATLAB scripts read `rstar_external.csv` and
`cond_*.mat` from the current folder and write `figures_paper/` into it. So run both Dynare and
MATLAB from `model/`:

```matlab
cd model
addpath('../scripts')
```

The scripts look for the estimation outputs in `results/` **and** `../results/`, so they find the
repository's `results/` folder whether you sit in `model/` or at the root.

`data/usdata_1959_2026.mat` and `model/usdata_1959_2026.mat` are currently byte-identical copies of
the same file. Only the one in `model/` is read by anything — `datafile=` resolves relative to the
run folder. Two tracked copies of one data file is how the two versions of `rstar_external.csv`
drifted apart earlier in this project, one of them silently corrupted, so consider deleting
`data/usdata_1959_2026.mat` and keeping `data/usdata_1959_2026.csv` as the readable export.

### model/

| File | Purpose |
|---|---|
| `frbstl_us_est_fastlr_lean_STEP1_MODE.mod` | **Step 1**, posterior mode: `mode_compute = 5` (newrat), `optim = ('MaxIter', 7)`, warm-started from `frbstl_us_est_fastlr_mode_ded_fixed.mat`. Writes `frbstl_us_est_fastlr_mode_ndc.mat`. |
| `frbstl_us_est_fastlr_lean_STEP2_PILOT.mod` | **Step 2**, short pilot: 8 chains × 5,000 draws with `mh_tune_jscale` and `MCMC_jumping_covariance = hessian`, to tune the proposal scale. |
| `frbstl_us_est_fastlr_lean_STEP3_MCMC.mod` | **Step 3**, full MCMC: 8 chains × 250,000 draws, `mh_drop = 0.5`, `mh_jscale = 0.23046` (the value step 2 returned for this model). |
| `frbstl_us_est_fastlr_lean_POSTEST.mod` | **Post-estimation**: adds the flexible-price / r\* block, reuses the draws (`load_mh_file`, `sub_draws = 1000`), and produces the smoother, the historical decompositions, `forecast = 40`, and the conditional scenarios. |
| `fevd_run.mod` | **Standalone fast run** (no MCMC): loads the mode and produces the FEVD arrays and the 40-quarter impulse responses in ~2 min. |
| `frbstl_us_est_fastlr_lean_STEP1b_CSMINWEL.mod` | Fallback mode search (`mode_compute = 4`). Not used in the shipped run. |
| `frbstl_us_est_fastlr_lean_STEP1c_BOUNDED.mod` | Bounded variant of step 1. Not used; see §5 for why bounds do not do what you would expect here. |
| `frbstl_us_est_fastlr_lean_STEP2b_PILOT_PRIORVAR.mod` | Pilot fallback using `MCMC_jumping_covariance = prior_variance`, for the case where the Hessian cannot be repaired. Not used. |
| `cluster.ini` | Local parallel-pool config, `CPUnbr = 8`. Must match `mh_nblocks`. |
| `frbstl_us_est_fastlr_mode_ndc.mat` | Raw step-1 mode. Its Hessian has one all-NaN row (§5). |
| `frbstl_us_est_fastlr_mode_ndc_fixed.mat` | **The mode every later step loads.** Step 1's mode with the broken Hessian row repaired by `scripts/repair_mode_hessian.m`. |
| `frbstl_us_est_fastlr_mode_ded_fixed.mat` | Mode from the earlier specification. Kept only because step 1 warm-starts from it; nothing else reads it. |
| `usdata_1959_2026.mat` | The 15 quarterly observables read by every `.mod` (`datafile=`). |
| `cond_constrate.mat`, `cond_infltarget.mat`, `cond_tighter.mat` | Conditional-forecast scenarios written by the post-estimation run; read by Fig 8 and the counterfactual. |
| `rstar_external.csv` | External r\* comparators for Fig 6: `date, LM_median, LW_rstar`. |
| `Laubach_Williams_current_estimates.xlsx` | The NY Fed workbook `add_lw_to_external.m` parses (§7). |

**Model specification.** The capital-utilization first-order condition carries the capital-tax
wedge, `(1 − τᵏ)Rᵏ = A′(ν)`, with the budget constraints as printed in the paper. This pins the
calibration `κ_a = (1 − τ̄ᵏ)R̄ᵏ = Γ̄/β − (1 − δ)`. Confirmed by the author (Faria-e-Castro, private
correspondence, August 2026). Every `.mod` here uses that branch, each one's header comment says so,
and the note's appendix derives it. `κ_a` is backed out of the steady state and targets no moment,
so the branch moves a calibration constant but not the estimated elasticity of the utilization cost.

### data/

| File | Purpose |
|---|---|
| `build_usdata.py` | Builds `usdata_1959_2026.mat` from FRED plus two non-FRED sources. Prints the `cobs_*` constants to paste into the `.mod`. |
| `load_fernald_spf.py` | Downloads and parses the Fernald utilization-adjusted TFP workbook (SF Fed) and the SPF 10-year CPI expectations file (Philadelphia Fed). |
| `check_vintage.py` | Read-only probe: prints the last available quarter for every input series. Run it before extending the sample; it writes nothing. |
| `usdata_1959_2026.csv` | Human-readable export of the observables. |

A rebuilt `.mat` has to be copied into `model/`, which is where `datafile=` looks.

Extending the sample is **not** an append. `obs_oil` and `obs_tfp` are demeaned over the sample and
the `cobs_*` constants are calibrated to sample means, so one extra quarter changes every historical
value and every printed constant, and the whole pipeline (§4) has to be re-run. `check_vintage.py`
exists so you can confirm the new quarter is complete in *every* series first: `build_usdata.py`
raises only when a column is entirely empty, so a single trailing NaN would pass through and be
treated as a missing observation by the Kalman filter, without warning.

The slowest input is the Flow of Funds transfer series (`BOGZ1FA3664040/30/20 05Q`), which arrives
with the quarterly Z.1 about ten weeks after the quarter ends. The Fernald TFP workbook and the SPF
file are both faster than that.

### scripts/  (MATLAB; run from `model/` with `addpath('../scripts')`)

| File | Produces | Reads |
|---|---|---|
| `check_postest.m` | Gate: verifies the post-estimation file has every variable and horizon the downstream scripts need. **Run this first.** | POSTEST results |
| `make_paper_figures.m` | Figs 4–8 + Table 5 (`table5_forecast.csv`, `table5_growth_decomp.csv`) | POSTEST results |
| `make_detailed_decomp.m` | Figs 9–10 (12-shock decompositions) | POSTEST results |
| `make_fevd.m` | `fevd_*.png` + FEVD table | `fevd_run` results |
| `make_irf_panel.m` | `fig_irf_page1/2.png` (12 shocks × 4 vars) | `fevd_run` results |
| `make_counterfactual.m` | `cf_mp_hike_100bp.png` | `fevd_run` results + `cond_*.mat` |
| `repair_mode_hessian.m` | `<mode>_fixed.mat` — surgical repair of all-NaN Hessian rows | a Dynare mode file |
| `mode_from_m1.m` | Rebuilds a mode file from `m1.mat` when a newrat search crashes mid-optimization | `m1.mat` + a reference mode |
| `slim_results.m` | Shrinks a post-estimation `*_results.mat` to ~20 MB | POSTEST results |
| `add_lw_to_external.m` | Fills the `LW_rstar` column of `rstar_external.csv` | the LW `.xlsx` |
| `peek_lw.m` | Column-by-column census of the LW workbook, for when the auto-detect picks the wrong series | the LW `.xlsx` |
| `diag_gap.m` | Output-gap sanity check and the gap-centering convention | POSTEST results |
| `diag_forecast.m` | Forecast / r\* / inflation-target diagnostic | POSTEST results |
| `diag_inflation_target.m` | Splits the longer-run inflation gap into the smoothed target `Pi_star` and the measurement error `me_inflexp` | POSTEST results |
| `diag_ogap_offset.m` | Tests the `ogap` vs computed-gap discrepancy against the predicted second-order linearisation error | POSTEST results |
| `diag_flexblock.m` | Does the flex block track? Amplitudes plus correlations in levels, first differences and the HP(1600) cycle | POSTEST results |
| `diag_fevd_split.m` | Per-shock-**type** FEVD for actual vs natural output, splitting technology from markup | `fevd_run` results |

Every script takes the results path as its first argument, so an explicit
`make_paper_figures('../results/frbstl_us_est_fastlr_lean_results.mat')` always works if the
auto-detect picks the wrong file.

### note/

| File | Purpose |
|---|---|
| `frbstl_replication_note.tex` | The companion note (methodology, estimation, results, assumptions, references). |
| `app_derivations.tex` | Full derivations, `\input` by the note as the appendix. |
| `figures/` | Where the note looks for its figures. Populate it from your own run — see §6. |

### results/

`results/` lets the figure scripts run with no Dynare step: they read only the post-estimation
`*_results.mat`, never the MH draws.

| File | Feeds |
|---|---|
| `frbstl_us_est_fastlr_lean_results.mat` | `make_paper_figures`, `check_postest`, the `diag_*` scripts, `make_detailed_decomp` |
| `fevd_run_results.mat` | `make_fevd`, `make_irf_panel`, `make_counterfactual`, `diag_fevd_split` |

The post-estimation file is slimmed by `scripts/slim_results.m`, which drops the unused smoother
containers (`oo_.Smoother/SmoothedShocks/Constant/UpdatedVariables`, ~50 MB each) and keeps
`oo_.SmoothedVariables.Mean`. Every field the scripts read survives. The MH draws
(`model/frbstl_us_est_fastlr_lean/metropolis/`, ~1 GB) are not shipped: they are needed only to
re-run the post-estimation step or to recover the full posterior distribution.

---

## 3. The one naming rule that matters

Dynare writes the MH draws into a folder named after the running `.mod`
(`frbstl_us_est_fastlr_lean/metropolis/`), and the post-estimation step finds them via
`load_mh_file`, which looks in a folder with the **same base name**. The step files carry distinct
names only to keep them apart on disk. At run time each is copied onto the single base name
`frbstl_us_est_fastlr_lean.mod`:

```bat
cd model
copy frbstl_us_est_fastlr_lean_STEP3_MCMC.mod  frbstl_us_est_fastlr_lean.mod
dynare frbstl_us_est_fastlr_lean parallel conffile=cluster.ini

copy frbstl_us_est_fastlr_lean_POSTEST.mod     frbstl_us_est_fastlr_lean.mod   :: overwrite
dynare frbstl_us_est_fastlr_lean                                               :: reuses the draws
```

`frbstl_us_est_fastlr_lean.mod` is therefore a scratch copy; it is in `.gitignore` and is not
tracked.

The other hard-coded names:

- **`mode_file`.** Steps 2 and 3, POSTEST and `fevd_run` all load
  `frbstl_us_est_fastlr_mode_ndc_fixed`. Rename the mode and you must update the `mode_file =`
  option in each of them.
- **`datafile = 'usdata_1959_2026.mat'`.** Appears in the `estimation(...)` call of every `.mod`,
  and resolves relative to `model/`.
- **Results path.** The scripts search, in order, `results/`, `../results/`, and the Dynare `Output/`
  folder. If you change the base name, pass the results path explicitly.
- **Figure filenames.** The scripts write PNGs into `model/figures_paper/` under fixed names and the
  note includes them by those names. Mapping in §6.
- **Conditional-scenario files.** POSTEST writes `cond_constrate.mat`, `cond_infltarget.mat` and
  `cond_tighter.mat` into the run folder; `make_paper_figures` (Fig 8) and `make_counterfactual`
  look for them by name in the current MATLAB folder, so stay in `model/`.

---

## 4. End-to-end run order

Runtimes are what this machine actually took (8 physical cores, `use_dll`). Everything below runs
with `model/` as the working folder.

```
(0) data      python data/check_vintage.py            confirm the sample is complete   [seconds]
              python data/build_usdata.py             -> usdata_1959_2026.mat          [optional]
              copy the .mat into model/

(1) mode      copy ..._STEP1_MODE.mod   -> ..._lean.mod
              dynare frbstl_us_est_fastlr_lean        -> ..._mode_ndc.mat

(2) repair    (MATLAB) repair_mode_hessian('frbstl_us_est_fastlr_mode_ndc.mat')
                                                      -> ..._mode_ndc_fixed.mat        [seconds]

(3) pilot     copy ..._STEP2_PILOT.mod  -> ..._lean.mod
              dynare frbstl_us_est_fastlr_lean parallel conffile=cluster.ini
                                                      -> tuned mh_jscale               [~10 min]

(4) MCMC      paste the tuned jscale into ..._STEP3_MCMC.mod
              copy ..._STEP3_MCMC.mod   -> ..._lean.mod
              dynare frbstl_us_est_fastlr_lean parallel conffile=cluster.ini
                                                      -> posterior draws               [~2 h 36 min]

(5) postest   copy ..._POSTEST.mod      -> ..._lean.mod
              dynare frbstl_us_est_fastlr_lean        -> smoother/forecast/decomp      [~18 min]

(6) fast run  dynare fevd_run                         -> FEVD + IRFs                   [~2 min]

(7) checks    (MATLAB) check_postest                  gate before anything downstream

(8) figures   (MATLAB) make_paper_figures; make_detailed_decomp; make_fevd;
                       make_irf_panel; make_counterfactual;
                       diag_gap; diag_forecast; diag_inflation_target;
                       diag_ogap_offset; diag_flexblock; diag_fevd_split

(9) external  add_lw_to_external           (§7, then re-run make_paper_figures for Fig 6)

(10) archive  (MATLAB) slim_results('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
                                    '../results/frbstl_us_est_fastlr_lean_results.mat')
              copy fevd_run\Output\fevd_run_results.mat ..\results\

(11) note     copy model\figures_paper\*.png -> note\figures\ ;  pdflatex ×2            (§6)
```

Steps 1, 3, 4 and 5 all run under the same base name (§3). Step 6 is independent and cheap: it is
the quickest way to regenerate the FEVD and IRFs without touching the 18-minute smoother.

Step 4 does not have to be one wave. `mh_replic` can be raised in stages and continued with
`load_mh_file`; the shipped run took the eight chains to 250,000 draws in a single pass.

---

## 5. The mode search, and what to do when it breaks

Two things go wrong on this model, both worth knowing before you burn a day on them.

**(a) The mode search walks into a complex steady state.** `alpha`, `nu_C` and `nu_Y` enter the
steady state through expressions that are undefined outside a bounded region. When newrat's line
search steps outside it, `static_resid` throws *"y must be a real dense numeric array"* and the run
dies mid-optimization.

Adding bounds in `estimated_params` does **not** fix this. Dynare's `newrat` (`mode_compute = 5`)
alternates Sims' **unconstrained** `csminit1` line search with the bounded univariate `mr_gstep`;
the `estimated_params` bounds constrain `mr_gstep` and the prior support, not the line search. A
bounded re-run is bit-identical to an unbounded one, which is what `STEP1c_BOUNDED.mod` is there to
demonstrate. The bounds are kept in all the files anyway: they are harmless and they document the
admissible region.

What does work is capping the iteration count. The newrat trajectory is deterministic and its
per-iteration improvements decay geometrically (0.104, 0.118, 0.123, 0.094, 0.043, 0.026, 0.013), so
`optim = ('MaxIter', 7)` stops one iteration short of the crash with the objective already flat.
That is what `STEP1_MODE.mod` ships with, and it returns a log posterior of **6807.574866**.

If a search does crash, `scripts/mode_from_m1.m` rebuilds a usable mode file from the `m1.mat`
newrat writes at the end of every completed iteration.

**(b) The Hessian comes back with an all-NaN row.** The near-unit-root persistences sit close enough
to the stationarity boundary that the numerical second difference is undefined on one side. In the
shipped run exactly one row broke, `rho_mup`. A single NaN row makes the whole matrix unusable as an
MH proposal even though every other direction is fine, so `scripts/repair_mode_hessian.m` keeps the
good block exactly as estimated and gives each broken parameter an independent proposal with a small
step (default implied sd 0.0010). That preserves all the posterior correlations, which
`MCMC_jumping_covariance = prior_variance` would throw away.

One consequence: the **Laplace approximation** to the marginal likelihood is computed from that
repaired matrix and is not a clean model-comparison statistic. Use the **modified harmonic mean**,
which the shipped run puts at **−7050.235409**.

---

## 6. Figures the note expects

| Note figure | Filename | Produced by |
|---|---|---|
| Fig 4a / 4b | `fig4a_core_pce_yoy.png`, `fig4b_headline_pce_yoy.png` | `make_paper_figures` |
| Fig 5a / 5b | `fig5a_output_gap.png`, `fig5b_rstar_spot.png` | `make_paper_figures` |
| Fig 6 (r\*) | `fig6_rstar.png` | `make_paper_figures` (after `add_lw_to_external`) |
| Fig 7 / 8 | `fig7_uncond_forecast.png`, `fig8_cond_scenarios.png` | `make_paper_figures` |
| Figs 9–10 | `fig9a_core_detailed.png`, `fig9b_headline_detailed.png`, `fig10a_gap_detailed.png`, `fig10b_rstar_detailed.png` | `make_detailed_decomp` |
| FEVD | `fevd_gdp.png`, `fevd_obs_pi.png`, `fevd_obs_ffr.png`, `fevd_rstar.png`, `fevd_ogap.png`, `fevd_iv.png` | `make_fevd` |
| IRFs | `fig_irf_page1.png`, `fig_irf_page2.png` | `make_irf_panel` |
| Counterfactual | `cf_mp_hike_100bp.png` | `make_counterfactual` |

The scripts save into `model/figures_paper/`. Copy those PNGs into `note/figures/` under the same
names, then:

```
cd note
pdflatex frbstl_replication_note.tex
pdflatex frbstl_replication_note.tex     :: second pass for the ToC and cross-references
```

The note uses `\graphicspath{{figures/}}`, so only the folder name matters, not where
`figures_paper/` lives.

---

## 7. External r\* data for Figure 6

`model/rstar_external.csv` ships with the **Lubik–Matthes** median filled (1985Q1–2026Q1) and the
**Laubach–Williams** one-sided estimate filled from the workbook vintage in `model/`. To refresh it:

1. Download `Laubach_Williams_current_estimates.xlsx` from the New York Fed's r\* page and put it in
   `model/`, replacing the shipped copy.
2. In MATLAB, from `model/`: `add_lw_to_external` — rewrites the `LW_rstar` column.
3. Re-run `make_paper_figures`; Fig 6 picks up the new series.

The workbook's layout changes between vintages, so the script prints what it detected and refuses to
guess: a date it cannot parse becomes a hole, never a value filed under the wrong quarter. If the
auto-detect picks the wrong column — the workbook carries the one-sided and two-sided estimates,
trend growth `g`, and the `z` component, and they do not all span the same rows — run `peek_lw` for a
column-by-column census, then pass the sheet and columns explicitly:

```matlab
peek_lw('Laubach_Williams_current_estimates.xlsx')
add_lw_to_external('Laubach_Williams_current_estimates.xlsx', 'data', 1, COL)
```

The figure legend and the note both label this series **Laubach–Williams (one-sided)**.

---

## 8. Reproducibility notes

- **Software:** Dynare 7.1, `use_dll`, an 8-worker local pool. Different Dynare or BLAS versions can
  shift low digits.
- **The shipped run:** mode log posterior 6807.574866 (7 newrat iterations); pilot 8 × 5,000 draws
  in 10 min 16 s, acceptance 31.4–33.9%, tuned `mh_jscale` 0.23046; full MCMC 8 × 250,000 draws in
  2 h 36 min 18 s, acceptance 32.10–32.46% in every chain, worst inefficiency factor 707 and a
  typical one near 230; modified harmonic mean −7050.235409; post-estimation 18 min 08 s.
- **Against the paper:** 59 of the 63 posterior means fall inside the paper's 10th-to-90th
  percentile bands (Tables 2 and 3). The four outside are `rho_me_inflexp`, `Gammabar_Z`, `eta_p`
  and `phi_Y`; the last two miss by a few hundredths of a band width. The note discusses the two
  substantive ones.
- **What is an assumption, not a model output:** the ≈0.4 pp population add-on in Table 5 and the
  growth decomposition (`POP_GROWTH` in `make_paper_figures.m`) is a flat placeholder for the Census
  population-growth trend. The flexible-price / r\* block and the gap-centering convention are set
  out in the note's assumptions section; read it before interpreting r\*, the gap or the forecast
  table.
- **Known data gap:** the SPF 10-year CPI series in `Mean_CPI10_Level.xlsx` begins 1991Q4, while the
  paper's footnote 13 describes the series as available from 1979Q4. The missing 48 quarters cover
  the Volcker disinflation, the most informative episode for identifying the time-varying inflation
  target, and this is the leading candidate for the difference between this replication's `Pi_star`
  path and the paper's. Documented in `data/build_usdata.py`.

---

## 9. Quick start

**No-Dynare path.** Once `results/` holds `frbstl_us_est_fastlr_lean_results.mat` and
`fevd_run_results.mat`:

```matlab
cd model
addpath('../scripts')
check_postest
make_paper_figures; make_detailed_decomp; make_fevd; make_irf_panel; make_counterfactual
diag_gap; diag_forecast; diag_inflation_target; diag_flexblock
```

then build the note (§6).

**From the posterior draws.** If you have `model/frbstl_us_est_fastlr_lean/metropolis/` and the
repaired mode:

```bat
cd model
copy frbstl_us_est_fastlr_lean_POSTEST.mod frbstl_us_est_fastlr_lean.mod
dynare frbstl_us_est_fastlr_lean
dynare fevd_run
```

then the MATLAB block above, then `copy figures_paper\*.png ..\note\figures\` and `pdflatex` twice.
