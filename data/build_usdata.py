"""
Build the 15-observable dataset for the St. Louis Fed DSGE model (US, 1959Q1-).
Implements the data transforms of Section 4.2 and produces a Dynare-ready file
(usdata_1959_2026.mat and .csv) whose column names match the .mod's `varobs`.

USAGE
    pip install fredapi pandas numpy statsmodels
    export FRED_API_KEY=<your key>          # https://fred.stlouisfed.org/docs/api/api_key.html
    python build_usdata.py

TWO SERIES ARE NOT ON FRED; they are downloaded and parsed by load_fernald_spf.py
(imported below) directly from the source Excel files:
    * Fernald (2012) business-sector TFP + utilization (SF Fed quarterly_tfp.xlsx)
      -> observable TFP growth (106): dtfp/(1-alpha) - alpha/(1-alpha)*dutil, /400, demeaned.
        (If the file has no capital-share column, a constant alpha=1/3 is used; the
         difference is absorbed by the estimated TFP measurement error.)
    * SPF 10y mean CPI expectations (Philadelphia Fed Mean_CPI10_Level.xlsx),
      mapped to HEADLINE PCE by subtracting a constant 0.5pp -- the paper's
      footnote 13, following Del Negro et al. (2017): 0.5 is "the average
      difference between CPI and PCE inflation over the sample".  The series is
      an input to eq. (109), which is an expectation of Pi_h (headline), NOT of
      core Pi; core PCE enters separately as obs_pi (eq. 105).
      SAMPLE-START DISCREPANCY, UNRESOLVED: footnote 13 states the series is
      "only available from 1979Q4 onwards", but Mean_CPI10_Level.xlsx (SPF
      CPI10) begins 1991Q4, which is what this build produces.  The paper may
      splice an earlier long-range survey (Blue Chip's 10y CPI forecast starts
      ~1979Q4) for 1979Q4-1991Q3.  The missing 48 quarters cover the Volcker
      disinflation, the most informative episode for identifying the
      time-varying target Pi_star, so this gap is a candidate explanation for
      the difference between our Pi_star path and the paper's.  Ask the author.

The observation equations these must match (see the .mod, eqs 97-111):
    obs_gdp  = 400*dlog(gdp_pc_real)
    obs_c/i/g/tr = 400*dlog(real pc)          (+ cobs constant, calibrated to sample mean)
    obs_w    = 400*dlog(real wage)            (+ cobs)
    obs_n    = 400*dlog(hours pc)
    obs_ffr  = FEDFUNDS (already annualized %)
    obs_pi   = 400*dlog(PCEPILFE)
    obs_tfp  = Fernald model-consistent TFP growth, demeaned
    obs_oil  = log(real WTI), demeaned
    obs_gas  = 400*dlog(DGOERG3Q086SBEA price index)
    obs_inflexp = SPF 10y (headline) - 0.5
    obs_tr10 = DGS10                          (annualized %)
    obs_cy   = AAA - DGS20
"""
import os, numpy as np, pandas as pd
from fredapi import Fred

fred = Fred(api_key=os.environ.get("FRED_API_KEY"))
START = "1959-01-01"

def q(series_id, how="mean"):
    """Fetch a FRED series and aggregate to quarterly (period-average or last)."""
    s = fred.get_series(series_id, observation_start="1958-01-01")
    s.index = pd.to_datetime(s.index)
    if how == "mean":
        s = s.resample("QS").mean()
    else:
        s = s.resample("QS").last()
    return s

# ---- raw series ----
GDP     = q("GDP")            # nominal GDP
GDPDEF  = q("GDPDEF")         # GDP deflator
POP     = q("CNP16OV")        # civilian noninstitutional pop 16+
PCE     = q("PCE")            # nominal personal consumption
GPDI    = q("GPDI")           # nominal gross private domestic investment
GCE     = q("GCE")            # nominal govt consumption+investment
COMPNFB = q("COMPNFB")        # nonfarm business hourly compensation index
AWH     = q("AWHNONAG")       # avg weekly hours, total private (nonag)
EMP     = q("CE16OV")         # civilian employment
FF      = q("FEDFUNDS")       # effective federal funds rate (annualized %)
PCEPILFE= q("PCEPILFE")       # core PCE price index
WTI     = q("WTISPLC")        # spot WTI
GASPI   = q("DGOERG3Q086SBEA")# PCE gasoline & other energy price index
DGS10   = q("DGS10")          # 10y Treasury (annualized %)
AAA     = q("AAA")            # Moody's Aaa corporate yield (annualized %)
DGS20   = q("DGS20")          # 20y Treasury (annualized %)
# transfers: sum of social benefits + other current transfers + subsidies (Flow of Funds)
TR = (q("BOGZ1FA366404005Q") + q("BOGZ1FA366403005Q") + q("BOGZ1FA366402005Q"))

# ---- HP-filter population (removes census discontinuities), smoothed pop for scaling ----
from statsmodels.tsa.filters.hp_filter import hpfilter
pop_log = np.log(POP.dropna())
_, pop_trend = hpfilter(pop_log, lamb=1600)
POP_S = np.exp(pop_trend).reindex(GDP.index).ffill()

def real_pc_growth(nom, deflator=GDPDEF, pop=POP_S):
    real = (nom / deflator) / pop
    return 400 * np.log(real).diff()

obs = pd.DataFrame(index=GDP.index)
obs["obs_gdp"] = real_pc_growth(GDP)
obs["obs_c"]   = real_pc_growth(PCE)
obs["obs_i"]   = real_pc_growth(GPDI)
obs["obs_g"]   = real_pc_growth(GCE)
obs["obs_tr"]  = real_pc_growth(TR)
obs["obs_w"]   = 400 * np.log(COMPNFB / GDPDEF).diff()          # real wage growth (no pop)
HOURS = (AWH * EMP) / POP_S
obs["obs_n"]   = 400 * np.log(HOURS).diff()
obs["obs_ffr"] = FF
obs["obs_pi"]  = 400 * np.log(PCEPILFE).diff()
obs["obs_gas"] = 400 * np.log(GASPI).diff()
obs["obs_oil"] = np.log(WTI / GDPDEF)                            # real oil, demeaned below
obs["obs_tr10"]= DGS10
obs["obs_cy"]  = AAA - DGS20

# restrict to the estimation sample, then finalize the level / expectation series
obs = obs.loc[START:]

# obs_oil (107): real oil price, demeaned over the sample
obs["obs_oil"] = obs["obs_oil"] - obs["obs_oil"].mean()

# pre-create so obs[order] never KeyErrors even if a download/parse fails
obs["obs_tfp"] = np.nan
obs["obs_inflexp"] = np.nan

# obs_tfp (106): Fernald model-consistent, utilization-adjusted TFP growth (quarterly), demeaned
from load_fernald_spf import load_fernald_tfp   # hard dependency: see data/README
tfp = load_fernald_tfp()                           # quarterly log growth (deduped index)
obs["obs_tfp"] = tfp.reindex(obs.index)
obs["obs_tfp"] = obs["obs_tfp"] - obs["obs_tfp"].mean()

# obs_inflexp (109): SPF 10y-ahead CPI expectation - 0.5 (PCE)
from load_fernald_spf import load_spf_cpi10    # hard dependency: see data/README
obs["obs_inflexp"] = load_spf_cpi10().reindex(obs.index)

# ---- COVID special-shock observables (Section 4.3): 0 outside the active windows,
#      NaN (missing) inside, so the Kalman smoother infers the shock only then ----
def _win(a, b):
    return (obs.index >= pd.Timestamp(a)) & (obs.index <= pd.Timestamp(b))
obs["ec_chi"] = 0.0; obs.loc[_win("2020-01-01", "2020-07-01"), "ec_chi"] = np.nan  # chi: 2020Q1-Q3
obs["ec_gn"]  = 0.0; obs.loc[_win("2020-01-01", "2020-10-01"), "ec_gn"]  = np.nan  # Gamma_N: 2020Q1-Q4
obs["ec_tr"]  = 0.0; obs.loc[_win("2020-04-01", "2022-01-01"), "ec_tr"]  = np.nan  # transfers: 2020Q2-2022Q1

# ---- calibrate the cobs_* constants to sample-mean growth differentials ----
# obs_x mean should equal model trend 400*log(Gammabar); cobs_x = mean(obs_x) - mean(obs_gdp)
g_gdp = obs["obs_gdp"].mean()
cobs = {k: obs[k].mean() - g_gdp for k in ["obs_c","obs_i","obs_g","obs_tr","obs_w"]}
print("Calibrated cobs constants (paste into the .mod):")
for k,v in cobs.items():
    print(f"  cobs_{k.split('_')[1]}y = {v:.4f};")

# ---- save for Dynare (column order = varobs order) ----
order = ["obs_gdp","obs_c","obs_i","obs_g","obs_tr","obs_w","obs_n","obs_ffr",
         "obs_pi","obs_tfp","obs_oil","obs_gas","obs_inflexp","obs_tr10","obs_cy",
         "ec_chi","ec_tr","ec_gn"]
obs = obs[order]

# ---- fail loudly rather than writing a silently unusable dataset ----
_empty = [c for c in order if c not in ("ec_chi", "ec_tr", "ec_gn") and obs[c].notna().sum() == 0]
if _empty:
    raise RuntimeError(
        "These observables are entirely missing: " + ", ".join(_empty) +
        ".  Dynare would run on NaN columns and the estimation would be wrong.  "
        "Check the FRED key and that data/load_fernald_spf.py is importable."
    )

obs.to_csv("usdata_1959_2026.csv")
try:
    from scipy.io import savemat
    savemat("usdata_1959_2026.mat", {c: obs[c].values.reshape(-1,1) for c in order})
    print("Wrote usdata_1959_2026.mat and .csv  (Dynare: datafile=usdata_1959_2026)")
except Exception as e:
    print("Wrote .csv; install scipy for .mat:", e)
print(f"Sample: {obs.index[0].date()} .. {obs.index[-1].date()}, {len(obs)} quarters")
print("Non-missing counts per observable:")
print(obs.notna().sum().to_string())
print("Leading NaNs (hours pre-1964, SPF pre-1979Q4, etc.) are handled by Dynare's "
      "Kalman filter as missing values.")