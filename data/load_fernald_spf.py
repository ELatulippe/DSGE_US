"""
Loaders for the two non-FRED observables of the St. Louis Fed DSGE model:
  * obs_tfp     -- Fernald (2012) utilization-adjusted, model-consistent TFP growth (eq. 106)
  * obs_inflexp -- SPF 10-year mean CPI inflation expectation, PCE-adjusted (eq. 109)

Both return a pandas Series indexed by quarterly period-start timestamps, ready to
drop into the `obs` DataFrame in build_usdata.py.  Column auto-detection is used
because the source files' headers change between vintages; each loader prints what
it detected so you can verify.

Data sources (public):
  Fernald: https://www.frbsf.org/wp-content/uploads/quarterly_tfp.xlsx
  SPF CPI10: https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/
             survey-of-professional-forecasters/data-files/files/Mean_CPI10_Level.xlsx
"""
import re, numpy as np, pandas as pd

# Some official .xlsx files (e.g. the SPF file) carry a core-properties "modified"
# date that openpyxl >=3.1 on Python 3.13/3.14 fails to parse, raising a TypeError
# before any data is read.  Those properties are irrelevant here, so make openpyxl
# fall back to empty properties instead of crashing.
try:
    from openpyxl.reader import excel as _oxl
    from openpyxl.packaging.core import DocumentProperties as _DP
    _orig_read_props = _oxl.ExcelReader.read_properties
    def _safe_read_properties(self):
        try:
            _orig_read_props(self)
        except Exception:
            self.wb.properties = _DP()
    _oxl.ExcelReader.read_properties = _safe_read_properties
except Exception:
    pass

FERNALD_URL = "https://www.frbsf.org/wp-content/uploads/quarterly_tfp.xlsx"
SPF_CPI10_URL = ("https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/"
                 "survey-of-professional-forecasters/data-files/files/Mean_CPI10_Level.xlsx")


def _yq_to_ts(year, quarter):
    """Year & quarter arrays -> (quarterly period-start timestamps, positional mask)."""
    year = pd.to_numeric(pd.Series(year).reset_index(drop=True), errors="coerce")
    quarter = pd.to_numeric(pd.Series(quarter).reset_index(drop=True), errors="coerce")
    ok = (year.notna() & quarter.notna()).to_numpy()
    pstr = [f"{int(y)}Q{int(q)}" for y, q in zip(year[ok], quarter[ok])]
    idx = pd.PeriodIndex(pstr, freq="Q").to_timestamp(how="start")
    return idx, ok


def _find_col(cols, patterns, exclude=()):
    """First column whose lowercased name contains any pattern and no exclude term."""
    for c in cols:
        name = str(c).lower()
        if any(p in name for p in patterns) and not any(x in name for x in exclude):
            return c
    return None


def _find_header_row(url, sheet, needle="dtfp", maxrows=8):
    """Some vintages have a title/description row above the code header.
    Return the row index whose cells contain `needle`, else 0."""
    probe = pd.read_excel(url, sheet_name=sheet, header=None, nrows=maxrows)
    for i in range(len(probe)):
        if probe.iloc[i].astype(str).str.lower().str.contains(needle).any():
            return i
    return 0


def load_fernald_tfp(url=FERNALD_URL, sheet="quarterly", alpha_const=1/3.0):
    """Model-consistent TFP growth, eq. (106):
         obs_tfp = dtfp/(1-alpha) - alpha/(1-alpha)*dutil   [annualized %]
       returned in QUARTERLY log-growth units (/400) and demeaned by the caller.
    """
    hdr = _find_header_row(url, sheet, needle="dtfp")
    raw = pd.read_excel(url, sheet_name=sheet, header=hdr).dropna(axis=1, how="all")

    # date column: values like "1947:Q2" / "1947Q2" / "1947:2"
    datecol = None
    for c in raw.columns:
        s = raw[c].astype(str)
        if s.str.match(r"^\s*\d{4}[:\sqQ]").mean() > 0.5:
            datecol = c; break
    if datecol is None:
        datecol = raw.columns[0]
    d = raw[datecol].astype(str).str.extract(r"(\d{4})\D*([1-4])")
    idx, ok = _yq_to_ts(d[0], d[1])

    # prefer exact dtfp (not utilization-adjusted) ; utilization ; capital share
    c_dtfp = _find_col(raw.columns, ["dtfp"], exclude=["util", "_util"]) \
             or _find_col(raw.columns, ["dtfp"])
    c_dutil = _find_col(raw.columns, ["dutil"]) or _find_col(raw.columns, ["util"], exclude=["dtfp"])
    c_alpha = _find_col(raw.columns, ["alpha", "capital share", "share of capital", "capital_share"])

    dtfp = pd.to_numeric(raw[c_dtfp], errors="coerce").to_numpy()[ok]
    dutil = pd.to_numeric(raw[c_dutil], errors="coerce").to_numpy()[ok]
    if c_alpha is not None:
        alpha = pd.to_numeric(raw[c_alpha], errors="coerce").to_numpy()[ok]
        if np.nanmean(alpha) > 1.0:            # header gave percent (e.g. 33 -> 0.33)
            alpha = alpha / 100.0
    else:
        alpha = alpha_const
    print(f"[Fernald] header row={hdr}; dtfp='{c_dtfp}', dutil='{c_dutil}', "
          f"alpha={'col '+repr(c_alpha) if c_alpha is not None else f'const {alpha_const:.3f}'}")

    tfp_annual = dtfp/(1-alpha) - alpha/(1-alpha)*dutil     # annualized %
    s = pd.Series(tfp_annual/400.0, index=idx).sort_index()      # quarterly log growth
    return s[~s.index.duplicated(keep="first")]                  # drop any repeated quarters


def load_spf_cpi10(url=SPF_CPI10_URL, pce_adjust=0.5):
    """SPF 10y-ahead mean CPI inflation (annualized %), mapped to headline PCE by
       subtracting `pce_adjust` (Del Negro et al. 2017).  Eq. (109) observable."""
    raw = pd.read_excel(url).dropna(axis=1, how="all")
    up = {str(c).strip().upper(): c for c in raw.columns}
    ycol = up.get("YEAR")   or _find_col(raw.columns, ["year"])
    qcol = up.get("QUARTER") or _find_col(raw.columns, ["quarter", "qtr"])
    vcol = up.get("CPI10")  or _find_col(raw.columns, ["cpi10", "cpi_10", "cpi 10"])
    print(f"[SPF] year='{ycol}', quarter='{qcol}', value='{vcol}'")
    idx, ok = _yq_to_ts(raw[ycol], raw[qcol])
    val = pd.to_numeric(raw[vcol], errors="coerce").to_numpy()[ok]
    s = pd.Series(val, index=idx).sort_index() - pce_adjust
    return s[~s.index.duplicated(keep="first")]


if __name__ == "__main__":
    # quick live smoke test (requires internet)
    tfp = load_fernald_tfp()
    exp = load_spf_cpi10()
    print("\nFernald obs_tfp  :", tfp.dropna().index.min().date(), "->",
          tfp.dropna().index.max().date(), "| sd(quarterly)=%.4f" % tfp.std())
    print("SPF obs_inflexp  :", exp.dropna().index.min().date(), "->",
          exp.dropna().index.max().date(), "| mean=%.2f%%" % exp.mean())