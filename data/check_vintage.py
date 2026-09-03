"""
check_vintage.py -- read-only availability probe.

Answers one question: for each input of build_usdata.py, what is the last
quarter that actually has a number in it today?  Writes nothing, touches
neither usdata_1959_2026.mat nor .csv, and does not re-demean anything.

    export FRED_API_KEY=<your key>
    python check_vintage.py

The binding constraint for a new row is always the SLOWEST series in the list.
If any line prints a quarter earlier than the target, the row would be built
with a hole in it, and Dynare would silently treat that cell as missing.
"""
import os, sys
import pandas as pd
from fredapi import Fred

TARGET = "2026Q2"          # the quarter you want to add

fred = Fred(api_key=os.environ.get("FRED_API_KEY"))

FRED_IDS = {
    "GDP": "GDP", "GDPDEF": "GDPDEF", "POP": "CNP16OV", "PCE": "PCE",
    "GPDI": "GPDI", "GCE": "GCE", "COMPNFB": "COMPNFB", "AWH": "AWHNONAG",
    "EMP": "CE16OV", "FF": "FEDFUNDS", "PCEPILFE": "PCEPILFE",
    "WTI": "WTISPLC", "GASPI": "DGOERG3Q086SBEA", "DGS10": "DGS10",
    "AAA": "AAA", "DGS20": "DGS20",
    "TR (soc. benefits)": "BOGZ1FA366404005Q",
    "TR (other transf.)": "BOGZ1FA366403005Q",
    "TR (subsidies)":     "BOGZ1FA366402005Q",
}

def last_q(s):
    s = s.dropna()
    if s.empty:
        return None
    s.index = pd.to_datetime(s.index)
    return pd.Period(s.resample("QS").mean().dropna().index[-1], freq="Q")

target = pd.Period(TARGET, freq="Q")
rows = []

print("FRED series")
for label, sid in FRED_IDS.items():
    try:
        lq = last_q(fred.get_series(sid, observation_start="2020-01-01"))
    except Exception as e:
        lq = None
        print(f"  {label:<20s} {sid:<22s} ERROR: {e}")
        continue
    ok = (lq is not None and lq >= target)
    rows.append((label, lq, ok))
    print(f"  {label:<20s} {sid:<22s} last = {lq}   {'OK' if ok else '<-- SHORT'}")

print("\nNon-FRED series")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from load_fernald_spf import load_fernald_tfp, load_spf_cpi10
    for label, fn in (("Fernald TFP", load_fernald_tfp), ("SPF CPI10", load_spf_cpi10)):
        try:
            lq = last_q(fn())
            ok = (lq is not None and lq >= target)
            rows.append((label, lq, ok))
            print(f"  {label:<20s} {'':<22s} last = {lq}   {'OK' if ok else '<-- SHORT'}")
        except Exception as e:
            print(f"  {label:<20s} ERROR: {e}")
except ImportError as e:
    print("  could not import load_fernald_spf.py:", e)

short = [r for r in rows if not r[2]]
print("\n" + "=" * 62)
if short:
    print(f"{TARGET} is NOT complete.  Missing / short:")
    for label, lq, _ in short:
        print(f"  {label:<20s} stops at {lq}")
    print("Building the row now would put NaNs in those observables for that")
    print("quarter.  build_usdata.py's guard only fires on an ENTIRELY empty")
    print("column, so it would not catch this.")
else:
    print(f"{TARGET} is complete in every observable.  Safe to rebuild.")
print("=" * 62)
