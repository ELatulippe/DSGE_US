% ============================================================================
%  check_postest.m
%
%  Gate between the POSTEST run and everything downstream.  The figure scripts
%  and the diagnostics all read oo_ out of the POSTEST results file, and each
%  of them fails in its own confusing way if a variable it wants was never put
%  in the smoother list, or if the forecast horizon stops short of the calendar
%  window Table 5 needs.  Cheaper to find that here than three scripts later.
%
%     >> check_postest
%     >> check_postest('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat')
% ============================================================================
function check_postest(RESULTS)

if nargin<1 || isempty(RESULTS)
  cand = {'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
          'results/frbstl_us_est_fastlr_lean_results.mat', ...
          'frbstl_us_est_fastlr_lean_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no POSTEST results .mat found -- run the postest first.');
end
S=load(RESULTS); oo_=S.oo_; M_=S.M_;
fprintf('\n================== POSTEST COMPLETENESS CHECK ==================\n');
fprintf('results : %s\n', RESULTS);
fprintf('model   : %d endogenous variables\n\n', M_.endo_nbr);

ok = true;

% ---- 1. did it reuse the STEP3 draws, or silently resample? ----------------
if isfield(oo_,'MarginalDensity') && isfield(oo_.MarginalDensity,'ModifiedHarmonicMean')
  mhm = oo_.MarginalDensity.ModifiedHarmonicMean;
  fprintf('modified harmonic mean : %.6f\n', mhm);
  fprintf('  STEP3 reported       : -7050.235409\n');
  if abs(mhm - (-7050.235409)) < 1e-3
    fprintf('  MATCH -> load_mh_file reused the STEP3 chains. No resampling.\n\n');
  else
    fprintf('  *** DIFFERENT. The postest did NOT reuse the STEP3 draws.\n');
    fprintf('  *** Check mh_replic=0 and load_mh_file before trusting anything below.\n\n');
    ok = false;
  end
end

% ---- 2. the smoother variables the downstream scripts index by name --------
SV = struct();
if isfield(oo_,'SmoothedVariables')
  if isfield(oo_.SmoothedVariables,'Mean'); SV=oo_.SmoothedVariables.Mean;
  else;                                     SV=oo_.SmoothedVariables; end
end
need = { 'obs_pi','obs_ffr','obs_gdp',           'core observables'
         'rstar_ann','ogap','gdp_f',             'natural-rate / gap block'
         'Pi_star','obs_inflexp','me_inflexp',   'inflation-target diagnostic' };
fprintf('SMOOTHED VARIABLES\n');
for r = 1:size(need,1)
  grp = need{r,4-1+1};
  names = need(r,1:3);
  miss = names(~cellfun(@(n) isfield(SV,n), names));
  if isempty(miss)
    fprintf('  OK      %-32s %s\n', grp, strjoin(names,', '));
  else
    fprintf('  MISSING %-32s %s\n', grp, strjoin(miss,', '));
    ok = false;
  end
end

% ---- 3. forecast horizon vs the calendar window Table 5 needs --------------
FC = struct();
if     isfield(oo_,'MeanForecast') && isfield(oo_.MeanForecast,'Mean'); FC=oo_.MeanForecast.Mean;
elseif isfield(oo_,'forecast')     && isfield(oo_.forecast,'Mean');     FC=oo_.forecast.Mean; end
fprintf('\nFORECAST\n');
if isempty(fieldnames(FC))
  fprintf('  MISSING  no forecast in oo_ -- the `forecast = N` option did not run.\n'); ok=false;
else
  Tsm  = numel(reshape(SV.obs_pi,[],1));
  JOFF = 1959 + (Tsm-1)/4;
  nfc  = numel(reshape(FC.obs_pi,[],1));
  last = JOFF + nfc/4;
  fprintf('  smoothed sample ends : %s   (%d quarters from 1959Q1)\n', qlab(JOFF), Tsm);
  fprintf('  forecast horizon     : %d quarters -> reaches %s\n', nfc, qlab(last));
  if last >= 2035.75 - 1e-9
    fprintf('  OK  covers the 2031-2035 longer-run window Table 5 reports.\n');
  else
    fprintf('  *** TOO SHORT for the 2031-2035 longer-run row. Raise `forecast =`.\n');
    ok = false;
  end
end

% ---- 4. the conditional-forecast scenario files ---------------------------
fprintf('\nCONDITIONAL FORECASTS (Figure 8 scenarios)\n');
for f = {'cond_constrate.mat','cond_infltarget.mat','cond_tighter.mat'}
  if exist(f{1},'file'); fprintf('  OK      %s\n', f{1});
  else;                  fprintf('  MISSING %s  -> fig8_cond_scenarios will fail\n', f{1}); ok=false; end
end

fprintf('\n---------------------------------------------------------------\n');
if ok; fprintf('ALL CHECKS PASSED -- safe to run the figure and diagnostic scripts.\n');
else;  fprintf('SOMETHING IS MISSING -- fix it before running the downstream scripts,\n');
       fprintf('otherwise they fail late and it is hard to see which step was at fault.\n'); end
fprintf('===============================================================\n');
end

function s=qlab(d)
  y=floor(d+1e-9); q=round((d-y)*4)+1; s=sprintf('%dQ%d',y,q);
end
