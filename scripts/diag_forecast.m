% ============================================================================
%  diag_forecast.m  --  localise WHY the forecast table differs from the paper.
%  Splits each gap into (i) the steady state the forecast reverts to, (ii) the
%  smoothed 2026Q1 jump-off state, and (iii) the forecast path. Reads the postest
%  results (needs rstar_ann, obs_pi, obs_ffr in the smoother/forecast var list;
%  Pi_star / obs_inflexp too if you added them).
%     >> diag_forecast
%     >> diag_forecast('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat')
% ============================================================================
function diag_forecast(RESULTS)
if nargin<1 || isempty(RESULTS)
  cand={'results/frbstl_us_est_fastlr_lean_results.mat', 'results/fevd_run_results.mat', 'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat','frbstl_us_est_fastlr_lean_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no results .mat -- run the postest first.');
end
S=load(RESULTS); oo_=S.oo_; M_=S.M_; endo=cellstr(M_.endo_names);
SV=struct(); if isfield(oo_,'SmoothedVariables'); if isfield(oo_.SmoothedVariables,'Mean'); SV=oo_.SmoothedVariables.Mean; else; SV=oo_.SmoothedVariables; end; end
FC=struct(); if isfield(oo_,'MeanForecast')&&isfield(oo_.MeanForecast,'Mean'); FC=oo_.MeanForecast.Mean;
elseif isfield(oo_,'forecast')&&isfield(oo_.forecast,'Mean'); FC=oo_.forecast.Mean; end
ss=@(nm) oo_.steady_state(find(strcmp(nm,endo),1));
sm=@(nm) (isfield(SV,nm))*0 + reshape(getfield_or(SV,nm),[],1);   %#ok
gv=@(st,nm) reshape(getfield_or(st,nm),[],1);

Tsm = numel(gv(SV,'obs_pi')); JOFF = 1959 + (Tsm-1)/4;   % last data quarter

% ---- "longer run" = CALENDAR 2031-2035, not a fixed offset from the jump-off.
% Forecast index h has date JOFF + h/4, so the window has to be recomputed each
% vintage; otherwise a 2026Q2 jump-off silently reports 2031Q2-2036Q2 and the
% row stops being comparable across vintages.  make_paper_figures.m already
% anchors on calendar years via yravg(...,2031,2035,...); this matches it.
LRY0 = 2031; LRY1 = 2035;
NFC  = numel(gv(FC,'obs_pi'));
[LRA, LRB] = lrwin(JOFF, NFC, LRY0, LRY1);

fprintf('\n===================== FORECAST DIAGNOSTIC (model vs paper) =====================\n');
fprintf('jump-off = %s.  Reads: SS the forecast reverts to | smoothed jump-off | fcst long run\n', qlabel(JOFF));
if isnan(LRA)
  fprintf('longer run = %d-%d : NOT REACHED by the current forecast horizon (%d quarters).\n\n', LRY0, LRY1, NFC);
else
  fprintf('longer run = %d-%d  ->  forecast quarters %d-%d (%s to %s)\n\n', ...
          LRY0, LRY1, LRA, LRB, qlabel(JOFF+LRA/4), qlabel(JOFF+LRB/4));
end

rows = { 'rstar_ann','r* (annualized)',        'lin',  1.1; ...
         'obs_pi',   'Core PCE inflation',     'lin',  2.5; ...
         'obs_ffr',  'Federal funds rate',     'lin',  3.8 };
% (also handle rstar gross if rstar_ann missing)
for i=1:size(rows,1)
  nm=rows{i,1}; lab=rows{i,2}; papLR=rows{i,4};
  if ~isfield(SV,nm); fprintf('  %-20s : not in var list (add %s to the estimation var list)\n',lab,nm); continue; end
  jo = last(gv(SV,nm));                         % smoothed jump-off
  f  = gv(FC,nm);                               % forecast path
  lr = mean_or(f, LRA, LRB);                    % calendar 2031-35 average
  s  = ss(nm);
  fprintf('  %-20s : SS %6.2f | jump-off %6.2f | fcst long-run %6.2f   (paper LR %.1f)\n', lab, s, jo, lr, papLR);
end

% ---- inflation target (if Pi_star was added to the var list) ---------------
fprintf('\n');
if isfield(SV,'Pi_star')
  pj = 400*log(last(gv(SV,'Pi_star')));  ps = 400*log(ss('Pi_star'));
  pf = gv(FC,'Pi_star'); plr = 400*log(mean_or(pf,LRA,LRB));
  fprintf('  Inflation target Pi* : SS %.2f | jump-off %.2f | fcst long-run %.2f  (paper ~2.5 elevated)\n', ps, pj, plr);
  dpi = pj - ps;
  if dpi >= 0.15
    fprintf('   -> jump-off Pi* is %+.2f pp ABOVE Pi-bar: elevated, same direction as the paper.\n', dpi);
  elseif dpi <= -0.15
    fprintf('   -> jump-off Pi* is %+.2f pp BELOW Pi-bar: DEPRESSED, the OPPOSITE direction from the\n', dpi);
    fprintf('      paper''s elevated target. See diag_inflation_target for the split against me_inflexp.\n');
  else
    fprintf('   -> jump-off Pi* is at Pi-bar (%+.2f pp): the estimation reads no target deviation.\n', dpi);
  end
else
  fprintf('  Inflation target Pi*: not in var list. Long-run inflation reverts to Pibar = %.2f%%.\n', 400*log(1.02^0.25));
  fprintf('   -> paper long-run inflation is 2.5%% (elevated target); ours reverts to 2%%. Add Pi_star obs_inflexp\n');
  fprintf('      to the smoother var list + re-run to see the smoothed target directly.\n');
end
% 10y inflation expectations (the observable that identifies the target)
if isfield(SV,'obs_inflexp')
  ie_jo = last(gv(SV,'obs_inflexp')); ie_ss = ss('obs_inflexp');
  fprintf('  10y infl. expectations (obs_inflexp): SS %.2f | smoothed jump-off %.2f\n', ie_ss, ie_jo);
  if ie_jo <= ie_ss + 0.15
    fprintf('   -> the expectations DATA is not elevated at the jump-off (%.2f vs SS %.2f), so nothing in\n', ie_jo, ie_ss);
    fprintf('      the observable pushes the target up. If the paper reads an elevated target off the same\n');
    fprintf('      series, the difference is in the series itself (sample start, source, or splice).\n');
  else
    fprintf('   -> expectations ARE elevated but the target did not rise -> the model dumps it into measurement\n');
    fprintf('      error: tighten sig_me_inflexp (or loosen sig_pistar) so the TARGET absorbs it, then re-estimate.\n');
  end
end

% ---- FFR gap decomposition -------------------------------------------------
if isfield(SV,'obs_ffr') && isfield(SV,'rstar_ann') && isfield(SV,'obs_pi')
  ffrLR = mean_or(gv(FC,'obs_ffr'),LRA,LRB);
  rsLR  = mean_or(gv(FC,'rstar_ann'),LRA,LRB);
  piLR  = mean_or(gv(FC,'obs_pi'),LRA,LRB);
  fprintf('\n  FFR long-run check:  r*(%.2f) + inflation(%.2f) approx = %.2f ; model FFR = %.2f\n', rsLR, piLR, rsLR+piLR, ffrLR);
  fprintf('  paper FFR 3.8 = r*(1.1) + infl(2.5) + ~0.2.  Our FFR gap is downstream of the r* and inflation gaps.\n');
end
fprintf('\n  For the inflation gap specifically, run  diag_inflation_target  -- it splits the\n');
fprintf('  2.0 vs 2.5 difference into the smoothed target Pi* and the measurement error on\n');
fprintf('  obs_inflexp, and scales the required wedge by the unconditional SD of me_inflexp.\n');
fprintf('==============================================================================\n');
end

% helpers
function v=getfield_or(s,nm); if isfield(s,nm); v=s.(nm); else; v=[]; end; end
function y=last(x); if isempty(x); y=NaN; else; y=x(end); end; end
function m=mean_or(x,a,b)
x=x(:);
if isnan(a) || isnan(b) || numel(x)<a; m=NaN; return; end
b=min(b,numel(x)); m=mean(x(a:b));
end
function [a,b]=lrwin(joff, nfc, Y0, Y1)
% Forecast-index window covering calendar Y0Q1..Y1Q4.  Index h has date joff+h/4.
a = round((Y0        - joff)*4);
b = round((Y1 + 0.75 - joff)*4);
a = max(a,1); b = min(b,nfc);
if b < a; a=NaN; b=NaN; end
end
function s=qlabel(d)
y = floor(d + 1e-9); q = round((d - y)*4) + 1;
s = sprintf('%dQ%d', y, q);
end
