% ============================================================================
%  diag_inflation_target.m
%
%  WHY does our longer-run Core PCE forecast revert to ~2.0% while the paper's
%  Table 5 reports 2.5%?  Both models share Pibar = 2%, both feed in the SAME
%  observable (SPF 10y CPI expectations minus 0.5pp, paper footnote 13), and
%  both have rho_pistar ~ 0.99.  So the difference has to show up in ONE of two
%  smoothed objects at the 2026Q1 jump-off:
%
%    (a) Pi_star     -- the time-varying inflation target, or
%    (b) me_inflexp  -- the measurement error on obs_inflexp (eq. 109).
%
%  Equation (109) is   obs_inflexp = 100*( exp(4*piq10) - 1 + me_inflexp ),
%  so 100*me_inflexp is the wedge, IN PERCENTAGE POINTS, between the model's
%  own 10y expected HEADLINE inflation and the survey series we hand it.  In
%  the long run headline converges to core, so a longer-run core forecast of F%
%  against an observable reading D% requires a persistent me_inflexp of
%  (D - F) pp.  This script computes that, and scales it by the unconditional
%  SD of me_inflexp implied by the estimated (rho, sigma), so you can see how
%  strained each reading is.
%
%  Reads the POSTEST results.  Pi_star, obs_inflexp and me_inflexp must be in
%  the smoother variable list (POSTEST already has Pi_star obs_inflexp; add
%  me_inflexp if it is missing and re-run the postest -- no MCMC needed, the
%  smoother reuses the stored draws via load_mh_file).
%
%     >> diag_inflation_target
%     >> diag_inflation_target('results/frbstl_us_est_fastlr_lean_results.mat')
% ============================================================================
function diag_inflation_target(RESULTS, OUTDIR)

if nargin<1 || isempty(RESULTS)
  cand={'results/frbstl_us_est_fastlr_lean_results.mat', ...
        'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
        'frbstl_us_est_fastlr_lean_results.mat', ...
        'results/fevd_run_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no results .mat found -- run the postest first.');
end
if nargin<2 || isempty(OUTDIR); OUTDIR='figures_paper'; end
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

S=load(RESULTS); oo_=S.oo_; M_=S.M_;
endo=cellstr(M_.endo_names); pnam=cellstr(M_.param_names);

SV=struct();
if isfield(oo_,'SmoothedVariables')
  if isfield(oo_.SmoothedVariables,'Mean'); SV=oo_.SmoothedVariables.Mean;
  else;                                     SV=oo_.SmoothedVariables; end
end
FC=struct();
if     isfield(oo_,'MeanForecast') && isfield(oo_.MeanForecast,'Mean'); FC=oo_.MeanForecast.Mean;
elseif isfield(oo_,'forecast')     && isfield(oo_.forecast,'Mean');     FC=oo_.forecast.Mean; end

gv = @(st,nm) reshape(getf(st,nm),[],1);
ss = @(nm) oo_.steady_state(find(strcmp(nm,endo),1));
pv = @(nm) par(M_,pnam,nm);

have_pi   = isfield(SV,'Pi_star');
have_obs  = isfield(SV,'obs_inflexp');
have_me   = isfield(SV,'me_inflexp');

Tsm  = numel(gv(SV,'obs_pi'));
dates = 1959 + (0:Tsm-1)'/4;          % 1959Q1 .. jump-off
JOFF  = dates(end);

fprintf('\n================= INFLATION-TARGET / EXPECTATIONS DIAGNOSTIC =================\n');
fprintf('results : %s\n', RESULTS);
fprintf('jump-off: %.2f  (%d smoothed quarters)\n\n', JOFF, Tsm);

% ---------------------------------------------------------------- 1. the wedge
rho = pv('rho_me_inflexp'); sig = pv('sig_me_inflexp');
if ~isnan(rho) && ~isnan(sig) && abs(rho)<1
  sd_pp = 100*sig/sqrt(1-rho^2);                 % unconditional SD, in pp
  hl    = log(0.5)/log(rho);
  fprintf('me_inflexp process:  rho = %.4f, sigma = %.4f\n', rho, sig);
  fprintf('  unconditional SD = %.3f pp on obs_inflexp;  half-life = %.1f quarters\n\n', sd_pp, hl);
else
  sd_pp = NaN;
  fprintf('me_inflexp process:  rho/sigma not found in M_.params.\n\n');
end

if have_obs
  D = gv(SV,'obs_inflexp'); D = D(end);
  fprintf('obs_inflexp at the jump-off (smoothed) : %.2f %%\n', D);
  fprintf('Persistent me_inflexp required to sustain a given longer-run core forecast:\n');
  for F = [2.0 2.5]
    need = D - F;
    if isnan(sd_pp)
      fprintf('   longer-run core %.1f%%  ->  me_inflexp = %+5.2f pp\n', F, need);
    else
      fprintf('   longer-run core %.1f%%  ->  me_inflexp = %+5.2f pp  = %.2f unconditional SD\n', ...
              F, need, abs(need)/sd_pp);
    end
  end
  fprintf('   (ours reverts to ~2.0; the paper''s Table 5 longer-run row is 2.5)\n\n');
end

% ------------------------------------------------- 2. the two smoothed series
fprintf('Smoothed path, 2019Q1 to the jump-off:\n');
fprintf('  %-8s %10s %12s %12s\n','quarter','Pi* (ann%)','me (pp)','obs_inflexp');
sel = find(dates >= 2019 - 1e-9);
for k = sel'
  q  = sprintf('%dQ%d', floor(dates(k)+1e-9), round((dates(k)-floor(dates(k)+1e-9))*4)+1);
  s1 = NaN; s2 = NaN; s3 = NaN;
  if have_pi;  v=gv(SV,'Pi_star');     s1 = 400*log(v(k)); end
  if have_me;  v=gv(SV,'me_inflexp');  s2 = 100*v(k);      end
  if have_obs; v=gv(SV,'obs_inflexp'); s3 = v(k);          end
  fprintf('  %-8s %10.2f %12.2f %12.2f\n', q, s1, s2, s3);
end
fprintf('\n');

% ------------------------------------------------------------- 3. the verdict
if have_pi
  v   = gv(SV,'Pi_star');
  pjo = 400*log(v(end));
  pss = 400*log(ss('Pi_star'));
  fprintf('Pi* : steady state %.2f  |  smoothed jump-off %.2f\n', pss, pjo);
  if have_me
    v = gv(SV,'me_inflexp'); mjo = 100*v(end);
    fprintf('me_inflexp at the jump-off: %+.2f pp\n\n', mjo);
    dgap = pjo - pss;   % SIGNED: >0 = above target, <0 = below target
    if abs(dgap) < 0.15 && abs(mjo) < 0.20
      fprintf('  VERDICT: our smoother reads NEITHER an elevated target NOR a large measurement\n');
      fprintf('  error -- it simply believes the survey.\n');
    elseif dgap >= 0.15
      fprintf('  VERDICT: our target is ELEVATED at the jump-off (%+.2f pp above Pi-bar), the same\n', dgap);
      fprintf('  direction as the paper.  The forecast difference is then about horizon and\n');
      fprintf('  convergence, not about identifying Pi*.\n');
    elseif dgap <= -0.15
      fprintf('  VERDICT: our target is DEPRESSED at the jump-off (%+.2f pp relative to Pi-bar),\n', dgap);
      fprintf('  which is the OPPOSITE direction from the paper''s elevated target.  Both smoothers\n');
      fprintf('  see the same ~1.9%% expectations series, so they are reading it in opposite ways.\n');
      fprintf('  With me_inflexp only %+.2f pp, ours is essentially believing the survey.\n', mjo);
      fprintf('  LEADING SUSPECT: the sample start of obs_inflexp.  Paper footnote 13 says the series\n');
      fprintf('  runs from 1979Q4; SPF CPI10 begins 1991Q4, which is what we feed in.  Those 48\n');
      fprintf('  quarters span the Volcker disinflation, the most informative episode for a\n');
      fprintf('  near-unit-root target.  Ask the author whether an earlier source is spliced in.\n');
    else
      fprintf('  VERDICT: the target sits at Pi-bar but the measurement error is large, i.e. our\n');
      fprintf('  smoother is dumping the expectations gap into me_inflexp rather than into Pi*.\n');
      fprintf('  Compare sig_me_inflexp / sig_pistar against the paper''s Table 2 posteriors.\n');
    end
  end
elseif ~have_pi
  fprintf('  Pi_star is not in the smoother variable list -- add it (and me_inflexp) to the\n');
  fprintf('  POSTEST var list and re-run the postest.  No MCMC needed.\n');
end

% ------------------------------------------------------------------ 4. figure
if have_pi || have_me
  figure('visible','off','position',[100 100 780 520],'Color','w');
  d = dates(sel);

  subplot(2,1,1); hold on; grid on; box on;
  if have_pi
    v=gv(SV,'Pi_star'); plot(d, 400*log(v(sel)), 'LineWidth', 1.6);
  end
  yline_compat(400*log(1.02^0.25), 'k:');
  yline_compat(2.5, 'r--');
  title('Smoothed inflation target \Pi^*_t (annualized %)');
  legend({'\Pi^*_t','\Pi-bar = 2%','paper long-run 2.5%'},'Location','best'); legend boxoff;
  xlim([d(1) d(end)]);

  subplot(2,1,2); hold on; grid on; box on;
  if have_me
    v=gv(SV,'me_inflexp'); plot(d, 100*v(sel), 'LineWidth', 1.6);
  end
  if ~isnan(sd_pp)
    yline_compat( sd_pp, 'k:'); yline_compat(-sd_pp, 'k:');
  end
  yline_compat(0,'k-');
  title('Smoothed measurement error on obs\_inflexp (pp; dotted = \pm1 uncond. SD)');
  xlabel('year'); xlim([d(1) d(end)]);

  set(gcf,'InvertHardcopy','off');
  fn=fullfile(OUTDIR,'diag_inflation_target.png');
  print(fn,'-dpng','-r120'); close;
  fprintf('\n  wrote %s\n', fn);
end
fprintf('==============================================================================\n');
end

% ------------------------------------------------------------------- helpers
function v=getf(s,nm); if isfield(s,nm); v=s.(nm); else; v=[]; end; end
function x=par(M_,pnam,nm)
  i=find(strcmp(nm,pnam),1);
  if isempty(i); x=NaN; else; x=M_.params(i); end
end
function yline_compat(y, sty)
% yline() only exists from R2018b; fall back to plot() on older releases.
  if exist('yline','builtin') || exist('yline','file')
    yline(y, sty);
  else
    xl=xlim; plot(xl, [y y], sty);
  end
end
