% ============================================================================
%  diag_ogap_offset.m
%
%  Is the gap-vs-ogap discrepancy linearisation error, or a real inconsistency?
%
%     computed gap = 100*(log(gdp) - log(gdp_f))   -- EXACT log difference of the
%                                                     smoothed LEVELS
%     model  ogap  = the smoothed value of eq. (384), which the first-order
%                    solution carries as the LINEARISED version of that same
%                    expression
%
%  So a discrepancy is not just expected, it is PREDICTABLE. Writing
%  xg = (gdp-gdp_ss)/gdp_ss and xf likewise, the exact-minus-linear difference is
%
%       100*[log(1+xg) - xg] - 100*[log(1+xf) - xf]  ~=  -50*(xg^2 - xf^2)
%
%  and that is a series this script can compute from the smoothed levels alone.
%  Comparing it against the OBSERVED discrepancy settles the question outright:
%  if they track, the two objects are the same thing to first order and nothing
%  is wrong. An arbitrary "std ratio < 5%" threshold cannot settle it and, on
%  this model, gives the wrong answer -- the honest test is the predicted series.
%
%  Presample: the estimation uses presample = 4 with lik_init = 2, so the first
%  four smoothed quarters carry an enormous prior variance and are not
%  identified. They are dropped throughout; 1959Q1 alone was responsible for the
%  whole apparent discrepancy on the first pass.
%
%     >> diag_ogap_offset
% ============================================================================
function diag_ogap_offset(RESULTS, NPRE)

if nargin<2 || isempty(NPRE); NPRE = 4; end
if nargin<1 || isempty(RESULTS)
  cand={'results/frbstl_us_est_fastlr_lean_results.mat', ...
        'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no results .mat found.');
end
S=load(RESULTS); oo_=S.oo_; M_=S.M_; endo=cellstr(M_.endo_names);
SV=oo_.SmoothedVariables; if isfield(SV,'Mean'); SV=SV.Mean; end
for f={'gdp','gdp_f','ogap'}
  assert(isfield(SV,f{1}), '%s not in the smoother var list', f{1});
end
ss=@(n) oo_.steady_state(find(strcmp(n,endo),1));

g=SV.gdp(:); gf=SV.gdp_f(:); om=SV.ogap(:);
T0=numel(g); yr=1959+(0:T0-1)'/4;
% ---- trim EVERYTHING together (getting this wrong is how the last version
% ---- crashed: gap was trimmed, ogap was not) --------------------------------
k = NPRE+1:T0;
g=g(k); gf=gf(k); om=om(k); yr=yr(k); T=numel(g);

gss=ss('gdp'); gfss=ss('gdp_f');
gap = 100*(log(g)-log(gf));
d   = gap - om;

fprintf('\n=============== ogap OFFSET DIAGNOSTIC ===============\n');
fprintf('sample: %d quarters, %s to %s   (first %d dropped: presample)\n\n', ...
        T, qlab(yr(1)), qlab(yr(end)), NPRE);
fprintf('                       mean      std       min      max\n');
fprintf('  computed gap    %9.3f %8.3f %9.3f %8.3f\n', mean(gap), std(gap), min(gap), max(gap));
fprintf('  model  ogap     %9.3f %8.3f %9.3f %8.3f\n', mean(om),  std(om),  min(om),  max(om));
fprintf('  difference      %9.3f %8.3f %9.3f %8.3f\n', mean(d), std(d), min(d), max(d));
fprintf('  corr(gap,ogap)  %9.4f\n\n', corr2v(gap,om));

% ---- THE TEST: predicted second-order error vs the observed discrepancy -----
xg = (g-gss)/gss;  xf = (gf-gfss)/gfss;
dhat = 100*((log(1+xg)-xg) - (log(1+xf)-xf));

fprintf('---- is the discrepancy just second-order error? ----\n');
fprintf('predicted  -50*(xg^2 - xf^2):  mean %+7.3f  std %6.3f\n', mean(dhat), std(dhat));
fprintf('observed   gap - ogap       :  mean %+7.3f  std %6.3f\n', mean(d),    std(d));
rho  = corr2v(dhat, d);
resid = d - dhat;
fprintf('correlation predicted vs observed = %.4f\n', rho);
fprintf('residual after removing the predicted error: std %.4f pp ', std(resid));
fprintf('(%.2f%% of the gap)\n\n', 100*std(resid)/std(gap));

if rho > 0.8 && std(resid) < 0.1*std(gap)
  fprintf('  VERDICT: LINEARISATION ERROR, nothing more.\n');
  fprintf('  The discrepancy is predicted by the second-order term to r = %.3f, and what\n', rho);
  fprintf('  remains is %.3f pp against a gap that moves %.2f pp. ogap and the computed\n', std(resid), std(gap));
  fprintf('  gap are the same object; Figure 5a and the Table 5 gap column are sound.\n');
elseif rho > 0.5
  fprintf('  VERDICT: MOSTLY linearisation error (r = %.3f), but %.3f pp of unexplained\n', rho, std(resid));
  fprintf('  residual remains. Worth a look if the gap column matters to a decimal.\n');
else
  fprintf('  VERDICT: NOT explained by the second-order term (r = %.3f). The smoothed\n', rho);
  fprintf('  ogap is genuinely not the log-difference of the smoothed levels -- inspect\n');
  fprintf('  eq. (384) and the gdp_f block before trusting Figure 5a.\n');
end

% ---- the level wedge, and why the correlation in levels looks low ----------
fprintf('\n---- the level wedge (constant-markup convention) ----\n');
pn=cellstr(M_.param_names); i=find(strcmp('rho_mup',pn),1);
if ~isempty(i)
  r=M_.params(i);
  fprintf('rho_mup = %.4f -> half-life %.0f quarters (%.0f years) vs a %.0f-year sample.\n', ...
          r, log(.5)/log(r), log(.5)/log(r)/4, T/4);
  fprintf('The flex block holds mu_p at steady state, so this near-random-walk shock moves\n');
  fprintf('gdp but not gdp_f. The gap absorbs the wandering: that is why its mean is %+.2f pp\n', mean(gap));
  fprintf('rather than zero, and why the LEVEL correlation understates the fit.\n\n');
end
cl = corr2v(g, gf);
dg = 400*diff(log(g)); dgf = 400*diff(log(gf));
cg = corr2v(dg, dgf);
fprintf('  corr(gdp, gdp_f)          levels        = %.3f\n', cl);
fprintf('  corr(gdp, gdp_f)          growth rates  = %.3f   <-- the honest statistic\n', cg);
fprintf('  (levels share a wandering component only gdp has; growth rates do not,\n');
fprintf('   so the growth correlation is what says whether the flex block tracks.)\n');

fprintf('\n  decade means of the computed gap:\n');
for y0=1960:10:2020
  sel = yr>=y0 & yr<y0+10;
  if any(sel); fprintf('    %ds  %+6.2f pp\n', y0, mean(gap(sel))); end
end
fprintf('=====================================================\n');
end

function s=qlab(d); y=floor(d+1e-9); q=round((d-y)*4)+1; s=sprintf('%dQ%d',y,q); end
function c=corr2v(a,b)
  a=a(:); b=b(:); n=min(numel(a),numel(b));
  m=corrcoef(a(1:n),b(1:n)); c=m(1,2);
end
