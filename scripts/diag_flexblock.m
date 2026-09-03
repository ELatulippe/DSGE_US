% ============================================================================
%  diag_flexblock.m  --  does the flexible-price block actually track?
%
%  WHY THIS EXISTS.  Two correlations have been used to answer that question and
%  both are poor:
%
%    corr in LEVELS  = 0.717.  Understated. A near-random-walk price-markup shock
%      (rho_mup = 0.9968, 54-year half-life) moves gdp but NOT gdp_f, because the
%      flex block holds mu_p at steady state. That puts a wandering component in
%      one series and not the other and drags the level correlation down for a
%      reason that says nothing about tracking.
%
%    corr in GROWTH RATES = 0.036.  Worse, and misleading. gdp and gdp_f are
%      already STATIONARY (both detrended by the same Gamma). Differencing a
%      stationary series OVER-differences it: it flattens the spectrum and hands
%      the correlation to the highest frequencies, which are exactly the
%      component least likely to be shared. A near-zero growth correlation here
%      is what over-differencing produces, not evidence of a broken block.
%
%  The right frequency band is the business cycle. This reports the correlation
%  of the HP-filtered cyclical components (lambda = 1600, quarterly) alongside
%  the other two, plus the amplitudes, so you can see whether gdp_f moves at all
%  and whether it moves with gdp where it should.
%
%     >> diag_flexblock
% ============================================================================
function diag_flexblock(RESULTS, NPRE)

if nargin<2 || isempty(NPRE); NPRE = 4; end
if nargin<1 || isempty(RESULTS)
  cand={'results/frbstl_us_est_fastlr_lean_results.mat', ...
        'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no results .mat found.');
end
S=load(RESULTS); oo_=S.oo_; M_=S.M_; endo=cellstr(M_.endo_names);
SV=oo_.SmoothedVariables; if isfield(SV,'Mean'); SV=SV.Mean; end
ss=@(n) oo_.steady_state(find(strcmp(n,endo),1));

g=SV.gdp(:); gf=SV.gdp_f(:); T0=numel(g); yr=(1959+(0:T0-1)'/4);
k=NPRE+1:T0; g=g(k); gf=gf(k); yr=yr(k); T=numel(g);
gd  = 100*log(g/ss('gdp'));      % sticky output, % dev from SS
gfd = 100*log(gf/ss('gdp_f'));   % flex   output, % dev from SS
gap = gd - gfd;

fprintf('\n============== FLEX-BLOCK TRACKING ==============\n');
fprintf('sample: %d quarters, %s to %s\n\n', T, qlab(yr(1)), qlab(yr(end)));

% ---- amplitudes: does gdp_f move at all? -----------------------------------
fprintf('AMPLITUDE (%% deviation from steady state)\n');
fprintf('  std(gdp)    = %6.2f\n', std(gd));
fprintf('  std(gdp_f)  = %6.2f   (%.0f%% of gdp)\n', std(gfd), 100*std(gfd)/std(gd));
fprintf('  std(gap)    = %6.2f\n', std(gap));
if std(gfd) < 0.15*std(gd)
  fprintf('  *** gdp_f barely moves. A natural-output series that is nearly flat means\n');
  fprintf('  *** the gap is just detrended output, and r*/gap lose their content.\n');
else
  fprintf('  OK  gdp_f has real variation, so the gap is not merely detrended output.\n');
end

% ---- correlation at three frequency bands ----------------------------------
cl = corr2v(gd, gfd);
dg = diff(gd); dgf = diff(gfd);
cg = corr2v(dg, dgf);
[cyc_g,  ~] = hp(gd, 1600);
[cyc_gf, ~] = hp(gfd, 1600);
cc = corr2v(cyc_g, cyc_gf);

fprintf('\nCORRELATION gdp vs gdp_f\n');
fprintf('  levels                     %+.3f   understated by the markup wandering\n', cl);
fprintf('  first differences          %+.3f   over-differenced, dominated by noise\n', cg);
fprintf('  HP cycle (lambda = 1600)   %+.3f   <-- business-cycle frequencies\n', cc);

fprintf('\n');
if cc > 0.6
  fprintf('  VERDICT: the flex block TRACKS at business-cycle frequencies (%.2f).\n', cc);
  fprintf('  The low level and growth correlations are artefacts of the wandering\n');
  fprintf('  markup component and of over-differencing, not of a broken block.\n');
elseif cc > 0.3
  fprintf('  VERDICT: MODERATE co-movement at business-cycle frequencies (%.2f).\n', cc);
  fprintf('  Defensible for a model where demand and markup shocks drive much of the\n');
  fprintf('  cycle and natural output responds mainly to supply, but worth checking\n');
  fprintf('  against the FEVD: gdp_f should be supply-driven, gdp should not be.\n');
else
  fprintf('  VERDICT: LITTLE co-movement even at business-cycle frequencies (%.2f).\n', cc);
  fprintf('  That is the case worth investigating. Check the FEVD for gdp_f: if it is\n');
  fprintf('  not dominated by the supply shocks (e_gz, e_gn, e_zeta), the flex\n');
  fprintf('  allocation is not doing what the natural-rate concept intends.\n');
end

% ---- where do the two diverge most? ----------------------------------------
fprintf('\nLARGEST GAPS (centred on the sample mean, %+.2f pp)\n', mean(gap));
gc = gap - mean(gap);
[~,ix] = sort(abs(gc),'descend');
for j = 1:min(6,numel(ix))
  i = ix(j);
  fprintf('  %-7s  gap %+6.2f   gdp %+6.2f   gdp_f %+6.2f\n', ...
          qlab(yr(i)), gc(i), gd(i), gfd(i));
end
fprintf('=================================================\n');
end

% ---- helpers ---------------------------------------------------------------
function s=qlab(d); y=floor(d+1e-9); q=round((d-y)*4)+1; s=sprintf('%dQ%d',y,q); end
function c=corr2v(a,b); a=a(:); b=b(:); n=min(numel(a),numel(b)); m=corrcoef(a(1:n),b(1:n)); c=m(1,2); end
function [cyc,trend]=hp(y,lam)
% Hodrick-Prescott, written out so this needs no toolbox.
  y=y(:); T=numel(y);
  I=speye(T);
  D=spdiags(ones(T-2,1)*[1 -2 1],0:2,T-2,T);
  trend=(I+lam*(D'*D))\y;
  cyc=y-trend;
end
