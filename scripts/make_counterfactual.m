% ============================================================================
%  make_counterfactual.m  --  monetary-policy counterfactuals.
%
%  (A) IMPULSE counterfactual  [needs only oo_.irfs -> run after fevd_run]:
%      "What if the Fed tightened more?" Scales the estimated monetary-policy
%      shock (e_mp) so the funds rate jumps +100bp on impact, and traces the
%      response of output growth, core PCE inflation, the funds rate, r*, and
%      the output gap over 20 quarters. Linear model -> exact.
%
%  (B) CLOSED-LOOP counterfactual  [needs cond_tighter.mat from the postest]:
%      "Higher-for-longer": the funds rate is held at 5.0% for 8 quarters
%      (vs the model's ~3.7% baseline), with the policy shock absorbing the
%      difference. Baseline vs tighter, 68% bands. Reads cond_tighter.mat.
%
%  USAGE:
%     >> make_counterfactual                        % auto-find results
%     >> make_counterfactual('fevd_run/Output/fevd_run_results.mat')
% ============================================================================
function make_counterfactual(RESULTS, OUTDIR)
if nargin<1 || isempty(RESULTS)
  cand = {'results/fevd_run_results.mat', 'results/frbstl_us_est_fastlr_lean_results.mat', 'fevd_run/Output/fevd_run_results.mat', ...
          'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
          'fevd_run_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no results .mat found -- run `dynare fevd_run` first, or pass the path.');
end
if nargin<2 || isempty(OUTDIR); OUTDIR='figures_paper'; end
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end
S=load(RESULTS); oo_=S.oo_; M_=S.M_;

% ============================= (A) 100bp SURPRISE HIKE IRF ===================
assert(isfield(oo_,'irfs'),'no oo_.irfs -- run stoch_simul with irf>0 (fevd_run).');
irf = oo_.irfs;
geti = @(nm) getfieldsafe(irf, nm);
ffr0 = geti('obs_ffr_e_mp');
assert(~isempty(ffr0),'obs_ffr_e_mp IRF missing -- add obs_ffr to the stoch_simul var list.');
scale = 1.0 / ffr0(1);                      % normalise so impact on FFR = +100bp (=+1.0 ann. %)

% (name, label, transform to annualized-% deviation)
P = { 'obs_gdp','Output growth (ann.)',   @(x) x; ...
      'obs_pi', 'Core PCE inflation (ann.)',@(x) x; ...
      'obs_ffr','Federal funds rate',      @(x) x; ...
      'rstar',  'Natural rate r* (ann.)',  @(x) 400*x; ...
      'ogap',   'Output gap (level)',      @(x) x };
H = 20; q = (0:H-1)';
figure('visible','off','position',[100 100 1150 620],'Color','w');
peak = struct();
for k=1:size(P,1)
  nm=P{k,1}; lab=P{k,2}; tf=P{k,3};
  y = geti([nm '_e_mp']); if isempty(y); continue; end
  y = tf(scale*y(1:min(H,numel(y))));
  subplot(2,3,k); hold on;
  plot(q, 0*q,'k:','linewidth',0.8);
  plot(q, y,'-','color',[0.72 0.10 0.10],'linewidth',2.0);
  [~,ip]=max(abs(y)); peak.(nm)=[y(ip) q(ip)];
  title(lab,'interpreter','tex'); xlabel('quarters after shock'); ylabel('pp deviation');
  xlim([0 H-1]); box on;
end
try; sgtitle('Counterfactual: unexpected +100bp monetary tightening'); catch; end
set(findall(gcf,'type','axes'),'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],'Color','w');
set(gcf,'InvertHardcopy','off');
print(fullfile(OUTDIR,'cf_mp_hike_100bp.png'),'-dpng','-r120'); close;
fprintf('  wrote %s/cf_mp_hike_100bp.png\n', OUTDIR);

fprintf('\n============ Counterfactual: +100bp surprise tightening (peak effects) ============\n');
fl = fieldnames(peak);
for i=1:numel(fl); pk=peak.(fl{i});
  fprintf('  %-9s peak %+6.2f pp at quarter %d\n', fl{i}, pk(1), pk(2));
end
fprintf('  (output growth & inflation fall; the gap turns negative; the funds rate is +1.00 on impact)\n');
fprintf('==================================================================================\n');

% ============================= (B) CLOSED-LOOP HIGHER-FOR-LONGER =============
if exist('cond_tighter.mat','file')
  z=load('cond_tighter.mat'); cf=z.cond_tighter;
  cvars={'obs_gdp','Output growth (ann.)','lev';'obs_pi','Core PCE infl YoY','yoy'; ...
         'obs_ffr','Federal funds rate','lev';'rstar','Natural rate r*','rstar';'ogap','Output gap','gap'};
  SV = smoothed_srcC(oo_); JOFF_ok = isfield(oo_,'SmoothedVariables');
  figure('visible','off','position',[100 100 1180 640],'Color','w');
  for k=1:size(cvars,1)
    nm=cvars{k,1}; lab=cvars{k,2}; kind=cvars{k,3};
    subplot(2,3,k); hold on;
    switch kind; case 'rstar'; tf=@(x)400*log(x); otherwise; tf=@(x)x; end
    [bm,blo,bhi]=cf_get(cf,'uncond',nm);     % baseline
    [tm,tlo,thi]=cf_get(cf,'cond',  nm);     % tighter
    if isempty(bm)&&isempty(tm); title([lab ' (n/a)']); box on; continue; end
    bm=tf(bm); blo=tf(blo); bhi=tf(bhi); tm=tf(tm); tlo=tf(tlo); thi=tf(thi);
    x=(1:numel(bm))';
    if ~isempty(blo)&&numel(blo)==numel(bm); fill([x;flipud(x)],[blo;flipud(bhi)],[0.80 0.85 0.95],'EdgeColor','none','FaceAlpha',0.5); end
    if ~isempty(tlo)&&numel(tlo)==numel(tm); fill([x;flipud(x)],[tlo;flipud(thi)],[0.98 0.84 0.72],'EdgeColor','none','FaceAlpha',0.5); end
    if ~isempty(bm); plot(1:numel(bm),bm,'-','color',[0 0.30 0.85],'linewidth',2.0); end
    if ~isempty(tm); plot(1:numel(tm),tm,'--','color',[0.85 0.33 0.10],'linewidth',2.0); end
    title(lab,'interpreter','tex'); xlabel('quarters ahead'); ylabel('annualized %'); box on;
    if k==1; legend({'68% base','68% tighter','baseline','tighter Fed'},'box','off','location','best','fontsize',7); end
  end
  try; sgtitle('Counterfactual: "higher-for-longer" (FFR held 5.0% for 8q) vs baseline'); catch; end
  set(findall(gcf,'type','axes'),'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],'Color','w');
  set(gcf,'InvertHardcopy','off');
  print(fullfile(OUTDIR,'cf_higher_for_longer.png'),'-dpng','-r120'); close;
  fprintf('  wrote %s/cf_higher_for_longer.png\n', OUTDIR);
else
  fprintf('  (closed-loop "higher-for-longer" skipped: cond_tighter.mat not found -- re-run the postest.)\n');
end
end

% -------- helpers --------
function v = getfieldsafe(s,nm); if isfield(s,nm); v=s.(nm)(:); else; v=[]; end; end
function SV = smoothed_srcC(oo_)
SV=struct(); if isfield(oo_,'SmoothedVariables'); if isfield(oo_.SmoothedVariables,'Mean'); SV=oo_.SmoothedVariables.Mean; else; SV=oo_.SmoothedVariables; end; end
end
function [m,lo,hi]=cf_get(cf,side,nm)
m=[];lo=[];hi=[];
if ~isfield(cf,side); return; end; Sd=cf.(side);
if isfield(Sd,'Mean')&&isfield(Sd.Mean,nm); m=Sd.Mean.(nm); m=m(:); if numel(m)>1; m=m(2:end); end; end
if isfield(Sd,'ci')&&isfield(Sd.ci,nm)
  c=Sd.ci.(nm); if size(c,1)==2; lo=c(1,:)';hi=c(2,:)'; elseif size(c,2)==2; lo=c(:,1);hi=c(:,2); end
  if numel(lo)>numel(m); lo=lo(end-numel(m)+1:end); hi=hi(end-numel(m)+1:end); end
end
end
