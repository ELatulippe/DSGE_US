% ============================================================================
%  diag_gap.m -- diagnostic: sticky GDP vs flexible-price GDP (gdp_f) + gap.
%  Run AFTER the postest (needs gdp and gdp_f in the smoother var list).
%     >> diag_gap
%     >> diag_gap('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat')
%  Tells you whether the flex block (gdp_f) is sensible:
%   * gdp_f should track gdp closely OUTSIDE 2020 (small gap in normal times),
%   * and diverge the RIGHT way in 2020 (actual falls more than flex -> gap<0),
%   * the computed gap must equal the model's own ogap (sanity of the definition).
% ============================================================================
function diag_gap(RESULTS)
if nargin<1 || isempty(RESULTS)
    cand = {'results/frbstl_us_est_fastlr_lean_results.mat', 'results/fevd_run_results.mat', 'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
            'frbstl_us_est_fastlr_lean_results.mat'};
    RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
    assert(~isempty(RESULTS),'no results .mat found -- pass the path explicitly.');
end
S=load(RESULTS); oo_=S.oo_; M_=S.M_;

% robust smoothed source (mode = flat, MH/load_mh_file = nested under .Mean)
SV=struct();
if isfield(oo_,'SmoothedVariables')
  if isfield(oo_.SmoothedVariables,'Mean'); SV=oo_.SmoothedVariables.Mean; else; SV=oo_.SmoothedVariables; end
end
assert(isfield(SV,'gdp') && isfield(SV,'gdp_f'), ...
  'need gdp AND gdp_f in the estimation() smoother var list (re-run the postest).');

g=SV.gdp(:); gf=SV.gdp_f(:); T=numel(g); yr=1959+(0:T-1)/4;
endo=cellstr(M_.endo_names);
if isfield(oo_,'steady_state') && ~isempty(oo_.steady_state)
  gss=oo_.steady_state(find(strcmp('gdp',endo),1)); gfss=oo_.steady_state(find(strcmp('gdp_f',endo),1));
else
  gss=exp(mean(log(g))); gfss=exp(mean(log(gf)));   % proxy if SS not stored
end
gd  = 100*log(g/gss);           % sticky output, % dev from SS
gfd = 100*log(gf/gfss);         % flex   output, % dev from SS
gap = 100*(log(g)-log(gf));     % output gap (SS cancels since gss=gfss)
ogap=[]; if isfield(SV,'ogap'); ogap=SV.ogap(:); end

% ---- printed diagnostics ---------------------------------------------------
cc=corrcoef(g,gf);
fprintf('\n================= flex-block diagnostic =================\n');
fprintf('SS: gdp = %.5f,  gdp_f = %.5f   (should be equal)\n', gss, gfss);
fprintf('corr(gdp, gdp_f) full sample     = %.3f   (want high, ~>0.9)\n', cc(1,2));
fprintf('output gap: mean %.2f%%, std %.2f%%\n', mean(gap), std(gap));
pre = yr<2020;
fprintf('pre-2020 gap: std %.2f%%, max|gap| %.2f%%   (want modest, single digits)\n', std(gap(pre)), max(abs(gap(pre))));
[mn,imn]=min(gap); q=round(mod(yr(imn),1)*4)+1;
fprintf('most negative gap: %.1f%% in %dQ%d   (paper Fig 5a: ~ -8 to -10 in 2020Q2)\n', mn, floor(yr(imn)), q);
[mx,imx]=max(gap); q2=round(mod(yr(imx),1)*4)+1;
fprintf('most positive gap: %+.1f%% in %dQ%d\n', mx, floor(yr(imx)), q2);
if ~isempty(ogap)
  fprintf('max|computed gap - model ogap|   = %.2e   (should be ~0 -> definition OK)\n', max(abs(gap-ogap)));
end
% quick 2020 read
i20=find(abs(yr-2020.25)<1e-6,1);
if ~isempty(i20)
  fprintf('2020Q2:  gdp %+.1f%%,  gdp_f %+.1f%%,  gap %+.1f%%\n', gd(i20), gfd(i20), gap(i20));
end
% ---- centered gap (remove constant efficient-benchmark wedge) --------------
off = mean(gap);
gapc = gap - off;
fprintf('---- CENTERED gap (offset %+.2f pp removed; this is what make_paper_figures reports) ----\n', off);
[mnc,imnc]=min(gapc); qc=round(mod(yr(imnc),1)*4)+1;
fprintf('centered most-negative gap: %.1f%% in %dQ%d   (paper Fig 5a: ~ -8 to -10 in 2020Q2)\n', mnc, floor(yr(imnc)), qc);
if ~isempty(i20); fprintf('centered 2020Q2 gap: %+.1f%%\n', gapc(i20)); end
fprintf('centered gap std: %.2f%% (level wedge is what differed from the paper, not the dynamics)\n', std(gapc));
fprintf('=========================================================\n');
fprintf('READS:\n');
fprintf(' * corr high & pre-2020 gap small  -> flex block tracks; gap definition fine.\n');
fprintf(' * gap wrong SIGN in 2020 (positive) -> flex output falling MORE than actual: flex-block issue.\n');
fprintf(' * corr low / gdp_f wild in normal times -> genuine bug in the flex allocation.\n');
fprintf(' * computed gap != model ogap -> the ogap equation, not the block.\n\n');

% ---- plots -----------------------------------------------------------------
figure('visible','off','position',[100 100 1150 760],'Color','w');
subplot(2,2,1); hold on;
plot(yr, gd, '-','color',[0 0.30 0.85],'linewidth',1.3);
plot(yr, gfd,'-','color',[0.85 0.33 0.10],'linewidth',1.3);
plot(yr, 0*yr,'k:'); legend({'GDP (sticky)','GDP flex (gdp\_f)'},'box','off','location','best');
title('Detrended output: sticky vs flexible-price'); xlabel('year'); ylabel('% dev. from SS');
xlim([min(yr) max(yr)]); box on;
subplot(2,2,2); hold on; w=yr>=2015;
plot(yr(w), gd(w), '-','color',[0 0.30 0.85],'linewidth',1.8);
plot(yr(w), gfd(w),'-','color',[0.85 0.33 0.10],'linewidth',1.8);
plot(yr(w), 0*yr(w),'k:'); legend({'sticky','flex'},'box','off','location','best');
title('Zoom 2015-2026'); xlabel('year'); ylabel('% dev. from SS'); xlim([2015 max(yr)]); box on;
subplot(2,2,3); hold on;
plot(yr, gap,'-','color',[0.20 0.55 0.25],'linewidth',1.3); plot(yr,0*yr,'k:');
title('Output gap = 100(log gdp - log gdp\_f)'); xlabel('year'); ylabel('%');
xlim([min(yr) max(yr)]); box on;
subplot(2,2,4); hold on;
plot(yr(w), gap(w),'-','color',[0.20 0.55 0.25],'linewidth',1.8); plot(yr(w),0*yr(w),'k:');
title('Output gap, zoom 2015-2026'); xlabel('year'); ylabel('%'); xlim([2015 max(yr)]); box on;
set(findall(gcf,'type','axes'),'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],'Color','w');
set(gcf,'InvertHardcopy','off');
if ~exist('figures','dir'); mkdir('figures'); end
print('figures/diag_gdp_vs_gdpf.png','-dpng','-r120'); close;
fprintf('wrote figures/diag_gdp_vs_gdpf.png\n');
end
