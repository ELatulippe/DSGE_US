% ============================================================================
%  make_irf_panel.m  --  impulse responses to the 12 structural shocks, laid out
%  as a grid with SHOCKS in columns and VARIABLES in rows (paper Appendix A, but
%  re-arranged). Two pages of 6 shocks each; 4 rows:
%     Output growth | Core inflation YoY | Headline inflation YoY | Federal funds rate
%  Levels (steady state + response) in pp, 40 quarters, dotted SS line -- exactly
%  the paper's units. Responses are to a +1 standard-deviation structural shock.
%
%  USAGE (after `dynare fevd_run`, which now runs irf=40):
%     >> make_irf_panel
%     >> make_irf_panel('fevd_run/Output/fevd_run_results.mat')
%  Outputs: figures_paper/fig_irf_page1.png, fig_irf_page2.png
% ============================================================================
function make_irf_panel(RESULTS, OUTDIR)
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
assert(isfield(oo_,'irfs') && ~isempty(fieldnames(oo_.irfs)), 'no oo_.irfs -- run stoch_simul with irf>0 (fevd_run).');
irf=oo_.irfs;

% ---- steady-state levels for the 4 display variables -----------------------
endo=cellstr(M_.endo_names); ss=@(nm) oo_.steady_state(find(strcmp(nm,endo),1));
PiH_ss = ss('Pi_h'); head_scale = 400/PiH_ss;         % Pi_h (gross) -> annualized pp
SSv = struct('obs_gdp',ss('obs_gdp'),'obs_pi',ss('obs_pi'),'obs_ffr',ss('obs_ffr'), ...
             'HEAD', 400*log(PiH_ss));

% ---- the 12 structural shocks (paper names) + 2 pages ----------------------
SH = { 'e_gz','Productivity'; 'e_gn','Labor disutility'; 'e_zeta','Marg. eff. of investment'; ...
       'e_mup','Price markup'; 'e_muw','Wage markup'; 'e_oS','Oil supply'; ...
       'e_chi','Marginal utility'; 'e_vth','Convenience yield'; 'e_g','Govt. consumption'; ...
       'e_tr','Transfers'; 'e_mp','Monetary policy'; 'e_pistar','Inflation target' };
% rows: {irf field or 'HEAD', label, kind}
RV = { 'obs_gdp','Output growth','lev'; 'obs_pi','Core inflation YoY','yoy'; ...
       'HEAD','Headline inflation YoY','head'; 'obs_ffr','Federal funds rate','lev' };

for pg = 1:2
  cols = (pg-1)*6 + (1:6);
  figure('visible','off','position',[60 60 1560 940],'Color','w');
  for r = 1:4
    vfield = RV{r,1}; vlab = RV{r,2}; kind = RV{r,3};
    for c = 1:6
      sidx = cols(c); sh = SH{sidx,1}; shlab = SH{sidx,2};
      ax = subplot(4,6,(r-1)*6 + c); hold on;
      y = irf_series(irf, vfield, sh, kind, head_scale);
      if isempty(y); title('n/a'); box on; continue; end
      switch kind
        case 'yoy';  dev = mov4pad(y);            base = SSv.obs_pi;
        case 'head'; dev = mov4pad(y);            base = SSv.HEAD;
        otherwise;   dev = y;                     base = SSv.(vfield);
      end
      H=numel(dev); q=1:H; lvl = base + dev(:)';
      plot(q, base+0*q, 'k:','linewidth',0.7);
      plot(q, lvl, '-','color',[0 0.35 0.75],'linewidth',1.6);
      xlim([1 H]); box on;
      if r==1; title(shlab,'interpreter','none','fontsize',9,'fontweight','bold'); end
      if c==1; ylabel(vlab,'fontsize',8); end
      set(ax,'FontSize',7);
    end
  end
  try; sgtitle(sprintf('Impulse responses to structural shocks (+1 s.d.)  --  page %d of 2  (levels, pp; dotted = steady state)', pg)); catch; end
  set(findall(gcf,'type','axes'),'XColor',[.25 .25 .25],'YColor',[.25 .25 .25],'Color','w');
  set(gcf,'InvertHardcopy','off');
  fn=fullfile(OUTDIR,sprintf('fig_irf_page%d.png',pg));
  print(fn,'-dpng','-r150'); close;
  fprintf('  wrote %s\n', fn);
end
fprintf('Done. IRF panels: %s/fig_irf_page1.png, fig_irf_page2.png\n', OUTDIR);
end

% ---- helpers ---------------------------------------------------------------
function y = irf_series(irf, vfield, sh, kind, head_scale)
y=[];
if strcmp(kind,'head')
  f=['Pi_h_' sh]; if isfield(irf,f); y = head_scale*irf.(f)(:); end
else
  f=[vfield '_' sh]; if isfield(irf,f); y = irf.(f)(:); end
end
end

% 4q trailing mean assuming zero response before the shock (proper YoY build-up)
function z = mov4pad(x)
x=x(:); n=numel(x); z=zeros(n,1);
for t=1:n
  acc=0; for k=0:3; if t-k>=1; acc=acc+x(t-k); end; end
  z(t)=acc/4;
end
end
