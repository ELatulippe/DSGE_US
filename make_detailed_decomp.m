% ============================================================================
%  make_detailed_decomp.m  --  the paper's Appendix-B "Detailed Historical
%  Decompositions" (Figs 9-10): the SAME four objects as the main text, but
%  decomposed into the 12 INDIVIDUAL structural shocks (not the 5 groups):
%     Fig 9a  YoY Core PCE Inflation      Fig 10a  Output Gap
%     Fig 9b  YoY Headline PCE Inflation  Fig 10b  Spot natural rate r*
%  Window 2019Q4-2026Q1, relative to steady state. "Other" = initial conditions
%  + measurement error. COVID special shocks fold into their parents
%  (chi_c->Marginal utility, tr_c->Transfers, gn_c->Labor disutility).
%
%  Reads the POST-ESTIMATION results (oo_.shock_decomposition at posterior mean).
%  USAGE:
%     >> make_detailed_decomp
%     >> make_detailed_decomp('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat')
%  Outputs: figures_paper/fig9a_core_detailed.png, fig9b_headline_detailed.png,
%           fig10a_gap_detailed.png, fig10b_rstar_detailed.png
% ============================================================================
function make_detailed_decomp(RESULTS, OUTDIR)
if nargin<1 || isempty(RESULTS)
  cand = {'results/frbstl_us_est_fastlr_lean_results.mat', 'results/fevd_run_results.mat', 'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
          'frbstl_us_est_fastlr_lean_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no results .mat found -- run the postest first, or pass the path.');
end
if nargin<2 || isempty(OUTDIR); OUTDIR='figures_paper'; end
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end
S=load(RESULTS); oo_=S.oo_; M_=S.M_;
assert(isfield(oo_,'shock_decomposition'),'no oo_.shock_decomposition -- run the postest shock_decomposition block.');
fprintf('Loaded %s\n', RESULTS);

exo=cellstr(M_.exo_names); endo=cellstr(M_.endo_names);
idx=@(nm,lst) find(strcmp(strtrim(nm),strtrim(lst)),1);
Pibar = 1.02^0.25;
T = size(oo_.shock_decomposition,3); yr = 1959+(0:T-1)/4;

% ---- gap centering (same convention as make_paper_figures: long-run gap = 0) --
gapoff = 0; SV=struct();
if isfield(oo_,'SmoothedVariables'); if isfield(oo_.SmoothedVariables,'Mean'); SV=oo_.SmoothedVariables.Mean; else; SV=oo_.SmoothedVariables; end; end
FC=struct();
if isfield(oo_,'MeanForecast')&&isfield(oo_.MeanForecast,'Mean'); FC=oo_.MeanForecast.Mean;
elseif isfield(oo_,'forecast')&&isfield(oo_.forecast,'Mean');      FC=oo_.forecast.Mean; end
if isfield(FC,'ogap'); ogf=FC.ogap(:); nT=min(8,numel(ogf)); gapoff=mean(ogf(end-nT+1:end));
elseif isfield(SV,'ogap'); gapoff=mean(SV.ogap(:)); end

% ---- the 12 detailed shock groups (COVID specials folded into parents) + Other
G = { 'TFP',{{'e_gz'}}; 'P markup',{{'e_mup'}}; 'W markup',{{'e_muw'}}; 'MEI',{{'e_zeta'}}; ...
      'Labor disut.',{{'e_gn','e_gn_c'}}; 'Oil supply',{{'e_oS'}}; 'Mg. util.',{{'e_chi','e_chi_c'}}; ...
      'Convenience',{{'e_vth'}}; 'Govt. cons.',{{'e_g'}}; 'Transfers',{{'e_tr','e_tr_c'}}; ...
      'Mon. pol.',{{'e_mp'}}; 'Infl. target',{{'e_pistar'}} };
gname = G(:,1);
col = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.49 0.18 0.56; 0.47 0.67 0.19; ...
       0.93 0.10 0.55; 0.64 0.08 0.18; 0.93 0.69 0.13; 0.30 0.75 0.93; ...
       0.55 0.35 0.15; 0.10 0.60 0.55; 0.80 0.45 0.20; 0.00 0.20 0.55];   % 12 group colors
othercol = [0.6 0.6 0.6];

spec = { 'obs_pi',   'Core PCE Inflation YoY',           true,  1.0,       'fig9a_core_detailed'; ...
         'Pi_h',     'Headline PCE Inflation YoY',       true,  400/Pibar, 'fig9b_headline_detailed'; ...
         'ogap',     'Output Gap',                       false, 1.0,       'fig10a_gap_detailed'; ...
         'rstar_ann','Spot natural rate of interest (r*)',false, 1.0,      'fig10b_rstar_detailed' };

w0=2019.75; w1=2026.0; wm=(yr>=w0-1e-6)&(yr<=w1+1e-6);
for d=1:size(spec,1)
  vname=spec{d,1}; vttl=spec{d,2}; doYoY=spec{d,3}; vscale=spec{d,4}; fn=spec{d,5};
  try
    ei=idx(vname,endo); assert(~isempty(ei),sprintf('%s not in shock_decomposition',vname));
    nexo=M_.exo_nbr; sd=squeeze(oo_.shock_decomposition(ei,:,:));   % (nexo+2) x T
    tot=sd(nexo+2,:);
    Gm=zeros(numel(gname),T);
    for g=1:numel(gname)
      shlist=G{g,2}{1};
      for s=1:numel(shlist); j=idx(shlist{s},exo); if ~isempty(j); Gm(g,:)=Gm(g,:)+sd(j,:); end; end
    end
    other = tot - sum(Gm,1);            % init cond + measurement error + residual
    Gall=[Gm; other];
    if strcmp(vname,'ogap') && gapoff~=0    % center the gap (attribute constant wedge to Other)
      tot=tot-gapoff; Gall(end,:)=Gall(end,:)-gapoff;
    end
    Gall=vscale*Gall; tot=vscale*tot;
    if doYoY; Gall=mov4(Gall); tot=mov4(tot); end
    Gw=Gall(:,wm); tw=tot(wm); yw=yr(wm);
    figure('visible','off','position',[100 100 1050 520],'Color','w');
    hpos=max(Gw,0); hneg=min(Gw,0); cols=[col; othercol];
    b1=bar(yw,hpos','stacked','EdgeColor','none'); hold on;
    for k=1:numel(b1); set(b1(k),'FaceColor',cols(k,:)); end
    b2=bar(yw,hneg','stacked','EdgeColor','none');
    for k=1:numel(b2); set(b2(k),'FaceColor',cols(k,:)); end
    hline=plot(yw,tw,'k-','linewidth',2.0);
    xlim([min(yw)-0.2 max(yw)+0.2]); set(gca,'XTick',ceil(min(yw)):1:floor(max(yw)));
    legend([b1(:)', hline], [gname; {'Other';'Data'}], 'location','eastoutside','box','off','fontsize',8);
    ttlnote=''; if strcmp(vname,'ogap')&&gapoff~=0; ttlnote=sprintf(', centered (-%.1fpp)',gapoff); end
    title(sprintf('%s   (2019Q4-2026Q1, dev. from steady state%s)', vttl, ttlnote));
    xlabel('year'); ylabel('percent');
    set(findall(gcf,'type','axes'),'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],'Color','w');
    set(gcf,'InvertHardcopy','off');
    print(fullfile(OUTDIR,[fn '.png']),'-dpng','-r120'); close;
    fprintf('  wrote %s.png\n', fn);
  catch err
    warning('detailed decomp %s failed: %s', vname, err.message);
  end
end
fprintf('Done. Detailed decompositions in ./%s/\n', OUTDIR);
end

function Y=mov4(X)   % 4q trailing mean along columns
[m,T]=size(X); Y=zeros(m,T);
for t=1:T; k=max(1,t-3):t; Y(:,t)=mean(X(:,k),2); end
end
