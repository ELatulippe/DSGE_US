% ============================================================================
%  make_fevd.m  --  Forecast Error Variance Decomposition (FEVD) at the posterior
%  point estimate, aggregated into the paper's 5 shock groups.
%
%  Computes the conditional FEVD directly from the model's order-1 decision rule
%  (oo_.dr) -> version-proof, and independent of Dynare's own array layout. For
%  each variable and horizon it reports the share of forecast-error variance due
%  to Supply / Demand / Fiscal / Monetary / Oil / Other(measurement) shocks.
%
%  USAGE (after `dynare fevd_run`, or after the postest which now also computes it):
%     >> make_fevd
%     >> make_fevd('fevd_run/Output/fevd_run_results.mat')
%
%  Outputs: a printed table + figures_paper/fevd_<var>.png stacked-bar charts.
% ============================================================================
function make_fevd(RESULTS, OUTDIR)
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
assert(isfield(oo_,'dr') && isfield(oo_.dr,'ghx'), 'no decision rule (oo_.dr.ghx) -- run stoch_simul (fevd_run).');

% ---- order-1 state space in DR order ---------------------------------------
% (Dynare 7.x keeps the DR block counts in M_, older versions in oo_.dr -- read robustly)
dr = oo_.dr;
ghx = dr.ghx; ghu = dr.ghu;              % rows: all endo in DR order; cols: states / shocks
getcnt = @(f) pickfield(M_, dr, f);
nstatic = getcnt('nstatic');
nsx = size(ghx, 2);                       % # state variables = # columns of ghx (robust)
sidx = nstatic + (1:nsx);                 % rows (DR order) that ARE the states
assert(max(sidx) <= size(ghx,1), 'state-row indexing off (nstatic=%d, nsx=%d, rows=%d)', nstatic, nsx, size(ghx,1));
T = ghx(sidx, :);                         % state transition (nsx x nsx)
R = ghu(sidx, :);                         % state loading on shocks (nsx x nexo)
nexo = M_.exo_nbr;
sig2 = diag(M_.Sigma_e);                  % shock variances (=1 here; sig_x are inside ghu)
name_dr = cellstr(M_.endo_names); name_dr = name_dr(dr.order_var);   % DR-order names
exo = cellstr(M_.exo_names);

% ---- horizons + variables to report ----------------------------------------
H   = [1 4 8 16 40 1000];                 % 1000 ~ unconditional
Hlab= {'1q','1y','2y','4y','10y','\infty'};
vars = {'gdp','Output';'gdp_f','Natural output (flex)';'obs_pi','Core PCE infl (ann.)';'obs_ffr','Fed funds rate'; ...
        'rstar','Natural rate r*';'ogap','Output gap';'iv','Investment'};

% ---- shock groups (COVID specials folded in, exactly as the decompositions) -
groups = struct('Supply',{{'e_gz','e_gn','e_zeta','e_mup','e_muw','e_gn_c'}}, ...
                'Demand',{{'e_chi','e_vth','e_chi_c'}}, ...
                'Fiscal',{{'e_g','e_tr','e_tr_c'}}, ...
                'Monetary',{{'e_mp','e_pistar'}}, ...
                'Oil',{{'e_oS'}});
gnames=fieldnames(groups);
gmap = zeros(nexo,1);                      % shock -> group index (0 = Other)
for g=1:numel(gnames)
  for s=groups.(gnames{g}); j=find(strcmp(s{1},exo),1); if ~isempty(j); gmap(j)=g; end; end
end
gcol=[0.00 0.45 0.74;0.85 0.33 0.10;0.49 0.18 0.56;0.93 0.69 0.13;0.47 0.67 0.19;0.30 0.75 0.93];
glabels=[gnames;{'Other'}];

% ---- FEVD per variable -----------------------------------------------------
maxH=max(H(H<1e6)); bigH=max(H);
fprintf('\n===================== FEVD (%% of forecast-error variance) =====================\n');
for vi=1:size(vars,1)
  vn=vars{vi,1}; jj=find(strcmp(vn,name_dr),1);
  if isempty(jj); warning('FEVD: %s not in decision rule (add it to the stoch_simul list).',vn); continue; end
  % MA coefficients of variable jj wrt each shock, accumulate variance by horizon
  Vg = zeros(numel(H), numel(gnames)+1);   % horizon x (groups+Other)
  % precompute cumulative variance contribution per shock at each requested horizon
  contrib = zeros(nexo,1); Sprev = R; h=0; Hset=sort(H); ptr=1; store=zeros(nexo,numel(H));
  % h=0 term
  MA = ghu(jj,:)';                          % nexo x1
  contrib = contrib + (MA.^2).*sig2;
  hcur=1;                                    % we've now covered horizon=1 (h=0 only)
  while ptr<=numel(Hset) && Hset(ptr)==1; store(:,ptr)=contrib; ptr=ptr+1; end
  for h=1:bigH-1
    MA = (ghx(jj,:)*Sprev)';                 % nexo x1
    contrib = contrib + (MA.^2).*sig2;
    Sprev = T*Sprev;
    hcur=h+1;
    while ptr<=numel(Hset) && Hset(ptr)==hcur; store(:,ptr)=contrib; ptr=ptr+1; end
  end
  % map Hset order back to H order
  [~,ord]=ismember(H,Hset);
  storeH=store(:,ord);
  for c=1:numel(H)
    tot=sum(storeH(:,c)); if tot<=0; tot=1; end
    for g=1:numel(gnames); Vg(c,g)=100*sum(storeH(gmap==g,c))/tot; end
    Vg(c,end)=100*sum(storeH(gmap==0,c))/tot;
  end
  % print
  fprintf('\n%-22s  |', vars{vi,2});
  fprintf(' %8s', glabels{:}); fprintf('\n'); fprintf('%s\n',repmat('-',1,22+2+9*numel(glabels)));
  for c=1:numel(H)
    fprintf('  horizon %-11s |', Hlab{c}); fprintf(' %7.1f ', Vg(c,:)); fprintf('\n');
  end
  % stacked-bar figure
  figure('visible','off','position',[100 100 720 420],'Color','w');
  b=bar(1:numel(H), Vg, 0.75, 'stacked','EdgeColor','none'); hold on;
  for k=1:numel(b); set(b(k),'FaceColor',gcol(k,:)); end
  set(gca,'XTick',1:numel(H),'XTickLabel',Hlab); xlim([0.4 numel(H)+0.6]); ylim([0 100]);
  legend(glabels,'location','eastoutside','box','off');
  title(sprintf('FEVD: %s', vars{vi,2}),'interpreter','tex');
  xlabel('forecast horizon'); ylabel('% of variance');
  set(findall(gcf,'type','axes'),'XColor',[.2 .2 .2],'YColor',[.2 .2 .2],'Color','w');
  set(gcf,'InvertHardcopy','off');
  fn=fullfile(OUTDIR,sprintf('fevd_%s.png',vn));
  print(fn,'-dpng','-r120'); close;
  fprintf('  wrote %s\n', fn);
end
fprintf('===============================================================================\n');
fprintf('Reads: which shock GROUP drives each variable''s forecast uncertainty, by horizon.\n');
fprintf('(point FEVD at the posterior estimate; Other = measurement error + unmatched shocks)\n');
end

% read a DR block-count from M_ (Dynare 7.x) or oo_.dr (older); error if neither
function v = pickfield(M_, dr, f)
if isfield(M_,f);        v = M_.(f);
elseif isfield(dr,f);    v = dr.(f);
else; error('field "%s" not found in M_ or oo_.dr -- cannot locate state block.', f);
end
end
