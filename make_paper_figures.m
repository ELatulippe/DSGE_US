% ============================================================================
%  make_paper_figures.m  --  reproduce the paper's results objects for comparison
%  St. Louis Fed DSGE model (Faria-e-Castro, 2026).  Reads the oo_/M_ from the
%  banded post-estimation run (frbstl_us_est_fastlr_lean_results.mat) and builds,
%  in the paper's exact format:
%     Fig 4  YoY Core & Headline PCE inflation historical decompositions
%     Fig 5  Output-gap & spot natural-rate (r*) historical decompositions
%            [Figs 4-5 window: 2019Q4 -> 2026Q1, deviations from steady state]
%     Fig 6  Model r* (1y & 5y forward) vs external HLW / Lubik-Matthes (slots)
%     Fig 7  Unconditional forecast (Output growth, Core PCE YoY, FFR, r*)
%     Fig 8  Conditional scenarios (constant-rate & inflation-target) vs baseline
%     Table 5  Unconditional forecast point estimates + paper values, side by side
%
%  USAGE (from the working folder, after the postest run):
%     >> make_paper_figures
%     >> make_paper_figures('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat')
%
%  SHOCK GROUPS (paper's 5, with the COVID special shocks folded into their
%  economic category, exactly as the paper does):
%     Supply   = e_gz e_gn e_zeta e_mup e_muw   (+ e_gn_c)
%     Demand   = e_chi e_vth                    (+ e_chi_c)
%     Fiscal   = e_g e_tr                        (+ e_tr_c)
%     Monetary = e_mp e_pistar
%     Oil      = e_oS
%     Other    = everything else (measurement error + initial conditions)
%
%  NOTES / CAVEATS (see the companion message):
%   * YoY inflation = 4-quarter trailing mean of the annualized quarterly rate;
%     the decomposition is aggregated the same way (4q mean of each contribution).
%   * Headline (Pi_h) is not a %-observable, so its decomposition is scaled from
%     the gross-Pi_h decomposition by 400/Pibar (exact for the linearized model).
%   * r* is shown continuously annualized (400*log gross) -- tamer than compound.
%   * Table 5 / Fig 7 need forecast=40 in the postest to reach 2029Q4 & the
%     2031-35 'longer run'; with forecast=12 you get 2026-2028 (rest -> n/a).
%   * Fig 8 conditional scenarios ship with oil-conditioning OFF (see postest).
% ============================================================================
function make_paper_figures(RESULTS, OUTDIR)
if nargin < 1 || isempty(RESULTS)
    cand = {'results/frbstl_us_est_fastlr_lean_results.mat', 'results/fevd_run_results.mat',  'frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
             'frbstl_us_est_fastlr_lean_results.mat' };
    RESULTS = '';
    for i = 1:numel(cand); if exist(cand{i},'file'); RESULTS = cand{i}; break; end; end
    assert(~isempty(RESULTS), 'No results .mat found -- pass the path explicitly.');
end
if nargin < 2 || isempty(OUTDIR); OUTDIR = 'figures_paper'; end
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end
S = load(RESULTS); oo_ = S.oo_; M_ = S.M_;
fprintf('Loaded %s\nSaving to ./%s/\n', RESULTS, OUTDIR);

% ---- output-gap centering --------------------------------------------------
% The flex/natural block holds markups at steady state (efficient-output concept),
% so the model gap carries a ~constant level wedge vs the actual economy (the gap
% averages ~+5pp). r* is a rate (a growth/difference object) and does NOT inherit
% this wedge -- its level already matches the paper. Removing the constant offset
% lands the 2020Q2 trough near the paper's -8..-10 and the forecast gap near 0.
% Set CENTER_GAP=false to see the raw (uncentered) efficient-benchmark gap.
CENTER_GAP = true;
SV0 = smoothed_src(oo_);
FC0 = forecast_src(oo_);
gapoff = 0;
if CENTER_GAP
  ogf = getfc(FC0,'ogap','Mean');           % forecast path of the model's ogap
  if ~isempty(ogf)
    nT = min(8, numel(ogf));
    gapoff = mean(ogf(end-nT+1:end));       % long-run (stochastic-SS) level -> gap=0 in the long run (paper's convention)
  elseif isfield(SV0,'ogap')
    gapoff = mean(SV0.ogap(:));             % fallback (no forecast): full-sample mean
  end
end
if gapoff~=0
  fprintf('Output gap centered on long-run forecast level: removing %+.2f pp so the gap is ~0 in the long run.\n', gapoff);
end

% ---- population add-on for AGGREGATE output growth --------------------------
% Paper Table 5 / footnote 3: reported "output growth" is AGGREGATE real GDP
% growth = model per-capita growth (obs_gdp) + a Census population-growth
% forecast. We keep the estimation DATA-ANCHORED (our Gammabar_Z ~1.0039 gives a
% ~1.44% per-capita trend, consistent with realized avg GDP growth), so the ONLY
% legitimate add-on here is the true population growth (~0.4pp Census).
% NOTE: this leaves our aggregate growth ~0.4pp below the paper, which is the
% Gammabar_Z trend-growth difference (paper's 1.0049 sits above realized data) --
% a deliberate data-anchored choice, not a residual to hide in POP_GROWTH.
POP_GROWTH = 0.4;

INTERP = 'tex';
endo = cellstr(M_.endo_names);  exo = cellstr(M_.exo_names);
idx  = @(nm,lst) find(strcmp(strtrim(nm), strtrim(lst)));
Pibar = 1.02^0.25;

% ---- sample dates ----------------------------------------------------------
T    = size(oo_.shock_decomposition, 3);
yr   = 1959 + (0:T-1)/4;                 % calendar year of each quarter
JOFF = 1959 + (T-1)/4;                   % last data quarter (= 2026Q1)

% ---- shock groups (COVID special shocks folded into economic categories) ---
groups = struct( ...
  'Supply',   {{'e_gz','e_gn','e_zeta','e_mup','e_muw','e_gn_c'}}, ...
  'Demand',   {{'e_chi','e_vth','e_chi_c'}}, ...
  'Fiscal',   {{'e_g','e_tr','e_tr_c'}}, ...
  'Monetary', {{'e_mp','e_pistar'}}, ...
  'Oil',      {{'e_oS'}} );
gnames = fieldnames(groups);
gcol = [0.00 0.45 0.74;   % Supply   blue
        0.85 0.33 0.10;   % Demand   orange
        0.49 0.18 0.56;   % Fiscal   purple
        0.93 0.69 0.13;   % Monetary gold
        0.47 0.67 0.19];  % Oil      green
othercol = [0.30 0.75 0.93];   % Other   light blue (paper)

% ============================ FIGS 4 & 5: HISTORICAL DECOMPOSITIONS ==========
% window 2019Q4 -> 2026Q1
w0 = 2019.75; w1 = 2026.0;
wm = (yr >= w0-1e-6) & (yr <= w1+1e-6);

% (variable, title, YoY?, scale) for the four panels
decspec = { ...
  'obs_pi', 'YoY Core PCE Inflation',      true,  1.0, 'fig4a_core_pce_yoy'; ...
  'Pi_h',   'YoY Headline PCE Inflation',  true,  400/Pibar, 'fig4b_headline_pce_yoy'; ...
  'ogap',   'Output Gap',                  false, 1.0, 'fig5a_output_gap'; ...
  'rstar_ann','Spot natural rate of interest (r*)', false, 1.0, 'fig5b_rstar_spot' };

for d = 1:size(decspec,1)
  vname = decspec{d,1}; vttl = decspec{d,2}; doYoY = decspec{d,3};
  vscale = decspec{d,4}; fn = decspec{d,5};
  try
    [G, tot] = group_decomp(oo_, M_, vname, groups, gnames, exo, idx);   % (nG+1)xT, 1xT
    if strcmp(vname,'ogap') && gapoff~=0        % center the gap: subtract constant wedge (attribute to 'Other')
      tot = tot - gapoff; G(end,:) = G(end,:) - gapoff;
    end
    G = vscale*G; tot = vscale*tot;
    if doYoY; G = mov4(G); tot = mov4(tot); end
    Gw = G(:, wm); tw = tot(wm); yw = yr(wm);
    figure('visible','off','position',[100 100 980 470],'Color','w');
    hpos = max(Gw,0); hneg = min(Gw,0);
    cols = [gcol; othercol];
    b1 = bar(yw, hpos', 'stacked', 'EdgeColor','none'); hold on;
    for k=1:numel(b1); set(b1(k),'FaceColor',cols(k,:)); end
    b2 = bar(yw, hneg', 'stacked', 'EdgeColor','none');
    for k=1:numel(b2); set(b2(k),'FaceColor',cols(k,:)); end
    hline = plot(yw, tw, 'k-', 'linewidth', 2.0);   % black line = actual (smoothed) series = sum of the bars
    xlim([min(yw)-0.2 max(yw)+0.2]); set(gca,'XTick', ceil(min(yw)):1:floor(max(yw)));
    % map labels to explicit handles so 'Data' is the BLACK line (not a blue bar series)
    legend([b1(:)', hline], [gnames; {'Other'; 'Data'}], 'location','eastoutside','box','off');
    ttlnote = '';
    if strcmp(vname,'ogap') && gapoff~=0; ttlnote = sprintf(', centered (-%.1fpp wedge)', gapoff); end
    title(sprintf('%s   (2019Q4-2026Q1, dev. from steady state%s)', vttl, ttlnote), 'interpreter',INTERP);
    xlabel('year'); ylabel('percent');
    lighttheme(gcf);
    print(fullfile(OUTDIR,[fn '.png']),'-dpng','-r120'); close;
    fprintf('  wrote %s.png\n', fn);
  catch err
    warning('decomp %s failed: %s', vname, err.message);
  end
end

% ============================ FIG 6: NATURAL RATE r* =========================
% Model spot + 1y/5y forward (continuous annualized). External comparators are
% data (FRB NY = Holston-Laubach-Williams; FRB Richmond = Lubik-Matthes): drop
% them into HLW / LM below, quarterly-aligned to yrr, to complete the figure.
try
  SVr = smoothed_src(oo_);
  comp2cont = @(pc) 100*log(1 + pc/100);          % 100*(g^4-1) -> 400*log(g)
  if isfield(SVr,'rstar'); spot = 400*log(SVr.rstar(:)); else; spot = comp2cont(SVr.rstar_ann(:)); end
  Tr = numel(spot); yrr = 1959 + (0:Tr-1)/4;
  figure('visible','off','position',[100 100 1000 460],'Color','w'); hold on;
  plot(yrr, 0*yrr, 'k:','linewidth',0.8);
  leg = {}; yref = [];
  if isfield(SVr,'rstar_1y'); f1=SVr.rstar_1y(:); plot(yrr,f1,'-','color',[0 0.30 0.85],'linewidth',1.4); leg{end+1}='DSGE, 1y forward'; yref=[yref f1(:)']; end
  if isfield(SVr,'rstar_5y'); f5=SVr.rstar_5y(:); plot(yrr,f5,'-','color',[0.85 0.10 0.10],'linewidth',2.2); leg{end+1}='DSGE, 5y forward'; yref=[yref f5(:)']; end
  % ---- external comparators: read rstar_external.csv (date,LM_median,HLW_rstar) ----
  %   LM = Lubik-Matthes (Richmond Fed) -- shipped filled 1985Q1-2026Q1.
  %   HLW = Holston-Laubach-Williams (NY Fed) -- blank until you run
  %         add_hlw_to_external.m on the official xlsx (see that file's header).
  ext = read_external_rstar(yrr);
  if ~isempty(ext.LM) && any(isfinite(ext.LM))
    plot(yrr, ext.LM, '--','color',[0.1 0.5 0.1],'linewidth',1.5);
    leg{end+1}='Lubik-Matthes'; yref=[yref ext.LM(isfinite(ext.LM))'];
  end
  if ~isempty(ext.HLW) && any(isfinite(ext.HLW))
    plot(yrr, ext.HLW, '--','color',[0.4 0.4 0.4],'linewidth',1.5);
    leg{end+1}='Holston-Laubach-Williams'; yref=[yref ext.HLW(isfinite(ext.HLW))'];
  end
  if isempty(yref); yref = spot(:)'; end
  xlim([1985 max(yrr)]); set(gca,'XTick',1985:5:2025);
  % show the full forward range (the 5y-forward plunges in 2009/2020 were being clipped);
  % floor at -8 / cap at +8 so a single COVID spike can't dominate the axis
  yv = yref(isfinite(yref)); ylo = max(min(yv),-8); yhi = min(max(yv),8);
  if yhi<=ylo; yhi=ylo+1; end; dd=yhi-ylo; ylim([ylo-0.05*dd, yhi+0.05*dd]);
  box on; legend(leg,'location','best','box','off');
  title('Natural rate of interest r^*: model forwards vs standard measures','interpreter',INTERP);
  xlabel('year'); ylabel('annualized %');
  lighttheme(gcf);
  print(fullfile(OUTDIR,'fig6_rstar.png'),'-dpng','-r120'); close;
  fprintf('  wrote fig6_rstar.png\n');
catch err
  warning('Fig 6 failed: %s', err.message);
end

% ============================ FIG 7: UNCONDITIONAL FORECAST ==================
% Panels: Output growth (obs_gdp, per-capita*), Core PCE inflation YoY, FFR, r*.
% * paper aggregate growth = per-capita + a Census population-growth add-on.
try
  SVr = smoothed_src(oo_); FC = forecast_src(oo_); PF = pointfc_src(oo_);
  if isempty(fieldnames(PF)); PF = FC; end   % fall back to MeanForecast HPD if no PointForecast
  % Output growth shown YoY (4q) so transitory shocks cancel -> readable band (paper uses Q4/Q4).
  % r* shown as the smoothed 1-year forward (rstar_1y), not the volatile spot rate.
  panels = { 'obs_gdp',  'Output growth (YoY)',      'yoy'; ...
             'obs_pi',   'Core PCE inflation (YoY)',  'yoy'; ...
             'obs_ffr',  'Federal funds rate',        'lev'; ...
             'rstar_1y', 'Natural rate r* (1y fwd)',  'lev' };
  figure('visible','off','position',[100 100 1050 700],'Color','w');
  for v = 1:4
    nm = panels{v,1}; ttl = panels{v,2}; kind = panels{v,3};
    subplot(2,2,v); hold on;
    h  = getf1(SVr, nm);                            % smoothed history
    fc = getfc(FC, nm, 'Mean');                     % central path
    ilo= getfc(FC,nm,'HPDinf'); ihi= getfc(FC,nm,'HPDsup');   % INNER band: mean-forecast (parameter) uncertainty -- tight
    olo= getfc(PF,nm,'HPDinf'); ohi= getfc(PF,nm,'HPDsup');   % OUTER band: predictive (param + future shocks) -- wide
    if isempty(fc); title([ttl ' (n/a)']); box on; continue; end
    switch kind; case 'rstar'; tf=@(x) 400*log(x); otherwise; tf=@(x) x; end
    hh=tf(h(:)); ff=tf(fc(:)); ilo=tf(ilo(:)); ihi=tf(ihi(:)); olo=tf(olo(:)); ohi=tf(ohi(:));
    if strcmp(kind,'yoy')                            % 4q trailing mean, stitched on history
      ff=yoystitch(hh,ff); ilo=yoystitch(hh,ilo); ihi=yoystitch(hh,ihi);
      olo=yoystitch(hh,olo); ohi=yoystitch(hh,ohi); hh=mov4v(hh);
    end
    if strcmp(nm,'obs_gdp')      % report AGGREGATE growth (per-capita + population add-on), as in the paper
      hh=hh+POP_GROWTH; ff=ff+POP_GROWTH; ilo=ilo+POP_GROWTH; ihi=ihi+POP_GROWTH; olo=olo+POP_GROWTH; ohi=ohi+POP_GROWTH;
    end
    Kh=20; hh=hh(max(1,end-Kh+1):end);
    % show the forecast only through ~2030 (the longer tail is flat and uninformative)
    Hshow = max(1, round((2030 - JOFF)*4));
    cutf = @(x) x(1:min(Hshow, numel(x)));
    ff=cutf(ff); ilo=cutf(ilo); ihi=cutf(ihi); olo=cutf(olo); ohi=cutf(ohi);
    xh=JOFF-(numel(hh)-1)/4:0.25:JOFF; xf=JOFF+(1:numel(ff))/4;
    % ONE band only: the tighter 68% interval (mean-forecast / parameter uncertainty).
    % (The wide predictive band is available as olo/ohi -- set SHOW_PREDICTIVE=true to overlay it.)
    SHOW_PREDICTIVE = false;
    if SHOW_PREDICTIVE && ~isempty(olo)&&numel(olo)==numel(ff)
      fill([xf fliplr(xf)],[olo' fliplr(ohi')],[0.90 0.93 0.98],'EdgeColor','none');
    end
    if ~isempty(ilo)&&numel(ilo)==numel(ff); fill([xf fliplr(xf)],[ilo' fliplr(ihi')],[0.62 0.74 0.90],'EdgeColor','none'); end
    plot(xh, hh, '-', 'color',[0.15 0.35 0.55], 'linewidth',1.6);
    plot(xf, ff, '-', 'color',[0 0.30 0.85], 'linewidth',2.2);
    scaleset = [hh(:);ff(:);ilo(:);ihi(:)];
    if SHOW_PREDICTIVE; scaleset = [scaleset; olo(:); ohi(:)]; end
    yl = niceylim(scaleset, 0.22, 2.0); ylim(yl);      % roomy range around the band + history
    plot([JOFF JOFF], yl, '--','color',[0.4 0.4 0.4]);
    xlim([2022 2030]); set(gca,'XTick',2022:2:2030);
    box on; title(ttl,'interpreter',INTERP); ylabel('annualized %');
    if v==1; legend({'68% band','history','forecast'},'box','off','location','best','fontsize',8); end
  end
  try; sgtitle('Unconditional model forecast  (dashed line = 2026Q1 jump-off)'); catch; end
  lighttheme(gcf);
  print(fullfile(OUTDIR,'fig7_uncond_forecast.png'),'-dpng','-r120'); close;
  fprintf('  wrote fig7_uncond_forecast.png\n');
catch err
  warning('Fig 7 failed: %s', err.message);
end

% ============================ FIG 8: CONDITIONAL SCENARIOS ===================
% Two rows = two scenarios; each overlays the UNCONDITIONAL baseline (blue solid)
% vs the CONDITIONAL scenario mean (orange dashed). Vars: output growth, core PCE
% YoY, FFR, r*.  (Oil-conditioning is OFF in the shipped postest -- see caveats.)
try
  FC = forecast_src(oo_); SVr = smoothed_src(oo_); PF = pointfc_src(oo_);
  if isempty(fieldnames(PF)); PF = FC; end
  sc = {};
  if exist('cond_constrate.mat','file');  z=load('cond_constrate.mat');  sc{end+1}={'Constant rate',    z.cond_constrate};  end
  if exist('cond_infltarget.mat','file'); z=load('cond_infltarget.mat'); sc{end+1}={'Inflation target', z.cond_infltarget}; end
  cvars = { 'obs_gdp','Output growth (YoY)','yoy'; 'obs_pi','Core PCE inflation (YoY)','yoy'; ...
            'obs_ffr','Federal funds rate','lev'; 'rstar_1y','Natural rate r* (1y fwd)','lev' };
  if ~isempty(sc)
    nRow = numel(sc);
    figure('visible','off','position',[100 100 1200 360*nRow],'Color','w'); p=0;
    for si = 1:nRow
      nmsc = sc{si}{1}; cf = sc{si}{2};
      for v = 1:4
        nm=cvars{v,1}; ttl=cvars{v,2}; kind=cvars{v,3}; p=p+1; subplot(nRow,4,p); hold on;
        switch kind; case 'rstar'; tf=@(x) 400*log(x); otherwise; tf=@(x) x; end
        uf = getfc(FC, nm, 'Mean');                       % unconditional mean
        ulo0 = getfc(FC,nm,'HPDinf'); uhi0 = getfc(FC,nm,'HPDsup'); % uncond band: tighter mean-forecast 68%
        cfm = cond_mean(cf, nm);                          % conditional mean (drop origin)
        [clo0,chi0] = cond_ci(cf, nm);                    % conditional 68% band
        h  = getf1(SVr, nm);
        if isempty(uf) && isempty(cfm); title([nmsc ': ' ttl ' (n/a)']); box on; continue; end
        hh = tf(h(:)); uu = tf(uf(:)); cc = tf(cfm(:));
        ulo=tf(ulo0(:)); uhi=tf(uhi0(:)); clo=tf(clo0(:)); chi=tf(chi0(:));
        % align band lengths to their mean paths (trim any origin element from the front)
        trimf = @(x,n) x(max(1,numel(x)-n+1):end);
        if ~isempty(uu); ulo=trimf(ulo,numel(uu)); uhi=trimf(uhi,numel(uu)); end
        if ~isempty(cc); clo=trimf(clo,numel(cc)); chi=trimf(chi,numel(cc)); end
        if strcmp(kind,'yoy')                              % YoY = 4q trailing mean, stitched on history
          uu=yoystitch(hh,uu); ulo=yoystitch(hh,ulo); uhi=yoystitch(hh,uhi);
          cc=yoystitch(hh,cc); clo=yoystitch(hh,clo); chi=yoystitch(hh,chi);
          hh = mov4v(hh);
        end
        % fig 8 spans the CONDITIONAL horizon (paper stops ~2031): cap the longer
        % unconditional path (forecast=40) to the conditional length so both align.
        if ~isempty(cc)
          Lc = numel(cc);
          if numel(uu)>Lc;  uu = uu(1:Lc);  end
          if numel(ulo)>Lc; ulo = ulo(1:Lc); uhi = uhi(1:Lc); end
        end
        Kh=20; hh=hh(max(1,end-Kh+1):end); xh=JOFF-(numel(hh)-1)/4:0.25:JOFF;
        % bands removed on request (too wide to be informative) -- lines only
        plot(xh,hh,'-','color',[0.4 0.4 0.4],'linewidth',1.2);
        if ~isempty(uu); xu=JOFF+(1:numel(uu))/4; plot(xu,uu,'-','color',[0 0.30 0.85],'linewidth',2.0); end
        if ~isempty(cc); xc=JOFF+(1:numel(cc))/4; plot(xc,cc,'--','color',[0.85 0.33 0.10],'linewidth',2.0); end
        % scale axis to the central paths + history
        yl = niceylim([hh(:);uu(:);cc(:)], 0.22, 2.0); ylim(yl);
        plot([JOFF JOFF],yl,':','color',[0.4 0.4 0.4]);
        xhi = JOFF + max([numel(uu) numel(cc) 1])/4;
        xlim([2022 xhi]); set(gca,'XTick',2022:3:ceil(xhi));
        box on; title(sprintf('%s: %s', nmsc, ttl),'interpreter',INTERP);
        if v==1; ylabel('annualized %'); end
        if p==1; legend({'history','Unconditional','Conditional'},'box','off','location','best','fontsize',7); end
      end
    end
    try; sgtitle('Conditional forecast scenarios  (baseline vs conditioned; oil-conditioning OFF)'); catch; end
    lighttheme(gcf);
    print(fullfile(OUTDIR,'fig8_cond_scenarios.png'),'-dpng','-r120'); close;
    fprintf('  wrote fig8_cond_scenarios.png\n');
  else
    warning('Fig 8: no cond_*.mat scenario files found in the folder.');
  end
catch err
  warning('Fig 8 failed: %s', err.message);
end

% ============================ TABLE 5: FORECAST POINT ESTIMATES ==============
% Q4/Q4 for growth & inflation; Q4 level for FFR, r*, output gap.
% Paper values (Table 5) hard-coded for a side-by-side comparison.
try
  SVr = smoothed_src(oo_); FC = forecast_src(oo_);
  % combined quarterly series: smoothed history (t<=T) then forecast (t>T)
  comb = @(nm) [colv(getf1(SVr,nm)); colv(getfc(FC,nm,'Mean'))];
  gdpc = comb('obs_gdp'); pic = comb('obs_pi'); ffrc = comb('obs_ffr');
  rsc  = 400*log(comb('rstar')); ogc = comb('ogap');
  if gapoff~=0; ogc = ogc - gapoff; end        % center the gap to match the paper's benchmark level
  Ntot = numel(gdpc); ycal = 1959 + (0:Ntot-1)/4;
  q4idx = @(Y) find(abs(ycal-(Y+0.75))<1e-6, 1);
  yr4   = @(v,Y) local_yoy_at(v, ycal, Y);        % Q4/Q4 annual = mean of the 4 quarters of year Y
  lvlQ4 = @(v,Y) local_q4(v, ycal, Y);
  paper = struct('yrs',[2026 2027 2028 2029], ...
     'g',[1.7 2.5 2.5 2.4], 'pi',[3.3 2.5 2.5 2.5], 'ffr',[4.0 3.9 3.9 3.9], ...
     'rs',[0.9 1.0 1.0 1.1], 'og',[-0.1 -0.1 0.1 0.2]);
  paperLR = struct('g',2.2,'pi',2.5,'ffr',3.8,'rs',1.1,'og',0.1);
  fid = fopen(fullfile(OUTDIR,'table5_forecast.csv'),'w');
  fprintf(fid,'Year,OutputGrowth_model,OutputGrowth_paper,CorePCE_model,CorePCE_paper,FFR_model,FFR_paper,rstar_model,rstar_paper,OutputGap_model,OutputGap_paper\n');
  fprintf('\n===================== TABLE 5: model vs paper =====================\n');
  fprintf('%-10s | %-13s | %-13s | %-13s | %-11s | %-12s\n','Year','OutputGrowth','CorePCE infl','Fed funds','r*','Output gap');
  fprintf('%-10s | %-13s | %-13s | %-13s | %-11s | %-12s\n','','model /paper','model /paper','model /paper','model/paper','model /paper');
  for i = 1:numel(paper.yrs)
    Y = paper.yrs(i);
    g = yr4(gdpc,Y) + POP_GROWTH; pv = yr4(pic,Y); fr = lvlQ4(ffrc,Y); rs = lvlQ4(rsc,Y); og = lvlQ4(ogc,Y);
    prow = @(x) sprintf('%s', nastr(x));
    fprintf('%-10d | %6s /%4.1f | %6s /%4.1f | %6s /%4.1f | %5s/%3.1f | %6s /%4.1f\n', ...
        Y, prow(g),paper.g(i), prow(pv),paper.pi(i), prow(fr),paper.ffr(i), ...
        prow(rs),paper.rs(i), prow(og),paper.og(i));
    fprintf(fid,'%d,%s,%.1f,%s,%.1f,%s,%.1f,%s,%.1f,%s,%.1f\n', Y, ...
        nastr(g),paper.g(i), nastr(pv),paper.pi(i), nastr(fr),paper.ffr(i), ...
        nastr(rs),paper.rs(i), nastr(og),paper.og(i));
  end
  % longer run 2031-35 (needs forecast>=40)
  glr = yravg(gdpc,ycal,2031,2035,'yoy') + POP_GROWTH; plr = yravg(pic,ycal,2031,2035,'yoy');
  frlr= yravg(ffrc,ycal,2031,2035,'q4');  rslr= yravg(rsc,ycal,2031,2035,'q4'); oglr=yravg(ogc,ycal,2031,2035,'q4');
  fprintf('%-10s | %6s /%4.1f | %6s /%4.1f | %6s /%4.1f | %5s/%3.1f | %6s /%4.1f\n', ...
    'Longer run', nastr(glr),paperLR.g, nastr(plr),paperLR.pi, nastr(frlr),paperLR.ffr, ...
    nastr(rslr),paperLR.rs, nastr(oglr),paperLR.og);
  fprintf(fid,'Longer run,%s,%.1f,%s,%.1f,%s,%.1f,%s,%.1f,%s,%.1f\n', ...
    nastr(glr),paperLR.g, nastr(plr),paperLR.pi, nastr(frlr),paperLR.ffr, nastr(rslr),paperLR.rs, nastr(oglr),paperLR.og);
  fclose(fid);
  % ---- output growth: per-capita vs aggregate breakdown (transparency; matches the note) ----
  fprintf('\n--------- Output growth: per-capita -> aggregate (Q4/Q4, %%) ---------\n');
  fprintf('%-10s | %-8s | %-6s | %-8s | %-6s\n','Year','per-cap','+pop','=aggr','paper');
  fid2 = fopen(fullfile(OUTDIR,'table5_growth_decomp.csv'),'w');
  fprintf(fid2,'Year,PerCapita_model,Population,Aggregate_model,Paper_aggregate\n');
  for i = 1:numel(paper.yrs)
    Y=paper.yrs(i); pc=yr4(gdpc,Y); ag=pc+POP_GROWTH;
    fprintf('%-10d | %8s | %6.1f | %8s | %6.1f\n', Y, nastr(pc), POP_GROWTH, nastr(ag), paper.g(i));
    fprintf(fid2,'%d,%s,%.2f,%s,%.1f\n', Y, nastr(pc), POP_GROWTH, nastr(ag), paper.g(i));
  end
  pclr = yravg(gdpc,ycal,2031,2035,'yoy'); aglr = pclr+POP_GROWTH;
  fprintf('%-10s | %8s | %6.1f | %8s | %6.1f\n','Longer run', nastr(pclr), POP_GROWTH, nastr(aglr), paperLR.g);
  fprintf(fid2,'Longer run,%s,%.2f,%s,%.1f\n', nastr(pclr), POP_GROWTH, nastr(aglr), paperLR.g);
  fclose(fid2);
  fprintf('  (per-capita is the model''s native object; aggregate = per-capita + %.1fpp population, the paper''s definition)\n', POP_GROWTH);
  fprintf('===================================================================\n');
  fprintf('(''n/a'' = beyond the current forecast horizon; set forecast=40 in the\n postest and re-run to fill 2029 and the 2031-35 longer run.)\n');
  fprintf('  wrote table5_forecast.csv, table5_growth_decomp.csv\n');
catch err
  warning('Table 5 failed: %s', err.message);
end

fprintf('Done. Outputs in ./%s/\n', OUTDIR);
end

% ======================= helpers ===========================================
% grouped shock decomposition of one endogenous variable -> (nG+1) x T + total
function [G, tot] = group_decomp(oo_, M_, vname, groups, gnames, exo, idx)
endo = cellstr(M_.endo_names);
ei = idx(vname, endo);
assert(~isempty(ei), sprintf('variable %s not in shock_decomposition', vname));
nexo = M_.exo_nbr;
sd = squeeze(oo_.shock_decomposition(ei, :, :));   % (nexo+2) x T
tot = sd(nexo+2, :);
T = size(sd,2);
G = zeros(numel(gnames), T);
for g = 1:numel(gnames)
  for s = groups.(gnames{g})
    j = idx(s{1}, exo);
    if ~isempty(j); G(g,:) = G(g,:) + sd(j,:); end
  end
end
other = tot - sum(G,1);            % measurement error + initial conditions + residual
G = [G; other];
end

% 4-quarter trailing mean along columns (time = dim 2), edge-safe
function Y = mov4(X)
[m,T] = size(X); Y = zeros(m,T);
for t = 1:T; k = max(1,t-3):t; Y(:,t) = mean(X(:,k),2); end
end
function y = mov4v(x)   % vector version
x = x(:); n = numel(x); y = zeros(n,1);
for t = 1:n; k = max(1,t-3):t; y(t) = mean(x(k)); end
end

% source structs robust to mode (flat) vs MH (nested under .Mean)
function SV = smoothed_src(oo_)
SV = struct();
if isfield(oo_,'SmoothedVariables')
  if isfield(oo_.SmoothedVariables,'Mean'); SV = oo_.SmoothedVariables.Mean; else; SV = oo_.SmoothedVariables; end
end
end
function FC = forecast_src(oo_)
FC = struct();
if isfield(oo_,'forecast') && isfield(oo_.forecast,'Mean');              FC = oo_.forecast;
elseif isfield(oo_,'MeanForecast') && isfield(oo_.MeanForecast,'Mean');  FC = oo_.MeanForecast;
elseif isfield(oo_,'PointForecast') && isfield(oo_.PointForecast,'Mean');FC = oo_.PointForecast;
end
end
function v = getf1(SV, nm)          % smoothed history vector or []
v = [];
if isfield(SV, nm); v = SV.(nm)(:); end
end
function v = getfc(FC, nm, fld)     % forecast Mean/HPDinf/HPDsup vector or []
v = [];
if isfield(FC, fld) && isfield(FC.(fld), nm); v = FC.(fld).(nm)(:); end
end
function m = cond_mean(cf, nm)      % conditional-forecast mean (drop origin element)
m = [];
if isfield(cf,'cond') && isfield(cf.cond,'Mean') && isfield(cf.cond.Mean, nm)
  mm = cf.cond.Mean.(nm); m = mm(2:end);
end
end
function y = colv(x); if isempty(x); y = []; else; y = x(:); end; end

% Q4/Q4 annual (= mean of the four quarters of year Y); NaN if unavailable
function v = local_yoy_at(series, ycal, Y)
qs = [Y+0.00 Y+0.25 Y+0.50 Y+0.75]; ii = zeros(1,4);
for k=1:4; f=find(abs(ycal-qs(k))<1e-6,1); if isempty(f); v=NaN; return; end; ii(k)=f; end
if max(ii)>numel(series); v=NaN; else; v = mean(series(ii)); end
end
function v = local_q4(series, ycal, Y)
f = find(abs(ycal-(Y+0.75))<1e-6,1);
if isempty(f) || f>numel(series); v=NaN; else; v = series(f); end
end
function v = yravg(series, ycal, Y0, Y1, mode)
acc = []; for Y=Y0:Y1
  if strcmp(mode,'yoy'); x=local_yoy_at(series,ycal,Y); else; x=local_q4(series,ycal,Y); end
  acc(end+1)=x; end
if any(isnan(acc)); v=NaN; else; v=mean(acc); end
end
function s = nastr(x); if isnan(x); s='  n/a'; else; s=sprintf('%5.1f',x); end; end

function yl = robustlim(v, p)
v = v(isfinite(v)); v = sort(v(:));
if isempty(v); yl=[-1 1]; return; end
n=numel(v); lo=v(max(1,round(p/100*n))); hi=v(min(n,max(1,round((1-p/100)*n))));
d=hi-lo; if d<=0; d=max(1,abs(hi)); end; yl=[lo-0.08*d, hi+0.08*d];
end

% PointForecast struct (param + shock uncertainty; wider bands than MeanForecast)
function PF = pointfc_src(oo_)
PF = struct();
if isfield(oo_,'PointForecast') && isfield(oo_.PointForecast,'Mean'); PF = oo_.PointForecast; end
end

% conditional-forecast 68% band for one var: [lo,hi] column vectors (or [])
function [lo,hi] = cond_ci(cf, nm)
lo=[]; hi=[];
if isfield(cf,'cond') && isfield(cf.cond,'ci') && isfield(cf.cond.ci, nm)
  c = cf.cond.ci.(nm);
  if size(c,1)==2;      lo=c(1,:)'; hi=c(2,:)';    % 2 x H  (Dynare default)
  elseif size(c,2)==2;  lo=c(:,1);  hi=c(:,2);     % H x 2  (fallback)
  end
end
end

% 4q trailing mean of a forecast vector, stitched onto history (for YoY panels)
function w = yoystitch(hist, v)
if isempty(v); w=v; return; end
tmp = mov4v([hist(:); v(:)]); w = tmp(numel(hist)+1:end);
end

% read rstar_external.csv (date,LM_median,HLW_rstar) -> vectors aligned to yrr (NaN where absent)
function ext = read_external_rstar(yrr)
ext = struct('LM',[],'HLW',[]);
fn = 'rstar_external.csv';
if ~exist(fn,'file'); return; end
LM = nan(numel(yrr),1); HLW = nan(numel(yrr),1);
fid = fopen(fn,'r'); hdr = fgetl(fid);  %#ok<NASGU>
while ~feof(fid)
  ln = fgetl(fid); if ~ischar(ln) || isempty(ln); continue; end
  p = strsplit(ln, ',');
  d = strtrim(p{1});                          % 'YYYYQq'
  q = str2double(d(6));                        % quarter digit
  Y = str2double(d(1:4));                       % year
  if isnan(Y)||isnan(q); continue; end
  yv = Y + (q-1)/4;
  k = find(abs(yrr - yv) < 1e-6, 1);
  if isempty(k); continue; end
  if numel(p)>=2 && ~isempty(strtrim(p{2})); LM(k)  = str2double(p{2}); end
  if numel(p)>=3 && ~isempty(strtrim(p{3})); HLW(k) = str2double(p{3}); end
end
fclose(fid);
ext.LM = LM; ext.HLW = HLW;
end

% y-limits from a set of series (min/max + padding, optional minimum span), NaN-safe
function yl = niceylim(vals, padfrac, minspan)
v = vals(isfinite(vals));
if isempty(v); yl=[-1 1]; return; end
lo=min(v); hi=max(v); d=hi-lo; if d<=0; d=max(1,abs(hi)); end
yl=[lo-padfrac*d, hi+padfrac*d];
if nargin>=3 && ~isempty(minspan) && (yl(2)-yl(1))<minspan
  c=mean(yl); yl=[c-minspan/2, c+minspan/2];
end
end

function lighttheme(fh)
set(fh,'Color','w','InvertHardcopy','off');
ax = findall(fh,'type','axes');
for k=1:numel(ax)
  set(ax(k),'Color','w','XColor',[0.2 0.2 0.2],'YColor',[0.2 0.2 0.2],'Layer','top');
end
end
