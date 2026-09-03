% ============================================================================
%  diag_fevd_split.m  --  per-shock-TYPE FEVD for actual vs natural output.
%
%  WHY.  make_fevd pools technology and markup shocks into one "Supply" group:
%        Supply = {e_gz, e_gn, e_zeta, e_mup, e_muw, e_gn_c}
%  and that makes the headline flex-block comparison unreadable, because the
%  constant-markup convention means e_mup and e_muw cannot reach gdp_f AT ALL.
%  "gdp 19.9% supply vs gdp_f 0.4% supply" is therefore mostly the statement
%  that the markup shocks were removed by construction, not a finding about
%  technology. Splitting Supply into TECH and MARKUP makes it honest.
%
%  This reuses make_fevd's engine -- the FEVD is built directly from the order-1
%  decision rule (ghx, ghu) with variables looked up BY NAME in DR order. An
%  earlier version of this script indexed oo_.conditional_variance_decomposition
%  by a hardcoded row order, which is fragile and silently returns the wrong
%  variable if the stoch_simul list ever changes.
%
%     >> diag_fevd_split
%     >> diag_fevd_split('fevd_run/Output/fevd_run_results.mat')
% ============================================================================
function diag_fevd_split(RESULTS)

if nargin<1 || isempty(RESULTS)
  cand={'fevd_run/Output/fevd_run_results.mat','fevd_run_results.mat', ...
        'results/fevd_run_results.mat'};
  RESULTS=''; for i=1:numel(cand); if exist(cand{i},'file'); RESULTS=cand{i}; break; end; end
  assert(~isempty(RESULTS),'no fevd_run results .mat -- run  dynare fevd_run.mod  first.');
end
S=load(RESULTS); oo_=S.oo_; M_=S.M_;

% ---- order-1 state space in DR order (same as make_fevd) -------------------
dr=oo_.dr; ghx=dr.ghx; ghu=dr.ghu;
nstatic = pickfield(M_,dr,'nstatic');
nsx  = size(ghx,2);
sidx = nstatic + (1:nsx);
assert(max(sidx)<=size(ghx,1),'state-row indexing off');
T = ghx(sidx,:); R = ghu(sidx,:);
nexo = M_.exo_nbr; sig2 = diag(M_.Sigma_e);
name_dr = cellstr(M_.endo_names); name_dr = name_dr(dr.order_var);
exo = cellstr(M_.exo_names);

% ---- the finer grouping ----------------------------------------------------
G = struct('Tech',    {{'e_gz','e_gn','e_zeta','e_gn_c'}}, ...
           'Markup',  {{'e_mup','e_muw'}}, ...
           'Demand',  {{'e_chi','e_vth','e_chi_c'}}, ...
           'Fiscal',  {{'e_g','e_tr','e_tr_c'}}, ...
           'Monetary',{{'e_mp','e_pistar'}}, ...
           'Oil',     {{'e_oS'}});
gn = fieldnames(G)';
gmap = zeros(nexo,1);
for g=1:numel(gn)
  for s = G.(gn{g}); j=find(strcmp(s{1},exo),1); if ~isempty(j); gmap(j)=g; end; end
end

H=[1 4 8 16 40 1000]; Hlab={'1q','1y','2y','4y','10y','\infty'};

fprintf('\n============ PER-SHOCK-TYPE FEVD: actual vs natural output ============\n');
fprintf('results: %s\n', RESULTS);
fprintf('Supply is split into TECHNOLOGY and MARKUP. The flex block holds mu_p and\n');
fprintf('mu_w at steady state, so the MARKUP column for gdp_f is ~0 BY CONSTRUCTION.\n\n');

res = struct();
for v = {'gdp','gdp_f'}
  vn=v{1}; jj=find(strcmp(vn,name_dr),1);
  if isempty(jj); fprintf('  %s not in the decision rule -- add it to the stoch_simul list.\n\n',vn); continue; end
  Vg = fevd_one(jj, ghx, ghu, T, R, sig2, nexo, gmap, numel(gn), H);
  res.(vn) = Vg;
  fprintf('%-22s |', upper(vn));
  fprintf('%9s', gn{:}); fprintf('\n');
  fprintf('%s\n', repmat('-',1,22+1+9*numel(gn)));
  for c=1:numel(H)
    fprintf('  horizon %-11s |', Hlab{c});
    fprintf('%9.1f', Vg(c,:)); fprintf('\n');
  end
  fprintf('\n');
end

% ---- the read --------------------------------------------------------------
if isfield(res,'gdp') && isfield(res,'gdp_f')
  iT=find(strcmp(gn,'Tech')); iM=find(strcmp(gn,'Markup'));
  iD=find(strcmp(gn,'Demand')); iF=find(strcmp(gn,'Fiscal'));
  fprintf('READ (1-quarter horizon)\n');
  fprintf('  technology : gdp %.1f%%   gdp_f %.1f%%\n', res.gdp(1,iT), res.gdp_f(1,iT));
  fprintf('  markup     : gdp %.1f%%   gdp_f %.1f%%   <- 0 by construction for gdp_f\n', ...
          res.gdp(1,iM), res.gdp_f(1,iM));
  fprintf('  demand     : gdp %.1f%%   gdp_f %.1f%%\n', res.gdp(1,iD), res.gdp_f(1,iD));
  fprintf('  fiscal     : gdp %.1f%%   gdp_f %.1f%%\n\n', res.gdp(1,iF), res.gdp_f(1,iF));
  if res.gdp(1,iM) > 2*res.gdp(1,iT)
    fprintf('  -> gdp''s "Supply" share is mostly MARKUP, and markups are absent from the\n');
    fprintf('     flex block by construction. That, not a broken flex allocation, is why\n');
    fprintf('     the two series barely co-move at business-cycle frequencies.\n');
  end
  if res.gdp_f(1,iT) < 15
    fprintf('  -> natural output has little technology content in ABSOLUTE terms too.\n');
    fprintf('     Both series are detrended by the same stochastic trend Gamma, which\n');
    fprintf('     removes the PERMANENT TFP component from both, and rho_gz near zero\n');
    fprintf('     leaves a negligible transitory remainder. Natural output here is a\n');
    fprintf('     preference-and-fiscal object, not a technology object -- worth saying\n');
    fprintf('     plainly in the note, since Figs 5b, 6 and 10b all rest on it.\n');
  end
end
fprintf('======================================================================\n');
end

% ---- helpers ---------------------------------------------------------------
function Vg = fevd_one(jj, ghx, ghu, T, R, sig2, nexo, gmap, ng, H)
  bigH=max(H); Hset=sort(H); ptr=1; store=zeros(nexo,numel(H));
  contrib=zeros(nexo,1); Sprev=R;
  MA = ghu(jj,:)';  contrib = contrib + (MA.^2).*sig2;
  while ptr<=numel(Hset) && Hset(ptr)==1; store(:,ptr)=contrib; ptr=ptr+1; end
  for h=1:bigH-1
    MA = (ghx(jj,:)*Sprev)';
    contrib = contrib + (MA.^2).*sig2;
    Sprev = T*Sprev;
    while ptr<=numel(Hset) && Hset(ptr)==h+1; store(:,ptr)=contrib; ptr=ptr+1; end
  end
  [~,ord]=ismember(H,Hset); storeH=store(:,ord);
  Vg=zeros(numel(H),ng);
  for c=1:numel(H)
    tot=sum(storeH(:,c)); if tot<=0; tot=1; end
    for g=1:ng; Vg(c,g)=100*sum(storeH(gmap==g,c))/tot; end
  end
end

function v = pickfield(M_, dr, f)
  if isfield(M_,f); v=M_.(f); elseif isfield(dr,f); v=dr.(f);
  else; error('cannot find %s in M_ or oo_.dr',f); end
end
