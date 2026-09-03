% ============================================================================
%  mode_from_m1.m
%
%  SALVAGE a mode search that crashed mid-optimization.
%
%  Dynare's newrat (mode_compute = 5) writes m1.mat in the working directory at
%  the end of EVERY completed iteration.  So when the line search later walks
%  into a region where the steady state goes complex and static_resid throws,
%  the last GOOD iterate is still on disk -- Dynare just never got to write the
%  final <fname>_mode.mat.
%
%  This turns that m1.mat into a proper Dynare mode file you can point
%  mode_file at.  Parameter names are copied from a reference mode file (the
%  estimated_params block must be the same length and order, which it is as
%  long as you have not added or removed an estimated parameter).
%
%     >> mode_from_m1
%     >> mode_from_m1('m1.mat', 'frbstl_us_est_fastlr_mode_ded_fixed.mat', ...
%                     'frbstl_us_est_fastlr_mode_ndc.mat')
%
%  NOTE ON THE HESSIAN.  m1.mat's hh is newrat's WORKING Hessian at that
%  iteration, not the final outer-product-gradient one Dynare would have
%  written.  It is fine as a restart point for another mode pass.  Do NOT feed
%  it straight to MCMC_jumping_covariance = hessian without checking it: run
%  repair_mode_hessian on the output first, or better, do one more short mode
%  pass from here so Dynare computes the Hessian properly.
% ============================================================================
function mode_from_m1(M1FILE, REFMODE, OUTFILE)

if nargin<1 || isempty(M1FILE);  M1FILE  = 'm1.mat'; end
if nargin<2 || isempty(REFMODE); REFMODE = 'frbstl_us_est_fastlr_mode_ded_fixed.mat'; end
if nargin<3 || isempty(OUTFILE); OUTFILE = 'frbstl_us_est_fastlr_mode_ndc.mat'; end

assert(exist(M1FILE,'file')==2, ...
  ['%s not found. newrat writes it in the folder you ran dynare from, and only ' ...
   'once iteration 1 has completed. If it is missing, just re-run the bounded ' ...
   'step 1 from scratch -- it is about 10 s per iteration.'], M1FILE);

S = load(M1FILE);
fprintf('\n=== mode_from_m1 ===\nm1 file : %s\ncontains: %s\n', M1FILE, strjoin(fieldnames(S)',', '));

% ---- the iterate -----------------------------------------------------------
if     isfield(S,'x');       x = S.x(:);
elseif isfield(S,'xparam1'); x = S.xparam1(:);
else;  error('no x / xparam1 in %s -- inspect it by hand with whos(''-file'',''%s'')', M1FILE, M1FILE);
end

% ---- the objective value ---------------------------------------------------
fval = NaN;
for f = {'fval0','fval','f'}
    if isfield(S,f{1}); v = S.(f{1}); fval = v(find(v~=0,1,'last')); break; end
end
if isempty(fval); fval = NaN; end

% ---- the Hessian -----------------------------------------------------------
if     isfield(S,'hh');  hh = S.hh;
elseif isfield(S,'H1');  hh = S.H1;
elseif isfield(S,'hhg'); hh = S.hhg;
else;  error('no hh / H1 / hhg in %s', M1FILE);
end
hh = (hh+hh.')/2;

% ---- names from the reference mode file ------------------------------------
parameter_names = [];
if exist(REFMODE,'file')==2
    R = load(REFMODE);
    if isfield(R,'parameter_names')
        parameter_names = R.parameter_names;
        assert(numel(cellstr(parameter_names))==numel(x), ...
          ['%s has %d parameters but m1.mat has %d. The estimated_params block ' ...
           'changed length -- do not reuse these names.'], ...
           REFMODE, numel(cellstr(parameter_names)), numel(x));
    end
else
    fprintf('reference mode file %s not found; writing without parameter_names.\n', REFMODE);
end

fprintf('recovered: %d parameters, fval %.6f\n', numel(x), fval);
e = eig(hh); [~,p] = chol(hh);
fprintf('hessian  : eig %.4e .. %.4e, negatives %d, cholesky %s\n', ...
        min(e), max(e), sum(e<0), tern(p==0,'OK','FAILS'));
nnf = sum(~isfinite(hh(:)));
if nnf>0; fprintf('           %d non-finite entries -- run repair_mode_hessian next\n', nnf); end

if ~isempty(parameter_names)
    nm = cellstr(parameter_names);
    fprintf('\nfirst few / last few:\n');
    idx = [1:min(4,numel(x)), max(1,numel(x)-3):numel(x)];
    for k = unique(idx); fprintf('   %2d  %-18s %.6f\n', k, nm{k}, x(k)); end
end

xparam1 = x;                                                  %#ok<NASGU>
if isempty(parameter_names); save(OUTFILE,'xparam1','hh','fval');
else;                        save(OUTFILE,'xparam1','hh','fval','parameter_names'); end
fprintf('\nwrote %s\n', OUTFILE);
fprintf('-> point mode_file at %s and re-run the mode step.\n', regexprep(OUTFILE,'\.mat$',''));
end

function s = tern(c,a,b); if c; s=a; else; s=b; end; end
