% ============================================================================
%  repair_mode_hessian.m
%
%  A newrat mode search can return a Hessian that is complete and positive
%  definite EXCEPT for a small number of rows/columns that come back all-NaN.
%  That happens on the near-unit-root persistences: the numerical differencing
%  step walks the parameter across the stationarity boundary and the likelihood
%  is not finite on both sides, so the second difference is undefined.
%
%  A NaN row makes the whole matrix unusable as an MH proposal, but the other
%  directions are fine.  So the sensible repair is surgical: keep the good block
%  exactly as estimated, and give each broken parameter an INDEPENDENT proposal
%  with a sensible step size (zero off-diagonals, a diagonal chosen to imply the
%  step you want).  That is far better than throwing the whole matrix away and
%  falling back to MCMC_jumping_covariance = prior_variance, which discards
%  every good direction and every posterior correlation.
%
%  IN THE SHIPPED RUN (non-deductible utilization branch) exactly one row was
%  broken, rho_mup; rho_tr, which had been broken under the earlier branch, came
%  back clean.  The script reports what it finds, so do not assume the count.
%
%  Choosing the step.  diag(hh) = 1/sd^2 for that parameter, so
%      sd 0.0005 -> diag 4.0e+06
%      sd 0.0010 -> diag 1.0e+06
%      sd 0.0020 -> diag 2.5e+05
%  The parameters that break sit ~0.001-0.002 below 1 and the model is near
%  non-stationary there, so a small step is right.  The working near-unit-root
%  persistences in this model imply sd ~0.004 (rho_pistar) to ~0.010 (rho_g);
%  0.001 is deliberately tighter because these are the ones the differencing
%  choked on.  mh_tune_jscale rescales the WHOLE proposal afterwards anyway, so
%  this sets the relative step, not the absolute one.
%
%     >> repair_mode_hessian('frbstl_us_est_fastlr_mode_ndc.mat')
%     >> repair_mode_hessian('frbstl_us_est_fastlr_mode_ndc.mat', 0.0005)
%
%  Writes <input>_fixed.mat.  Point mode_file at it (without the .mat) and keep
%  MCMC_jumping_covariance = hessian.
% ============================================================================
function repair_mode_hessian(MODEFILE, STEP_SD)

if nargin<1 || isempty(MODEFILE); MODEFILE='frbstl_us_est_fastlr_mode_ndc.mat'; end
if nargin<2 || isempty(STEP_SD);  STEP_SD = 0.0010; end

S = load(MODEFILE);
assert(isfield(S,'hh') && isfield(S,'xparam1'), 'not a Dynare mode file (need hh and xparam1)');
hh = S.hh; x = S.xparam1(:); n = numel(x);
if isfield(S,'parameter_names'); nm = cellstr(S.parameter_names); else; nm = repmat({''},n,1); end

% ENTIRELY non-finite rows/cols only.  Using any() instead of all() here is a
% trap: every row has 2 non-finite entries (the ones in the two broken COLUMNS),
% so any() flags all 63 rows and the repair would blank the whole matrix.
bad = find(all(~isfinite(hh),2) | all(~isfinite(hh),1)');
fprintf('\n=== repair_mode_hessian ===\n');
fprintf('file      : %s\n', MODEFILE);
fprintf('fval      : %.6f\n', S.fval);
fprintf('size      : %dx%d,  non-finite entries %d\n', n, n, sum(~isfinite(hh(:))));
if isempty(bad)
  fprintf('nothing to repair -- the Hessian is already finite.\n'); return
end
fprintf('broken rows/cols (%d):\n', numel(bad));
for k=1:numel(bad)
  fprintf('   %2d  %-18s mode %.6f\n', bad(k), nm{bad(k)}, x(bad(k)));
end

good = setdiff((1:n)', bad);
Hg = hh(good,good);
assert(all(isfinite(Hg(:))), 'the surviving block still has non-finite entries -- inspect by hand');
Hg = (Hg+Hg.')/2;
eg = eig(Hg);
fprintf('\nsurviving %dx%d block: eig %.4e .. %.4e, negatives %d\n', ...
        numel(good), numel(good), min(eg), max(eg), sum(eg<0));

% ---- repair: zero the coupling, set an independent diagonal ----------------
dval = 1/STEP_SD^2;
H = hh;
H(bad,:) = 0; H(:,bad) = 0;
for k=1:numel(bad); H(bad(k),bad(k)) = dval; end
H = (H+H.')/2;

e = eig(H);
fprintf('repaired  : diagonal %.4e for the two (implied step sd %.4g)\n', dval, STEP_SD);
fprintf('            eig %.4e .. %.4e, negatives %d\n', min(e), max(e), sum(e<0));
[~,p] = chol(H);
if p==0
  fprintf('            CHOLESKY OK -> usable as MCMC_jumping_covariance = hessian\n');
else
  error('repaired matrix is still not positive definite -- raise STEP_SD or inspect');
end

hh = H;                                                     %#ok<NASGU>
xparam1 = S.xparam1; fval = S.fval;                          %#ok<NASGU>
out = regexprep(MODEFILE, '\.mat$', '_fixed.mat');
if isfield(S,'parameter_names')
  parameter_names = S.parameter_names;                       %#ok<NASGU>
  save(out, 'xparam1', 'hh', 'fval', 'parameter_names');
else
  save(out, 'xparam1', 'hh', 'fval');
end
fprintf('\nwrote %s\n', out);
fprintf('-> set  mode_file = %s  in the pilot/MCMC .mod\n', regexprep(out,'\.mat$',''));
fprintf('   and keep MCMC_jumping_covariance = hessian.\n');
fprintf('   NOTE: this changes only the PROPOSAL. The posterior being sampled is\n');
fprintf('   unchanged, so the draws remain valid; a poor proposal costs mixing, not\n');
fprintf('   correctness. Check the acceptance rate in the pilot.\n');
end
