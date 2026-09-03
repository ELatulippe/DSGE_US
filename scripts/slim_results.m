function slim_results(INFILE, OUTFILE)
% slim_results  --  shrink a post-estimation *_results.mat so it can be uploaded.
%
% Keeps only what the post-estimation figure scripts actually read
% (make_paper_figures, make_detailed_decomp, diag_forecast, diag_gap):
%   forecasts, smoother, shock_decomposition, steady state, decision rule.
% Drops the large theoretical-moment / IRF / variance-decomposition blobs,
% which those scripts do not use (FEVD & IRF come from fevd_run_results.mat).
%
% USAGE:
%   >> slim_results                         % in/out default names
%   >> slim_results('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat')
%
% Output: frbstl_us_est_fastlr_lean_results.mat (v7) in the current folder,
% ready to drop into  results/.

if nargin<1 || isempty(INFILE);  INFILE  = 'frbstl_us_est_fastlr_lean_results.mat'; end
if nargin<2 || isempty(OUTFILE); OUTFILE = 'frbstl_us_est_fastlr_lean_results.mat'; end

S = load(INFILE);
oo_ = S.oo_;

% ---- report per-field sizes so you can see what is big ---------------------
fprintf('\n%-40s %10s\n','oo_ field','MB');
f = fieldnames(oo_);
for i=1:numel(f)
  x = oo_.(f{i}); w = whos('x');           %#ok<NASGU>
  fprintf('  oo_.%-34s %8.1f\n', f{i}, w.bytes/1e6);
end

% ---- drop large fields that the post-estimation figure scripts never read --
%  The four smoother containers below (~50 MB each) are unused: the scripts only
%  ever read oo_.SmoothedVariables.Mean, oo_.MeanForecast/PointForecast/
%  conditional_forecast, oo_.shock_decomposition, oo_.dr, oo_.steady_state.
drop = {'Smoother','SmoothedShocks','Constant','UpdatedVariables', ...   % ~200 MB, unused
        'convergence','posterior_density', ...                          % diagnostics, unused
        'gamma_y','autocorr','variance_decomposition', ...
        'conditional_variance_decomposition','conditional_variance_decomposition_ME', ...
        'irfs','mean','var','PosteriorTheoreticalMoments','endo_simul','exo_simul'};
for i=1:numel(drop)
  if isfield(oo_,drop{i}); oo_ = rmfield(oo_,drop{i}); end
end

% ---- collapse SmoothedVariables to just the posterior-mean sub-struct -------
%  (drops the HPD / deciles / density arrays the scripts never read: 51 MB -> ~0.5 MB)
if isfield(oo_,'SmoothedVariables') && isfield(oo_.SmoothedVariables,'Mean')
  sv = struct('Mean', oo_.SmoothedVariables.Mean);
  oo_.SmoothedVariables = sv;
end

M_            = S.M_;
options_      = S.options_;
if isfield(S,'estim_params_'); estim_params_ = S.estim_params_; else; estim_params_ = struct(); end
if isfield(S,'bayestopt_');    bayestopt_    = S.bayestopt_;    else; bayestopt_    = struct(); end

save(OUTFILE,'oo_','M_','options_','estim_params_','bayestopt_','-v7');
w = dir(OUTFILE);
fprintf('\nwrote %s  (%.1f MB)\n', OUTFILE, w.bytes/1e6);
fprintf('If still too large, also rmfield oo_.Smoother subfields you do not need,\n');
fprintf('or paste the size table above and I will target the exact culprit.\n');
end
