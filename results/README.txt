results/ -- estimation outputs, so the figures regenerate WITHOUT re-running Dynare.

THIS FOLDER IS EMPTY IN THE REPOSITORY AND MUST BE POPULATED FROM YOUR OWN RUN.
The two .mat files below are 20 MB and 0.9 MB; drop them in after the post-estimation
and fevd_run steps and every make_* / diag_* script runs with no Dynare step.

WHAT GOES HERE:
  frbstl_us_est_fastlr_lean_results.mat  POST-ESTIMATION output, slimmed by scripts/slim_results.m
                                         -> make_paper_figures, diag_*, make_detailed_decomp
  fevd_run_results.mat                   stoch_simul(order=1, irf=40, cond. var. decomp.)
                                         -> make_fevd, make_irf_panel, make_counterfactual

The conditional-scenario files (cond_constrate.mat, cond_infltarget.mat, cond_tighter.mat)
live in the RUN FOLDER (the repository root), not here: make_paper_figures and
make_counterfactual look for them by name in the current MATLAB folder.

SLIMMING.  scripts/slim_results.m drops the unused smoother containers
(oo_.Smoother / SmoothedShocks / Constant / UpdatedVariables, ~50 MB each) and reduces
oo_.SmoothedVariables to its .Mean sub-struct.  Every field the scripts read survives:
dr, steady_state, shock_decomposition, MeanForecast, PointForecast, conditional_forecast,
SmoothedVariables.Mean.  Re-run the estimation only if you need the full posterior
distribution.

    >> slim_results('frbstl_us_est_fastlr_lean/Output/frbstl_us_est_fastlr_lean_results.mat', ...
                    'results/frbstl_us_est_fastlr_lean_results.mat')

WHICH SCRIPT READS WHICH FILE (each loads exactly one *_results.mat):
  frbstl_us_est_fastlr_lean_results.mat  (smoother, forecast=40, shock_decomposition)
       -> make_paper_figures, diag_forecast, diag_gap, diag_flexblock, diag_ogap_offset,
          diag_inflation_target, make_detailed_decomp, check_postest
  fevd_run_results.mat                   (stoch_simul order=1, irf=40, cond. var. decomp.)
       -> make_fevd, make_irf_panel, make_counterfactual, diag_fevd_split
The scripts' search lists are ordered so each grabs its correct file first when both are here.

NOT SHIPPED: the MH draws (frbstl_us_est_fastlr_lean/metropolis/, ~1 GB).  They are needed
only to re-run the post-estimation step or to recover the full posterior distribution, which
the shipped objects already summarize.

MAKE SURE THESE COME FROM THE CURRENT SPECIFICATION.  Results produced before the
capital-utilization branch was corrected are NOT interchangeable with the current .mod files:
the posterior, every figure and Table 5 all move.  If in doubt, re-run rather than reuse.
