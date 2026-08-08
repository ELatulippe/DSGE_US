results/ -- estimation outputs, so the figures can be regenerated WITHOUT re-running Dynare.
This folder is now COMPLETE: every figure/table script runs with no Dynare run.

SHIPPED:
  frbstl_us_est_fastlr_lean_results.mat  POST-ESTIMATION output, slimmed (~20 MB) by scripts/slim_results.m
                                         -> make_paper_figures, diag_forecast, diag_gap, make_detailed_decomp
  fevd_run_results.mat                   stoch_simul(order=1, irf=40, cond. var. decomp.)
                                         -> make_fevd, make_irf_panel, make_counterfactual
  cond_constrate.mat                     conditional-forecast scenario 1 (constant rate)   -> Fig 8, counterfactual baseline
  cond_infltarget.mat                    conditional-forecast scenario 2 (inflation target)-> Fig 8

NOTE on the slimmed post-estimation file: scripts/slim_results.m dropped the unused smoother
containers (oo_.Smoother/SmoothedShocks/Constant/UpdatedVariables, ~50 MB each) and reduced
oo_.SmoothedVariables to its .Mean sub-struct. All fields the scripts read are intact
(dr, steady_state, shock_decomposition, MeanForecast, PointForecast, conditional_forecast,
SmoothedVariables.Mean). Re-run the estimation only if you need the full posterior distribution.

OPTIONAL: cond_tighter.mat (higher-for-longer counterfactual) -- drop in if you generated it.

WHICH SCRIPT READS WHICH FILE (each loads exactly one *_results.mat):
  frbstl_us_est_fastlr_lean_results.mat  (post-estimation: smoother, forecast=40, shock_decomposition)
       -> make_paper_figures, diag_forecast, diag_gap, make_detailed_decomp
  fevd_run_results.mat                   (stoch_simul order=1, irf=40, conditional_variance_decomposition)
       -> make_fevd, make_irf_panel, make_counterfactual
  The scripts' search lists are ordered so each grabs its correct file first when both are in results/.

With those two *_results.mat present here, run make_paper_figures / make_detailed_decomp /
make_fevd / make_irf_panel / make_counterfactual / diag_gap / diag_forecast with NO Dynare run.

NOT shipped: the MH draws (frbstl_us_est_fastlr_lean/metropolis/, ~1 GB). They are only needed
to *re-run* the post-estimation step; the *_results.mat above already contains everything the
figure scripts read.
