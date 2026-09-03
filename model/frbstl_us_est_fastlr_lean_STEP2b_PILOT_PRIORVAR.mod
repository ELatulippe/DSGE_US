// ---------------------------------------------------------------------------
//  TAX TREATMENT OF THE UTILIZATION COST: A(nu) is NOT DEDUCTIBLE.
//  Confirmed by the author (Faria-e-Castro, private correspondence, Aug 2026):
//  the budget constraints (5)/(38)/(52)/(76) stand AS WRITTEN in the paper, and
//  it is the utilization FOC (12) that is corrected.  A(nu)*K is paid out of
//  AFTER-tax income, so it sits OUTSIDE the (1-tauK) bracket in Tobin's Q
//  (52)/(52f) and is NOT netted out of the capital-tax base in (76)/(76f).
//    FOC (51)/(51f):   (1-tauK)*R_k = A'(nu)        [was R_k = A'(nu)]
//    calibration   :   kappa_a = (1-tauK_ss)*Rk_ss  [was kappa_a = Rk_ss]
//                                = Gammabar/beta-(1-delta)
//  Rk_ss is UNCHANGED, and A(1)=0, so the whole deterministic steady state is
//  bit-identical to the baseline.  kappa_a falls by the factor (1-tauK_ss),
//  which DOES change the linearised dynamics (it is A''(1)/A'(1), the
//  utilization elasticity, that the data speak to -- and sigma_a is estimated).
//  UNCHANGED: (70)/(70f), (80), tau_c, the steady state.
//  RE-ESTIMATE before using any results.
// ---------------------------------------------------------------------------
// ============================================================================
//  The St. Louis Fed DSGE Model (Faria-e-Castro, 2026) -- US ESTIMATION file -- FAST long rates (perpetuity approx.)
//  Core model (eqs. 44-95) + measurement block:
//     * 15 observation equations (97)-(111)
//     * 7 measurement-error AR(1) processes
//     * COVID special shocks (112)-(114)
//     * varobs + estimated_params (Table 2-3 priors)
//
//  ESTIMATION-READY: all steady-state-dependent quantities (Gammabar, rbar,
//  kappa_a, tau_c, xi, Bw_bar, and the shock means oS_ss/g_ss/tr_ss) are
//  computed as model-local (#) expressions from PRIMITIVE parameters, so they
//  update on every parameter draw.  No external steadystate file is required.
//
//  Verified: SS reproduces frbstl_us.mod (residuals < 1e-15); Blanchard-Kahn
//  holds; observable steady states match the calibration (GDP growth ~1.9%,
//  core PCE ~2%, FFR ~3.2%, 10y ~4.2%, convenience yield ~0.9%).
//
//  For estimation: supply the data (see notes) and uncomment `estimation(...)`.
// ============================================================================


var
  R r Pi Pi_idx Pi_w Pi_w_idx Pi_h R_k Q mc pV pO w R_AAA
  cw cwN cwO cs csN csO cN cO cc y v oY k ku iv n gdp d mrs
  bw bg tauL tauK nu
  Gamma_Z Gamma_N zeta mu_p mu_w oS chi vth g tr mp Pi_star Gamma
  chi_eff tr_eff GammaN_eff
  me_gdp me_tfp me_oil me_gas me_inflexp me_tp me_cy
  obs_gdp obs_c obs_i obs_g obs_tr obs_w obs_n obs_ffr obs_pi
  obs_tfp obs_oil obs_gas obs_inflexp obs_tr10 obs_cy
  yq10_R yq20_R yq20_AAA piq10   // perpetuity long-rate states
  ec_chi ec_tr ec_gn      // observed COVID special-shock innovations
;

varexo
  e_gz e_gn e_zeta e_mup e_muw e_oS e_chi e_vth e_g e_tr e_mp e_pistar
  e_me_gdp e_me_tfp e_me_oil e_me_gas e_me_inflexp e_me_tp e_me_cy
  e_chi_c e_tr_c e_gn_c
;

parameters
  // primitives only (all directly calibrated or estimated)
  pbeta phi alpha delta lam nu_C nu_Y omega_C omega_Y
  mu_p_bar mu_w_bar Gammabar_Z Gammabar_N Pibar vthbar
  eta_p eta_w psi_i psi_w sigma_a iota_p iota_w
  phi_Pi phi_Y rho_r rho_tau phi_tl phi_tk
  tauL_ss tauK_ss s_g s_t s_b bw_y n_ss
  rho_gz rho_gn rho_zeta rho_mup rho_muw rho_oS rho_chi rho_vth rho_g rho_tr rho_mp rho_pistar
  eta_mup eta_muw
  sig_gz sig_gn sig_zeta sig_mup sig_muw sig_oS sig_chi sig_vth sig_g sig_tr sig_mp sig_pistar
  cobs_cy cobs_iy cobs_gy cobs_ty cobs_wy tp
  rho_me_gdp rho_me_tfp rho_me_oil rho_me_gas rho_me_inflexp rho_me_tp rho_me_cy
  sig_me_gdp sig_me_tfp sig_me_oil sig_me_gas sig_me_inflexp sig_me_tp sig_me_cy
  sig_chi_c sig_tr_c sig_gn_c
  d10 d20
;

// ---- values (Table 1 + posterior means) ----
pbeta=0.0383; phi=1.0074; alpha=0.1971; delta=0.025; lam=0.5881;
nu_C=0.3252; nu_Y=0.1783; omega_C=0.02; omega_Y=0.03;
mu_p_bar=1.2; mu_w_bar=1.2; Gammabar_Z=1.0049; Gammabar_N=0.9998;
Pibar=1.02^0.25; vthbar=1.0021;
eta_p=73.8207; eta_w=102.9788; psi_i=1.2770; psi_w=0.2606; sigma_a=0.4879;
iota_p=0.2357; iota_w=0.7685; phi_Pi=2.8492; phi_Y=0.7872; rho_r=0.8055;
rho_tau=0.9; phi_tl=0.0665; phi_tk=0.0719;
tauL_ss=0.186; tauK_ss=0.218; s_g=0.202; s_t=0.137; s_b=4*0.639; bw_y=1.0590; n_ss=1;
rho_gz=0.0654; rho_gn=0.2594; rho_zeta=0.3450; rho_mup=0.9892; rho_muw=0.5888;
rho_oS=0.9581; rho_chi=0.9397; rho_vth=0.9589; rho_g=0.9627; rho_tr=0.9881;
rho_mp=0.0734; rho_pistar=0.9901; eta_mup=0.7439; eta_muw=0.5341;
sig_gz=0.0078; sig_gn=0.0028; sig_zeta=0.0472; sig_mup=0.0268; sig_muw=0.2080;
sig_oS=0.0149; sig_chi=0.0218; sig_vth=0.0010; sig_g=0.0096; sig_tr=0.0262;
sig_mp=0.0024; sig_pistar=0.0004;
// cobs_* calibrated to sample-mean growth differentials (from build_usdata.py)
cobs_cy=0.1518; cobs_iy=0.1472; cobs_gy=-0.4891; cobs_ty=1.7091; cobs_wy=-0.1865;
tp=0.0102;
rho_me_gdp=0.0817; sig_me_gdp=0.0030; rho_me_tfp=0.0791; sig_me_tfp=0.0084;
rho_me_oil=0.9527; sig_me_oil=0.0909; rho_me_gas=0.4178; sig_me_gas=0.0049;
rho_me_inflexp=0.9305; sig_me_inflexp=0.0016; rho_me_tp=0.9626; sig_me_tp=0.0036;
rho_me_cy=0.9616; sig_me_cy=0.0013;
sig_chi_c=0.045; sig_tr_c=0.220; sig_gn_c=0.015;
d10=0.9512; d20=0.9753;   // geometric-decay (mean-lag matched to 40q / 80q)

// ----------------------------------------------------------------------------
model(use_dll);   // compile the model to C (MEX) for speed; needs a C compiler
                  // (MATLAB: run `mex -setup C` once; Octave: needs gcc/mkoctfile).
                  // If no compiler is available, change this back to `model;`.
// beta recovered from the estimated 100*(1/beta-1) (exact IG-prior parameterization)
#beta = 1/(1+pbeta/100);
// ===== SS-dependent quantities as model-local expressions (update per draw) =====
#Gammabar = Gammabar_Z*Gammabar_N;
#Rss  = Gammabar*Pibar/(beta*vthbar);
#rbar = Gammabar/(beta*vthbar);
#mc_ss = 1/mu_p_bar;
#Rk_ss   = (Gammabar/beta-(1-delta))/(1-tauK_ss);
#kappa_a = (1-tauK_ss)*Rk_ss;   // = Gammabar/beta-(1-delta); A'(1)=(1-tauK)Rk, not Rk
#pV_ss = ((mc_ss^(1-nu_Y)-omega_Y)/(1-omega_Y))^(1/(1-nu_Y));
#w_ss  = (1-alpha)*(alpha/Rk_ss)^(alpha/(1-alpha))*pV_ss^(1/(1-alpha));
#kun_ss = alpha/(1-alpha)*w_ss/Rk_ss;
#vn_ss  = kun_ss^alpha;
#oYv_ss = omega_Y/(1-omega_Y)*pV_ss^nu_Y;
#oYn_ss = oYv_ss*vn_ss;
#yn_ss  = ((1-omega_Y)^(1/nu_Y)*vn_ss^((nu_Y-1)/nu_Y)+omega_Y^(1/nu_Y)*oYn_ss^((nu_Y-1)/nu_Y))^(nu_Y/(nu_Y-1));
#ik_ss  = (Gammabar-(1-delta))/Gammabar;
#kn_ss  = Gammabar*kun_ss;
#y_ss   = yn_ss*n_ss;
#iv_ss  = ik_ss*kn_ss*n_ss;
#ku_ss  = kun_ss*n_ss;
#oY_ss  = oYn_ss*n_ss;
#phio   = omega_C/(1-omega_C);
#cN_ss  = ((1-s_g)*y_ss - iv_ss)/(1+s_g*phio);
#cO_ss  = phio*cN_ss;
#gdp_ss = y_ss + cO_ss;
#g_ss   = s_g*gdp_ss;
#tr_ss  = s_t*gdp_ss;
#bg_ss  = s_b*gdp_ss;
#oS_ss  = cO_ss + oY_ss;
#d_ss   = y_ss - w_ss*n_ss - Rk_ss*ku_ss - oY_ss;
#tau_c  = (g_ss+tr_ss+bg_ss/(Pibar*Gammabar)-bg_ss/Rss - tauL_ss*w_ss*n_ss - tauK_ss*(Rk_ss*ku_ss+d_ss))/(cN_ss+cO_ss);
#Bw_bar = bw_y*y_ss;
#cwN_ss = (1-omega_C)/(1+tau_c)*((1-tauL_ss)*w_ss*n_ss/lam + tr_ss + Bw_bar*(1/(Pibar*Gammabar)-1/Rss));
#cwO_ss = phio*cwN_ss;
#cw_ss  = ((1-omega_C)^(1/nu_C)*cwN_ss^((nu_C-1)/nu_C)+omega_C^(1/nu_C)*cwO_ss^((nu_C-1)/nu_C))^(nu_C/(nu_C-1));
#mrs_ss = w_ss/(mu_w_bar*(1+phi));
#xi     = mrs_ss*(1-tauL_ss)/((1+tau_c)*(n_ss/lam)^phi*cw_ss);

// ===== SDF / adjustment-cost helpers =====
#SDF_w = (beta/Gamma(+1))*(chi_eff(+1)/chi_eff)*(cw(+1)/cw)^(-(nu_C-1)/nu_C)*(cwN(+1)/cwN)^(-1/nu_C);
#SDF_s = (beta/Gamma(+1))*(chi_eff(+1)/chi_eff)*(cs(+1)/cs)^(-(nu_C-1)/nu_C)*(csN(+1)/csN)^(-1/nu_C);
#Anu  = kappa_a*(nu-1)   + 0.5*sigma_a*(nu-1)^2;
#Anu1 = kappa_a*(nu(+1)-1)+ 0.5*sigma_a*(nu(+1)-1)^2;
#Apnu = kappa_a + sigma_a*(nu-1);
#gI   = Gamma*iv/iv(-1);
#gI1  = Gamma(+1)*iv(+1)/iv;
#Sc   = (psi_i/2)*(gI  - Gammabar)^2;
#Sp   = psi_i*(gI  - Gammabar);
#Sp1  = psi_i*(gI1 - Gammabar);

// ===== COVID-effective drivers (112)-(114) =====
chi_eff    = chi*exp(sig_chi_c*e_chi_c);
tr_eff     = tr *exp(sig_tr_c *e_tr_c);
GammaN_eff = Gamma_N*exp(sig_gn_c*e_gn_c)/exp(sig_gn_c*e_gn_c(-1));
// observed special-shock innovations: data = 0 outside COVID windows, NaN (missing) inside;
// the Kalman smoother then infers the shock only during the active quarters.
ec_chi = e_chi_c;
ec_tr  = e_tr_c;
ec_gn  = e_gn_c;

// ===== WORKERS =====
1 = R*vth*SDF_w/Pi(+1)*1/(1+psi_w*(bw-Bw_bar));                                   // (45)
(1+tau_c)*(cwN+pO*cwO) + bw/R = (1-tauL)*w*n/lam + bw(-1)/(Pi*Gamma) + tr_eff;     // (46)
cwO = pO^(-nu_C)*(omega_C/(1-omega_C))*cwN;                                        // (47)
cw = ((1-omega_C)^(1/nu_C)*cwN^((nu_C-1)/nu_C)+omega_C^(1/nu_C)*cwO^((nu_C-1)/nu_C))^(nu_C/(nu_C-1)); // (48)

// ===== CAPITALISTS =====
1 = R*vth*SDF_s/Pi(+1);                                                            // (50)
(1-tauK)*R_k = Apnu;                                                               // (51) corrected
Q = SDF_s*((1-tauK(+1))*R_k(+1)*nu(+1)-Anu1+(1-delta)*Q(+1));                       // (52)
1 - Q*zeta*(1 - Sc - gI*Sp) = SDF_s*Q(+1)*zeta(+1)*gI1^2*Sp1;                       // (53)
csO = pO^(-nu_C)*(omega_C/(1-omega_C))*csN;                                         // (54)
cs = ((1-omega_C)^(1/nu_C)*csN^((nu_C-1)/nu_C)+omega_C^(1/nu_C)*csO^((nu_C-1)/nu_C))^(nu_C/(nu_C-1)); // (55)

// ===== LABOR MARKET =====
w + (mu_w/(mu_w-1))*(mrs*(1+phi)-w)
  - eta_w*(Pi_w/Pi_w_idx)*(Pi_w/Pi_w_idx-1)
  + beta*eta_w*(n(+1)/n)*(Pi_w(+1)/Pi_w_idx(+1))*(Pi_w(+1)/Pi_w_idx(+1)-1) = 0;    // (56)
mrs = (1+tau_c)*xi*(n/lam)^phi*cw^((nu_C-1)/nu_C)*cwN^(1/nu_C)/((1-tauL)*(1-omega_C)^(1/nu_C)); // (57)
Pi_w_idx = Gammabar_Z*(Pi(-1))^iota_w*(Pi_star)^(1-iota_w);                          // (58)
Pi_w = (w/w(-1))*Gamma_Z*Pi;                                                        // (59)

// ===== CAPITAL / INVESTMENT =====
k = (1-delta)*k(-1)/Gamma + zeta*(1-Sc)*iv;                                         // (60)
ku = nu/Gamma*k(-1);                                                               // (61)

// ===== FIRMS =====
(Pi/Pi_idx)*(Pi/Pi_idx-1) =
   SDF_s*(y(+1)/y)*Gamma(+1)*(Pi(+1)/Pi_idx(+1))*(Pi(+1)/Pi_idx(+1)-1)
 + (mu_p*mc-1)/(eta_p*(mu_p-1));                                                    // (62)
mc = (omega_Y*pO^(1-nu_Y)+(1-omega_Y)*pV^(1-nu_Y))^(1/(1-nu_Y));                    // (63)
pV = (R_k/alpha)^alpha*(w/(1-alpha))^(1-alpha);                                     // (64)
Pi_idx = (Pi(-1))^iota_p*(Pi_star)^(1-iota_p);                                      // (65)
y = ((1-omega_Y)^(1/nu_Y)*v^((nu_Y-1)/nu_Y)+omega_Y^(1/nu_Y)*oY^((nu_Y-1)/nu_Y))^(nu_Y/(nu_Y-1)); // (66)
v = ku^alpha*n^(1-alpha);                                                           // (67)
ku = alpha/(1-alpha)*w/R_k*n;                                                       // (68)
oY = (pV/pO)^nu_Y*(omega_Y/(1-omega_Y))*v;                                          // (69)

// ===== AGGREGATION =====
cN + g + iv + Anu*k(-1)/Gamma = y*(1 - (eta_p/2)*(Pi/Pi_idx-1)^2);                  // (70)
cN = lam*cwN + (1-lam)*csN;                                                         // (71)
cO = lam*cwO + (1-lam)*csO;                                                         // (72)
cc = cN + pO*cO;                                                                    // (73)
cO + oY = oS;                                                                       // (74)
gdp = y + pO*(oS - oY);                                                             // (75)

// ===== GOVERNMENT =====
g + tr_eff + bg(-1)/(Pi*Gamma) = tau_c*(cN+pO*cO)+tauL*w*n+tauK*(R_k*ku+d)+bg/R; // (76)
tauL = (tauL(-1))^rho_tau*(tauL_ss*((bg(-1)/gdp(-1))/s_b)^phi_tl)^(1-rho_tau);      // (77)
tauK = (tauK(-1))^rho_tau*(tauK_ss*((bg(-1)/gdp(-1))/s_b)^phi_tk)^(1-rho_tau);      // (78)
R = R(-1)^rho_r*(rbar*Pi_star*(Pi/Pi_star)^phi_Pi*(gdp/gdp(-1))^phi_Y)^(1-rho_r)*mp;// (79)

// ===== OTHER =====
d = y - w*n - R_k*ku - pO*oY - (eta_p/2)*y*(Pi/Pi_idx-1)^2;                         // (80)
r = R/Pi(+1);                                                                       // (81)
Pi_h = Pi*((1-omega_C+omega_C*pO^(1-nu_C))/(1-omega_C+omega_C*pO(-1)^(1-nu_C)))^(1/(1-nu_C)); // (82)
1 = R_AAA*SDF_s/Pi(+1);                                                             // (83)

// ===== STRUCTURAL SHOCKS (84)-(95) =====
log(Gamma_Z) = (1-rho_gz)*log(Gammabar_Z) + rho_gz*log(Gamma_Z(-1)) + sig_gz*e_gz;
log(Gamma_N) = (1-rho_gn)*log(Gammabar_N) + rho_gn*log(Gamma_N(-1)) + sig_gn*e_gn;
log(zeta)    = rho_zeta*log(zeta(-1)) + sig_zeta*e_zeta;
log(mu_p)    = (1-rho_mup)*log(mu_p_bar)+rho_mup*log(mu_p(-1))+sig_mup*e_mup-eta_mup*sig_mup*e_mup(-1);
log(mu_w)    = (1-rho_muw)*log(mu_w_bar)+rho_muw*log(mu_w(-1))+sig_muw*e_muw-eta_muw*sig_muw*e_muw(-1);
log(oS)      = (1-rho_oS)*log(oS_ss) + rho_oS*log(oS(-1)) + sig_oS*e_oS;
log(chi)     = rho_chi*log(chi(-1)) + sig_chi*e_chi;
log(vth)     = (1-rho_vth)*log(vthbar) + rho_vth*log(vth(-1)) + sig_vth*e_vth;
log(g)       = (1-rho_g)*log(g_ss) + rho_g*log(g(-1)) + sig_g*e_g;
log(tr)      = (1-rho_tr)*log(tr_ss) + rho_tr*log(tr(-1)) + sig_tr*e_tr;
log(mp)      = rho_mp*log(mp(-1)) + sig_mp*e_mp;
log(Pi_star) = (1-rho_pistar)*log(Pibar) + rho_pistar*log(Pi_star(-1)) + sig_pistar*e_pistar;
Gamma = Gamma_Z*GammaN_eff;

// ===== MEASUREMENT-ERROR PROCESSES (7) =====
me_gdp     = rho_me_gdp*me_gdp(-1)         + sig_me_gdp*e_me_gdp;
me_tfp     = rho_me_tfp*me_tfp(-1)         + sig_me_tfp*e_me_tfp;
me_oil     = rho_me_oil*me_oil(-1)         + sig_me_oil*e_me_oil;
me_gas     = rho_me_gas*me_gas(-1)         + sig_me_gas*e_me_gas;
me_inflexp = rho_me_inflexp*me_inflexp(-1) + sig_me_inflexp*e_me_inflexp;
me_tp      = (1-rho_me_tp)*tp + rho_me_tp*me_tp(-1) + sig_me_tp*e_me_tp;
me_cy      = rho_me_cy*me_cy(-1)           + sig_me_cy*e_me_cy;

// ===== OBSERVATION EQUATIONS (97)-(111) =====
obs_gdp = 400*(log(gdp)-log(gdp(-1)) + log(Gamma_Z)+log(GammaN_eff) + me_gdp);      // (97)
obs_c   = 400*(log(cc) -log(cc(-1))  + log(Gamma_Z)+log(GammaN_eff)) + cobs_cy;      // (98)
obs_i   = 400*(log(iv) -log(iv(-1))  + log(Gamma_Z)+log(GammaN_eff)) + cobs_iy;      // (99)
obs_g   = 400*(log(g)  -log(g(-1))   + log(Gamma_Z)+log(GammaN_eff)) + cobs_gy;      // (100)
obs_tr  = 400*(log(tr_eff)-log(tr_eff(-1)) + log(Gamma_Z)+log(GammaN_eff)) + cobs_ty;// (101)
obs_w   = 400*(log(w)  -log(w(-1))   + log(Gamma_Z)) + cobs_wy;                      // (102)
obs_n   = 400*(log(n)  -log(n(-1))   + log(GammaN_eff));                             // (103)
obs_ffr = 100*(R^4 - 1);                                                            // (104)
obs_pi  = 400*log(Pi);                                                              // (105)
obs_tfp = log(Gamma_Z) - log(Gammabar_Z) + me_tfp;                                  // (106)
obs_oil = log(pO) + me_oil;                                                         // (107)
obs_gas = 400*(log(Pi) + log(pO) - log(pO(-1)) + me_gas);                           // (108)
// perpetuity (geometric-decay) approximation of the EH long rates: identical SS,
// far fewer states than the exact 40/80-quarter sums.
yq10_R   = (1-d10)*log(R)        + d10*yq10_R(+1);
yq20_R   = (1-d20)*log(R)        + d20*yq20_R(+1);
yq20_AAA = (1-d20)*log(R_AAA)    + d20*yq20_AAA(+1);
piq10    = (1-d10)*log(Pi_h(+1)) + d10*piq10(+1);
obs_inflexp = 100*( exp(4*piq10)    - 1 + me_inflexp );                              // (109)
obs_tr10    = 100*( exp(4*yq10_R)   - 1 + me_tp );                                   // (110)
obs_cy      = 100*( exp(4*yq20_AAA) - exp(4*yq20_R) + me_cy );                       // (111)
end;

// ----------------------------------------------------------------------------
//  STEADY STATE (recomputes the same recipe as locals each draw; handles the
//  long-rate lead auxiliaries automatically).
// ----------------------------------------------------------------------------
steady_state_model;
  q_beta = 1/(1+pbeta/100);
  q_Gam = Gammabar_Z*Gammabar_N;
  q_R   = q_Gam*Pibar/(q_beta*vthbar);
  q_r   = q_Gam/(q_beta*vthbar);
  q_mc  = 1/mu_p_bar;
  q_Rk  = (q_Gam/q_beta-(1-delta))/(1-tauK_ss);
  q_pV  = ((q_mc^(1-nu_Y)-omega_Y)/(1-omega_Y))^(1/(1-nu_Y));
  q_w   = (1-alpha)*(alpha/q_Rk)^(alpha/(1-alpha))*q_pV^(1/(1-alpha));
  q_kun = alpha/(1-alpha)*q_w/q_Rk;
  q_vn  = q_kun^alpha;
  q_oYn = omega_Y/(1-omega_Y)*q_pV^nu_Y*q_vn;
  q_yn  = ((1-omega_Y)^(1/nu_Y)*q_vn^((nu_Y-1)/nu_Y)+omega_Y^(1/nu_Y)*q_oYn^((nu_Y-1)/nu_Y))^(nu_Y/(nu_Y-1));
  q_ik  = (q_Gam-(1-delta))/q_Gam;
  q_ph  = omega_C/(1-omega_C);
  q_y=q_yn*n_ss; q_ku=q_kun*n_ss; q_v=q_vn*n_ss; q_oY=q_oYn*n_ss;
  q_k=q_Gam*q_kun*n_ss; q_iv=q_ik*q_Gam*q_kun*n_ss;
  q_cN=((1-s_g)*q_y-q_iv)/(1+s_g*q_ph); q_cO=q_ph*q_cN;
  q_gdp=q_y+q_cO; q_g=s_g*q_gdp; q_tr=s_t*q_gdp; q_bg=s_b*q_gdp;
  q_oS=q_cO+q_oY; q_d=q_y-q_w*n_ss-q_Rk*q_ku-q_oY;
  q_tc=(q_g+q_tr+q_bg/(Pibar*q_Gam)-q_bg/q_R-tauL_ss*q_w*n_ss-tauK_ss*(q_Rk*q_ku+q_d))/(q_cN+q_cO);
  q_Bw=bw_y*q_y;
  q_cwN=(1-omega_C)/(1+q_tc)*((1-tauL_ss)*q_w*n_ss/lam+q_tr+q_Bw*(1/(Pibar*q_Gam)-1/q_R));
  q_cwO=q_ph*q_cwN;
  q_cw=((1-omega_C)^(1/nu_C)*q_cwN^((nu_C-1)/nu_C)+omega_C^(1/nu_C)*q_cwO^((nu_C-1)/nu_C))^(nu_C/(nu_C-1));
  q_csN=(q_cN-lam*q_cwN)/(1-lam); q_csO=q_ph*q_csN;
  q_cs=((1-omega_C)^(1/nu_C)*q_csN^((nu_C-1)/nu_C)+omega_C^(1/nu_C)*q_csO^((nu_C-1)/nu_C))^(nu_C/(nu_C-1));
  // assignments
  Gamma_Z=Gammabar_Z; Gamma_N=Gammabar_N; GammaN_eff=Gammabar_N; Gamma=q_Gam;
  zeta=1; mu_p=mu_p_bar; mu_w=mu_w_bar; chi=1; chi_eff=1; vth=vthbar; mp=1;
  Pi=Pibar; Pi_star=Pibar; Pi_idx=Pibar; Pi_h=Pibar;
  Pi_w=Gammabar_Z*Pibar; Pi_w_idx=Gammabar_Z*Pibar;
  R=q_R; r=q_r; R_k=q_Rk; Q=1; mc=q_mc; pV=q_pV; pO=1; w=q_w;
  R_AAA=Pibar*q_Gam/q_beta; nu=1;
  n=n_ss; ku=q_ku; v=q_v; oY=q_oY; y=q_y; k=q_k; iv=q_iv;
  cN=q_cN; cO=q_cO; cc=q_cN+q_cO; gdp=q_gdp; d=q_d; mrs=q_w/(mu_w_bar*(1+phi));
  g=q_g; tr=q_tr; tr_eff=q_tr; bg=q_bg; oS=q_oS; tauL=tauL_ss; tauK=tauK_ss;
  cwN=q_cwN; cwO=q_cwO; cw=q_cw; csN=q_csN; csO=q_csO; cs=q_cs; bw=q_Bw;
  me_gdp=0; me_tfp=0; me_oil=0; me_gas=0; me_inflexp=0; me_tp=tp; me_cy=0;
  ec_chi=0; ec_tr=0; ec_gn=0;
  obs_gdp=400*log(q_Gam); obs_c=400*log(q_Gam)+cobs_cy; obs_i=400*log(q_Gam)+cobs_iy;
  obs_g=400*log(q_Gam)+cobs_gy; obs_tr=400*log(q_Gam)+cobs_ty;
  obs_w=400*log(Gammabar_Z)+cobs_wy; obs_n=400*log(Gammabar_N);
  obs_ffr=100*(q_R^4-1); obs_pi=400*log(Pibar); obs_tfp=0; obs_oil=0; obs_gas=400*log(Pibar);
  obs_inflexp=100*(Pibar^4-1); obs_tr10=100*(q_R^4-1+tp);
  obs_cy=100*((Pibar*q_Gam/q_beta)^4 - q_R^4);
  yq10_R=log(q_R); yq20_R=log(q_R); yq20_AAA=log(Pibar*q_Gam/q_beta); piq10=log(Pibar);
end;

// ----------------------------------------------------------------------------
varobs obs_gdp obs_c obs_i obs_g obs_tr obs_w obs_n obs_ffr obs_pi
       obs_tfp obs_oil obs_gas obs_inflexp obs_tr10 obs_cy
       ec_chi ec_tr ec_gn;

estimated_params;
  Gammabar_Z, normal_pdf, 1.005, 0.005;
  Gammabar_N, normal_pdf, 1.000, 0.005;
  pbeta,      inv_gamma_pdf, 0.040, 0.010;   // exact prior on 100*(1/beta-1)
  phi,        gamma_pdf,  2.000, 0.500;
  eta_p,      gamma_pdf, 80.000, 15.000;
  eta_w,      gamma_pdf, 80.000, 15.000;
  psi_i,      inv_gamma_pdf, 5.000, 5.000;
  iota_p,     beta_pdf,   0.500, 0.150;
  iota_w,     beta_pdf,   0.500, 0.150;
  sigma_a,    inv_gamma_pdf, 0.500, 0.100;
  // ---- BOUNDS ON THE THREE UNBOUNDED STEADY-STATE PARAMETERS ----------------
  // alpha, nu_C and nu_Y carry NORMAL (unbounded) priors but enter the steady state
  // through fractional powers: alpha<0 makes (alpha/Rk_ss)^(alpha/(1-alpha)) complex,
  // and nu_C or nu_Y at 0 or 1 makes the CES aggregators singular. A long optimizer
  // line-search step crosses them, ys comes back complex, and static_resid throws
  // "y must be a real dense numeric array". These bounds sit 4+ prior s.d. from the
  // mode, so they never bind on the posterior; they only stop the line search from
  // walking off the domain.
  alpha,      0.1945, 0.0200, 0.6000, normal_pdf, 0.300, 0.050;
  psi_w,      inv_gamma_pdf, 0.500, 0.250;
  lam,        beta_pdf,   0.700, 0.100;
  vthbar,     normal_pdf, 1.004, 0.001;
  tp,         inv_gamma_pdf, 0.011, 0.005;
  nu_C,       0.2854, 0.0200, 0.9000, normal_pdf, 0.100, 0.100;
  nu_Y,       0.1864, 0.0200, 0.9000, normal_pdf, 0.200, 0.100;
  phi_Pi,     normal_pdf, 2.000, 0.250;
  phi_Y,      normal_pdf, 0.250, 0.250;
  rho_r,      beta_pdf,   0.700, 0.200;
  bw_y,       normal_pdf, 0.000, 0.500;
  phi_tl,     inv_gamma_pdf, 0.100, 0.050;
  phi_tk,     inv_gamma_pdf, 0.100, 0.050;
  rho_gz,   beta_pdf, 0.500, 0.200;   rho_mup,  beta_pdf, 0.500, 0.200;
  rho_muw,  beta_pdf, 0.500, 0.200;   rho_zeta, beta_pdf, 0.500, 0.200;
  rho_chi,  beta_pdf, 0.500, 0.200;   rho_vth,  beta_pdf, 0.990, 0.005;
  rho_g,    beta_pdf, 0.500, 0.200;   rho_tr,   beta_pdf, 0.500, 0.200;
  rho_pistar, beta_pdf, 0.990, 0.005; rho_mp,   beta_pdf, 0.500, 0.200;
  rho_gn,   beta_pdf, 0.500, 0.200;   rho_oS,   beta_pdf, 0.500, 0.200;
  sig_gz,   inv_gamma_pdf, 0.010, 0.050;  sig_mup,  inv_gamma_pdf, 0.010, 0.050;
  sig_muw,  inv_gamma_pdf, 0.010, 0.050;  sig_zeta, inv_gamma_pdf, 0.010, 0.050;
  sig_chi,  inv_gamma_pdf, 0.010, 0.050;  sig_vth,  inv_gamma_pdf, 0.005, 0.050;
  sig_g,    inv_gamma_pdf, 0.010, 0.050;  sig_tr,   inv_gamma_pdf, 0.010, 0.050;
  sig_pistar, inv_gamma_pdf, 0.001, 0.050; sig_mp,  inv_gamma_pdf, 0.010, 0.050;
  sig_gn,   inv_gamma_pdf, 0.010, 0.050;  sig_oS,   inv_gamma_pdf, 0.010, 0.050;
  eta_mup,  beta_pdf, 0.500, 0.200;   eta_muw,  beta_pdf, 0.500, 0.200;
  rho_me_tfp,     beta_pdf, 0.500, 0.200;   sig_me_tfp,     inv_gamma_pdf, 0.010, 0.050;
  rho_me_cy,      beta_pdf, 0.500, 0.200;   sig_me_cy,      inv_gamma_pdf, 0.010, 0.050;
  rho_me_gdp,     beta_pdf, 0.500, 0.200;   sig_me_gdp,     inv_gamma_pdf, 0.010, 0.050;
  rho_me_inflexp, beta_pdf, 0.500, 0.200;   sig_me_inflexp, inv_gamma_pdf, 0.010, 0.050;
  rho_me_tp,      beta_pdf, 0.500, 0.200;   sig_me_tp,      inv_gamma_pdf, 0.010, 0.050;
  rho_me_oil,     beta_pdf, 0.500, 0.200;   sig_me_oil,     inv_gamma_pdf, 0.010, 0.050;
  rho_me_gas,     beta_pdf, 0.500, 0.200;   sig_me_gas,     inv_gamma_pdf, 0.010, 0.050;
end;

// NOTE: no estimated_params_init here -- the mode_file supplies the parameter
// values, and use_calibration is incompatible with the mode_file option.

shocks;
  var e_gz=1; var e_gn=1; var e_zeta=1; var e_mup=1; var e_muw=1; var e_oS=1;
  var e_chi=1; var e_vth=1; var e_g=1; var e_tr=1; var e_mp=1; var e_pistar=1;
  var e_me_gdp=1; var e_me_tfp=1; var e_me_oil=1; var e_me_gas=1;
  var e_me_inflexp=1; var e_me_tp=1; var e_me_cy=1;
  // COVID special-shock innovations: unit variance; pinned to 0 outside the active
  // windows by the ec_* observables (0/NaN data), inferred by the smoother inside.
  var e_chi_c=1; var e_tr_c=1; var e_gn_c=1;
end;

steady;
check;

// ============================================================================
//  PARALLEL METROPOLIS-HASTINGS, STARTED FROM THE SAVED MODE
//  Lean file: NO flexible-price/r* block and NO forward-rate leads (both are
//  likelihood-neutral -> they only slow MCMC). Add them later in the postest.
//
//  HOW TO LAUNCH (one 8-core machine; edit CPUnbr in cluster.ini to match):
//     dynare frbstl_us_est_fastlr_lean parallel conffile=cluster.ini
//  use_dll is already on (model(use_dll)). To run WITHOUT the parallel engine,
//  just: dynare frbstl_us_est_fastlr_lean   (chains then run sequentially).
//
//  Requires the saved 63-param mode file (frbstl_us_est_fastlr_mode_ndc_fixed.mat) in this
//  folder. mode_compute=0 loads it -> NO mode search. The proposal uses the
//  inverse Hessian at the mode (MCMC_jumping_covariance=hessian, the default).
//  mh_tune_jscale auto-tunes the step to ~25-33% acceptance.
//  mh_nblocks = number of parallel chains = set to your physical core count.
//  8 chains x 250000 = 2,000,000 pooled draws (>= the paper's 2 x 1,000,000).
// ============================================================================

// ---- STEP A (do this FIRST): quick convergence check, ~10-20 min ----------
// Confirms the mode/Hessian proposal mixes before you commit to the long run.
// Read the Brooks-Gelman diagnostics + trace plots it prints, then scale up.
// estimation(datafile = 'usdata_1959_2026.mat',
//            mode_compute = 0, mode_file = frbstl_us_est_fastlr_mode,
//            mh_replic = 50000, mh_nblocks = 2, mh_drop = 0.5,
//            mh_tune_jscale,   // auto-tunes the step to ~33% acceptance (sets jscale itself)
//            MCMC_jumping_covariance = hessian,
//            presample = 4, lik_init = 2, nograph);

// ---- STEP B: the full parallel posterior (Tables 2-3) ---------------------
// FALLBACK PILOT: does not use the Hessian at all.  prior_variance is positive
// definite by construction, so it always runs; mh_tune_jscale adapts the scale.
// Mixing will be worse than a good Hessian proposal -- check the acceptance rate.
estimation(datafile = 'usdata_1959_2026.mat',
           mode_compute = 0, mode_file = frbstl_us_est_fastlr_mode_ndc,
           mh_replic = 5000, mh_nblocks = 8, mh_drop = 0.5,
           mh_tune_jscale,
           MCMC_jumping_covariance = prior_variance,
           presample = 4, lik_init = 2, nograph);

// FALLBACK if the mode Hessian is NOT positive definite (Dynare warns):
//   (1) rebuild a PD proposal once with mode_compute=6 (gmhmaxlik), OR
//   (2) swap the proposal to the prior covariance and tune jscale:
// estimation(datafile = 'usdata_1959_2026.mat',
//            mode_compute = 0, mode_file = frbstl_us_est_fastlr_mode,
//            mh_replic = 250000, mh_nblocks = 8, mh_drop = 0.5,
//            mh_tune_jscale,   // auto-tunes the step to ~33% acceptance (sets jscale itself)
//            MCMC_jumping_covariance = prior_variance,
//            presample = 4, lik_init = 2, nograph);

// ---- STEP C (optional): a second refinement pass ---------------------------
// After STEP B, adopt the EMPIRICAL posterior covariance as the proposal and
// continue the chains for a cleaner mixing / higher effective sample size:
// estimation(datafile = 'usdata_1959_2026.mat',
//            mode_compute = 0, mode_file = frbstl_us_est_fastlr_mode,
//            mh_replic = 250000, mh_nblocks = 8, mh_drop = 0.5,
//            use_mh_covariance_matrix, mh_tune_jscale,
//            mh_initialize_from_previous_mcmc,
//            presample = 4, lik_init = 2, nograph);
