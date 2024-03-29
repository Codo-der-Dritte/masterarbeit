# Title: Model Zezza

# Purpose : This is the model developed in the paper called
# U.S. growth, the housing market, and the distribution of
# income. 
#
# Equations that appear are in the paper are includes with the
# same number they have in the paper. They are marked as such [xx].
# Equations that do not appear in the paper but are present in
# Gennaros code are marked as such [xxCx]. They belong to
# a equation that is in the paper they have the same number and
# only a C as a suffix. If they only appear in the code and are
# not helper equations they are marked as such [00xxC].

# Author: Alexander Toplitsch
# Contact details: alexander.toplitsch@s.wu.ac.at

# Date script created: Fr 02.12.2023 21:03 -------------
# Date script last modified:  

# Check necessary packages for installation and call them

libs <- c("tidyverse", "leaflet", "lubridate", "readxl", "sfcr", "scales", "ggraph", "networkD3")
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
  install.packages(libs[!installed_libs])
}
invisible(lapply(libs, library, character.only = T)) 

#############################################################
##########################Model##############################
#############################################################

# Set up the model equations
model_eqs <- sfcr_set(
  #----------------------------------#
  # Capitalists ("rich households")
  #----------------------------------#
  # [1] - Disposable income - Sum of wage income, distributed profits from firms and banks,
  # interest income from deposits and treasuries, and rents, net of direct taxes
  # paid to the government.
  Yc ~ wc * Nc + Rents + FD + rm[-1] * Mc[-1] + FB + rb[-1] * Bh[-1] - Tdc,
  # [2] Saving - augments the stock of wealth.
  Shc ~ Yc - Cc, 
  # [3] Consumption inflator.
  Cc ~ cc * p,
  # [4] Consumption - depend on expected disposable income, the open stock of wealth  and
  # expected capital gains, all measured in real terms.
  cc ~ alpha_1c * yc_e + alpha_2c * vc[-1] + alpha_3c * 
    (cge_e + cghc_e - p_e * vc[-1]/(1 + p_e)),
  # [5] Stock of Wealth - includes saving and capitals gains which are given from changes in
  # market price for equities [7] and homes [8].
  Vc ~ Vc[-1] + Shc +  CGE + CGHc,
  # [5C1] Wealth deflator.
  vc ~ Vc / p,
  # [6] Disposable income deflator.
  yc ~ Yc / p,
  # [7] Changes in market price for equities -  Capital gains.
  CGE ~ (pe - pe[-1]) * E[-1],
  # [8] Changes in market price for homes - Capital gains.
  CGHc ~ (ph - ph[-1]) * Hc[-1],
  #-----Portfolio Choice-----#
  # [9] Cash - depend on current consumption.
  HPc ~ eta * Cc,
  # [10] Bank deposits.
  Mc ~ Vc - HPc - Bh - E * pe - ph * Hc,
  # [11] Bonds.
  Bh ~ (Vc_e - HPc) * 
    (lambda_10 - lambda_11 * rrm  - lambda_12 * 
       rre_e - lambda_13 * (Yc_e/Vc_e) + lambda_14 *
       rrb - lambda_15 * rrh_e),
  # [12] Price of Equities.
  pe ~ ((Vc_e - HPc) * 
         (lambda_20 - lambda_21 * rrm + lambda_22 * 
            rre_e - lambda_23 * (Yc_e/Vc_e) - lambda_24 *
            rrb - lambda_25 * rrh_e) - xi * (I-FU)
  ) / E[-1],
  # [13] Homes.
  Hc ~ ((Vc_e - HPc) * 
          (lambda_30 - lambda_31 * rrm - lambda_32 * 
             rre_e - lambda_33 * (Yc_e/Vc_e) - lambda_34 *
             rrb + lambda_35 * rrh_e)
  ) / ph,
  # [14] Return on equities - Given by distributed profits and expected capital gains.
  re ~ (FD + CGE_e)/(pe[-1] * E[-1]),
  # [15] Return on Housing - Given by rents and expected capital gains.
  rh ~ (Rents + CGHc_e)/(ph[-1] * Hc[-1]),
  #----------------------------------#
  # Workers ("other households")
  #----------------------------------#
  # [16] Disposable income - Sum of wage income, rent to pay,
  # interest payments from bank deposits, interest payments for mortgages,
  # and rents, net of direct taxes paid to the government.
  Yo ~ wo * No - Rents + rm[-1] * Mo[-1] - rmo[-1] * MOo[-1] - Tdo,
  # [17] Saving
  Sho ~ Yo - Co,
  # [18] Consumption inflator.
  Co ~ co * p,
  # [19] Consumption - depends on expected real disposable income, past 
  # real wealth, and expected real capital gains on homes minus past
  # wealth normalized by expected inflation.
  co ~ alpha_0o + alpha_1o * yo_e + alpha_2o * vo[-1] + 
    alpha_3o * (cgho_e - p_e * vo[-1]/(1 + p_e)) +
    iec - alpha_4o * morp * MOo[-1]/Yo[-1],
  # [20] Wealth - Past wealth plus saving plus capital gains from homes.
  Vo ~ Vo[-1] + Sho + CGHo,
  # [20C1] Wealth deflator.
  vo ~ Vo / p,
  # [21] Disposable income deflator.
  yo ~ Yo/p,
  # [22] Capital gains from homes.
  CGHo ~ (ph - ph[-1]) * Ho[-1],
  # [23] Imitation paramter.
  iec ~ im * No * (cc[-1]/Nc - co[-1]/No),
  # # [MEMO] Real income of households combined.
  # yc_yo ~ yc + yo,
  # # [MEMO] Nominal income of households combined.
  # Yc_Yo ~ yc_yo * p,
  # # [MEMO] Nominal consumption of households combined.
  # C ~ Cc + Co,
  #-----Portfolio Choice-----#
  # [24] Cash - depends on current consumption.
  HPo ~ eta * Co,
  # [25] Bank deposits - residual.
  Mo ~ Vo - HPo - ph * Ho + MOo,
  # [26] Demand for Homes - depends on population growth,
  # expected real income and lagged debt repayment ratios.
  Ho_d ~ (((No - No[-1])/No[-1] + mu_1 * ((yo_e - yo_e[-1])/yo_e) -
    mu_2 * delta_debt_rep[-1]) * Ho[-1]),
  # [26C1] Number of Homes.
  Ho ~ Ho[-1] + Ho_d,
  # [0026C1] Debt repayment ratio.
  debt_rep ~ (rmo[-1] + morp[-1]) * MOo_1[-1] / Yo[-1],
  # [0026C2] Change in debt repayment.
  delta_debt_rep ~ debt_rep - debt_rep[-1],
  # [0026C3] Mo lag1
  MOo_1 ~ MOo[-1],
  # [27] Change Mortgages.
  MOo ~ (MOo[-1] * (1 - morp) + 
           MOo_cond * (ph * (Ho - Ho[-1]) - Sho)),
  # [0027C1] Condition change Mortgages.
  MOo_cond ~ ifelse(ph * (Ho - Ho[-1]) - Sho > 0, 1, 0),
  # [28] Share of rented homes owned by capitalist.
  Rents ~ rent * Hc[-1],
  # [29] Rent increases
  rent ~ rent[-1] * (1 + y_e),
  #----------------------------------#
  # Nonfinancial firms
  #----------------------------------#
  # [30] Investment decision - Growth of k depends on actual profits, Tobin's q,
  # borrowing costs from banks and utilization rate. This the the growth rate.
  kgr ~ (iota_0 + iota_1 * FU[-1]/K_1[-1] - 
          iota_2 * rll[-1] * lev[-1] +
          iota_3 * q[-1] +
          iota_4 * (ut[-1] - unorm)),
  # [30C1] Leverage.
  lev ~ L/K,
  # [30C2] Tobin's Q.
  q ~ pe * E/K,
  # [30C3] Real Capital.
  k ~ (1 + kgr) * k[-1],
  # [30C4] Capital in non real terms.
  K ~ k * p,
  # [30C5] Real Investment,
  i ~ k - k[-1],
  # [30C6] Nominal Investment.
  I ~ i * p,
  # [30C7] Real rate on loans.
  rll ~ (1 + rl) / (1 + p_e)-1,
  # [0030C1] Capital from one periods ago.
  K_1 ~ K[-1],
  # [31] Prices - are set with a mark-up on wages.
  p ~ (1 + rho) * w / (prod * (1 - tau)),
  # [32] Total Profits - are determined relative to the wage bill.
  FT ~ rho * WB,
  # [33] Distributed Profits - fixed share net of taxes and interest payments
  # is distributed to capitalists.
  FD ~ (1 - beta) * (FT - rl[-1] * L[-1] - TF),
  # [34] Mark-up - depends on relative strength of workers and capitalists.
  rho ~ ((rho_1 * (prodg - wo_g) + 1)*(1 + rho[-1])) - 1,
  # [35] Utilization rate - ratio of real sales to "normal" sales, which in
  # turn are in a fixed ratio (lambda) to the stock of real capital.
  ut ~ s/sfc,
  # [35C1] Normal Sales - fixed ratio of real capital.
  sfc ~ (lambda * k[-1]),
  # [36] Retained Profits.
  FU ~ FT - rl[-1] * L[-1] - FD - TF,
  # [37] New equities issued.
  E ~ E[-1] + xi * (I - FU) / pe,
  # [38] Loans - rewritten to display capital.
  L ~ L[-1] + I - FU - pe * (E - E[-1]),
  # # [Memo] Firms wealth.
  Vf ~ p * K + HU * ph - L - pe * E,
  #----------------------------------#
  # Banks
  #----------------------------------#
  # [39] Demand for Bonds.
  Bb ~ chi_1 * M,
  # [40] Reserve requirement.
  HPb ~ chi_2 * M,
  # [39C1] Deposits.
  M ~ Mo + Mc,
  # [41] Advancements from Central Bank - if internal funds are
  # not sufficient to cover for demand for loans the banks gets
  # advances from the central bank.
  A ~ L + HPb + Bb + MOo - M,
  # [42] Interest rates on loans.
  rl ~ ra + spread_1,
  # [43] Interest rates on mortgages.
  rmo ~ ra + spread_2,
  # [44] Interest rates on deposits.
  rm ~ ra + spread_3,
  # [45] Banks profits, all distributed.
  FB ~ rl[-1] * L[-1] + rb[-1] * Bb[-1] + rmo[-1] * MOo[-1] -
    (rm[-1] * M[-1] + ra[-1] * A[-1]),
  #----------------------------------#
  # Central bank
  #----------------------------------#
  # [46] Accommodates demand for advances and buys bonds that are
  # not absorbed by Households and Banks.
  Bcb ~ B - Bh - Bb,
  # [47] Interest income is redistributed to the government.
  FC ~ ra[-1] * A[-1] + rb[-1] * Bcb[-1],
  # [47C1] Rate on Advances.
  ra ~ ( 1 + rra) * (1 + p_e) - 1,
  # [MEMO] Total High powered money.
  HP ~ HPb + HPo + HPc,
  #----------------------------------#
  # Government
  #----------------------------------#
  # [48] Deficit -  Collects taxes production, wages, and profits,
  # and any deficit is financed by issuing bonds.
  GD ~ (G + rb[-1] * Bh[-1] + rb[-1] * Bcb[-1] + rb[-1] * Bb[-1]) - 
    (IT + Td + TF + FC),
  # [48C1] Rate on bonds.
  rb ~ (1 + rrb) * (1 + p_e)-1,
  # [49] Taxes on production.
  IT ~ tau * S,
  # [50] Tax on capitalists.
  Tdc ~ tau_d * wc * Nc,
  # [51] Tax on workers.
  Tdo ~ tau_d * wo * No,
  # [52] Total tax from wages.
  Td ~ Tdc + Tdo,
  # [53] Tax on profits.
  TF ~ tau_f * FT,
  # [54] New bonds issued.
  B ~ B[-1] + GD,
  # [55] Government spending inflator.
  G ~ g * p,
  # [56] Real Government spending.
  g ~ g[-1] * (1 + y_e),
  # [MEMO] Taxes.
  TT ~ IT + Td + TF,
  #----------------------------------#
  # The housing market
  #----------------------------------#
  # Demand for houses is laid down for capitalists and
  # workers in equations [13] and [26].
  
  # [57] Number of unsold homes - changes when the number of
  # newly built homes exceed the demand for homes. Can't be negative.
  HU ~ HU_cond * ((HU[-1] + HN - (Hc - Hc[-1])-(Ho - Ho[-1]))) + 0,
  # [0057C1] Condition helper.
  HU_cond ~ ifelse(((HU[-1] + HN - (Hc - Hc[-1])-(Ho - Ho[-1]))) > 0, 1, 0),
  # [0057C2] HU lag 1.
  HU_1 ~ HU[-1],
  # [58] Check for HND.
  HN ~ HN_cond * HND + 0,
  # [0058C1] Condition for HND.
  HN_cond ~ ifelse(HND - HU[-1] > 0, 1, 0) ,
  # [58C1] Supply of new homes - is a function of expected demand 
  # and past capital gains.
  HND ~ HND_cond * (v_1 * (Hc[-1] * y_e + (Ho - Ho[-1])) +
                      v_2 * (ph_3 - ph_3[-1])) + 0,
  # [0058C1] Condition for Supply of new homes.
  HND_cond ~ ifelse(v_1 * (Hc[-1] * y_e + (Ho - Ho[-1])) +
                v_2 * (ph_3 - ph_3[-1]) > 0, 1, 0),
  # [0058C2] ph lag 1,
  ph_1 ~ ph[-1],
  # [0058C3] ph lag 2.
  ph_2 ~ ph_1[-1],
  # [0058C4] ph lag 3.
  ph_3 ~ ph_2[-1],
  # [59] Market price of Homes.
  ph ~ ph_cond1 * ph1 + ph_cond2 * ph2,
  # [0059C1] Condition for market price of Homes.
  ph_cond1 ~ ifelse(HU < 10, 1, 0),
  # [0059C2] Condition for market price of Homes.
  ph_cond2 ~ ifelse(HU >= 10, 1, 0),
  # [59C1] Total supply of Houses.
  HNS ~ HN + HU[-1],
  # [59C2] Helper for House prices ph2.
  ph2 ~ ph[-1] - v_3 * (HU_1 - HU_1[-1]),
  # [59C3] Helper for House prices ph1.
  ph1 ~ ph[-1] + ph[-1] * v_3 * (Hc[-1] * y_e + (Ho_1 - Ho_1[-1]) - HNS[-1]),
  # [59C4] Lag 1 of Ho
  Ho_1 ~ Ho[-1],
  
  #----------------------------------#
  # Aggregate demand, un-/employment, wages
  #----------------------------------#
  # [60] Aggregate Demand -  Sales
  s ~ c + i + ih + g,
  # [60C1] Helper for Aggregate Demand - total consumption.
  c ~ cc + co,
  # [60C2] Real value of houses.
  ih ~ HN * p,
  # [60C3] Nominal value of houses
  IH ~ ih * p,
  # [61] Aggregate Demand inflator.
  S ~ s * p,
  # [62] Number of workers
  N ~ s/prod,
  # [63] Share of capitalists
  Nc ~ omega_c * N,
  # [64] Share of workers
  No ~ N - Nc,
  # [65] growth in real income
  y ~ s/s[-1]-1,
  # [66] Unemployment rate - Follows some sort of Okun's laws
  ur ~ ur_cond * ((ur[-1]) - (y - ny) / okun),
  # [0066C1] Condition Unemployment rate
  ur_cond ~ ifelse((ur[-1] - (y - ny) / okun) > 0, 1, 0),
  # [66C1] helper for ur.
  ur1 ~ ur1_cond1 * 0 + ur1_cond2 * 0.2 + ur1_cond3 * ur,
  # [0066C1] helper for ur.
  ur1_cond1 ~ ifelse(ur <= 0, 1, 0),
  # [0066C2] helper for ur.
  ur1_cond2 ~ ifelse(ur >= 0.2, 1, 0),
  # [0066C3] helper for ur.
  ur1_cond3 ~ ifelse(ur > 0 & ur < 0.2, 1,0),
  
  
  # [67] Wage Bill - All wages paid out.
  WB ~ w * N,
  # [67C1] Wage.
  w ~ (wc * omega_c + wo * (1 - omega_c)),
  # [68] Wage for capitalists.
  wc ~ wc[-1] * (1 + (wc_g)),
  # [69] Wage for workers.
  wo ~ wo[-1] * (1 + (wo_g)),
  # [70] Wage growth/inflation for capitalists.
  wc_g ~ p_e + omega * prodg_e + shockwc_g,
  # [71] Wage growth/inflation for workers.
  wo_g ~ p_e + omega * prodg_e + shockwo_g,
  # [72] Omega is the wage share,
  omega ~ o_0 - o_2 * sqrt(ur1/o_1),
  # [73] Productivity gains.
  prodg ~ pi_0 - pi_1 * dut + shockprodg,
  # [73C1] Production of Capitalists.
  prodc ~ prodc[-1] * (1 + prodg),
  # [73C1] Production of Workers.
  prodo ~ prodo[-1] * (1 + prodg),
  # [73C3] Total production.
  prod ~ omega_c * prodc + (1 - omega_c) * prodo,
  # [73C4] Change in utilization rate.
  dut ~ dut_cond * ((ut - unorm) / 100) + 0,
  # [0073C4] Change in utilization rate.
  dut_cond ~ ifelse((ut - unorm) > 0, 1, 0),
  # Accounting MEMO
  WBo ~ wo * No,
  WBc ~ WB - WBo,
  #----------------------------------#
  # Expectations X_e = X[-1] + sigma * (X_e[-1] - X[-1])
  #----------------------------------#
  # [74] Expected inflation
  p_e ~ pgr[-1] + sigma_p_e * (p_e[-1] - pgr[-1]),
  # [74C1] Inflation.
  pgr ~ p / p[-1] - 1,
  # [75] Expected productivity growth
  prodg_e ~ prodg[-1] + sigma_pg * (prodg_e[-1] - prodg[-1]),
  # [76] Expected income growth
  y_e ~ y[-1] + sigma_yg * (y_e[-1] - y[-1]),
  # [77] Expected worker real income
  yo_e ~ yo_e[-1] * (1 + y_e),
  # [78] Expected capitalist real income.
  yc_e ~ yc_e[-1] * (1 + y_e),
  # [78] Expected capitalist nominal income.
  Yc_e ~ Yc[-1] * (1 + y_e) * (1 + p_e),
  # [79] Expected capital gains for homes (capitalists)
  CGHc_e ~ (ph_e - ph[-1]) * Hc[-1],
  # [80] Helper for Expected growth on capital gains for homes (capitalists).
  phg_e ~ ph[-1] / ph_1[-1] -1 + sigma_pe *
    (phg_e - (ph[-1] / ph_1[-1] -1)) * shockphg_e,
  # [81] Expected price of houses.
  ph_e ~ ph[-1] * (1 + phg_e),
  # [82] Expected real capital gains for homes (capitalists)
  cghc_e ~ (ph_e - ph[-1]) * Hc[-1] / p[-1] * (1 + p_e),
  # [83] Expected capital gains for homes (workers)
  CGHo_e ~ (ph_e - ph[-1]) * Ho[-1],
  # [84] Expected real capital gains for homes (workers)
  cgho_e ~ (ph_e - ph[-1]) * Ho[-1] / p[-1] * (1 + p_e),
  # [85] Expected capital gains on equities
  CGE_e ~ (pe_e - pe[-1]) * E[-1],
  # [86] Helper for Expected growth of capital gains on equities.
  peg_e ~ pe[-1] / pe_1[-1] -1 + sigma_pe *
    (peg_e[-1] - (pe[-1] / pe_1[-1] -1)) * shockpeg_e,
  # [87] pe lag 1,
  pe_1 ~ pe[-1],
  # [88] Expected price of equities.
  pe_e ~ pe[-1] * (1 + peg_e),
  # [89] Expected real capital gains on equities
  cge_e ~ (pe_e - pe[-1]) * E[-1] / p[-1] * (1 + p_e),
  # [90] Expected wealth for capitalists.
  Vc_e ~ Vc[-1] * (1 + y_e) - Cc + CGE_e + CGHc_e,
  # [91] Expected return on equities
  re_e ~ re[-1] + sigma * (re_e[-1] - re[-1]),
  # [92] Expected real return on equities
  rre_e ~ (1 + re_e) / (1 + p_e) -1,
  # [93] Expected return on houses.
  rh_e ~ rh[-1] + sigma * (rh_e[-1] - rh[-1]),
  # [94] Expected real return on houses.
  rrh_e ~ (1 + rh_e) / (1 + p_e) - 1,
)



# Set up parameter and initial conditions

model_para <- sfcr_set(
  mu_2 ~ 20, # parameter
  mu_1 ~ 0.001, # parameter
  alpha_1c	~	0.7, # parameter
  alpha_2c	~	0.025, # parameter
  alpha_3c	~	0.08, # parameter
  alpha_0o	~	0, # parameter
  alpha_1o	~	0.8, # parameter
  alpha_2o	~	0.025, # parameter
  alpha_3o	~	0.08, # parameter
  alpha_4o	~	10, # parameter
  lambda_10	~	0.5, # parameter
  lambda_11	~	0.45, # parameter
  lambda_12	~	0.25, # parameter
  lambda_13	~	0.01, # parameter
  lambda_14	~	0.25, # parameter
  lambda_15	~	0.05, # parameter
  lambda_20	~	0.18, # parameter
  lambda_21	~	0.3, # parameter
  lambda_22	~	0.25, # parameter
  lambda_23	~	0.01, # parameter
  lambda_24	~	0.1, # parameter
  lambda_25	~	0.05, # parameter
  lambda_30	~	0.18, # parameter
  lambda_31	~	0.45, # parameter
  lambda_32	~	0.25, # parameter
  lambda_33	~	0.01, # parameter
  lambda_34	~	0.25, # parameter
  lambda_35	~	0.1, # parameter
  eta	~	0.2, # parameter
  morp	~	0.01, # exogenous
  im	~	0, #parameter
  xi	~	0.25, # exogenous
  spread_1	~	0.005, # exogenous
  spread_2	~	0.0025, # exogenous
  spread_3	~	-0.025, # exogenous
  iota_0	~	0.005, #parameter
  iota_1	~	2, #parameter
  iota_2	~	1, #parameter
  iota_3	~	0.2, #parameter
  iota_4	~	0.4, #parameter
  unorm	~	0.75, 
  beta	~	0.1,
  lambda	~	1.3,
  omega_c	~	0.05,
  rho_1	~	0.01,
  chi_2	~	0.25,
  chi_1	~	0.6,
  tau	~	0.1,
  tau_f	~	0.4,
  tau_d	~	0.2,
  v_1	~	0.5,
  v_2	~	1,
  v_3	~	0.0005,
  o_1	~	0.05,
  o_2	~	1,
  o_0	~	2,
  pi_0	~	0.02,
  pi_1	~	1,
  ny	~	0.0338,
  okun	~	3,
  sigma_p_e	~	0.75,
  sigma_pg	~	0.75,
  sigma_yg	~	0.75,
  sigma_pe	~	0.75,
  sigma	~	0.75,
  phge	~	0,
  rra	~	0.025,
  rrb	~	0.03,
  rrm	~	0.02,
  shockpeg_e	~	0,
  shockphg_e	~	0,
  shockwc_g	~	0,
  shockwo_g	~	0,
  shockprodg	~	0,
)

model_init <- sfcr_set(
  wc	~	4.16666666666667,
  Nc	~	106.4,
  Rents	~	375,
  FD	~	400.65984,
  Tdc	~	450.190222222222,
  Yc	~	798.449084444444,
  yc	~	1218.80851008,
  Cc	~	1820.62222222222,
  cc	~	1170.4,
  CGE	~	0,
  pe	~	5,
  E	~	410.467555555555,
  CGHc	~	0,
  ph	~	3	,
  Hc	~	37,
  re	~	0.195221198156682	,
  CGE_e	~	0,
  CGHc	~	0,
  rh	~	0.01875	,
  vc	~	4256	,
  Vc	~	6620.44444444444	,
  Mc	~	1383	,
  M	~	4501.90222222222	,
  Vc_e	~	6620.44444444444	,
  i	~	157.6,
  I	~	245.155555555555	,
  Co	~	1820.62222222222	,
  co	~	1170.4	,
  Yo	~	355.026762097778	,
  yo	~	228.23148992	,
  wo	~	0.833333333333333	,
  No	~	2021.6,
  Mo	~	2000,
  Vo	~	0,
  vo	~	0,
  rent	~	0.5,
  CGHo	~	0,
  morp	~	0.01	,
  MOo	~	1000	,
  Ho	~	1430	,
  y_e	~	0.02	,
  unorm	~	0.75	,
  WB	~	2128	,
  p	~	1.55555555555556	,
  s	~	2128	,
  prodg	~	0.02	,
  wo_g	~	0.02	,
  prod	~	1	,
  w	~	1	,
  Td	~	450.190222222222	,
  TF	~	340.48	,
  g	~	800	,
  HU	~	100	,
  HN	~	0	,
  IH	~	0	,
  ih	~	0	,
  HND	~	0	,
  HNS	~	0	,
  c	~	1170.4	,
  S	~	3310.22222222222	,
  y	~	0.02	,
  pgr	~	0	,
  prodc	~	4.16666666666667	,
  prodo	~	0.833333333333333	,
  dut	~	0	,
  wc_g	~	0.02	,
  shockprodg	~	0	,
  N	~	2128	,
  omega	~	1	,
  prodg_e	~	0.02	,
  re_e	~	0.195221198156682	,
  rh_e	~	0.01875	,
  pe_e	~	5	,
  phge	~	0	,
  peg_e	~	0	,
  cgho_e	~	0	,
  ur	~	0.05	,
  FT	~	851.2	,
  rrm	~	0.02	,
  Tdo	~	0
)

index <- sfcr_set_index(model_eqs)

# Create the Baseline (Policy Steady state)

model_zezza <- sfcr_baseline(
  equations = model_eqs,
  external = model_para,
  initial = model_init,
  periods = 100,
  max_iter = 350,
  method = "Broyden"
)

# Take a look at the Directed acyclic graph

sfcr_dag_blocks_plot(model_eqs, size = 5)
sfcr_dag_cycles_plot(model_eqs, size = 5)

bs_zezza <- sfcr_matrix(
  columns = c("Household (Top 5%)" ,"Household (Bottom 95%)", "Firms", "Banks", "Central Bank", "Government", "Sum"),
  codes = c("h5", "h95", "f", "b", "cb", "g", "s"),
  r1 = c("Productive Capital", f = "+p * K", s = "+p * K"),
  r2 = c("Homes", h5 = "+ph * Hc", h95 = "+ph * Ho", f = "+ph * HU", s = "+ph * HNS"),
  r3 = c("Cash", h5 = "+HPc", h95 = "+HPo", b = "+HPb", cb = "-HP"),
  r4 = c("CB Advances", b = "-A", cb = "+A"),
  r5 = c("Bank deposits", h5 = "+Mc", h95 = "+Mo", b = "-M"),
  r6 = c("Loans to firms", f = "-L", b = "+L"),
  r7 = c("Mortgages", h95 = "-MOo", b = "+MOo"),
  r8 = c("Treasuries", h5 = "+Bh", b = "+Bb", cb = "+Bcb", g = "-B"),
  r9 = c("Equities", h5 = "+pe * E", f = "-pe * E"),
  r10 = c("Balance", h5 = "-Vc", h95 = "-Vo", f = "-Vf", g = "+GD", s = "+(p * K + ph * HNS)")
)

sfcr_matrix_display(bs_zezza, "bs")

tfm_zezza <- sfcr_matrix(
  columns = c("Household (Top 5%)" ,"Household (Bottom 95%)", "Firms Cur.", "Firms Cap.", "Banks", "Central Bank", "Government", "Production"),
  codes = c("h5", "h95", "fcu", "fca", "b", "cb", "g", "prod"),
  r1 = c("Wages", h5 = "+WBc", h95 = "+WBo", prod = "-WB"),
  r2 = c("Consumption", h5 = "-p * cc", h95 = "-p * co", prod = "+p * c"),
  r3 = c("Profit Firms", h5 = "+FD", fcu = "+FT", fca = "+FU", prod = "-FT"),
  r4 = c("Profit Banks", h5 = "+FB", b = "-FB"),
  r5 = c("Profit Central Bank", cb = "-FC", g = "+FC"),
  r6 = c("Rents", h5 = "+Rents", h95 = "-Rents"),
  r7 = c("Government Spending", g = "-p * g", prod = "+p * g"),
  r8 = c("Taxes", h5 = "-Tdc", h95 = "-Tdo", fcu = "-TF", g = "+TT", prod = "-IT"),
  r9 = c("Investment in productive capital", fca = "+(k - k[-1]) * p", prod = "-(k - k[-1]) * p"),
  r10 = c("Investment in Housing", fca = "+(HN - HN[-1]) * ph", prod = "-(HN - HN[-1]) * ph"),
  r11 = c("Interest on Deposits", h5 = "+rm[-1]*Mc[-1]", h95 = "+rm[-1]*Mo[-1]", b = "-rm[-1]*M[-1]"),
  r12 = c("Interest on Advances", b = "-ra[-1]*A[-1]", cb = "+ra[-1]*A[-1]"),
  r13 = c("Interest on Loans", fcu = "-rl[-1]*L[-1]", b = "+rl[-1]*L[-1]"),
  r14 = c("Interest on Mortgages", h95 = "-rmo[-1]*MOo[-1]", b = "+rmo[-1]*MOo[-1]"),
  r15 = c("Interest on Bills", h5 = "+rb[-1]*Bh[-1]", b = "+rb[-1]*Bb[-1]", cb = "+rb[-1]*Bcb[-1]", g = "-rb[-1]*B[-1]"),
  r16 = c("Change in Cash", h5 = "-(HPc - HPc[-1])", h95 = "-(HPo - HPo[-1])", b = "-(HPb - HPb[-1])", cb = "+(HP - HP[-1])"),
  r17 = c("Change in Deposits", h5 = "-(Mc - Mc[-1])", h95 = "-(Mo - Mo[-1])", b = "+(M - M[-1])"),
  r18 = c("Change in Loans", fca = "+(L - L[-1])", b = "-(L - L[-1])"),
  r19 = c("Change in Mortgages", h95 = "+(MOo - MOo[-1])", b = "-(MOo - MOo[-1])"),
  r20 = c("Change in Bills", h5 = "-(Bh - Bh[-1])", b = "-(Bb - Bb[-1])", cb = "-(Bcb - Bcb[-1])", g = "+(B - B[-1])"),
  r21 = c("Change in Advances", b = "+(A - A[-1])", cb = "-(A - A[-1])"),
  r22 = c("Change in Equities", h5 = "-(E - E[-1]) * pe", fca = "+(E - E[-1]) * pe")
)
sfcr_matrix_display(tfm_zezza, "tfm")

#Check for model consistency
sfcr_validate(bs_zezza, model_zezza, which = "bs")
sfcr_validate(tfm_zezza, model_zezza, which = "tfm")



# Check Model for consistent

sfcr_validate(bs_pc, model_sim, which = "bs")

sfcr_validate(tfm_pc, model_sim, which = "tfm")

all.equal(model_sim$Hh, model_sim$Hh1)

sfcr_sankey(tfm_pc, model_sim, when = "end")
# Equilibrium values of variables

model_var <- model_sim %>% 
  select(period,Y, TX, r, Bh) %>% 
  mutate("YD" = Y - TX + dplyr::lag(r) * dplyr::lag(Bh)) %>% 
  
  pivot_longer(names_to = "vars", values_to = "val", -1)

ggplot(data = model_var, aes(x = period, y = val, color = vars)) +
  geom_line()

# Add shock to system

shocka0 <- sfcr_shock(
  variables = list(
    a0 ~ 60
  ),
  start = 10,
  end = 50
)

model_shock <- sfcr_scenario(
  baseline = model_sim,
  scenario = shocka0,
  periods = 50)


# Plot scenario

dtt <- model_shock %>% 
  pivot_longer(names_to = "names", values_to =  "values", - period)

dtt %>%
  filter(names %in% c("Yt", "YDt", "Lt", "Ct", "Kt")) %>%
  ggplot(aes(x = period, y = values, color = names)) +
  labs(x = "Periods", y = "Value") +
  scale_color_brewer(palette = "Paired") +
  geom_line(linewidth = 1.5) +
  theme_minimal()












