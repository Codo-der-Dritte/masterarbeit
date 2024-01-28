# Title: Model Zezza

# Purpose : This is the model developed in the paper called
# U.S. growth, the housing market, and the distibution of
# income.
#           

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
  # [3p] Consumption deflator.
  Cc ~ cc * p,
  # [4] Consumption - depend on expected disposable income, the open stock of wealth  and
  # expected capital gains, all measured in real terms.
  cc ~ alpha_1c * yc_e + alpha_2c * vc[-1] + alpha_3c * 
    (cge_e + cghc_e - p_e * vc[-1]/(1 + p_e)),
  # [5] Stock of Wealth - includes saving and capitals gains which are given from changes in
  # market price for equities [7] and homes [8].
  Vc ~ Vc[-1] + Shc +  CGE + CGHc,
  # [5p] Wealth deflator.
  vc ~ Vc / p,
  # [6p] Disposable income deflator.
  yc ~ Yc / p,
  # [7] Changes in market price for equities -  Capital gains.
  CGE ~ (pe - pe[-1]) * E[-1],
  # [8] Changes in market price for homes - Capital gains.
  CGHc ~ (ph - ph[-1]) * Hc[-1],
  #-----Portfolio Choice-----#
  # [9] Cash - depend on current consumption.
  HPc ~ eta * Cc,
  # [10] Bank deposits.
  Mc ~ Vc - HPc - Bc - E * pe - ph * Hc,
  # [11] Bonds.
  Bh ~ (Vc_e - HPc) * 
    (lambda_10 - lambda_11 * rrm  - lambda_12 * 
       rre_e - lambda_13 * (Yc_e/Vc_e) + lambda_14 *
       rrb - lambda_15 * rrh_e),
  # [12] Equities.
  pe ~ ((Vc_e - HPc) * 
         (lambda_20 - lambda_21 * rrm + lambda_22 * 
            rre_e - lambda_23 * (Yc_e/V_e) - lambda_24 *
            rrb - lambda_25 * rrh_e)
  ) / E,
  # [13] Homes.
  Hc ~ ((Vc_e - HPc) * 
          (lambda_30 - lambda_31 * rrm - lambda_32 * 
             rre_e - lambda_33 * (Yc_e/V_e) - lambda_34 *
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
  # [18p] Consumption deflator.
  Co ~ co * p,
  # [19] Consumption - depends on expected real disposable income, past 
  # real wealth, and expected real capital gains on homes minus past
  # wealth normalized by expected inflation.
  co ~ alpha_1o * yo_e + alpha_2o * vo[-1] + 
    alpha_3o * (cgho_e - p_e * vo[-1]/(1 + p_e)) +
    iec - alpha_4o * morp * MOo[-1]/Yo[-1],
  # [20] Wealth - Past wealth plus saving plus capital gains from homes.
  Vo ~ Vo[-1] + Sho + CGHo,
  # [20p] Wealth deflator.
  vo ~ Vo / p,
  # [21] Disposable income deflator.
  yo ~ Yo/p,
  # [22] Capital gains from homes.
  CGHo ~ (ph - ph[-1]) * Ho[-1],
  # [23] Imitation paramter.
  iec ~ alpha_4o * No * (cc[-1]/Nc - 
                         co[-1]/No),
  #-----Portfolio Choice-----#
  # [24] Cash - depends on current consumption.
  HPo ~ eta * Co,
  # [25] Bank deposits - residual.
  Mo ~ Vo - HPo - ph * Ho + (MOo - MOo[-1]),
  # [26] Demand for Homes - depends on population growth,
  # expected real income and lagged debt repayment ratios.
  Ho ~ (((No - No[-1])/No[-1] + mu_1 * ((yo_e - yo_e[-1])/yo_e) -
    mu_2 * delta_debt_rep[-1]) * Ho[-1]) + Ho[-1],
  # [26h] Debt repayment ratio.
  debt_rep ~ (rmo[-1] + morp) * Mo_1[-1] / Yo,
  # [26h] Change in debt repayment.
  delta_debt_rep ~ debt_rep - debt_rep[-1],
  # [26h] Mo lag1
  Mo_1 ~ Mo[-1],
  # [27] Mortgages.
  MOo ~ (MOo[-1] * (1 - morp) + 
          ifelse((ph * (Ho - Ho[-1]) - Sho) > 0, 1, 0) * 
          (ph * (Ho - Ho[-1]) - Sho)) + MOo[-1],
  # [28] Share of rented homes owned by capitalist.
  Rents ~ rent * Hc[-1],
  # [29] Rent increases
  rent ~ rent[-1] * (1 + y_e),
  #----------------------------------#
  # Nonfinancial firms
  #----------------------------------#
  # [30] Investment decision - Growth of k depends on actual profits, Tobin's q,
  # borrowing costs from banks and utilization rate.
  delta_k ~ (iota_0 + iota_1 * FU[-1]/K_1[-1] - 
          iota_2 * rll[-1] * lev[-1] +
          iota_3 * q[-1] +
          iota_4 * ut[-1] - unorm),
  # [30h] Leverage.
  lev ~ L/K,
  # [30h] Tobin's Q.
  q ~ pe * E/K,
  # [30h] Real Capital.
  k ~ (1 + delta_k) * k[-1],
  # [30h] Real Investment,
  i ~ k - k[-1],
  # [31] Prices - are set with a mark-up on wages.
  p ~ (1 + rho) * wage/(prod * (1 - tau)),
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
  # [35h] Normal Sales - fixed ratio of real capital.
  sfc ~ (lambda * k[-1]),
  # [36] Retained Profits.
  FU ~ FT - rl[-1] * L[-1] - FD - TF,
  # [37] New equities issued.
  E ~ E[-1] + xi * (i - FU) / pe,
  # [38] Loan changes - rewritten to display capital.
  L ~ L[-1] + i - FU - pe * (E - E[-1]),
  # [30h] Capital in non real terms.
  K ~ k * p,
  # [30h] Capital from one periods ago.
  K_1 ~ K[-1],
  # [Memo] Firms wealth.
  Vf ~ p * K + HU * ph - L - pe * E,
  #----------------------------------#
  # Banks
  #----------------------------------#
  # [39h] Deposits.
  M ~ Mo + Mc,
  # [39] Demand for Bonds.
  Bb ~ chi_1 * M,
  # [40] Reserve requirement.
  HPb ~ chi_2 * M,
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
    (rm[-1] * M + ra[-1] * A[-1]),
  #----------------------------------#
  # Central bank
  #----------------------------------#
  # [46] Accommodates demand for advances and buys bonds that are
  # not absorbed by Households and Banks.
  Bc ~ B - Bh - Bb,
  # [47] Interest income is redistributed to the government.
  FC ~ ra[-1] * A[-1] + rb[-1] * Bc[-1],
  #----------------------------------#
  # Government
  #----------------------------------#
  # [48] Deficit -  Collects taxes production, wages, and profits,
  # and any deficit is financed by issuing bonds.
  GD ~ (G + rb[-1] * Bh[-1] + rb[-1] * Bc[-1] + rb[-1] * Bb[-1]) - (IT + DT + TF + FC),
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
  # [55] Government spending deflator.
  G ~ g * p,
  # [56] Real Government spending.
  g ~ g[-1] * (1 + y_e),
  #----------------------------------#
  # The housing market
  #----------------------------------#
  # Demand for houses is laid down for capitalists and
  # workers in equations [13] and [26].
  # [57] Number of unsold homes - changes when the number of
  # newly built homes exceed the demand for homes. Can't be negative.
  HU ~ ifelse(((HU[-1] + HN - (Hc - Hc[-1])-(Ho - Ho[-1]))) > 0, 1, 0) *
    ((HU[-1] + HN - (Hc - Hc[-1])-(Ho - Ho[-1]))) + 0,
  # [57] HU lag 1.
  HU_1 ~ HU[-1],
  # [58] Supply of new homes - is a function of expected demand 
  # and past capital gains.
  HND ~ ifelse(v_1 * (Hc[-1] * y_e + (Ho - Ho[-1])) + v_2 * (ph_3 - ph_3[-1]) > 0, 1, 0)
  * v_1 * (Hc[-1] * y_e + (Ho - Ho[-1])) + v_2 * (ph_3 - ph_3[-1]) + 0,
  # [58h] Check for HND.
  HN ~ ifelse(HND - HU[-1] > 0, 1, 0) * HND + 0,
  # [58h] ph lag 1.
  ph_1 ~ ph[-1],
  # [58h] ph lag 2.
  ph_2 ~ ph_1[-1],
  # [58h] ph lag 3.
  ph_3 ~ ph_2[-1],
  # [59] Market price of Homes.
  ph ~ ifelse(HU < 10, 1, 0) * ph1 +
    ifelse(HU >= 10, 1, 0) * ph2,
  # [59h] Total supply of Houses.
  HNS ~ HN + HU[-1],
  # [59h] Helper for House prices ph2.
  ph2 ~ ph(-1) - v_3 * (HU_1 - HU_1(-1)),
  # [59h] Helper for House prices ph1.
  ph1 ~ ph(-1) + ph(-1) * v_3 * (Hc(-1) * y_e + (Ho - Ho[-1]) - NHS(-1)),
  #----------------------------------#
  # Aggregate demand, un-/employment, wages
  #----------------------------------#
  # [60] Aggregate Demand -  Sales
  s ~ c + i + ih  + g,
  # [60h] Helper for Aggregate Demand - total consumption.
  c ~ cc + co,
  # [60h] Real value of houses.
  ih ~ NH * p,
  # [60h] Nominal value of houses
  IH ~ ih * p,
  # [61] Aggregate Demand deflator.
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
  ur ~ -psi * (((y-y[-1])/y[-1]) - y_n) + ur[-1],
  # [67] Wage Bill - All wages paid out.
  WB ~ w * N,
  # [68] Wage for capitalists.
  wc ~ wc[-1] * (1 + (wc_g)),
  # [69] Wage for workers.
  wo ~ wo[-1] * (1 + (wo_g)),
  # [68 + 69] Wage.
  w ~ (wc * omega_c + wo * (1 - omega_c)),
  # [70] Wage growth/inflation for capitalists.
  wc_g ~ p_e + omega * prodg_e + shockwc_g,
  # [71] Wage growth/inflation for workers.
  wo_g ~ p_e + omega * prodg_e + shockwo_g,
  # [xx] Inflation.
  p_ ~ p / p[-1] - 1,
  # [72] Omega,
  omega ~ o0 - o2 * sqrt(ur1/o1),
  # [73] Productivity gains.
  prodg ~ pi_0 - pi_1 * u + shockprodg,
  # [xx] Unemployment rate.
  u ~ ifelse(ut - unorm > 0, 1, 0) * ((ut - unorm) / 100) + 0,
  # [xx] Okuns law.
  ur ~ ifelse(((ur[-1]) - (y - ny) / okun) < 0, 1, 0) * 0 +
    ifelse(((ur[-1]) - (y - ny) / okun) > 0, 1, 0) * ((ur[-1]) - (y - ny) / okun),
  ur1 ~ ifelse(ur <=0, 1, 0) * 0 +
    ifelse(ur >= 0.2) * 0.2 +
    ifelse(ur > 0 & ur < 0.2, 1, 0) * ur,
  # Accounting MEMO
  WBo ~ wo * No,
  WBc ~ WB - WBo,
  #----------------------------------#
  # Expectations X_e = X[-1] + sigma * (X_e[-1] - X[-1])
  #----------------------------------#
  # [74] Expected inflation
  p_e ~ p_[-1] + sigma_p_e * (p_e[-1] - p_[-1]),
  # [75] Expected productivity growth
  prodg_e ~ prodg[-1] + sigma * (prodg_e[-1] - prodg[-1]),
  # [76] Expected income growth
  #missing,
  # [77] Expected worker income
  Yo_e ~ Yo[-1] + sigma * (Yo_e[-1] - Yo[-1]),
  # [78] Expected capitalist income
  Yc_e ~ Yc[-1] + sigma * (Yc_e[-1] - Yc[-1]),
  # [780] Expected capitalist income
  yc_e ~ Yc / p,
  # [79] Expected capital gains for homes (capitalists)
  CGHc_e ~ CGHc[-1] + sigma * (CGHc_e[-1] - CGHc[-1]),
  # [79p] Expected capital gains for homes (capitalists)
  cghc_e ~ CGHc_e / p,
  # [80] Expected capital gains for homes (workers)
  CGHo_e ~ CGHo[-1] + sigma * (CGHo_e[-1] - CGHo[-1]),
  # [80] Expected capital gains for homes (workers)
  cgho_e ~ CGHo_e / p,
  # [81] Expected capital gains on equities
  CGE_e ~ CGE[-1] + sigma * (CGE_e[-1] - CGE[-1]),
  # [81p] Expected real capital gains on equities
  cge_e ~ CGE_e / p,
  # [82] Expected wealth for capitalists.
  Vc_e ~ Vc[-1] + sigma * (Vc_e[-1] - Vc[-1]),
  # [83] Expected return on equities
  re_e ~ re[-1] + sigma * (re_e[-1] - re[-1]),
  # [83p] Expected real return on equities
  rre_e ~ re_e / p,
  # [84] Expected return on equities
  rh_e ~ rh[-1] + sigma * (rh_e[-1] - rh[-1]),
  # [84p] Expected real return on equities
  rrh_e ~ rh_e / p,
)

# Take a look at the Directed acyclic graph

sfcr_dag_blocks_plot(model_eqs)
sfcr_dag_cycles_plot(model_eqs)
# Set up parameter and initial conditions

bs_zezza <- sfcr_matrix(
  columns = c("Household (Top 5%)" ,"Household (Bottom 95%)", "Firms", "Banks", "Central Bank", "Government", "Sum"),
  codes = c("h5", "h95", "f", "b", "cb", "g", "s"),
  r1 = c("Productive Capital", f = "+p * K", s = "+p * K"),
  r2 = c("Homes", h5 = "+ph * Hc", h95 = "+ph * Ho", f = "+ph * HU", s = "+ph * H"),
  r3 = c("Cash", h5 = "+HPhc", h95 = "+HPho", b = "+HPb", cb = "-HP"),
  r4 = c("CB Advances", b = "-A", cb = "+A"),
  r5 = c("Bank deposits", h5 = "+Mc", h95 = "+Mo", b = "-M"),
  r6 = c("Loans to firms", f = "-L", b = "+L"),
  r7 = c("Mortgages", h95 = "-MO", b = "+MO"),
  r8 = c("Treasuries", h5 = "+Bh", b = "+Bb", cb = "+Bcb", g = "-B"),
  r9 = c("Equities", h5 = "+pe * E", f = "-pe * E"),
  r10 = c("Balance", h5 = "-Vc", h95 = "-Vo", f = "-Vf", g = "+GD", s = "+(p * K + ph * H)")
)

sfcr_matrix_display(bs_zezza, "bs")

tfm_zezza <- sfcr_matrix(
  columns = c("Household (Top 5%)" ,"Household (Bottom 95%)", "Firms Cur.", "Firms Cap.", "Banks", "Central Bank", "Government", "Production"),
  codes = c("h5", "h95", "fcu", "fca", "b", "cb", "g", "prod"),
  r1 = c("Wages", h5 = "+WBc", h95 = "+WBo", prod = "-WB"),
  r2 = c("Consumption", h5 = "-p * Cc", h95 = "-p * Co", prod = "+p * C"),
  r3 = c("Profit Firms", h5 = "+FD", fcu = "+FT", fca = "+FU", prod = "-FT"),
  r4 = c("Profit Banks", h5 = "+FB", b = "-FB"),
  r5 = c("Profit Central Bank", cb = "-FC", g = "+FC"),
  r6 = c("Rents", h5 = "+Rents", h95 = "-Rents"),
  r7 = c("Government Spending", g = "-p * G", prod = "+p * G"),
  r8 = c("Taxes", h5 = "-TDc", h95 = "-TDo", fcu = "-TF", g = "+T", prod = "-IT"),
  r9 = c("Investment in productive capital", fca = "+(K - K[-1]) * p", prod = "-(K - K[-1]) * p"),
  r10 = c("Investment in Housing", fca = "+(HN - HN[-1]) * ph", prod = "-(HN - HN[-1]) * ph"),
  r11 = c("Interest on Deposits", h5 = "+rm[-1]*Mc[-1]", h95 = "+rm[-1]*Mo[-1]", b = "-rm[-1]*M[-1]"),
  r12 = c("Interest on Advances", b = "-ra[-1]*A[-1]", cb = "+ra[-1]*A[-1]"),
  r13 = c("Interest on Loans", fcu = "-rl[-1]*L[-1]", b = "+rl[-1]*L[-1]"),
  r14 = c("Interest on Mortgages", h95 = "-rmo[-1]*MO[-1]", b = "+rmo[-1]*MO[-1]"),
  r15 = c("Interest on Bills", h5 = "+rb[-1]*Bh[-1]", b = "+rb[-1]*Bb[-1]", cb = "+rb[-1]*Bcb[-1]", g = "-rb[-1]*B[-1]"),
  r16 = c("Change in Cash", h5 = "-(HPhc - HPhc[-1])", h95 = "-(HPho - HPho[-1])", b = "-(HPb - HPb[-1])", cb = "+(HP - HP[-1])"),
  r17 = c("Change in Deposits", h5 = "-(Mc - Mc[-1])", h95 = "-(Mo - Mo[-1])", b = "+(M - M[-1])"),
  r18 = c("Change in Loans", fca = "+(L - L[-1])", b = "-(L - L[-1])"),
  r19 = c("Change in Mortgages", h95 = "+(MO - MO[-1])", b = "-(MO - MO[-1])"),
  r20 = c("Change in Bills", h5 = "-(Bh - Bh[-1])", b = "-(Bb - Bb[-1])", cb = "-(Bcb - Bcb[-1])", g = "+(B - B[-1])"),
  r21 = c("Change in Advances", b = "+(A - A[-1])", cb = "-(A - A[-1])"),
  r22 = c("Change in Equities", h5 = "-(E - E[-1]) * pe", fca = "+(E - E[-1]) * pe")
)
sfcr_matrix_display(tfm_zezza, "tfm")

index <- sfcr_set_index(model_eqs)
sfcr::sfcr_get_matrix(bs_zezza)

model_para <- sfcr_set(
  
  # Parameters

  # Exogenous variables

)

# Create the Baseline (Policy Steady state)

model_zezza <- sfcr_baseline(
  equations = model_eqs,
  external = model_para,
  periods = 100,
  
)

# Check Model for consisteny

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












