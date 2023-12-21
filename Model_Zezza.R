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
  cc ~ alpha_1c * y_e + alpha_2c * vc[-1] + alpha_3c * 
    (cge_e + cghc_e - p_e * vc[-1]/(1 + p_e)),
  # [5] Stock of Wealth - includes saving and capitals gains which are given from changes in
  # market price for equities [7] and homes [8].
  Vc ~ Vc[-1] + Shc +  CGE + CGHc,
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
  E ~ ((Vc_e - HPc) * 
         (lambda_20 - lambda_21 * rrm + lambda_22 * 
            rre_e - lambda_23 * (Yc_e/V_e) - lambda_24 *
            rrb - lambda_25 * rrh_e)
  ) / pe,
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
  Mo ~ Vo - HPo - ph * Ho,
  # [26] Demand for Homes - depends on population growth,
  # expected real income and lagged debt repayment ratios.
  Ho ~ (No - No[-1])/No[-1] - mu_1 * ((yo_e - yo_e[-1])/yo_e) -
    mu_2 * delta_debt_rep[-1],
  # [26h] Debt repayment ratio.
  debt_rep ~ ((rmo[-1]+morp) * Mo[-1]/Yo),
  # [26h] Change debt repayment ratio.
  delta_debt_rep ~ debt_rep - debt_rep[-1],
  # [27] Change in Mortgages
  MO ~ (ph * (Ho - Ho[-1]) - Sho - morp * MO[-1]) + MO[-1],
  # [28] Share of rented homes owned by capitalist.
  Rents ~ rent * Hc[-1],
  # [29] Rent increases
  rent ~ rent[-1] * (1 + y_e),
  #----------------------------------#
  # Nonfinancial firms
  #----------------------------------#
  # [30] Investment decision - depends on actual profits, Tobin's q,
  # borrowing costs from banks and utilization rate.
  k ~ ((iota_0 + iota_1 * FU[-1]/K_1[-1] - iota_2 * rll[-1] * 
    (L[-1]/K[-1]) + iota_3 * (pe[-1] * E[-1]/K[-1]) +
    iota_4 * u[-1]) * k[-1]) + k[-1],
  # [31] Prices - are set with a mark-up on wages.
  p ~ (1 + rho) * wage/(prod * (1 - tau)),
  # [32] Total Profits - are determined relative to the wage bill.
  FT ~ rho * WB,
  # [33] Distributed Profits - fixed share net of taxes and interest payments
  # is distributed to capitalists.
  FD ~ (1 - beta) * (FT - rl[-1] * L[-1] - TF),
  # [34] Mark-up - depends on relative strength of workers and capitalists.
  rho ~ ((rho_1 * (prodg - wo_) + 1)/(1 + rho[-1])) - 1,
  # [35] Utilization rate - ratio of real sales to "normal" sales, which in
  # turn are in a fixed ratio (lambda) to the stock of real capital.
  u ~ s/(lambda * k[-1]),
  # [36] Retained Profits.
  FU ~ FT -rl[-1] * L[-1] - FD - TF,
  # [37] New equities issued.
  pe ~ ((xi * ((K - K[-1]) - FU)) / (E - E[-1])),
  # [38] Loan changes - rewritten to display capital.
  L ~ ((K - K[-1]) - FU - pe * (E - E[-1])) + L[-1],
  # [30h] Capital in non real terms.
  K ~ k * p,
  # [30h] Capital from one periods ago.
  K_1 ~ K[-1],
  # [Memo] Firms wealth.
  Vf ~ p * K + HU * ph - L - pe * E,
  #----------------------------------#
  # Banks
  #----------------------------------#
  # [39] Demand for Bonds.
  Bb ~ chi_1 * (Mc + Mo),
  # [40] Reserve requirement.
  HPb ~ chi_2 * (Mc + Mo),
  # [41] Advancements from Central Bank - if internal funds are
  # not sufficient to cover for demand for loans the banks gets
  # advances from the central bank.
  A ~ L + HPb + Bb + MO - (Mc + Mo),
  # [42] Interest rates on loans.
  rl ~ ra + spread_1,
  # [43] Interest rates on mortgages.
  rmo ~ ra + spread_2,
  # [44] Interest rates on deposits.
  rm ~ ra + spread_3,
  # [45] Banks profits, all distributed.
  FB ~ rl[-1] * L[-1] + rb[-1] * Bb[-1] + rmo[-1] * MO[-1] -
    (rm[-1] * (Mc[-1] + Mo[-1]) + ra[-1] * A[-1]),
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
  GD ~ (G + rb[-1] * B[-1]) - (IT + DT + TF + FC),
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
  B ~ GD + B[-1],
  # [55] Government spending deflator.
  G ~ g * p,
  # [56] Government spending.
  g ~ g[-1] * (1 + y_e),
  #----------------------------------#
  # The housing market
  #----------------------------------#
  # Demand for houses is laid down for capitalists and
  # workers in equations [13] and [26].
  # [57] Number of unsold homes - changes when the number of
  # newly built homes exceed the demand for homes.
  HU ~ (HN - (Hc - Hc[-1])-(Ho - Ho[-1])) + HU[-1],
  # [58] Supply of new homes - is a function of expected demand 
  # and past capital gains.
  HN ~ v_1 * (Hc[-1] * y_e + (Ho - Ho[-1])) + v_2 * (ph_1 - ph_1[-1]),
  # [58h] Helper for ph.
  ph_1 ~ ph[-1],
  # [59] Market price of Homes.
  ph ~ (-v_3 * (HU - HU[-1]) * ph[-1]) + ph[-1],
  #----------------------------------#
  # Aggregate demand, un-/employment, wages
  #----------------------------------#
  # [60] Aggregate Demand -  Sales
  s ~ cc + co + (k - k[-1]) * p + HN * p + g,
  # [61] Aggregate Demand deflator.
  S ~ s * p,
  # [62] Number of workers
  N ~ s/prod,
  # [63] Share of capitalists
  Nc ~ omega_c * N,
  # [64] Share of workers
  No ~ N - Nc,
  # [65] growth in real income
  y ~ ((s/s[-1]-1)*y[-1]) + y[-1],
  # [66] Unemployment rate - Follows some sort of Okun's laws
  ur ~ -psi * (((y-y[-1])/y[-1]) - y_n) + ur[-1],
  # [67] Wage Bill - All wages paid out.
  WB ~ (wc * omega_c + wo * (1 - omega_c)) * N,
  # [68] Wage for capitalists.
  wc ~ wc[-1] * (1 + (wc_g)),
  # [69] Wage for workers.
  wo ~ wo[-1] * (1 + (wo_g )),
  # [70] Wage growth for capitalists.
  wc_g ~   omega * prodg_e,
  # [71] Wage growth for workers.
  wo_g ~   omega * prodg_e,
  # [72] Omega,
  omega ~ -2 * ur,
  # [73] Productivity gains.
  prodg ~ pi_0 - pi_1 * u,
  #----------------------------------#
  # Expectations X_e = X[-1] + sigma * (X_e[-1] - X[-1])
  #----------------------------------#
  # [74] Expected inflation
  p_e ~ p[-1] + sigma * (p_e[-1] - p[-1]),
  # [75] Expected productivity growth
  prodg_e ~ prodg[-1] + sigma * (prodg_e[-1] - prodg[-1]),
  # [76] Expected income growth
  #missing
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
  r3 = c("Cash", h5 = "+HP_hc", h95 = "+HP_ho", b = "+HP_b", cb = "-HP"),
  r4 = c("CB Advances", b = "-A", cb = "+A"),
  r5 = c("Bank deposits", h5 = "+M_c", h95 = "+M_o", b = "-M"),
  r6 = c("Loans to firms", f = "-L", b = "+L"),
  r7 = c("Mortgages", h95 = "-MO", b = "+MO"),
  r8 = c("Treasuries", h5 = "+B_h", b = "+B_b", cb = "+B_cb", g = "-B"),
  r9 = c("Equities", h5 = "+pe * E", f = "-pe * E"),
  r10 = c("Balance", h5 = "-Vc", h95 = "-Vo", f = "-Vf", g = "GD", s = "-(p * K + ph * H)")
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
  r8 = c("Taxes", h5 = "-TDc", h95 = "-TDo", fcu = "-TF", g = "T", prod = "-IT"),
  r9 = c("Investment in productive capital", fca = "+(K - K[-1]) * p", prod = "-(K - K[-1]) * p"),
  r10 = c("Investment in Housing", fca = "+(HN - HN[-1]) * ph", prod = "-(HN - HN[-1]) * ph"),
  r11 = c("Interest on Deposits", h5 = "+rm[-1]*Mc[-1]", h95 = "+rm[-1]*Mo[-1]", b = "-rm[-1]*M[-1]"),
  r12 = c("Interest on Advances", b = "-ra[-1]*A[-1]", cb = "+ra[-1]*A[-1]"),
  r13 = c("Interest on Loans", fcu = "-rl[-1]*L[-1]", b = "+rl[-1]*L[-1]"),
  r14 = c("Interest on Mortgages", h95 = "-rmo[-1]*MO[-1]", b = "+rmo[-1]*MO[-1]"),
  r15 = c("Interest on Bills", h5 = "+rb[-1]*Bh[-1]", b = "+rb[-1]*Bb[-1]", cb = "+rb[-1]*Bcb[-1]", g = "-rb[-1]*B[-1]"),
  r16 = c("Change in Cash", h5 = "-(HP_hc - HP_hc[-1])", h95 = "-(HP_ho - HP_ho[-1])", b = "-(HP_b - HP_b[-1])", cb = "+(HP - HP[-1])"),
  r17 = c("Change in Deposits", h5 = "-(Mc - Mc[-1])", h95 = "-(Mo - Mo[-1])", b = "+(M - M[-1])"),
  r18 = c("Change in Loans", fca = "+(L - L[-1])", b = "-(L - L[-1])"),
  r19 = c("Change in Mortgages", h95 = "+(MO - MO[-1])", b = "-(MO - MO[-1])"),
  r20 = c("Change in Bills", h5 = "-(Bh - Bh[-1])", b = "-(Bb - Bb[-1])", cb = "-(Bcb - Bcb[-1])", g = "+(B - B[-1])"),
  r21 = c("Change in Advances", b = "+(A - A[-1])", cb = "-(A - A[-1])"),
  r22 = c("Change in Equities", h5 = "-(E - E[-1]) * pe", fca = "+(E - E[-1]) * pe")
)
sfcr_matrix_display(tfm_zezza, "tfm")


model_para <- sfcr_set(
  
  # Parameters

  
  # Exogenous variables

)

# Create the Baseline (Policy Steady state)

model_sim <- sfcr_baseline(
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












