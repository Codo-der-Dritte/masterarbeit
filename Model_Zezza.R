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
  # [3p] v
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
  Mc ~ Vc - HPc - Bc - E * pe - ph *Hc,
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
  ) * pe,
  # [13] Homes.
  Hc ~ ((Vc_e - HPc) * 
          (lambda_30 - lambda_31 * rrm - lambda_32 * 
             rre_e - lambda_33 * (Yc_e/V_e) - lambda_34 *
             rrb + lambda_35 * rrh_e)
  ) * ph,
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
  # [24] Cash - depends on current conmsumption.
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
  MOo ~ (ph * (Ho - Ho[-1]) - Sho - morp * MOo[-1]) + MOo[-1],
  # [28] Share of rented homes owned by capitalist.
  Rents ~ rent * Hcr[-1],
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
  # [30h] Capital from one periods ago.
  K_1 ~ K[-1],
  # [31] Prices - are set with a mark-up on wages.
  p ~ (1 + rho) * wage/(prod * (1 - tau)),
  # [32] Total Profits - are determined relative to the wage bill.
  FT ~ rho * WB,
  # [33] Distributed Profits - fixed share net of taxes and interest payments
  # is ditributed to capitalists.
  FD ~ (1 - beta) * (FT - rl[-1] * L[-1] - TF),
  # [34] Mark-up - depends on relative strength of workers and capitalists.
  rho ~ ((rho_1 * (prodg - wo_) + 1)/(1 + rho[-1])) - 1,
  # [35] Utilization rate - ratio of real sales to "normal" sales, which in
  # turn are in a fixed ratio (lambda) to the stock of real capital.
  u ~ s/(lambda * k[-1]),
  # [36] Retained Profits.
  FU ~ FT -rl[-1] * L[-1] - FD - TF,
  # [37] New equities issued.
  pe ~ (xi * ((K - K[-1]) - FU)) / (E - E[-1]),
  # [38] Loan changes - rewritten to display capital.
  K ~ (L - L[-1]) + pe * (E - E[-1]) + FU + K[-1]
  #----------------------------------#
  # Banks and the central bank
  #----------------------------------#
  
  
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
  r10 = c("Balance", h5 = "+Vc", h95 = "+Vo", f = "+Vf", b = 0, cb = 0, g = "-B", s = "+p * K")
)

sfcr_matrix_display(bs_zezza, "bs")

tfm_zezza <- sfcr_matrix(

)

sfcr_matrix_display(tfm_zezza)



model_para <- sfcr_set(
  
  # Parameters
  alpha1 ~ 0.6, # Consumption YD
  alpha2 ~ 0.4, # Consumption V
  theta ~ 0.2, # Taxes
  lambda0 ~ 0.635, # Ratio of Bonds in wealth of households
  lambda1 ~ 0.05, # Parameter to interest rate
  lambda2 ~ 0.01, # Parameter to YD to V ratio
  
  # Exogenous variables
  r ~ 0.025,
  G ~ 20
)

# Create the Baseline (Policy Steady state)

model_sim <- sfcr_baseline(
  equations = model_eqs,
  external = model_para,
  periods = 100,
  hidden = c("Hh" = "Hs")
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












