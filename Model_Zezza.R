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
  # interest income from deposits and treasuries, and rents, net of direct taxes paid
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












