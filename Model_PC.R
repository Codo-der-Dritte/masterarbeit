# Title: Model PC

# Purpose : This is a training model of the Textbook
#           Monetary economics
#           

# Author: Alexander Toplitsch
# Contact details: alexander.toplitsch@s.wu.ac.at

# Date script created: Fr 24.11.2023 17:03 -------------
# Date script last modified:  Sa 28.10.2023 15:47

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

bs_pc <- sfcr_matrix(
  columns = c("Household" , "Production", " Government", "Central Bank", "Sum"),
  codes = c("h", "f", "g", "cb", "s"),
  r1 = c("Money", h = "+Hh", cb = "-Hs"),
  r2 = c("Bills", h = "+Bh", g = "-Bs", cb = "+Bcb"),
  r3 = c("Balance", h = "-V", g = "+V")
)

sfcr_matrix_display(bs_pc)

tfm_pc <- sfcr_matrix(
  columns = c("Households" , "Production", " Government", "CB current", "CB capital"),
  codes = c("h", "f", "g", "cbc", "cbk"),
  c("Consumption", h = "-C", f = "+C"),
  c("Govt. Expenditure", f = "+G", g = "-G"),
  c("Income", h = "+Y", f = "-Y"),
  c("Int. Payments", h = "+r[-1]*Bh[-1]", g = "-r[-1]*Bs[-1]", cbc = "+r[-1]*Bcb[-1]"),
  c("CB profits", g = "+r[-1]*Bcb[-1]", cbc = "-r[-1]*Bcb[-1]"),
  c("Taxes", h = "-TX", g = "+TX"),
  c("Ch. Money", h = "-(Hh - Hh[-1])", cbk = "+(Hh - Hh[-1])"),
  c("Ch. Bills", h = "-(Bh - Bh[-1])", g = "+(Bs - Bs[-1])", cbk = "-(Bcb - Bcb[-1])")
)

sfcr_matrix_display(tfm_pc)

# Set up the model equations

model_eqs <- sfcr_set(
  Y ~ C + G,
  YD ~ Y - TX + r[-1] * Bh[-1],
  TX ~ theta * (Y + r[-1] * Bh[-1]),
  V ~ V[-1] + (YD - C),
  C ~ alpha1 * YD + alpha2 * V[-1],
  Hh ~ V - Bh,
  Bh ~ V * (lambda0 + lambda1 * r - lambda2 * (YD / V)), 
  Hh1 ~ V * ((1 - lambda0) - lambda1 * r + lambda2 * (YD / V)), # EQ 4.6A
  Bs ~ Bs[-1] + (G + r[-1] * Bs[-1]) - (TX + r[-1] * Bcb[-1]),
  Hs ~ Bcb - Bcb[-1] + Hs[-1],
  Bcb ~ Bs - Bh
)

# Take a look at the Directed acyclic graph

sfcr_dag_blocks_plot(model_eqs)
sfcr_dag_cycles_plot(model_eqs)
# Set up parameter and initial conditions

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












