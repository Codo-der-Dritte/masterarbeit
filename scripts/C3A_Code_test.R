# Title: C3A Code test

# Purpose : This script sets up a model that was proposed in
#           the C3A Code test of Asjad Naqvi. The script in-
#           cludes the model set up, the model baseline, a 
#           shock to the model and the plotted outcome of that
#           shock.
#         

# Author: Alexander Toplitsch
# Contact details: alexander.toplitsch@s.wu.ac.at

# Date script created: Fr 27.10.2023 15:36 -------------
# Date script last modified:  Sa 28.10.2023 15:47

# Check necessary packages for installation and call them

libs <- c("tidyverse", "leaflet", "lubridate", "readxl", "sfcr", "scales", "ggraph")
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
  Yt ~ Ct + It,
  Lt ~ Lt[-1] * (1 - p) + It,
  YDt ~ WBt + Pt +rm * Mt[-1],
  Mt ~ Mt[-1] + YDt -Ct,
  WBt ~ w * Nt,
  Pt ~ Yt - WBt - (rl + p) * Lt[-1],
  Nt ~ Yt / el,
  Ct ~ a0 + a1 * YDt + a2 * Mt[-1],
  Kt ~ Kt[-1] * (1 - delta) + It,
  KTt ~ Yt / ek,
  It ~ gamma * (KTt - Kt[-1]) + delta * Kt[-1]
)

# Take a look at the Directed acyclic graph

sfcr_dag_blocks_plot(model_eqs)

# Set up parameter and initial conditions

model_para <- sfcr_set(
  a0 ~ 50,
  a1 ~ 0.75,
  a2 ~ 0.1,
  delta ~ 0.1,
  p ~ 0.1,
  gamma ~ 0.2,
  el ~ 1,
  ek ~ 1,
  rl ~ 0.04,
  rm ~ 0.04,
  w ~ 1
)

model_initial <- sfcr_set(
Lt ~ 400,
Mt ~ 400,
Kt ~ 400,
Yt ~ 400
)

# Create the Baseline (Policy Steady state)

model_sim <- sfcr_baseline(
  equations = model_eqs,
  external = model_para,
  initial = model_initial,
  periods = 50
)

# Equilibrium values of variables

endo <- model_sim %>% 
  slice(n()) %>% 
  select(c(2:12))

endo


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












