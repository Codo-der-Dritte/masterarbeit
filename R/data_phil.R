# This script is created so check if the phillips curve holds in Austria


# load data
dtt_unemp <- read_excel("data/data_raw/statat_unemp.xlsx", sheet = "form")
dtt_wage <- read_excel("data/data_raw/eurostat_wage.xlsx", sheet = "form")

# calculate wage growth rates and transform unemployment to non-percent

dtt_wage <- dtt_wage %>% 
  mutate(wert = wert/lag(wert)-1)

dtt_unemp <- dtt_unemp %>% 
  mutate(wert = wert/100)

dtt_phil <- left_join(dtt_unemp, dtt_wage, by = "jahr", suffix = c("unemp", "wage"))

# plot the relationship

plot_phil <- ggplot(data = dtt_phil) +
  geom_point(aes(x = wertunemp, y = wertwage)) +
  geom_line(aes(x = wertunemp, y = wertwage)) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(aes(x = wertunemp, y = wertwage, label = jahr)) +
  geom_smooth(aes(x = wertunemp, y = wertwage), color = "red") +
  labs(title = "The Phillips Curve in Austria",
       subtitle = "from 2000 to 2021",
       x = "Unemployment Rate", y = "Wage Growth")+
  theme_minimal()

