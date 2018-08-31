# stock sol vol 2

# stock_sol_conc is the concentration of the stock solution - current OD
# new_sol_conc is the concentration of the new solution - the desired OD
# diluent vol is the volume of the diluent you want to use (i.e. the amount of water/M9)


stock_sol_vol2 <- function (stock_sol_conc, new_sol_conc, diluent_vol)
{
  return((new_sol_conc * diluent_vol)/(stock_sol_conc - new_sol_conc))
}

stock_sol_vol2(0.2, 0.1, 200)
stock_sol_vol2(0.5, 0.1, 200)
stock_sol_vol2(1.2, 0.1, 200)
stock_sol_vol2(c(1.2, 0.8, 0.5), 0.1, 200)

