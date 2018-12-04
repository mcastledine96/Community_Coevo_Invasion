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

## expanding code to transform concentrations from OD to CFU/ul then normalise:
# new_sol_CFU is the desired number of CFUs/ul per species
# dilution_vol is the volume (ul) of M9/water added 

CFUs <- function (Pc_OD, Sr_OD, Aa_OD, Vg_OD, Od_OD, new_sol_CFU, dilution_vol) {
  CFU_Pc = (546 * Pc_OD - 20) * 16667
  CFU_Sr = (56625 * Sr_OD - 2246) * 16667
  CFU_Aa = (11676 * Aa_OD - 291) * 16667
  CFU_Vg = (955 * Vg_OD - 30) * 16667
  CFU_Od = (5548 * Od_OD - 67) * 16667
  CFUs <- c(CFU_Pc, CFU_Sr, CFU_Aa, CFU_Vg, CFU_Od)
  return((new_sol_CFU * dilution_vol)/(CFUs - new_sol_CFU))
}

Pc <- c(0.175, 0.222, 0.209, 0.17, 0.139, 0.148, 0.214, 0.144, 0.189)
Sr <- c(0.147, 0.133, 0.146, 0.148, 0.143, 0.138, 0.143, 0.152, 0.148)
Aa <- c(0.102, 0.157, 0.119, 0.1, 0.116, 0.109, 0.116, 0.109, 0.116, 0.133, 0.123)
Vg <- c(0.093, 0.086, 0.084, 0.091, 0.088, 0.103, 0.086, 0.097, 0.089)
Od <- c(0.143, 0.143, 0.125, 0.132, 0.128, 0.133, 0.141, 0.142, 0.12)

CFUs(Pc, Sr, Aa, Vg, Od, 10^5, 500)


CFUs2 <- function (Pc_OD, Sr_OD, Aa_OD, Vg_OD, Od_OD) {
  CFU_Pc = (546 * Pc_OD - 20) * 16667
  CFU_Sr = (56625 * Sr_OD - 2246) * 16667
  CFU_Aa = (11676 * Aa_OD - 291) * 16667
  CFU_Vg = (955 * Vg_OD - 30) * 16667
  CFU_Od = (5548 * Od_OD - 67) * 16667
  CFUs <- c(CFU_Pc, CFU_Sr, CFU_Aa, CFU_Vg, CFU_Od)
  return(CFUs)
}

CFUs2(0.1, 0.1, 0.1, 0.1, 0.1)

#CFU/ul values for 0.1 OD:
#PC = 576678 CFU/ul
#SR = 56942806
#AA = 14610292
#VG = 1091689
#OD = 8130163

1000000 / 576678
#1.7ul

1000000 / 56942806
#0.018ul

1000000 / 14610292
#0.068ul

1000000 / 1091689
#0.916ul

1000000 / 8130163
#0.123ul

#example data
CFUs(0.09, 0.05, 0.05, 0.06, 0.06, 10^5, 10000)

CFUs(0.447, 0.137, 0.097, 0.114, 0.138, 100000, 1000)

CFUs(0.296, 0.132, 0.152, 0.155, 0.137, 100000, 1000)

CFUs2 <- function (Pc_OD, Sr_OD, Aa_OD, Vg_OD, Od_OD) {
  CFU_Pc = (546 * Pc_OD - 20) * 16667
  CFU_Sr = (56625 * Sr_OD - 2246) * 16667
  CFU_Aa = (11676 * Aa_OD - 291) * 16667
  CFU_Vg = (955 * Vg_OD - 30) * 16667
  CFU_Od = (5548 * Od_OD - 67) * 16667
  return(c(CFU_Pc, CFU_Sr, CFU_Aa, CFU_Vg, CFU_Od)) 
}