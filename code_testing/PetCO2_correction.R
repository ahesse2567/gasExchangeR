# PetCO2 correction

# solubility of CO2 in blood (37°C): 0.0308 mmol / L / mmHg
# source https://academic.oup.com/bjaed/article/5/6/207/331369

# dVCO2 = [B] * Q * dPetCO2
# = 0.0308 mmol/L * L/min * mmHg = mmol/min
# convert mmol/min to mL/min
# convert mmol to mol
0.0308 * 10^-3
# multiply by molecular mass: 44.01 g/mol
0.0308 * 10^-3 * 44.01 # gCO2/L/mmHg
# divide by density of CO2 at 37°C and 1 atm: 1.713 kg/m3 = 1.713 g/L
0.0308 * 10^-3 * 44.01 / 1.713 # L CO2/L/mmHg
# multiple by 1000 mL/L
0.0308 * 10^-3 * 44.01 / 1.713 * 1000
# 0.7913065 mLCO2/L/mmHg


# use a golden section search or another algorithm to find the minimum
