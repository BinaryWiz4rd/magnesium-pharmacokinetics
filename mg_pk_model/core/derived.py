import math
from .constants import LITERATURE_CONSTANTS, ASSUMPTIONS

BW = ASSUMPTIONS["Body weight (kg)"]
T_HALF = ASSUMPTIONS["t½ elimination (h)"]

K_EL = math.log(2) / T_HALF
VD = LITERATURE_CONSTANTS["Vd (L/kg)"] * BW
CL = K_EL * VD

MW_MG = LITERATURE_CONSTANTS["MW_Mg (g/mol)"]

D_ORAL = ASSUMPTIONS["Oral dose – elem. Mg (mg)"] / MW_MG
D_IV = ASSUMPTIONS["IV dose – elem. Mg (mg)"] / MW_MG
T_INF = ASSUMPTIONS["IV infusion duration (h)"]

C_BASE = LITERATURE_CONSTANTS["C_baseline (mmol/L)"]
TW_LO = LITERATURE_CONSTANTS["C_therapeutic_low (mmol/L)"]
TW_HI = LITERATURE_CONSTANTS["C_therapeutic_high (mmol/L)"]

F_CIT = LITERATURE_CONSTANTS["F_citrate"]
F_OX = LITERATURE_CONSTANTS["F_oxide"]
KA_CIT = LITERATURE_CONSTANTS["ka_citrate (h⁻¹)"]
KA_OX = LITERATURE_CONSTANTS["ka_oxide (h⁻¹)"]