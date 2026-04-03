import math
from .constants import LITERATURE_CONSTANTS, ASSUMPTIONS

# basic patient inputs
BW     = ASSUMPTIONS["Body weight (kg)"]
T_HALF = ASSUMPTIONS["t½ elimination (h)"]

# pk parameters
K_EL = math.log(2) / T_HALF
VD   = LITERATURE_CONSTANTS["Vd (L/kg)"] * BW
CL   = K_EL * VD

# molecular weight
MW_MG = LITERATURE_CONSTANTS["MW_Mg (g/mol)"]

# doses converted to mmol
D_ORAL = ASSUMPTIONS["Oral dose – elem. Mg (mg)"] / MW_MG
D_IV   = ASSUMPTIONS["IV dose – elem. Mg (mg)"]   / MW_MG
T_INF  = ASSUMPTIONS["IV infusion duration (h)"]

# concentration
C_BASE = LITERATURE_CONSTANTS["C_baseline (mmol/L)"]

# oral parameters
F_CIT  = LITERATURE_CONSTANTS["F_citrate"]
F_OX   = LITERATURE_CONSTANTS["F_oxide"]
KA_CIT = LITERATURE_CONSTANTS["ka_citrate (h⁻¹)"]
KA_OX  = LITERATURE_CONSTANTS["ka_oxide (h⁻¹)"]