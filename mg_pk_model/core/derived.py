import math
from .constants import LITERATURE_CONSTANTS, ASSUMPTIONS, TWO_COMP_CONSTANTS

BW     = ASSUMPTIONS["Body weight (kg)"]
T_HALF = ASSUMPTIONS["t½ elimination (h)"]

K_EL = math.log(2) / T_HALF
VD   = LITERATURE_CONSTANTS["Vd (L/kg)"] * BW
CL   = K_EL * VD

MW_MG = LITERATURE_CONSTANTS["MW_Mg (g/mol)"]

D_ORAL = ASSUMPTIONS["Oral dose – elem. Mg (mg)"] / MW_MG
D_IV   = ASSUMPTIONS["IV dose – elem. Mg (mg)"]   / MW_MG
T_INF  = ASSUMPTIONS["IV infusion duration (h)"]

C_BASE = LITERATURE_CONSTANTS["C_baseline (mmol/L)"]

F_CIT  = LITERATURE_CONSTANTS["F_citrate"]
F_OX   = LITERATURE_CONSTANTS["F_oxide"]
KA_CIT = LITERATURE_CONSTANTS["ka_citrate (h⁻¹)"]
KA_OX  = LITERATURE_CONSTANTS["ka_oxide (h⁻¹)"]

V1 = TWO_COMP_CONSTANTS["Vc (L/kg)"] * BW # central volume (L)
V2 = TWO_COMP_CONSTANTS["Vp (L/kg)"] * BW # peripheral volume (L)
Q  = TWO_COMP_CONSTANTS["Q (L/h/kg)"] * BW # intercompartmental clearance (L/h)

K12   = Q / V1           # central to peripheral
K21   = Q / V2           # peripheral to central
K10   = CL / V1          # elimination from central compartment

# hybrid (macro) rate constants alpha > beta
_sum  = K10 + K12 + K21
_disc = math.sqrt(_sum ** 2 - 4.0 * K21 * K10)
ALPHA = (_sum + _disc) / 2.0   # fast-distribution rate constant
BETA  = (_sum - _disc) / 2.0   # slow-elimination rate constant

# bolus coefficients  (A + B = 1 when expressed as fractions of D/V1)
A_2C = (ALPHA - K21) / (ALPHA - BETA)
B_2C = (K21  - BETA) / (ALPHA - BETA)