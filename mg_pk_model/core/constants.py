import math

LITERATURE_CONSTANTS = {
    "C_baseline (mmol/L)": 0.85,
    "C_therapeutic_low (mmol/L)": 0.85,
    "C_therapeutic_high (mmol/L)": 1.10,
    "Vd (L/kg)": 0.50,
    "F_citrate": 0.31,
    "F_oxide": 0.055,
    "ka_citrate (h⁻¹)": 0.90,
    "ka_oxide (h⁻¹)": 0.35,
    "MW_Mg (g/mol)": 24.305,
}

ASSUMPTIONS = {
    "Body weight (kg)": 70.0,
    "t½ elimination (h)": 30.0,
    "Oral dose – elem. Mg (mg)": 300.0,
    "IV dose – elem. Mg (mg)": 500.0,
    "IV infusion duration (h)": 4.0,
    "F_IV": 1.0,
}