LITERATURE_CONSTANTS = {
    "C_baseline (mmol/L)":          0.85,
    "C_therapeutic_low (mmol/L)":   0.85,
    "C_therapeutic_high (mmol/L)":  1.10,
    "Vd (L/kg)":                    0.50,
    "F_citrate":                    0.31,
    "F_oxide":                      0.055,
    "ka_citrate (h⁻¹)":             0.90,
    "ka_oxide (h⁻¹)":               0.35,
    "MW_Mg (g/mol)":               24.305,
}

ASSUMPTIONS = {
    "Body weight (kg)":             70.0,
    "t½ elimination (h)":           30.0,
    "Oral dose – elem. Mg (mg)":   300.0,
    "IV dose – elem. Mg (mg)":     500.0,
    "IV infusion duration (h)":      4.0,
    "Model type":                  "One-compartment, first-order",
    "F_IV":                         1.0,
}

TWO_COMP_CONSTANTS = {
    "Vc (L/kg)":   0.20,   # central volume of distribution
    "Vp (L/kg)":   0.30,   # peripheral volume of distribution
    "Q (L/h/kg)":  0.12,   # intercompartmental clearance
}