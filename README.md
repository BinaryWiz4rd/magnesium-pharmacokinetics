# Magnesium PK Simulator

A modular Python tool for simulating magnesium concentration-time profiles using a **one-compartment model** with **first-order kinetics**.

This project computes systemic magnesium levels over time based on literature-derived pharmacokinetic (PK) parameters, accounting for baseline endogenous magnesium concentrations. 

## Supported Models
* **Intravenous (IV):** * **Bolus:** Instantaneous distribution.
  * **Infusion:** Constant rate infusion over a specified duration (default: 4 hours).
* **Oral:** First-order absorption modeling comparing different formulations:
  * **Magnesium Citrate:** Higher bioavailability (F = 0.31, ka = 0.90 h⁻¹).
  * **Magnesium Oxide:** Lower bioavailability (F = 0.055, ka = 0.35 h⁻¹).

## Project Structure
The repository is structured to separate configuration, mathematical models, and future analysis/reporting features.

```text
Pharmacokinetics/
├── .venv/                      # Virtual environment (recommended)
├── mg_pk_model/                # Core simulation package
│   ├── core/                   # Parameters and state configuration
│   │   ├── __init__.py
│   │   ├── constants.py        # Literature values & patient assumptions
│   │   └── derived.py          # Derived PK values (Kel, Vd, CL, Doses in mmol)
│   ├── metrics/                # (WIP) PK metrics calculation (Cmax, AUC, Tmax)
│   ├── models/                 # Mathematical PK equations
│   │   ├── iv.py               # IV bolus and infusion math
│   │   └── oral.py             # Oral absorption math
│   ├── reporting/              # (WIP) Exporting results to CSV/PDF
│   ├── utils/                  # (WIP) Helper functions
│   ├── visualization/          # (WIP) Matplotlib/Plotly plotting scripts
│   └── main.py                 # Simulation runner and CLI entry point
├── mg_real_analysis/           # (WIP) Scripts for analyzing real-world/clinical data
├── .gitignore                  
└── README.md                   
```

## Requirements
This project requires **Python 3.x** and relies on `numpy` for efficient array-based time-series calculations.

Install the required dependency via pip:
```bash
pip install numpy
```

## How to Run It

To execute the default 24-hour simulation for IV Bolus, IV Infusion, and Oral (Magnesium Citrate) administration, run the main script from the root directory:

```bash
python mg_pk_model/main.py
```

### Expected Output
The script will output a text-based report to your console, detailing the early time-series data and key pharmacokinetic interpretations:

```text
============================================================
  MAGNESIUM PK SIMULATION
============================================================

Simulating concentration–time profiles...
Time range: 0.0 → 24.0 hours (100 points)

------------------------------------------------------------
IV BOLUS (instant injection)
First 5 values: [1.43724391 1.42436853 1.41177699 1.39946467 1.38742709]
...

============================================================
INTERPRETATION
============================================================
IV bolus starts at high concentration: 1.437 mmol/L
Oral starts at baseline:              0.850 mmol/L
Oral peak occurs at ~3.15 hours
Oral peak concentration: 0.975 mmol/L

conclusionss:
- IV administration produces immediate systemic exposure
- Oral administration shows delayed absorption (Tmax > 0)
- Model behavior matches expected one-compartment PK
```

## Under the Hood: Key Parameters & Assumptions
The model runs on default assumptions defined in `core/constants.py`, which can be easily modified for different patient profiles:

* **Patient Weight:** 70 kg
* **Baseline Mg Concentration:** 0.85 mmol/L
* **Volume of Distribution (Vd):** 0.50 L/kg
* **Elimination Half-life (t½):** 30.0 hours
* **Default Doses:** 300 mg elemental Mg (Oral), 500 mg elemental Mg (IV)
* *Note: All doses are automatically converted to mmol inside `derived.py` using the molecular weight of Magnesium (24.305 g/mol).*