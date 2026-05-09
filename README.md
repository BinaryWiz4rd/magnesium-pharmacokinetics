# 🧸 Magnesium PK Simulator 🧸

A beautifully designed, interactive Python tool for simulating magnesium concentration-time profiles using a one-compartment model with first-order kinetics. 

<img width="1908" height="858" alt="obraz" src="https://github.com/user-attachments/assets/fba76075-a898-41ec-95de-5d29e82f3755" />

This project computes systemic magnesium levels over time based on literature-derived pharmacokinetic (PK) parameters, accounting for baseline endogenous magnesium concentrations. It features a fully interactive web dashboard built with Streamlit, allowing users to adjust patient parameters and dosing regimens in real-time.

### Key Features
* **Interactive Dashboard**: Real-time sliders for Body Weight, Elimination Half-Life, IV Dosing, and Oral Dosing.
* **Instant Visualizations**: Dynamic concentration-time curves mapping drug levels against the therapeutic window (0.85 - 1.10 mmol/L).
* **Automated PK Metrics**: Instantly calculates Cmax, Tmax, and AUC for all administration routes.
* **Clinical Insights**: Generates quick clinical takeaways based on the simulated curves.

##  Supported Models
* **Intravenous (IV)**: 
  * **Bolus**: Instantaneous distribution.
  * **Infusion**: Constant rate infusion over a dynamically adjustable duration.
* **Oral**: First-order absorption modeling comparing different formulations:
  * **Magnesium Citrate**: Higher bioavailability (Default F = 0.31, ka = 0.90 h⁻¹).
  * **Magnesium Oxide**: Lower bioavailability (Default F = 0.055, ka = 0.35 h⁻¹).

## Project Structure
The repository is structured to separate configuration, mathematical models, and the interactive frontend.

    Pharmacokinetics/
    ├── .venv/                      # Virtual environment (recommended)
    ├── mg_pk_model/                # Core simulation package
    │   ├── core/                   # Parameters and state configuration
    │   │   ├── __init__.py
    │   │   ├── constants.py        # Literature values & patient assumptions
    │   │   └── derived.py          # Derived PK values (Kel, Vd, CL, Doses in mmol)
    │   ├── metrics/                # PK metrics calculations (Cmax, AUC, Tmax)
    │   ├── models/                 # Mathematical PK equations
    │   │   ├── iv.py               # IV bolus and infusion math
    │   │   └── oral.py             # Oral absorption math
    │   ├── visualization/          # Matplotlib plotting scripts
    │   ├── main.py                 # Static CLI simulation runner 
    │   └── dashboard.py            # INTERACTIVE STREAMLIT WEB APP
    ├── mg_real_analysis/           # (WIP) Scripts for analyzing real-world/clinical data
    ├── .gitignore                  
    └── README.md                   

## Requirements
This project requires Python 3.8+ and relies on NumPy (>= 2.0.0), Matplotlib, and Streamlit. 

Install the required dependencies via pip:

    pip install numpy matplotlib streamlit

*(Note: NumPy 2.0+ is required for the np.trapezoid function used in AUC calculations).*

## How to Run It

**Option 1: The Interactive Web Dashboard (Recommended)**
To launch the interactive GUI with sliders, real-time graphs, and metric calculation, navigate to the mg_pk_model folder and run:

    streamlit run mg_pk_model/dashboard.py

This will open a local web server (usually at http://localhost:8501) in your default web browser.

**Option 2: The Static CLI Script**
If you prefer to run a standard static simulation that outputs metrics to the console and saves a .png plot, run:

    python main.py

## Under the Hood: Key Parameters & Assumptions
The math engine runs on default literature constants defined in core/constants.py, but you can easily modify these dynamically using the dashboard sliders:

* **Patient Weight:** 40 - 120 kg (Default: 70 kg)
* **Baseline Mg Concentration:** 0.85 mmol/L
* **Volume of Distribution (Vd):** 0.50 L/kg
* **Elimination Half-life (t½):** 10 - 50 hours (Default: 30 hours)
* **Doses:** Fully adjustable via sliders (converted automatically to mmol using the MW of Mg: 24.305 g/mol).
