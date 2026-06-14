import numpy as np

from core.derived import F_CIT, KA_CIT, F_OX, KA_OX
from core.constants import LITERATURE_CONSTANTS
from models.iv import cp_iv_bolus, cp_iv_infusion
from models.oral import cp_oral
from models.two_compartment import (
    cp_2c_iv_bolus, cp_2c_iv_infusion,
    cp_2c_oral_citrate, cp_2c_oral_oxide,
)
from metrics.calculate import calculate_pk_metrics, time_to_therapeutic
from visualization.plots import plot_pk_profiles

def _print_section(title: str) -> None:
    print(f"\n{'-' * 60}")
    print(f"  {title}")
    print(f"{'-' * 60}")

def _print_metrics(name: str, t: np.ndarray, cp: np.ndarray, thresh: float) -> None:
    metrics  = calculate_pk_metrics(t, cp)
    t_thera  = time_to_therapeutic(t, cp, thresh)
    t_str    = f"{t_thera:.2f} h" if t_thera is not None else "Never reached"
    print(f"\n  {name}:")
    print(f"    Cmax : {metrics['Cmax']:.4f} mmol/L")
    print(f"    Tmax : {metrics['Tmax']:.2f} h")
    print(f"    AUC  : {metrics['AUC']:.2f} mmol*h/L")
    print(f"    Time to >={thresh} mmol/L : {t_str}")


def main() -> None:
    print("MAGNESIUM PHARMACOKINETICS SIMULATION")
    print("1-Compartment  vs  2-Compartment models")

    t = np.linspace(0, 48, 500)
    thresh = LITERATURE_CONSTANTS["C_therapeutic_low (mmol/L)"]

    c1_iv_bolus   = cp_iv_bolus(t)
    c1_iv_inf     = cp_iv_infusion(t)
    c1_oral_cit   = cp_oral(t, F_CIT, KA_CIT)
    c1_oral_ox    = cp_oral(t, F_OX,  KA_OX)

    _print_section("1C metrics")
    scenarios_1c = {
        "IV bolus (1C)":        c1_iv_bolus,
        "IV infusion (1C)":     c1_iv_inf,
        "Oral citrate (1C)":  c1_oral_cit,
        "Oral citrate (1C)":    c1_oral_ox,
    }
    for name, cp in scenarios_1c.items():
        _print_metrics(name, t, cp, thresh)

    c2_iv_bolus  = cp_2c_iv_bolus(t)
    c2_iv_inf    = cp_2c_iv_infusion(t)
    c2_oral_cit  = cp_2c_oral_citrate(t)
    c2_oral_ox   = cp_2c_oral_oxide(t)

    _print_section("2C metrics")
    scenarios_2c = {
        "IV bolus (2C)":        c2_iv_bolus,
        "IV infusion (2C)":     c2_iv_inf,
        "Oral citrate (2C)":    c2_oral_cit,
        "Oral oxide (2C)":      c2_oral_ox,
    }
    for name, cp in scenarios_2c.items():
        _print_metrics(name, t, cp, thresh)

    _print_section("Analysis:")
    c_low  = LITERATURE_CONSTANTS['C_therapeutic_low (mmol/L)']
    c_high = LITERATURE_CONSTANTS['C_therapeutic_high (mmol/L)']
    print(
        "\n  The two-compartment model shows a biphasic decline:\n"
        "  a rapid alpha-phase (redistribution into peripheral tissues)\n"
        "  followed by the slower beta-phase (terminal elimination).\n"
        "  The 2C IV bolus peaks HIGHER than the 1C model because dose\n"
        "  initially occupies the smaller central volume V1 before\n"
        "  redistributing into V2.\n"
        "\n  IV Infusion remains the safest route -- both models agree that\n"
        f"  a 4-hour infusion keeps plasma Mg within the therapeutic band\n"
        f"  ({c_low}--{c_high} mmol/L).\n"
        f"\n  Mg Citrate (F={F_CIT}) is the preferred oral formulation;\n"
        f"  Mg Oxide (F={F_OX}) barely shifts baseline and is not\n"
        "  recommended for acute repletion.\n"
    )
    plot_pk_profiles(
        t,
        c1_iv_bolus, c1_iv_inf, c1_oral_cit, c1_oral_ox,
        c2_iv_bolus, c2_iv_inf, c2_oral_cit, c2_oral_ox,
    )

if __name__ == "__main__":
    main()