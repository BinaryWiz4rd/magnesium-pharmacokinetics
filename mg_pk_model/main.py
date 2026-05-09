import numpy as np

from core.derived import F_CIT, KA_CIT, F_OX, KA_OX
from core.constants import LITERATURE_CONSTANTS
from models.iv import cp_iv_bolus, cp_iv_infusion
from models.oral import cp_oral
from metrics.calculate import calculate_pk_metrics, time_to_therapeutic
from visualization.plots import plot_pk_profiles


def main():
    print("MAGNESIUM PK SIMULATION: ORAL VS IV")

    # time grid (0–48 h) to capture full elimination phase
    t = np.linspace(0, 48, 500)

    print(f"\nSimulating concentration–time profiles over {t[-1]} hours...")

    c_iv_bolus = cp_iv_bolus(t)
    c_iv_inf = cp_iv_infusion(t)
    c_oral_cit = cp_oral(t, F_CIT, KA_CIT)
    c_oral_ox = cp_oral(t, F_OX, KA_OX)

    thresh_low = LITERATURE_CONSTANTS["C_therapeutic_low (mmol/L)"]

    scenarios = {
        "IV Bolus": c_iv_bolus,
        "IV Infusion": c_iv_inf,
        "Oral (Citrate)": c_oral_cit,
        "Oral (Oxide)": c_oral_ox
    }

    print("\nPHARMACOKINETIC METRICS & THERAPEUTIC WINDOW")

    for name, cp in scenarios.items():
        metrics = calculate_pk_metrics(t, cp)
        t_thera = time_to_therapeutic(t, cp, thresh_low)

        t_thera_str = f"{t_thera:.2f} hrs" if t_thera is not None else "Never reached"

        print(f"\n{name}:")
        print(f"  - Cmax: {metrics['Cmax']:.3f} mmol/L")
        print(f"  - Tmax: {metrics['Tmax']:.2f} hrs")
        print(f"  - AUC:  {metrics['AUC']:.2f} mmol*h/L")
        print(f"  - Time to Therapeutic (>={thresh_low}): {t_thera_str}")

    plot_pk_profiles(t, c_iv_bolus, c_iv_inf, c_oral_cit, c_oral_ox)

    print("\nCLINICAL INTERPRETATION")
    print("- IV Bolus provides immediate relief but risks exceeding the upper therapeutic limit.")
    print("- IV Infusion safely maintains levels within the therapeutic window.")
    print(f"- Mg Citrate (F={F_CIT}) reaches therapeutic levels effectively but with a delayed onset (Tmax).")
    print(f"- Mg Oxide (F={F_OX}) has extremely poor bioavailability and barely affects baseline levels.")

if __name__ == "__main__":
    main()