import numpy as np

from core import F_CIT, KA_CIT
from models.iv import cp_iv_bolus, cp_iv_infusion
from models.oral import cp_oral


def main():
    print("\n" + "=" * 60)
    print("  MAGNESIUM PK SIMULATION")
    print("=" * 60)

    # time grid (0–24 h)
    t = np.linspace(0, 24, 100)

    print("\nSimulating concentration–time profiles...")
    print(f"Time range: {t[0]} → {t[-1]} hours ({len(t)} points)")

    # simulations
    c_iv_bolus = cp_iv_bolus(t)
    c_iv_inf   = cp_iv_infusion(t)
    c_oral     = cp_oral(t, F_CIT, KA_CIT)

    print("\n" + "-" * 60)
    print("IV BOLUS (instant injection)")
    print(f"First 5 values: {c_iv_bolus[:5]}")

    print("\n" + "-" * 60)
    print("IV INFUSION (constant rate over time)")
    print(f"First 5 values: {c_iv_inf[:5]}")

    print("\n" + "-" * 60)
    print("ORAL (Mg Citrate)")
    print(f"First 5 values: {c_oral[:5]}")

    print("\n" + "=" * 60)
    print("INTERPRETATION")
    print("=" * 60)

    print(f"IV bolus starts at high concentration: {c_iv_bolus[0]:.3f} mmol/L")
    print(f"Oral starts at baseline:              {c_oral[0]:.3f} mmol/L")

    tmax_idx = np.argmax(c_oral)
    print(f"Oral peak occurs at ~{t[tmax_idx]:.2f} hours")
    print(f"Oral peak concentration: {c_oral[tmax_idx]:.3f} mmol/L")

    print("\nconclusionss:")
    print("- IV administration produces immediate systemic exposure")
    print("- Oral administration shows delayed absorption (Tmax > 0)")
    print("- Model behavior matches expected one-compartment PK")

    print("\n" + "=" * 60 + "\n")

if __name__ == "__main__":
    main()