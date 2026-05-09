import matplotlib.pyplot as plt
import numpy as np
from core.constants import LITERATURE_CONSTANTS


def plot_pk_profiles(t, c_iv_bolus, c_iv_inf, c_oral_cit, c_oral_ox):
    """
    Plots the PK profiles for all administration routes (with reference therapeutic window)
    """
    plt.figure(figsize=(10, 6))

    # dosage simulations
    plt.plot(t, c_iv_bolus, label='IV Bolus', linestyle='--', color='red')
    plt.plot(t, c_iv_inf, label='IV Infusion (4h)', linestyle='-.', color='orange')
    plt.plot(t, c_oral_cit, label='Oral (Mg Citrate)', linewidth=2, color='blue')
    plt.plot(t, c_oral_ox, label='Oral (Mg Oxide)', linewidth=2, color='purple')

    c_low = LITERATURE_CONSTANTS["C_therapeutic_low (mmol/L)"]
    c_high = LITERATURE_CONSTANTS["C_therapeutic_high (mmol/L)"]
    plt.axhspan(c_low, c_high, color='green', alpha=0.15,
                label=f'Therapeutic Window ({c_low}-{c_high} mmol/L)')

    # Formatting
    plt.title('Magnesium concentration-time profiles: Oral vs IV', fontsize=14)
    plt.xlabel('time (hours)', fontsize=12)
    plt.ylabel('plasma concentration (mmol/L)', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig("pk_simulation_plot.png", dpi=300)
    print("\n[Visualization] Plot saved as 'pk_simulation_plot.png'")