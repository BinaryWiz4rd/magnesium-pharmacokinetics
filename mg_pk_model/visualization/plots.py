import matplotlib.pyplot as plt
import numpy as np
from core.constants import LITERATURE_CONSTANTS


def plot_pk_profiles(
    t,
    c1_iv_bolus, c1_iv_inf, c1_oral_cit, c1_oral_ox,
    c2_iv_bolus=None, c2_iv_inf=None, c2_oral_cit=None, c2_oral_ox=None,
):
    c_low  = LITERATURE_CONSTANTS["C_therapeutic_low (mmol/L)"]
    c_high = LITERATURE_CONSTANTS["C_therapeutic_high (mmol/L)"]
    two_comp = c2_iv_bolus is not None

    if two_comp:
        fig, (ax_iv, ax_or) = plt.subplots(1, 2, figsize=(14, 6), sharey=False)
        fig.suptitle(
            "Magnesium PK: 1-Compartment vs 2-Compartment",
            fontsize=15, fontweight="bold",
        )

        ax_iv.plot(t, c1_iv_bolus, label="Bolus – 1C",     ls="--", color="#E53935", lw=2.0, alpha=0.75)
        ax_iv.plot(t, c2_iv_bolus, label="Bolus – 2C",     ls="--", color="#B71C1C", lw=2.5)
        ax_iv.plot(t, c1_iv_inf,   label="Infusion – 1C",  ls="-.", color="#FB8C00", lw=2.0, alpha=0.75)
        ax_iv.plot(t, c2_iv_inf,   label="Infusion – 2C",  ls="-.", color="#E65100", lw=2.5)
        ax_iv.axhspan(c_low, c_high, color="green", alpha=0.12,
                      label=f"Therapeutic ({c_low}–{c_high})")
        ax_iv.set_title("IV routes", fontsize=13)
        ax_iv.set_xlabel("Time (hours)", fontsize=11)
        ax_iv.set_ylabel("Plasma Mg (mmol/L)", fontsize=11)
        ax_iv.legend(fontsize=9)
        ax_iv.grid(True, alpha=0.3)

        ax_or.plot(t, c1_oral_cit, label="Citrate – 1C",  color="#1E88E5", lw=2.0, alpha=0.75)
        ax_or.plot(t, c2_oral_cit, label="Citrate – 2C",  color="#0D47A1", lw=2.5)
        ax_or.plot(t, c1_oral_ox,  label="Oxide – 1C",    color="#8E24AA", lw=2.0, alpha=0.75)
        ax_or.plot(t, c2_oral_ox,  label="Oxide – 2C",    color="#4A148C", lw=2.5)
        ax_or.axhspan(c_low, c_high, color="green", alpha=0.12,
                      label=f"Therapeutic ({c_low}–{c_high})")
        ax_or.set_title("Oral routes", fontsize=13)
        ax_or.set_xlabel("Time (hours)", fontsize=11)
        ax_or.set_ylabel("Plasma Mg (mmol/L)", fontsize=11)
        ax_or.legend(fontsize=9)
        ax_or.grid(True, alpha=0.3)

    else:
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.suptitle(
            "Magnesium concentration-time profiles: Oral vs IV",
            fontsize=14,
        )
        ax.plot(t, c1_iv_bolus,  label="IV bolus",          ls="--", color="red")
        ax.plot(t, c1_iv_inf,    label="IV infusion (4h)",   ls="-.", color="orange")
        ax.plot(t, c1_oral_cit,  label="Oral (Mg citrate)",  lw=2,    color="blue")
        ax.plot(t, c1_oral_ox,   label="Oral (Mg oxide)",    lw=2,    color="purple")
        ax.axhspan(c_low, c_high, color="green", alpha=0.15,
                   label=f"Therapeutic window ({c_low}–{c_high} mmol/L)")
        ax.set_xlabel("Time (hours)", fontsize=12)
        ax.set_ylabel("Plasma concentration (mmol/L)", fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("pk_simulation_plot.png", dpi=300)
    print("\nplot saved as 'pk_simulation_plot.png'")
