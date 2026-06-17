import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

sys.path.insert(0, os.path.dirname(__file__))
from renal_wasting import simulate_renal_wasting_patient
from simulation.monte_carlo import run_monte_carlo


MW_MG = 24.305
C_BASE = 0.85
C_LOW = 0.85
C_HIGH = 1.10
KA_CIT = 0.90
KA_OX = 0.35


st.set_page_config(page_title="Magnesium PK Simulation", layout="wide")

st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Nunito:wght@400;700;900&display=swap');

    html, body, [class*="css"] {
        font-family: 'Nunito', sans-serif !important;
    }

    .stApp {
        background-color: #FFFAFA;
    }

    [data-testid="stSidebar"] {
        background-color: #FFFDE7;
        border-right: 3px dashed #FFCDD2;
    }

    h1, h2, h3, h4 {
        color: #FF8A80 !important;
        font-weight: 800 !important;
    }

    p, label, .stMarkdown {
        color: #8D6E63 !important;
        font-weight: 700 !important;
    }

    .title-box {
        background-color: #FFE082;
        border-radius: 20px;
        padding: 20px;
        text-align: center;
        box-shadow: 0px 5px 15px rgba(255, 182, 193, 0.4);
        margin-bottom: 25px;
    }

    .summary-box {
        background-color: #E1F5FE;
        border: 2px dashed #81D4FA;
        border-radius: 20px;
        padding: 20px;
        margin-top: 25px;
        box-shadow: 0px 5px 15px rgba(129, 212, 250, 0.2);
    }

    .mc-box {
        background-color: #F3E5F5;
        border: 2px dashed #CE93D8;
        border-radius: 20px;
        padding: 20px;
        margin-top: 20px;
    }

    div[data-testid="column"] {
        background-color: #FFFFFF;
        border: 2px solid #FFCDD2;
        border-radius: 20px;
        padding: 15px;
        box-shadow: 0px 8px 15px rgba(255, 182, 193, 0.2);
        text-align: center;
        transition: transform 0.2s ease-in-out;
    }

    div[data-testid="column"]:hover {
        transform: translateY(-5px);
    }
    </style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="title-box">
    <h1 style="color: #FFFFFF !important; margin-bottom: 0px; font-size: 2.5em;">Magnesium PK Simulator</h1>
    <p style="color: #8D6E63 !important; font-size: 1.1em; margin-top: 5px;">
        1-Compartment | 2-Compartment | Monte Carlo | Genetic Renal Wasting
    </p>
</div>
""", unsafe_allow_html=True)

st.sidebar.header("Patient Data")
bw = st.sidebar.slider("Body weight (kg)", 40.0, 120.0, 70.0, 1.0)
t_half = st.sidebar.slider("Elimination half-life (h)", 10.0, 50.0, 30.0, 1.0)

st.sidebar.markdown("---")
st.sidebar.header("IV Dosing")
dose_iv = st.sidebar.slider("IV dose - elemental Mg (mg)", 100.0, 1000.0, 500.0, 50.0)
t_inf = st.sidebar.slider("IV infusion duration (h)", 0.5, 12.0, 4.0, 0.5)

st.sidebar.markdown("---")
st.sidebar.header("Oral Dosing")
dose_oral_cit = st.sidebar.slider("Mg citrate dose (mg)", 100.0, 1000.0, 300.0, 50.0)
f_cit = st.sidebar.slider("Mg citrate bioavailability", 0.10, 0.50, 0.31, 0.01)
st.sidebar.markdown("<br>", unsafe_allow_html=True)
dose_oral_ox = st.sidebar.slider("Mg oxide dose (mg)", 100.0, 1000.0, 300.0, 50.0)

st.sidebar.markdown("---")
st.sidebar.header("2-Compartment Params")
vc_kg = st.sidebar.slider("Central V (L/kg)", 0.10, 0.40, 0.20, 0.01)
vp_kg = st.sidebar.slider("Peripheral V (L/kg)", 0.10, 0.60, 0.30, 0.01)
q_kg = st.sidebar.slider("Intercompartmental Q (L/h/kg)", 0.02, 0.50, 0.12, 0.01)

k_el = math.log(2) / t_half
vd = 0.50 * bw
d_iv = dose_iv / MW_MG
d_cit = dose_oral_cit / MW_MG
d_ox = dose_oral_ox / MW_MG

t = np.linspace(0, 48, 500)


def sim_1c_iv_bolus():
    return C_BASE + (d_iv / vd) * np.exp(-k_el * t)


def sim_1c_iv_infusion():
    r = d_iv / t_inf
    css = r / (k_el * vd)
    c_during = C_BASE + css * (1.0 - np.exp(-k_el * t))
    c_end = C_BASE + css * (1.0 - math.exp(-k_el * t_inf))
    c_after = C_BASE + (c_end - C_BASE) * np.exp(-k_el * (t - t_inf))
    return np.where(t <= t_inf, c_during, c_after)


def sim_1c_oral(d_mmol, f, ka):
    adjusted_ka = ka
    if abs(adjusted_ka - k_el) < 1e-6:
        adjusted_ka = k_el + 1e-5
    coef = (f * d_mmol * adjusted_ka) / (vd * (adjusted_ka - k_el))
    return C_BASE + coef * (np.exp(-k_el * t) - np.exp(-adjusted_ka * t))


v1 = vc_kg * bw
v2 = vp_kg * bw
q = q_kg * bw
cl = k_el * vd
k10 = cl / v1
k12 = q / v1
k21 = q / v2
sum_k = k10 + k12 + k21
disc = max(sum_k**2 - 4.0 * k21 * k10, 0.0)
alpha = (sum_k + math.sqrt(disc)) / 2.0
beta = (sum_k - math.sqrt(disc)) / 2.0
den = alpha - beta if abs(alpha - beta) > 1e-12 else 1e-12
A2 = (alpha - k21) / den
B2 = (k21 - beta) / den


def sim_2c_iv_bolus():
    scale = d_iv / v1
    return C_BASE + scale * (A2 * np.exp(-alpha * t) + B2 * np.exp(-beta * t))


def sim_2c_iv_infusion():
    r0 = d_iv / t_inf
    scale = r0 / v1
    q_alpha = (A2 / alpha) * (1.0 - math.exp(-alpha * t_inf))
    q_beta = (B2 / beta) * (1.0 - math.exp(-beta * t_inf))
    c_during = C_BASE + scale * (
        (A2 / alpha) * (1.0 - np.exp(-alpha * t))
        + (B2 / beta) * (1.0 - np.exp(-beta * t))
    )
    c_after = C_BASE + scale * (
        q_alpha * np.exp(-alpha * (t - t_inf))
        + q_beta * np.exp(-beta * (t - t_inf))
    )
    return np.where(t <= t_inf, c_during, c_after)


def sim_2c_oral(d_mmol, f, ka):
    adjusted_ka = ka
    if abs(adjusted_ka - alpha) < 1e-6:
        adjusted_ka = alpha + 1e-5
    if abs(adjusted_ka - beta) < 1e-6:
        adjusted_ka = beta + 1e-5

    p_coeff = (k21 - alpha) / ((adjusted_ka - alpha) * (beta - alpha))
    q_coeff = (k21 - beta) / ((adjusted_ka - beta) * (alpha - beta))
    r_coeff = -(p_coeff + q_coeff)
    scale = f * d_mmol * adjusted_ka / v1
    return C_BASE + scale * (
        p_coeff * np.exp(-alpha * t)
        + q_coeff * np.exp(-beta * t)
        + r_coeff * np.exp(-adjusted_ka * t)
    )


def trapezoid_auc(concentration, time, axis=None):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(concentration, time, axis=axis)
    return np.trapz(concentration, time, axis=axis)


def get_metrics(cp):
    idx = int(np.argmax(cp))
    return cp[idx], t[idx], float(trapezoid_auc(cp, t))


def style_ax(ax):
    ax.set_facecolor("#FFFFFF")
    ax.spines[["top", "right"]].set_color("none")
    ax.spines[["bottom", "left"]].set_linewidth(2)
    ax.spines[["bottom", "left"]].set_color("#D7CCC8")
    ax.tick_params(colors="#8D6E63", width=2)
    ax.set_xlabel("Time (hours)", fontsize=12, color="#8D6E63", fontweight="bold")
    ax.set_ylabel("Plasma Mg concentration (mmol/L)", fontsize=12, color="#8D6E63", fontweight="bold")
    ax.axhspan(C_LOW, C_HIGH, color="#FFF59D", alpha=0.4, label=f"Therapeutic window ({C_LOW}-{C_HIGH} mmol/L)")
    ax.grid(True, alpha=0.5, color="#EFEBE9", linestyle="--")
    ax.legend(frameon=True, facecolor="#FFFFFF", edgecolor="#FFE082", labelcolor="#8D6E63", prop={"weight": "bold"})


def metric_cols(label_list, cp_list):
    cols = st.columns(len(label_list))
    for col, label, cp in zip(cols, label_list, cp_list):
        cmax, tmax, auc = get_metrics(cp)
        with col:
            st.markdown(f"#### **{label}**")
            st.write(f"**Cmax:** {cmax:.3f} mmol/L")
            st.write(f"**Tmax:** {tmax:.2f} h")
            st.write(f"**AUC:** {auc:.2f} mmol*h/L")


c1_iv_bolus = sim_1c_iv_bolus()
c1_iv_inf = sim_1c_iv_infusion()
c1_oral_cit = sim_1c_oral(d_cit, f_cit, KA_CIT)
c1_oral_ox = sim_1c_oral(d_ox, 0.055, KA_OX)

c2_iv_bolus = sim_2c_iv_bolus()
c2_iv_inf = sim_2c_iv_infusion()
c2_oral_cit = sim_2c_oral(d_cit, f_cit, KA_CIT)
c2_oral_ox = sim_2c_oral(d_ox, 0.055, KA_OX)

tab1, tab2, tab3, tab4 = st.tabs([
    "1-Compartment",
    "2-Compartment",
    "Monte Carlo",
    "Genetic Renal Wasting Analysis",
])

with tab1:
    st.markdown("<h3 style='text-align:center;'>1-Compartment Concentration-Time Profiles</h3>", unsafe_allow_html=True)

    fig1, ax1 = plt.subplots(figsize=(10, 4.5))
    fig1.patch.set_alpha(0.0)
    ax1.plot(t, c1_iv_bolus, label="IV bolus", linestyle="--", color="#FF8A80", linewidth=3.5)
    ax1.plot(t, c1_iv_inf, label=f"IV infusion ({t_inf:.0f} h)", linestyle="-.", color="#FFD180", linewidth=3.5)
    ax1.plot(t, c1_oral_cit, label=f"Oral citrate ({dose_oral_cit:.0f} mg)", color="#80D8FF", linewidth=3.5)
    ax1.plot(t, c1_oral_ox, label=f"Oral oxide ({dose_oral_ox:.0f} mg)", color="#B388FF", linewidth=3.5)
    style_ax(ax1)
    st.pyplot(fig1)

    st.markdown("<br><h3 style='text-align:center;'>Key Metrics</h3>", unsafe_allow_html=True)
    metric_cols(
        ["IV Bolus", f"IV Infusion ({t_inf:.0f} h)", "Mg Citrate", "Mg Oxide"],
        [c1_iv_bolus, c1_iv_inf, c1_oral_cit, c1_oral_ox],
    )

    st.markdown(f"""
    <div class="summary-box">
        <h4 style="margin-top: 0px; color: #0288D1 !important;">Clinical Analysis - 1-Compartment</h4>
        <ul style="margin-bottom: 0px;">
            <li><b>IV bolus</b> peaks instantly at <b>{np.max(c1_iv_bolus):.3f} mmol/L</b>. If this exceeds the yellow band, transient hypermagnesemia is possible.</li>
            <li><b>IV infusion</b> ({t_inf:.0f} h) creates a smoother rise and gives the dose more time to distribute.</li>
            <li><b>Mg citrate</b> (F = {f_cit:.2f}) produces a stronger oral response than <b>Mg oxide</b> (F = 0.055).</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

with tab2:
    st.markdown("<h3 style='text-align:center;'>2-Compartment Concentration-Time Profiles</h3>", unsafe_allow_html=True)

    beta_half = math.log(2) / beta if beta > 0 else float("inf")
    col_info1, col_info2, col_info3 = st.columns(3)
    with col_info1:
        st.metric("Alpha distribution rate", f"{alpha:.4f} h^-1")
    with col_info2:
        st.metric("Beta elimination rate", f"{beta:.4f} h^-1")
    with col_info3:
        st.metric("Beta-phase half-life", f"{beta_half:.1f} h")

    st.markdown("<br>", unsafe_allow_html=True)

    fig2, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    fig2.patch.set_alpha(0.0)

    ax_iv = axes[0]
    ax_iv.set_title("IV routes: 1C vs 2C", color="#8D6E63", fontweight="bold")
    ax_iv.plot(t, c1_iv_bolus, label="Bolus 1-comp", linestyle="--", color="#FF8A80", linewidth=2.5, alpha=0.7)
    ax_iv.plot(t, c2_iv_bolus, label="Bolus 2-comp", linestyle="--", color="#D50000", linewidth=2.5)
    ax_iv.plot(t, c1_iv_inf, label="Infusion 1-comp", linestyle="-.", color="#FFD180", linewidth=2.5, alpha=0.7)
    ax_iv.plot(t, c2_iv_inf, label="Infusion 2-comp", linestyle="-.", color="#FF6D00", linewidth=2.5)
    style_ax(ax_iv)

    ax_oral = axes[1]
    ax_oral.set_title("Oral routes: 1C vs 2C", color="#8D6E63", fontweight="bold")
    ax_oral.plot(t, c1_oral_cit, label="Citrate 1-comp", color="#80D8FF", linewidth=2.5, alpha=0.7)
    ax_oral.plot(t, c2_oral_cit, label="Citrate 2-comp", color="#0091EA", linewidth=2.5)
    ax_oral.plot(t, c1_oral_ox, label="Oxide 1-comp", color="#B388FF", linewidth=2.5, alpha=0.7)
    ax_oral.plot(t, c2_oral_ox, label="Oxide 2-comp", color="#6200EA", linewidth=2.5)
    style_ax(ax_oral)
    ax_oral.set_ylabel("")

    plt.tight_layout()
    st.pyplot(fig2)

    st.markdown("<br><h3 style='text-align:center;'>2-Compartment Metrics</h3>", unsafe_allow_html=True)
    metric_cols(
        ["IV Bolus (2C)", "IV Infusion (2C)", "Citrate (2C)", "Oxide (2C)"],
        [c2_iv_bolus, c2_iv_inf, c2_oral_cit, c2_oral_ox],
    )

    iv_bolus_1c_peak = float(np.max(c1_iv_bolus))
    iv_bolus_2c_peak = float(np.max(c2_iv_bolus))

    st.markdown(f"""
    <div class="summary-box">
        <h4 style="margin-top: 0px; color: #0288D1 !important;">Clinical Analysis - 2-Compartment</h4>
        <ul style="margin-bottom: 0px;">
            <li><b>Biphasic decline:</b> alpha captures rapid redistribution into peripheral tissue, while beta reflects slower elimination.</li>
            <li><b>IV bolus peak</b> rises from <b>{iv_bolus_1c_peak:.3f}</b> (1C) to <b>{iv_bolus_2c_peak:.3f} mmol/L</b> (2C), because the initial central volume is smaller.</li>
            <li><b>Adjust Vc, Vp, and Q</b> in the sidebar to explore how distribution changes the curve shape.</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

with tab3:
    st.markdown("<h3 style='text-align:center;'>Monte Carlo Population Variability</h3>", unsafe_allow_html=True)

    st.markdown("""
    Simulates a virtual patient population by sampling PK parameters from log-normal distributions.
    The shaded band shows the 5th-95th percentile concentration range; the solid line is the population median.
    """)

    mc_col1, mc_col2 = st.columns([1, 2])

    with mc_col1:
        st.markdown("#### Simulation Settings")
        mc_route = st.selectbox(
            "Route / model",
            options=[
                "iv_infusion",
                "iv_bolus",
                "oral_cit",
                "oral_ox",
                "2c_iv_bolus",
                "2c_iv_infusion",
                "2c_oral_cit",
            ],
            index=0,
            format_func=lambda route: {
                "iv_infusion": "IV infusion (1-comp)",
                "iv_bolus": "IV bolus (1-comp)",
                "oral_cit": "Oral citrate (1-comp)",
                "oral_ox": "Oral oxide (1-comp)",
                "2c_iv_bolus": "IV bolus (2-comp)",
                "2c_iv_infusion": "IV infusion (2-comp)",
                "2c_oral_cit": "Oral citrate (2-comp)",
            }[route],
        )
        n_sims = st.slider("Number of simulations", 100, 2000, 500, 100)
        mc_seed = st.number_input("Random seed", 0, 9999, 42, 1)

        st.markdown("#### Inter-Individual Variability (CV)")
        cv_bw = st.slider("Body weight CV", 0.05, 0.40, 0.15, 0.01)
        cv_t_half = st.slider("Half-life CV", 0.10, 0.60, 0.35, 0.01)
        cv_vd = st.slider("Vd CV", 0.05, 0.40, 0.20, 0.01)
        cv_f = st.slider("Bioavailability CV", 0.05, 0.50, 0.25, 0.01) if "oral" in mc_route else 0.25
        cv_ka = st.slider("ka CV", 0.05, 0.50, 0.30, 0.01) if "oral" in mc_route else 0.30

        run_mc = st.button("Run simulation", type="primary")

    with mc_col2:
        if run_mc or "mc_result" not in st.session_state:
            with st.spinner(f"Running {n_sims} simulations..."):
                result = run_monte_carlo(
                    n_sims=n_sims,
                    route=mc_route,
                    t=t,
                    bw_mean=bw,
                    t_half_mean=t_half,
                    dose_iv_mg=dose_iv,
                    dose_oral_mg=dose_oral_cit if "cit" in mc_route else dose_oral_ox,
                    t_inf=t_inf,
                    f_cit=f_cit,
                    vc_kg_mean=vc_kg,
                    vp_kg_mean=vp_kg,
                    q_kg_mean=q_kg,
                    cv_bw=cv_bw,
                    cv_t_half=cv_t_half,
                    cv_vd=cv_vd,
                    cv_f=cv_f,
                    cv_ka=cv_ka,
                    seed=int(mc_seed),
                )
                st.session_state["mc_result"] = result
        else:
            result = st.session_state["mc_result"]

        p5 = result["percentiles"][5]
        p50 = result["percentiles"][50]
        p95 = result["percentiles"][95]

        fig3, ax3 = plt.subplots(figsize=(9, 4.5))
        fig3.patch.set_alpha(0.0)

        for curve in result["curves"][::max(1, n_sims // 80)]:
            ax3.plot(t, curve, color="#B39DDB", linewidth=0.4, alpha=0.25)

        ax3.fill_between(t, p5, p95, color="#CE93D8", alpha=0.30, label="5th-95th percentile")
        ax3.plot(t, p50, color="#6A1B9A", linewidth=3.0, label="Median (p50)")
        style_ax(ax3)
        ax3.set_title(f"Monte Carlo | n = {n_sims} | route: {mc_route}", color="#8D6E63", fontweight="bold")
        st.pyplot(fig3)

        cmax_all = result["curves"].max(axis=1)
        auc_all = trapezoid_auc(result["curves"], t, axis=1)
        frac_in_window = np.mean((result["curves"] >= C_LOW) & (result["curves"] <= C_HIGH)) * 100

        st.markdown("<br>", unsafe_allow_html=True)
        sc1, sc2, sc3, sc4 = st.columns(4)
        with sc1:
            st.metric("Median Cmax", f"{np.median(cmax_all):.3f} mmol/L")
        with sc2:
            st.metric("Cmax 90% CI", f"{np.percentile(cmax_all, 5):.3f} - {np.percentile(cmax_all, 95):.3f}")
        with sc3:
            st.metric("Median AUC", f"{np.median(auc_all):.1f} mmol*h/L")
        with sc4:
            st.metric("% Time in window", f"{frac_in_window:.1f}%")

    st.markdown(f"""
    <div class="mc-box">
        <h4 style="margin-top: 0px; color: #7B1FA2 !important;">Monte Carlo Results</h4>
        <ul style="margin-bottom: 0px;">
            <li>The <b>shaded band</b> captures 90% of the simulated population.</li>
            <li><b>% Time in Window</b> estimates how often the simulated profiles fall within {C_LOW}-{C_HIGH} mmol/L.</li>
            <li>Increase CV sliders to model a more variable population.</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

with tab4:
    st.header("Genetic Tubulopathies & Renal Magnesium Wasting")
    st.markdown("""
    This section explores how inherited tubulopathies can change magnesium pharmacokinetics.
    The model compares a healthy patient with Gitelman syndrome, Bartter syndrome type 4,
    EAST/SeSAME syndrome, and CNNM2-related renal magnesium wasting.
    """)

    renal_col1, renal_col2 = st.columns(2)
    with renal_col1:
        chosen_route = st.selectbox(
            "Administration route",
            ["IV Infusion", "IV Bolus", "Oral Citrate", "Oral Oxide"],
            key="renal_route",
        )
        chosen_dose = st.slider(
            "Elemental magnesium dose (mg)",
            min_value=100,
            max_value=1000,
            value=300,
            step=50,
            key="renal_dose",
        )

    with renal_col2:
        patient_weight = st.slider(
            "Patient weight (kg)",
            min_value=40,
            max_value=120,
            value=int(bw),
            key="renal_weight",
        )
        inf_duration = st.slider(
            "IV infusion duration (h)",
            min_value=1,
            max_value=12,
            value=max(1, int(round(t_inf))),
            key="renal_infusion_duration",
        )

    fig4, ax4 = plt.subplots(figsize=(10, 5.5))
    fig4.patch.set_alpha(0.0)
    ax4.set_facecolor("#FFFFFF")

    conditions = [
        "Healthy",
        "Gitelman Syndrome",
        "Bartter Syndrome Type 4",
        "EAST/SeSAME Syndrome",
        "CNNM2-Related",
    ]
    colors = ["#80D8FF", "#FF8A80", "#B388FF", "#FFD180", "#81C784"]

    for condition, color in zip(conditions, colors):
        renal_t, c_plasma = simulate_renal_wasting_patient(
            condition=condition,
            route=chosen_route,
            dose_mg=chosen_dose,
            duration_h=inf_duration,
            weight_kg=patient_weight,
        )
        is_healthy = condition == "Healthy"
        ax4.plot(
            renal_t,
            c_plasma,
            label=condition,
            color=color,
            linewidth=2.0 if is_healthy else 3.0,
            linestyle="--" if is_healthy else "-",
        )

    ax4.set_title(
        f"Dynamic Profiling: {chosen_route} ({chosen_dose} mg elemental Mg)",
        fontsize=13,
        color="#8D6E63",
        fontweight="bold",
        pad=10,
    )
    ax4.set_xlim(0, 48)
    ax4.set_xticks(range(0, 49, 4))
    style_ax(ax4)
    st.pyplot(fig4)

    st.markdown("""
    <div class="summary-box">
        <h4 style="margin-top: 0px; color: #0288D1 !important;">Clinical Interpretation</h4>
        <p>
            Patients with genetic renal magnesium wasting return toward their lower disease-specific baseline
            faster than the healthy comparison curve. In this simulation, oral magnesium oxide has the weakest
            effect because its assumed bioavailability is only about 5.5%.
        </p>
    </div>
    """, unsafe_allow_html=True)
