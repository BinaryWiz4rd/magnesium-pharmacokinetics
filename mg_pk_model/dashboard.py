import math

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from renal_wasting import simulate_renal_wasting_patient


MW_MG = 24.305
C_BASE = 0.85
VD_L_KG = 0.50
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
    <p style="color: #8D6E63 !important; font-size: 1.1em; margin-top: 5px;">Use the controls to compare magnesium dosing strategies in real time.</p>
</div>
""", unsafe_allow_html=True)

st.sidebar.header("Patient Basics")
bw = st.sidebar.slider("Body weight (kg)", min_value=40.0, max_value=120.0, value=70.0, step=1.0)
t_half = st.sidebar.slider("Elimination half-life (h)", min_value=10.0, max_value=50.0, value=30.0, step=1.0)

st.sidebar.markdown("---")

st.sidebar.header("IV Dosing")
dose_iv = st.sidebar.slider("IV dose - elemental Mg (mg)", min_value=100.0, max_value=1000.0, value=500.0, step=50.0)
t_inf = st.sidebar.slider("IV infusion duration (h)", min_value=0.5, max_value=12.0, value=4.0, step=0.5)

st.sidebar.markdown("---")

st.sidebar.header("Oral Dosing")
dose_oral_cit = st.sidebar.slider("Mg citrate dose (mg)", min_value=100.0, max_value=1000.0, value=300.0, step=50.0)
f_cit = st.sidebar.slider("Mg citrate bioavailability", min_value=0.10, max_value=0.50, value=0.31, step=0.01)

st.sidebar.markdown("<br>", unsafe_allow_html=True)

dose_oral_ox = st.sidebar.slider("Mg oxide dose (mg)", min_value=100.0, max_value=1000.0, value=300.0, step=50.0)

vd = VD_L_KG * bw
k_el = math.log(2) / t_half
d_iv_mmol = dose_iv / MW_MG
d_oral_cit_mmol = dose_oral_cit / MW_MG
d_oral_ox_mmol = dose_oral_ox / MW_MG

t = np.linspace(0, 48, 500)


def sim_iv_bolus():
    dc0 = d_iv_mmol / vd
    return C_BASE + dc0 * np.exp(-k_el * t)


def sim_iv_infusion():
    r = d_iv_mmol / t_inf
    css = r / (k_el * vd)
    c_during = C_BASE + css * (1 - np.exp(-k_el * t))
    c_end = C_BASE + css * (1 - math.exp(-k_el * t_inf))
    c_after = C_BASE + (c_end - C_BASE) * np.exp(-k_el * (t - t_inf))
    return np.where(t <= t_inf, c_during, c_after)


def sim_oral(d_mmol, f, ka):
    coef = (f * d_mmol * ka) / (vd * (ka - k_el))
    return C_BASE + coef * (np.exp(-k_el * t) - np.exp(-ka * t))


def trapezoid_auc(concentration, time):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(concentration, time)
    return np.trapz(concentration, time)


def get_metrics(cp):
    idx = np.argmax(cp)
    auc = trapezoid_auc(cp, t)
    return cp[idx], t[idx], auc


c_iv_bolus = sim_iv_bolus()
c_iv_inf = sim_iv_infusion()
c_oral_cit = sim_oral(d_oral_cit_mmol, f_cit, KA_CIT)
c_oral_ox = sim_oral(d_oral_ox_mmol, 0.055, KA_OX)

standard_tab, renal_tab = st.tabs(["Standard PK Profiles", "Genetic Renal Wasting Analysis"])

with standard_tab:
    with st.container():
        st.markdown("<h3 style='text-align: center;'>Concentration-Time Profiles</h3>", unsafe_allow_html=True)

        fig, ax = plt.subplots(figsize=(10, 4.5))
        fig.patch.set_alpha(0.0)
        ax.set_facecolor("#FFFFFF")

        ax.plot(t, c_iv_bolus, label="IV bolus", linestyle="--", color="#FF8A80", linewidth=3.5)
        ax.plot(t, c_iv_inf, label=f"IV infusion ({t_inf} h)", linestyle="-.", color="#FFD180", linewidth=3.5)
        ax.plot(t, c_oral_cit, label=f"Oral Mg citrate ({dose_oral_cit:.0f} mg)", linewidth=3.5, color="#80D8FF")
        ax.plot(t, c_oral_ox, label=f"Oral Mg oxide ({dose_oral_ox:.0f} mg)", linewidth=3.5, color="#B388FF")

        ax.axhspan(C_LOW, C_HIGH, color="#FFF59D", alpha=0.4, label=f"Therapeutic window ({C_LOW}-{C_HIGH} mmol/L)")

        ax.set_xlabel("Time (hours)", fontsize=12, color="#8D6E63", fontweight="bold")
        ax.set_ylabel("Plasma concentration (mmol/L)", fontsize=12, color="#8D6E63", fontweight="bold")
        ax.tick_params(colors="#8D6E63", width=2)

        ax.spines["bottom"].set_linewidth(2)
        ax.spines["left"].set_linewidth(2)
        ax.spines["bottom"].set_color("#D7CCC8")
        ax.spines["left"].set_color("#D7CCC8")
        ax.spines["top"].set_color("none")
        ax.spines["right"].set_color("none")

        ax.legend(frameon=True, facecolor="#FFFFFF", edgecolor="#FFE082", labelcolor="#8D6E63", prop={"weight": "bold"})
        ax.grid(True, alpha=0.6, color="#EFEBE9", linestyle="--")

        st.pyplot(fig)

    st.markdown("<br><h3 style='text-align: center;'>Key Metrics</h3>", unsafe_allow_html=True)

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.markdown("#### **IV Bolus**")
        cmax, _, auc = get_metrics(c_iv_bolus)
        st.write(f"**Cmax:** {cmax:.3f} mmol/L")
        st.write(f"**AUC:** {auc:.2f} mmol*h/L")

    with col2:
        st.markdown("#### **IV Infusion**")
        cmax, tmax, auc = get_metrics(c_iv_inf)
        st.write(f"**Cmax:** {cmax:.3f} mmol/L")
        st.write(f"**Tmax:** {tmax:.2f} h")
        st.write(f"**AUC:** {auc:.2f} mmol*h/L")

    with col3:
        st.markdown("#### **Mg Citrate**")
        cmax, tmax, auc = get_metrics(c_oral_cit)
        st.write(f"**Cmax:** {cmax:.3f} mmol/L")
        st.write(f"**Tmax:** {tmax:.2f} h")
        st.write(f"**AUC:** {auc:.2f} mmol*h/L")

    with col4:
        st.markdown("#### **Mg Oxide**")
        cmax, tmax, auc = get_metrics(c_oral_ox)
        st.write(f"**Cmax:** {cmax:.3f} mmol/L")
        st.write(f"**Tmax:** {tmax:.2f} h")
        st.write(f"**AUC:** {auc:.2f} mmol*h/L")

    st.markdown(f"""
    <div class="summary-box">
        <h4 style="margin-top: 0px; color: #0288D1 !important;">Quick Clinical Takeaways</h4>
        <ul style="margin-bottom: 0px;">
            <li><b>Immediate risk:</b> IV bolus peaks instantly at <b>{np.max(c_iv_bolus):.2f} mmol/L</b>. If this line crosses the top of the yellow band, there may be a transient hypermagnesemia risk.</li>
            <li><b>Sustained care:</b> The {t_inf}-hour IV infusion creates a smoother curve, giving the dose more time to distribute.</li>
            <li><b>Bioavailability matters:</b> Magnesium oxide has very low absorption, which explains why the curve moves much less than magnesium citrate at similar elemental doses.</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

with renal_tab:
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

    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_alpha(0.0)
    ax.set_facecolor("#FFFFFF")

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
        ax.plot(
            renal_t,
            c_plasma,
            label=condition,
            color=color,
            linewidth=2.0 if is_healthy else 3.0,
            linestyle="--" if is_healthy else "-",
        )

    ax.axhspan(C_LOW, C_HIGH, color="#FFF59D", alpha=0.4, label=f"Therapeutic window ({C_LOW}-{C_HIGH} mmol/L)")
    ax.set_title(
        f"Dynamic Profiling: {chosen_route} ({chosen_dose} mg elemental Mg)",
        fontsize=13,
        color="#8D6E63",
        fontweight="bold",
        pad=10,
    )
    ax.set_xlabel("Time after administration (hours)", fontsize=12, color="#8D6E63", fontweight="bold")
    ax.set_ylabel("Plasma Mg concentration (mmol/L)", fontsize=12, color="#8D6E63", fontweight="bold")
    ax.set_xlim(0, 48)
    ax.set_xticks(range(0, 49, 4))
    ax.tick_params(colors="#8D6E63", width=2)

    ax.spines["bottom"].set_linewidth(2)
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_color("#D7CCC8")
    ax.spines["left"].set_color("#D7CCC8")
    ax.spines["top"].set_color("none")
    ax.spines["right"].set_color("none")

    ax.legend(frameon=True, facecolor="#FFFFFF", edgecolor="#FFE082", labelcolor="#8D6E63", prop={"weight": "bold"})
    ax.grid(True, alpha=0.6, color="#EFEBE9", linestyle="--")

    st.pyplot(fig)

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
