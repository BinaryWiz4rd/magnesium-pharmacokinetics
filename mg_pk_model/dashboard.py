import sys
import os
import math
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from simulation.monte_carlo import run_monte_carlo

MW_MG  = 24.305
C_BASE = 0.85
C_LOW  = 0.85
C_HIGH = 1.10
KA_CIT = 0.90
KA_OX  = 0.35

st.set_page_config(page_title="Magnesium PK Simulation", layout="wide")

st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Nunito:wght@400;700;900&display=swap');
    html, body, [class*="css"] { font-family: 'Nunito', sans-serif !important; }
    .stApp { background-color: #FFFAFA; }
    [data-testid="stSidebar"] {
        background-color: #FFFDE7;
        border-right: 3px dashed #FFCDD2;
    }
    h1, h2, h3, h4 { color: #FF8A80 !important; font-weight: 800 !important; }
    p, label, .stMarkdown { color: #8D6E63 !important; font-weight: 700 !important; }
    .title-box {
        background-color: #FFE082;
        border-radius: 20px; padding: 20px; text-align: center;
        box-shadow: 0px 5px 15px rgba(255,182,193,0.4); margin-bottom: 25px;
    }
    .summary-box {
        background-color: #E1F5FE; border: 2px dashed #81D4FA;
        border-radius: 20px; padding: 20px; margin-top: 25px;
        box-shadow: 0px 5px 15px rgba(129,212,250,0.2);
    }
    .mc-box {
        background-color: #F3E5F5; border: 2px dashed #CE93D8;
        border-radius: 20px; padding: 20px; margin-top: 20px;
    }
    div[data-testid="column"] {
        background-color: #FFFFFF; border: 2px solid #FFCDD2;
        border-radius: 20px; padding: 15px;
        box-shadow: 0px 8px 15px rgba(255,182,193,0.2);
        text-align: center; transition: transform 0.2s ease-in-out;
    }
    div[data-testid="column"]:hover { transform: translateY(-5px); }
    </style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="title-box">
    <h1 style="color:#FFFFFF !important; margin-bottom:0px; font-size:2.5em;">
        Magnesium PK Simulator
    </h1>
    <p style="color:#8D6E63 !important; font-size:1.1em; margin-top:5px;">
        1-Compartment · 2-Compartment · Monte Carlo
    </p>
</div>
""", unsafe_allow_html=True)

# ── Shared sidebar controls ───────────────────────────────────────────────────
st.sidebar.header("Patient data")
bw     = st.sidebar.slider("Body weight (kg)",      40.0, 120.0, 70.0, 1.0)
t_half = st.sidebar.slider("Elimination t½ (h)",    10.0,  50.0, 30.0, 1.0)

st.sidebar.markdown("---")
st.sidebar.header("IV dosing")
dose_iv = st.sidebar.slider("IV dose – elem. Mg (mg)", 100.0, 1000.0, 500.0, 50.0)
t_inf   = st.sidebar.slider("IV infusion duration (h)", 0.5, 12.0, 4.0, 0.5)

st.sidebar.markdown("---")
st.sidebar.header("oral dosing")
dose_oral_cit = st.sidebar.slider("Mg Citrate dose (mg)", 100.0, 1000.0, 300.0, 50.0)
f_cit         = st.sidebar.slider("Mg Citrate bioavailability", 0.10, 0.50, 0.31, 0.01)
st.sidebar.markdown("<br>", unsafe_allow_html=True)
dose_oral_ox  = st.sidebar.slider("Mg Oxide dose (mg)", 100.0, 1000.0, 300.0, 50.0)

st.sidebar.markdown("---")
st.sidebar.header("2-Compartment params")
vc_kg = st.sidebar.slider("central V (L/kg)",        0.10, 0.40, 0.20, 0.01)
vp_kg = st.sidebar.slider("peripheral V (L/kg)",     0.10, 0.60, 0.30, 0.01)
q_kg  = st.sidebar.slider("intercompart. Q (L/h/kg)", 0.02, 0.50, 0.12, 0.01)

k_el  = math.log(2) / t_half
vd    = 0.50 * bw
d_iv  = dose_iv       / MW_MG
d_cit = dose_oral_cit / MW_MG
d_ox  = dose_oral_ox  / MW_MG

t = np.linspace(0, 48, 500)


def sim_1c_iv_bolus():
    return C_BASE + (d_iv / vd) * np.exp(-k_el * t)


def sim_1c_iv_infusion():
    r   = d_iv / t_inf
    css = r / (k_el * vd)
    c_d = C_BASE + css * (1.0 - np.exp(-k_el * t))
    c_e = C_BASE + css * (1.0 - math.exp(-k_el * t_inf))
    c_a = C_BASE + (c_e - C_BASE) * np.exp(-k_el * (t - t_inf))
    return np.where(t <= t_inf, c_d, c_a)


def sim_1c_oral(d_mmol, f, ka):
    _ka = ka
    if abs(_ka - k_el) < 1e-6:
        _ka = k_el + 1e-5
    coef = (f * d_mmol * _ka) / (vd * (_ka - k_el))
    return C_BASE + coef * (np.exp(-k_el * t) - np.exp(-_ka * t))


v1   = vc_kg * bw
v2   = vp_kg * bw
q    = q_kg  * bw
cl   = k_el  * vd
k10  = cl / v1
k12  = q  / v1
k21  = q  / v2
_s   = k10 + k12 + k21
_d   = max(_s ** 2 - 4.0 * k21 * k10, 0.0)
alpha = (_s + math.sqrt(_d)) / 2.0
beta  = (_s - math.sqrt(_d)) / 2.0
_den  = alpha - beta if abs(alpha - beta) > 1e-12 else 1e-12
A2    = (alpha - k21) / _den
B2    = (k21  - beta) / _den

def sim_2c_iv_bolus():
    scale = d_iv / v1
    return C_BASE + scale * (A2 * np.exp(-alpha * t) + B2 * np.exp(-beta * t))

def sim_2c_iv_infusion():
    R0    = d_iv / t_inf
    scale = R0 / v1
    qa = (A2 / alpha) * (1.0 - math.exp(-alpha * t_inf))
    qb = (B2 / beta)  * (1.0 - math.exp(-beta  * t_inf))
    c_d = C_BASE + scale * (
        (A2 / alpha) * (1.0 - np.exp(-alpha * t)) +
        (B2 / beta)  * (1.0 - np.exp(-beta  * t))
    )
    c_a = C_BASE + scale * (
        qa * np.exp(-alpha * (t - t_inf)) +
        qb * np.exp(-beta  * (t - t_inf))
    )
    return np.where(t <= t_inf, c_d, c_a)

def sim_2c_oral(d_mmol, f, ka):
    _ka = ka
    if abs(_ka - alpha) < 1e-6: _ka = alpha + 1e-5
    if abs(_ka - beta)  < 1e-6: _ka = beta  + 1e-5
    P  = (k21 - alpha) / ((_ka - alpha) * (beta  - alpha))
    Qc = (k21 - beta)  / ((_ka - beta)  * (alpha - beta))
    R  = -(P + Qc)
    scale = f * d_mmol * _ka / v1
    return C_BASE + scale * (
        P  * np.exp(-alpha * t) +
        Qc * np.exp(-beta  * t) +
        R  * np.exp(-_ka   * t)
    )

c1_iv_bolus  = sim_1c_iv_bolus()
c1_iv_inf    = sim_1c_iv_infusion()
c1_oral_cit  = sim_1c_oral(d_cit, f_cit, KA_CIT)
c1_oral_ox   = sim_1c_oral(d_ox,  0.055, KA_OX)

c2_iv_bolus  = sim_2c_iv_bolus()
c2_iv_inf    = sim_2c_iv_infusion()
c2_oral_cit  = sim_2c_oral(d_cit, f_cit, KA_CIT)
c2_oral_ox   = sim_2c_oral(d_ox,  0.055, KA_OX)

def get_metrics(cp):
    idx = int(np.argmax(cp))
    return cp[idx], t[idx], float(np.trapezoid(cp, t))

def _style_ax(ax):
    ax.set_facecolor("#FFFFFF")
    ax.spines[["top", "right"]].set_color("none")
    ax.spines[["bottom", "left"]].set_linewidth(2)
    ax.spines[["bottom", "left"]].set_color("#D7CCC8")
    ax.tick_params(colors="#8D6E63", width=2)
    ax.set_xlabel("Time (hours)", fontsize=12, color="#8D6E63", fontweight="bold")
    ax.set_ylabel("Plasma Mg (mmol/L)", fontsize=12, color="#8D6E63", fontweight="bold")
    ax.axhspan(C_LOW, C_HIGH, color="#FFF59D", alpha=0.4, label=f"Therapeutic ({C_LOW}–{C_HIGH})")
    ax.grid(True, alpha=0.5, color="#EFEBE9", linestyle="--")
    ax.legend(frameon=True, facecolor="#FFFFFF", edgecolor="#FFE082",
              labelcolor="#8D6E63", prop={"weight": "bold"})


def _metric_cols(label_list, cp_list):
    cols = st.columns(len(label_list))
    for col, label, cp in zip(cols, label_list, cp_list):
        cmax, tmax, auc = get_metrics(cp)
        with col:
            st.markdown(f"#### **{label}**")
            st.write(f"**Cmax:** {cmax:.3f} mmol/L")
            st.write(f"**Tmax:** {tmax:.2f} h")
            st.write(f"**AUC:** {auc:.2f} mmol·h/L")
tab1, tab2, tab3 = st.tabs([
    "1-Compartment",
    "2-Compartment",
    "🎲 Monte Carlo",
])

# ═══════════════════════════════════════════════════════════════════════════════
# TAB 1 — one-compartment
# ═══════════════════════════════════════════════════════════════════════════════
with tab1:
    st.markdown("<h3 style='text-align:center;'>1-Compartment concentration–time profiles</h3>",
                unsafe_allow_html=True)

    fig1, ax1 = plt.subplots(figsize=(10, 4.5))
    fig1.patch.set_alpha(0.0)
    ax1.plot(t, c1_iv_bolus,  label="IV bolus",
             ls="--",  color="#FF8A80", lw=3.5)
    ax1.plot(t, c1_iv_inf,    label=f"IV infusion ({t_inf:.0f}h)",
             ls="-.",  color="#FFD180", lw=3.5)
    ax1.plot(t, c1_oral_cit,  label=f"Oral citrate ({dose_oral_cit:.0f}mg)",
             color="#80D8FF", lw=3.5)
    ax1.plot(t, c1_oral_ox,   label=f"Oral oxide ({dose_oral_ox:.0f}mg)",
             color="#B388FF", lw=3.5)
    _style_ax(ax1)
    st.pyplot(fig1)

    st.markdown("<br><h3 style='text-align:center;'>key metrics</h3>", unsafe_allow_html=True)
    _metric_cols(
        ["IV Bolus", f"IV Infusion ({t_inf:.0f}h)", "Mg Citrate", "Mg Oxide"],
        [c1_iv_bolus, c1_iv_inf, c1_oral_cit, c1_oral_ox],
    )

    st.markdown(f"""
    <div class="summary-box">
        <h4 style="margin-top:0px; color:#0288D1 !important;">Clinical analysis - 1-Compartment</h4>
        <ul style="margin-bottom:0px;">
            <li><b>IV bolus</b> peaks instantly at
                <b>{np.max(c1_iv_bolus):.3f} mmol/L</b>.
                If this exceeds the yellow band, transient hypermagnesemia is possible.</li>
            <li><b>IV infusion</b> ({t_inf:.0f} h) gives a smooth rise - the body has time
                to distribute the dose before peak is reached.</li>
            <li><b>Mg Citrate</b> (F = {f_cit:.2f}) is the best oral option for repletion.
                <b>Mg Oxide</b> (F = 0.055) barely shifts the curve - try "cranking"
                its dose to 1000 mg and see!</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

with tab2:
    st.markdown("<h3 style='text-align:center;'>2-Compartment concentration–time profiles</h3>",
                unsafe_allow_html=True)

    beta_half = math.log(2) / beta if beta > 0 else float("inf")

    col_info1, col_info2, col_info3 = st.columns(3)
    with col_info1:
        st.metric("α (dist. rate)", f"{alpha:.4f} h⁻¹")
    with col_info2:
        st.metric("β (elim. rate)", f"{beta:.4f} h⁻¹")
    with col_info3:
        st.metric("β-phase t½", f"{beta_half:.1f} h")

    st.markdown("<br>", unsafe_allow_html=True)

    fig2, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    fig2.patch.set_alpha(0.0)

    ax_iv = axes[0]
    ax_iv.set_title("IV routes: 1C vs 2C", color="#8D6E63", fontweight="bold")
    ax_iv.plot(t, c1_iv_bolus,  label="Bolus 1-comp",     ls="--",  color="#FF8A80", lw=2.5, alpha=0.7)
    ax_iv.plot(t, c2_iv_bolus,  label="Bolus 2-comp",     ls="--",  color="#D50000", lw=2.5)
    ax_iv.plot(t, c1_iv_inf,    label="Infusion 1-comp",  ls="-.",  color="#FFD180", lw=2.5, alpha=0.7)
    ax_iv.plot(t, c2_iv_inf,    label="Infusion 2-comp",  ls="-.",  color="#FF6D00", lw=2.5)
    _style_ax(ax_iv)

    ax_or = axes[1]
    ax_or.set_title("Oral routes: 1C vs 2C", color="#8D6E63", fontweight="bold")
    ax_or.plot(t, c1_oral_cit,  label="Citrate 1-comp",  color="#80D8FF", lw=2.5, alpha=0.7)
    ax_or.plot(t, c2_oral_cit,  label="Citrate 2-comp",  color="#0091EA", lw=2.5)
    ax_or.plot(t, c1_oral_ox,   label="Oxide 1-comp",    color="#B388FF", lw=2.5, alpha=0.7)
    ax_or.plot(t, c2_oral_ox,   label="Oxide 2-comp",    color="#6200EA", lw=2.5)
    _style_ax(ax_or)
    ax_or.set_ylabel("")

    plt.tight_layout()
    st.pyplot(fig2)

    st.markdown("<br><h3 style='text-align:center;'>2-Compartment metrics</h3>", unsafe_allow_html=True)
    _metric_cols(
        ["IV bolus (2C)", f"IV infusion (2C)", "Citrate (2C)", "Oxide (2C)"],
        [c2_iv_bolus, c2_iv_inf, c2_oral_cit, c2_oral_ox],
    )

    iv_bolus_1c_peak = float(np.max(c1_iv_bolus))
    iv_bolus_2c_peak = float(np.max(c2_iv_bolus))

    st.markdown(f"""
    <div class="summary-box">
        <h4 style="margin-top:0px; color:#0288D1 !important;">Clinical analysis - 2-Compartment</h4>
        <ul style="margin-bottom:0px;">
            <li><b>Biphasic decline:</b> the α-phase (α = {alpha:.3f} h⁻¹) captures rapid
                redistribution from plasma into peripheral tissues; the slower β-phase
                (t½ ≈ {beta_half:.1f} h) reflects true elimination.</li>
            <li><b>IV Bolus peak</b> rises from <b>{iv_bolus_1c_peak:.3f}</b> (1C) to
                <b>{iv_bolus_2c_peak:.3f} mmol/L</b> (2C) — because dose initially
                enters the smaller central volume V1 = {vc_kg:.2f} × {bw:.0f} = {v1:.0f} L,
                creating a higher transient peak before redistribution into V2.
                The 1-comp model underestimates this early spike.</li>
            <li><b>IV infusion</b> is still the safest strategy; the 2C model confirms
                that the 4-hour infusion keeps plasma Mg safely within the therapeutic
                window while avoiding the sharp bolus spike.</li>
            <li><b>Adjust the sidebar sliders</b> (Vc, Vp, Q) to explore how tissue
                binding and perfusion alter the biphasic shape.</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

with tab3:
    st.markdown("<h3 style='text-align:center;'>🎲 Monte Carlo (for population variability)</h3>",
                unsafe_allow_html=True)

    st.markdown("""
    Simulates a **virtual patient population** by sampling PK parameters from
    log-normal distributions.  The shaded band shows the 5th–95th percentile
    concentration range; the solid line is the population median.
    """)

    mc_col1, mc_col2 = st.columns([1, 2])

    with mc_col1:
        st.markdown("#### Simulation settings")
        mc_route = st.selectbox(
            "route / model",
            options=[
                "iv_infusion", "iv_bolus",
                "oral_cit",    "oral_ox",
                "2c_iv_bolus", "2c_iv_infusion", "2c_oral_cit",
            ],
            index=0,
            format_func=lambda x: {
                "iv_infusion":    "IV infusion (1-comp)",
                "iv_bolus":       "IV bolus (1-comp)",
                "oral_cit":       "Oral citrate (1-comp)",
                "oral_ox":        "Oral oxide (1-comp)",
                "2c_iv_bolus":    "IV bolus (2-comp)",
                "2c_iv_infusion": "IV infusion (2-comp)",
                "2c_oral_cit":    "Oral citrate (2-comp)",
            }[x],
        )
        n_sims  = st.slider("number of simulations", 100, 2000, 500, 100)
        mc_seed = st.number_input("random seed", 0, 9999, 42, 1)

        st.markdown("#### Inter-individual variability (CV)")
        cv_bw     = st.slider("body weight CV",     0.05, 0.40, 0.15, 0.01)
        cv_t_half = st.slider("half-life CV",       0.10, 0.60, 0.35, 0.01)
        cv_vd     = st.slider("Vd CV",              0.05, 0.40, 0.20, 0.01)
        cv_f      = st.slider("bioavailability CV", 0.05, 0.50, 0.25, 0.01) \
                    if "oral" in mc_route else 0.25
        cv_ka     = st.slider("ka CV",              0.05, 0.50, 0.30, 0.01) \
                    if "oral" in mc_route else 0.30

        run_mc = st.button("▶ run simulation", type="primary")

    with mc_col2:
        if run_mc or "mc_result" not in st.session_state:
            with st.spinner(f"running {n_sims} simulations…"):
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

        p5  = result["percentiles"][5]
        p50 = result["percentiles"][50]
        p95 = result["percentiles"][95]

        fig3, ax3 = plt.subplots(figsize=(9, 4.5))
        fig3.patch.set_alpha(0.0)

        for curve in result["curves"][::max(1, n_sims // 80)]:
            ax3.plot(t, curve, color="#B39DDB", lw=0.4, alpha=0.25)

        ax3.fill_between(t, p5, p95, color="#CE93D8", alpha=0.30,
                         label="5th–95th percentile")
        ax3.plot(t, p50, color="#6A1B9A", lw=3.0, label="Median (p50)")

        _style_ax(ax3)
        ax3.set_title(
            f"Monte Carlo  |  n = {n_sims}  |  route: {mc_route}",
            color="#8D6E63", fontweight="bold",
        )
        st.pyplot(fig3)

        cmax_all = result["curves"].max(axis=1)
        auc_all  = np.trapezoid(result["curves"], t, axis=1)
        frac_in_window = np.mean(
            (result["curves"] >= C_LOW) & (result["curves"] <= C_HIGH)
        ) * 100

        st.markdown("<br>", unsafe_allow_html=True)
        sc1, sc2, sc3, sc4 = st.columns(4)
        with sc1:
            st.metric("Median Cmax",   f"{np.median(cmax_all):.3f} mmol/L")
        with sc2:
            st.metric("Cmax 90% CI",
                      f"{np.percentile(cmax_all, 5):.3f} – {np.percentile(cmax_all, 95):.3f}")
        with sc3:
            st.metric("Median AUC",    f"{np.median(auc_all):.1f} mmol·h/L")
        with sc4:
            st.metric("% Time in window", f"{frac_in_window:.1f}%")

    st.markdown(f"""
    <div class="mc-box">
        <h4 style="margin-top:0px; color:#7B1FA2 !important;">Monte Carlo results</h4>
        <ul style="margin-bottom:0px;">
            <li>The <b>shaded band</b> captures 90% of the simulated population.
                A wide band means high sensitivity to PK variability - consider
                therapeutic drug monitoring for that route.</li>
            <li>The <b>% Time in Window</b> metric estimates what fraction of the
                simulated plasma-time profile falls within the therapeutic range
                ({C_LOW}–{C_HIGH} mmol/L).</li>
            <li>Increase the CV sliders to model sicker or more variable patients
                (e.g. renal impairment → higher t½ CV).</li>
            <li>Compare <b>1-comp vs 2-comp</b> routes to see how distributional
                kinetics widen the uncertainty band for IV bolus dosing.</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)