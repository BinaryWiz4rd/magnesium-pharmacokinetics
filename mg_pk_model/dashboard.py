import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math

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

    /* Apply rounded font everywhere */
    html, body, [class*="css"] {
        font-family: 'Nunito', sans-serif !important;
    }

    /* Main application background - Soft cream */
    .stApp {
        background-color: #FFFAFA;
    }

    /* Sidebar styling */
    [data-testid="stSidebar"] {
        background-color: #FFFDE7; 
        border-right: 3px dashed #FFCDD2;
    }

    /* Text styling */
    h1, h2, h3, h4 {
        color: #FF8A80 !important;
        font-weight: 800 !important;
    }
    p, label, .stMarkdown {
        color: #8D6E63 !important; 
        font-weight: 700 !important;
    }

    /* --- CUSTOM BOXES (CARDS) --- */

    /* Title Banner */
    .title-box {
        background-color: #FFE082;
        border-radius: 20px;
        padding: 20px;
        text-align: center;
        box-shadow: 0px 5px 15px rgba(255, 182, 193, 0.4);
        margin-bottom: 25px;
    }

    /* Summary / Takeaway Box */
    .summary-box {
        background-color: #E1F5FE; /* Baby blue */
        border: 2px dashed #81D4FA;
        border-radius: 20px;
        padding: 20px;
        margin-top: 25px;
        box-shadow: 0px 5px 15px rgba(129, 212, 250, 0.2);
    }

    /* Streamlit's native columns (used for metrics) */
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
    <h1 style="color: #FFFFFF !important; margin-bottom: 0px; font-size: 2.5em;">🧸 Magnesium PK Simulator 🧸</h1>
    <p style="color: #8D6E63 !important; font-size: 1.1em; margin-top: 5px;">Slide the handles in the menu to see how different doses change the curves! ✨</p>
</div>
""", unsafe_allow_html=True)

st.sidebar.header("🎀 Patient Basics")
bw = st.sidebar.slider("Body weight (kg)", min_value=40.0, max_value=120.0, value=70.0, step=1.0)
t_half = st.sidebar.slider("Elimination t½ (h)", min_value=10.0, max_value=50.0, value=30.0, step=1.0)

st.sidebar.markdown("---")

st.sidebar.header("💉 IV Dosing")
dose_iv = st.sidebar.slider("IV Dose - elem. Mg (mg)", min_value=100.0, max_value=1000.0, value=500.0, step=50.0)
t_inf = st.sidebar.slider("IV Infusion duration (h)", min_value=0.5, max_value=12.0, value=4.0, step=0.5)

st.sidebar.markdown("---")

st.sidebar.header("💊 Oral Dosing")
dose_oral_cit = st.sidebar.slider("Mg Citrate Dose (mg)", min_value=100.0, max_value=1000.0, value=300.0, step=50.0)
f_cit = st.sidebar.slider("Mg Citrate bioavailability", min_value=0.10, max_value=0.50, value=0.31, step=0.01)

st.sidebar.markdown("<br>", unsafe_allow_html=True) # little spacer

dose_oral_ox = st.sidebar.slider("Mg Oxide Dose (mg)", min_value=100.0, max_value=1000.0, value=300.0, step=50.0)

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


c_iv_bolus = sim_iv_bolus()
c_iv_inf = sim_iv_infusion()
c_oral_cit = sim_oral(d_oral_cit_mmol, f_cit, KA_CIT)
c_oral_ox = sim_oral(d_oral_ox_mmol, 0.055, KA_OX)

with st.container():
    st.markdown("<h3 style='text-align: center;'>☁️ Concentration-Time Profiles ☁️</h3>", unsafe_allow_html=True)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    fig.patch.set_alpha(0.0)
    ax.set_facecolor('#FFFFFF')

    ax.plot(t, c_iv_bolus, label='IV Bolus', linestyle='--', color='#FF8A80', linewidth=3.5)
    ax.plot(t, c_iv_inf, label=f'IV Infusion ({t_inf}h)', linestyle='-.', color='#FFD180', linewidth=3.5)
    ax.plot(t, c_oral_cit, label=f'Oral Mg Citrate ({dose_oral_cit}mg)', linewidth=3.5, color='#80D8FF')
    ax.plot(t, c_oral_ox, label=f'Oral Mg Oxide ({dose_oral_ox}mg)', linewidth=3.5, color='#B388FF')

    ax.axhspan(C_LOW, C_HIGH, color='#FFF59D', alpha=0.4, label=f'Therapeutic Window ({C_LOW}-{C_HIGH})')

    ax.set_xlabel('Time (hours)', fontsize=12, color='#8D6E63', fontweight='bold')
    ax.set_ylabel('Plasma Concentration (mmol/L)', fontsize=12, color='#8D6E63', fontweight='bold')
    ax.tick_params(colors='#8D6E63', width=2)

    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_color('#D7CCC8')
    ax.spines['left'].set_color('#D7CCC8')
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')

    ax.legend(frameon=True, facecolor='#FFFFFF', edgecolor='#FFE082', labelcolor='#8D6E63', prop={'weight': 'bold'})
    ax.grid(True, alpha=0.6, color='#EFEBE9', linestyle='--')

    st.pyplot(fig)

st.markdown("<br><h3 style='text-align: center;'>✨ Key Metrics ✨</h3>", unsafe_allow_html=True)


def get_metrics(cp):
    idx = np.argmax(cp)
    auc = np.trapezoid(cp, t)
    return cp[idx], t[idx], auc


col1, col2, col3, col4 = st.columns(4)

with col1:
    st.markdown("#### **IV Bolus**")
    cmax, tmax, auc = get_metrics(c_iv_bolus)
    st.write(f"**Cmax:** {cmax:.3f} mmol/L")
    st.write(f"**AUC:** {auc:.2f} mmol·h/L")

with col2:
    st.markdown(f"#### **IV Infusion**")
    cmax, tmax, auc = get_metrics(c_iv_inf)
    st.write(f"**Cmax:** {cmax:.3f} mmol/L")
    st.write(f"**Tmax:** {tmax:.2f} hrs")
    st.write(f"**AUC:** {auc:.2f} mmol·h/L")

with col3:
    st.markdown("#### **Mg Citrate**")
    cmax, tmax, auc = get_metrics(c_oral_cit)
    st.write(f"**Cmax:** {cmax:.3f} mmol/L")
    st.write(f"**Tmax:** {tmax:.2f} hrs")
    st.write(f"**AUC:** {auc:.2f} mmol·h/L")

with col4:
    st.markdown("#### **Mg Oxide**")
    cmax, tmax, auc = get_metrics(c_oral_ox)
    st.write(f"**Cmax:** {cmax:.3f} mmol/L")
    st.write(f"**Tmax:** {tmax:.2f} hrs")
    st.write(f"**AUC:** {auc:.2f} mmol·h/L")

st.markdown(f"""
<div class="summary-box">
    <h4 style="margin-top: 0px; color: #0288D1 !important;">🩺 Quick Clinical Takeaways</h4>
    <ul style="margin-bottom: 0px;">
        <li><b>Immediate Risk:</b> The IV Bolus peaks instantly at <b>{np.max(c_iv_bolus):.2f} mmol/L</b>. If this line crosses the top of the yellow box, there is a risk of transient hypermagnesemia!</li>
        <li><b>Sustained Care:</b> The {t_inf}-hour IV Infusion provides a much smoother curve, allowing the body time to distribute the dose safely.</li>
        <li><b>Bioavailability Matters:</b> Look at the purple line! Magnesium Oxide's extremely low absorption shows exactly why it is a poor choice for quickly correcting a deficiency compared to Citrate. Try cranking the Oxide dose up to 1000mg and see how little the curve actually moves!</li>
    </ul>
</div>
""", unsafe_allow_html=True)