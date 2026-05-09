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

st.set_page_config(page_title="Magnesium PK simulation", layout="wide")
st.title("INTERACTIVE TOOL FOR MAGNESIUM PK")
st.markdown("Adjust the patient and dosing parameters in the sidebar to see real-time updates to the PK profile!!")

st.sidebar.header("patient parameters")
bw = st.sidebar.slider("body weight (kg)", min_value=40.0, max_value=120.0, value=70.0, step=1.0)
t_half = st.sidebar.slider("elimination t½ (h)", min_value=10.0, max_value=50.0, value=30.0, step=1.0)

st.sidebar.header("IV Dosing")
dose_iv = st.sidebar.slider("IV Dose - elem. Mg (mg)", min_value=100.0, max_value=1000.0, value=500.0, step=50.0)
t_inf = st.sidebar.slider("IV Infusion duration (h)", min_value=0.5, max_value=12.0, value=4.0, step=0.5)

st.sidebar.header("Oral dosing")
dose_oral = st.sidebar.slider("Oral Dose - elem. Mg (mg)", min_value=100.0, max_value=1000.0, value=300.0, step=50.0)
f_cit = st.sidebar.slider("Mg Citrate bioavailability (F)", min_value=0.10, max_value=0.50, value=0.31, step=0.01)

vd = VD_L_KG * bw
k_el = math.log(2) / t_half
d_iv_mmol = dose_iv / MW_MG
d_oral_mmol = dose_oral / MW_MG

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

def sim_oral(f, ka):
    coef = (f * d_oral_mmol * ka) / (vd * (ka - k_el))
    return C_BASE + coef * (np.exp(-k_el * t) - np.exp(-ka * t))

c_iv_bolus = sim_iv_bolus()
c_iv_inf = sim_iv_infusion()
c_oral_cit = sim_oral(f_cit, KA_CIT)
c_oral_ox = sim_oral(0.055, KA_OX) # keeping Oxide F fixed at 0.055 for comparison

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(t, c_iv_bolus, label='IV Bolus', linestyle='--', color='red')
ax.plot(t, c_iv_inf, label=f'IV Infusion ({t_inf}h)', linestyle='-.', color='orange')
ax.plot(t, c_oral_cit, label=f'Oral Mg Citrate (F={f_cit})', linewidth=2, color='blue')
ax.plot(t, c_oral_ox, label='Oral Mg Oxide (F=0.055)', linewidth=2, color='purple')

# Therapeutic Window
ax.axhspan(C_LOW, C_HIGH, color='green', alpha=0.15, label=f'Therapeutic Window ({C_LOW}-{C_HIGH} mmol/L)')

ax.set_title('Magnesium concentration-time profiles', fontsize=14)
ax.set_xlabel('time (hours)', fontsize=12)
ax.set_ylabel('plasma concentration (mmol/L)', fontsize=12)
ax.legend()
ax.grid(True, alpha=0.3)

st.pyplot(fig)

st.markdown("### Pharmacokinetic metrics")

def get_metrics(cp):
    idx = np.argmax(cp)
    auc = np.trapezoid(cp, t)
    return cp[idx], t[idx], auc

col1, col2, col3 = st.columns(3)

with col1:
    st.info("**IV Bolus**")
    cmax, tmax, auc = get_metrics(c_iv_bolus)
    st.write(f"**Cmax:** {cmax:.3f} mmol/L")
    st.write(f"**AUC:** {auc:.2f} mmol*h/L")

with col2:
    st.warning(f"**IV Infusion ({t_inf}h)**")
    cmax, tmax, auc = get_metrics(c_iv_inf)
    st.write(f"**Cmax:** {cmax:.3f} mmol/L")
    st.write(f"**Tmax:** {tmax:.2f} hrs")
    st.write(f"**AUC:** {auc:.2f} mmol*h/L")

with col3:
    st.success("**Oral (Mg Citrate)**")
    cmax, tmax, auc = get_metrics(c_oral_cit)
    st.write(f"**Cmax:** {cmax:.3f} mmol/L")
    st.write(f"**Tmax:** {tmax:.2f} hrs")
    st.write(f"**AUC:** {auc:.2f} mmol*h/L")