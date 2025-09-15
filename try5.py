import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

def boyd_kleinman(xi):
    return np.exp(-((xi - 2) / 0.8) ** 2)

def calc_waist_opt(m2, f, lam, din):
    return 2 * m2 * f * lam / (np.pi * din)

def calc_xi(w, f, lam, l):
    zr = np.pi * w**2 / lam
    b = 2 * zr / 1000
    return l / b

def calc_rgen(w_arr, w_opt, l, lam, f, pwr, bright, det_eff, coup_eff, normalized):
    xis = np.array([calc_xi(w, f, lam, l) for w in w_arr])
    vals = boyd_kleinman(xis)
    val_opt = boyd_kleinman(calc_xi(w_opt, f, lam, l))
    if normalized:
        pwr = 1.0
        bright = 1.0
    return bright * pwr * det_eff * coup_eff * vals / val_opt

st.title("SPDC Type - II Coincidence Rate Vs Beam Waist")

with st.sidebar:
    st.header("Input Parameters")
    m2 = st.number_input("Beam M²", min_value=1.0, max_value=5.0, value=1.1, step=0.01, format="%.2f")
    wavelength_nm = st.number_input("Wavelength (nm)", min_value=350.0, max_value=800.0, value=405.0, step=1.0, format="%.1f")
    focal_length = st.number_input("Focal Length (mm)", min_value=10.0, max_value=500.0, value=200.0, step=1.0, format="%.1f")
    din = st.number_input("Input Diameter (mm)", min_value=0.1, max_value=10.0, value=2.0, step=0.1, format="%.2f")
    l = st.number_input("Crystal Length (mm)", min_value=1.0, max_value=100.0, value=25.0, step=1.0, format="%.1f")
    brightness = st.number_input("Brightness (pairs/s/mW)", min_value=0.0, max_value=1e6, value=84200.0, step=1000.0, format="%.1f")
    pwr = st.number_input("Pump Power (mW)", min_value=0.01, max_value=1000.0, value=4.8, step=0.01, format="%.2f")
    det_eff = st.number_input("Detector efficiency", min_value=0.0, max_value=1.0, value=0.65, step=0.01, format="%.2f")
    coup_eff = st.number_input("Coupling efficiency", min_value=0.0, max_value=1.0, value=0.11, step=0.01, format="%.2f")
    normalized = st.checkbox("Normalized calculation", True)

lam = wavelength_nm / 1000.0
w_opt = calc_waist_opt(m2, focal_length, lam, din)
w_arr = np.linspace(max(0.1, w_opt / 2), w_opt * 2, 300)
rates = calc_rgen(w_arr, w_opt, l, lam, focal_length, pwr, brightness, det_eff, coup_eff, normalized)

idx_max = np.argmax(rates)
w_max = w_arr[idx_max]
r_max = rates[idx_max]

st.markdown(f"### Theortical optimum waist: {w_opt:.2f} µm")
st.markdown(f"### Peak waist from simulation: {w_max:.2f} µm")


fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(w_arr, rates, color='blue', linewidth=2)
ax.axvline(w_opt, color='red', linestyle='--', linewidth=2, label=f'Optimum waist {w_opt:.2f} µm')
ax.axvline(w_max, color='green', linestyle=':', linewidth=2, label=f'Peak waist {w_max:.2f} µm')
ax.set_xlabel('Waist (µm)', fontsize=14)
ax.set_ylabel('Norm. Coincidence rate', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True)
ax.text(0.5, -0.15, f'Peak at waist {w_max:.2f} µm',
        ha='center', va='top', fontsize=14, transform=ax.transAxes)
st.pyplot(fig)
