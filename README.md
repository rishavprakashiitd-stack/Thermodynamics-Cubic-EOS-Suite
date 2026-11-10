# Thermodynamics Project: Pure-Gas Properties from Cubic Equations of State

A standalone desktop engineering application built with Python and Tkinter to evaluate thermodynamic properties, residual departure functions, and phase equilibria across four cubic equations of state:
- **vdW** (van der Waals)
- **RK** (Redlich-Kwong)
- **SRK** (Soave-Redlich-Kwong)
- **PR** (Peng-Robinson)

---

## 👥 Contributors (Group 2)
- **Rishav Prakash** (Entry No.: 2025CH11031)
- **Prince Raj** (Entry No.: 2025CH11863)

---

## 📌 Features
- **10 Industrial Gases:** Methane, Ethane, Propane, n-Butane, Ethylene, Nitrogen, Oxygen, Carbon Dioxide, Hydrogen, and Water.
- **Thermodynamic Properties:** Compressibility factor ($Z$), molar volume ($V$), residual enthalpy ($H^R$), residual entropy ($S^R$), fugacity coefficient ($\phi$), and fugacity ($f$).
- **VLE Saturation Solver:** Nested Wilson-Picard equifugacity loop for vapor pressure ($P^{sat}$) and latent heat of vaporization ($\Delta H_{vap}$).
- **State Change Computations:** $\Delta H$ and $\Delta S$ evaluation anchored to ideal-gas reference state ($298.15\text{ K}$, $1\text{ bar}$).
- **Real-Time Data Visualization:** Zero-dependency Tkinter vector canvas plotting isotherms ($Z$ vs $P$), $P\text{--}T$ saturation boundaries, and EOS model benchmark bar charts.
- **CSV Data Export:** Generate and export pressure scan datasets directly.

---

## 🚀 How to Run

No external scientific libraries (`matplotlib`, `numpy`, or `scipy`) are required.

```bash
thermodyanmics simulator.py
