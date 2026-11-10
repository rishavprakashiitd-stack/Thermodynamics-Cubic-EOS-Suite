# ============================================================================
#  Pure-Gas Properties from Cubic Equations of State (vdW, RK, SRK, PR)
# ============================================================================

import math
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

# ----------------------------------------------------------------------------
#  PART 1 : THERMODYNAMIC ENGINE & 10 CORE GASES
# ----------------------------------------------------------------------------
R_J = 8.314       # J/(mol K)
R_CC = 83.14      # cm3 bar/(mol K)
T_REF = 298.15    # Reference state: 298.15 K, 1 bar
P_REF = 1.0

# 10 Selected Key Gases:
# Name, Tc (K), Pc (bar), omega, A, B, C, D (Cp/R = A + B*T + C*T^2 + D/T^2)
gases = [
    ("Methane",          190.6,  45.99,  0.012, 1.702,  9.081e-3,  -2.164e-6,  0.0),
    ("Ethane",           305.3,  48.72,  0.100, 1.131,  19.225e-3, -5.561e-6,  0.0),
    ("Propane",          369.8,  42.48,  0.152, 1.213,  28.785e-3, -8.824e-6,  0.0),
    ("n-Butane",         425.1,  37.96,  0.200, 1.935,  36.915e-3, -11.402e-6, 0.0),
    ("Ethylene",         282.3,  50.40,  0.087, 1.424,  14.394e-3, -4.392e-6,  0.0),
    ("Nitrogen",         126.2,  34.00,  0.038, 3.280,  0.593e-3,  0.0,        0.040e5),
    ("Oxygen",           154.6,  50.43,  0.022, 3.639,  0.506e-3,  0.0,       -0.227e5),
    ("Carbon dioxide",   304.2,  73.83,  0.224, 5.457,  1.045e-3,  0.0,       -1.157e5),
    ("Hydrogen",          33.19, 13.13, -0.216, 3.249,  0.422e-3,  0.0,        0.083e5),
    ("Water",            647.1, 220.55,  0.345, 3.470,  1.450e-3,  0.0,        0.121e5),
]

NAME, TC, PC, W, A, B, C, D = range(8)
eos_names = ["van der Waals", "Redlich-Kwong", "Soave-Redlich-Kwong", "Peng-Robinson"]


def sigma_of(e):
    if e == 0:
        return 0.0
    if e == 3:
        return 1.0 + math.sqrt(2.0)
    return 1.0


def eps_of(e):
    return 1.0 - math.sqrt(2.0) if e == 3 else 0.0


def omega_b(e):
    if e == 0:
        return 1.0 / 8.0
    return 0.07780 if e == 3 else 0.08664


def psi_a(e):
    if e == 0:
        return 27.0 / 64.0
    return 0.45724 if e == 3 else 0.42748


def kappa_of(e, w):
    if e == 2:
        return 0.480 + 1.574 * w - 0.176 * w * w
    return 0.37464 + 1.54226 * w - 0.26992 * w * w


def alpha_of(e, Tr, w):
    if e == 0:
        return 1.0
    if e == 1:
        return 1.0 / math.sqrt(Tr)
    s = 1.0 + kappa_of(e, w) * (1.0 - math.sqrt(Tr))
    return s * s


def dln_alpha(e, Tr, w):
    if e == 0:
        return 0.0
    if e == 1:
        return -0.5
    k = kappa_of(e, w)
    return -k * math.sqrt(Tr) / (1.0 + k * (1.0 - math.sqrt(Tr)))


def get_params(e, g, T, P):
    Tr = T / g[TC]
    Pr = P / g[PC]
    beta = omega_b(e) * Pr / Tr
    q = psi_a(e) * alpha_of(e, Tr, g[W]) / (omega_b(e) * Tr)
    return beta, q


def solve_cubic(c2, c1, c0):
    Q = (c2 * c2 - 3.0 * c1) / 9.0
    R = (2.0 * c2 ** 3 - 9.0 * c2 * c1 + 27.0 * c0) / 54.0
    if R * R < Q ** 3:
        theta = math.acos(max(-1.0, min(1.0, R / math.sqrt(Q ** 3))))
        m = -2.0 * math.sqrt(Q)
        return [m * math.cos((theta + 2.0 * math.pi * k) / 3.0) - c2 / 3.0 for k in range(3)]
    s = 1.0 if R >= 0 else -1.0
    val = abs(R) + math.sqrt(max(0.0, R * R - Q ** 3))
    A_ = -s * (val ** (1.0 / 3.0))
    B_ = 0.0 if abs(A_) < 1e-12 else Q / A_
    return [A_ + B_ - c2 / 3.0]


def get_roots(e, g, T, P):
    beta, q = get_params(e, g, T, P)
    s = sigma_of(e) + eps_of(e)
    p = sigma_of(e) * eps_of(e)
    c2 = s * beta - (1.0 + beta)
    c1 = p * beta ** 2 - (1.0 + beta) * s * beta + q * beta
    c0 = -(1.0 + beta) * p * beta ** 2 - q * beta ** 2
    good = [z for z in solve_cubic(c2, c1, c0) if z > beta + 1e-9]
    good.sort()
    return good


def evaluate(e, g, T, P, Z):
    beta, q = get_params(e, g, T, P)
    sg = sigma_of(e)
    ep = eps_of(e)
    if abs(sg - ep) < 1e-12:
        I = beta / (Z + ep * beta)
    else:
        I = math.log(max(1e-15, (Z + sg * beta) / (Z + ep * beta))) / (sg - ep)
    dl = dln_alpha(e, T / g[TC], g[W])
    st = {"Z": Z, "V": Z * R_CC * T / P}
    st["HR"] = R_J * T * (Z - 1.0 + (dl - 1.0) * q * I)
    arg_log = max(1e-15, Z - beta)
    st["SR"] = R_J * (math.log(arg_log) + dl * q * I)
    st["lnphi"] = Z - 1.0 - math.log(arg_log) - q * I
    st["phi"] = math.exp(min(50.0, max(-50.0, st["lnphi"])))
    st["f"] = st["phi"] * P
    return st


def saturation(e, g, T):
    if T >= g[TC]:
        return False, 0.0, None, None
    P = g[PC] * math.exp(5.373 * (1.0 + g[W]) * (1.0 - g[TC] / T))
    for _ in range(2000):
        Z = get_roots(e, g, T, P)
        if len(Z) == 0:
            return False, 0.0, None, None
        if len(Z) < 3:
            P *= 0.9 if Z[0] < 0.3 else 1.1
            continue
        vap = evaluate(e, g, T, P, Z[2])
        liq = evaluate(e, g, T, P, Z[0])
        ratio = math.exp(liq["lnphi"] - vap["lnphi"])
        if abs(ratio - 1.0) < 1e-9:
            return True, P, vap, liq
        P *= ratio
    return False, 0.0, None, None


def stable_state(e, g, T, P):
    Z = get_roots(e, g, T, P)
    if len(Z) == 0:
        return evaluate(e, g, T, P, 1.0), "?"
    if len(Z) == 3:
        v = evaluate(e, g, T, P, Z[2])
        l = evaluate(e, g, T, P, Z[0])
        if l["lnphi"] < v["lnphi"]:
            return l, "Liquid"
        return v, "Vapor"
    st = evaluate(e, g, T, P, Z[0])
    if T >= g[TC]:
        return st, "Supercritical Gas"
    found, Ps, _, _ = saturation(e, g, T)
    if found:
        return st, ("Liquid" if P > Ps else "Vapor")
    return st, ("Liquid" if st["Z"] < 0.3 else "Vapor")


def ideal_h(g, T):
    T0 = T_REF
    return R_J * (g[A] * (T - T0) + g[B] / 2.0 * (T ** 2 - T0 ** 2)
                 + g[C] / 3.0 * (T ** 3 - T0 ** 3) - g[D] * (1.0 / T - 1.0 / T0))


def ideal_s(g, T, P):
    T0 = T_REF
    return R_J * (g[A] * math.log(T / T0) + g[B] * (T - T0) + g[C] / 2.0 * (T ** 2 - T0 ** 2)
                 - g[D] / 2.0 * (1.0 / T ** 2 - 1.0 / T0 ** 2)) - R_J * math.log(P / P_REF)


# ----------------------------------------------------------------------------
#  PART 2 : TEXT FORMATTERS
# ----------------------------------------------------------------------------
def state_block(title, s):
    return (f"--- {title} ---\n"
            f"  Z (Compressibility)       = {s['Z']:.5f}\n"
            f"  Molar Volume V            = {s['V']:.3f} cm3/mol\n"
            f"  Residual Enthalpy H^R     = {s['HR']:.2f} J/mol\n"
            f"  Residual Entropy S^R      = {s['SR']:.4f} J/(mol K)\n"
            f"  Fugacity Coefficient phi  = {s['phi']:.5f}\n"
            f"  Fugacity f = phi*P        = {s['f']:.4f} bar\n\n")


def state_text(g, e, T, P):
    Z = get_roots(e, g, T, P)
    st, phase = stable_state(e, g, T, P)
    t = f"Compound: {g[NAME]}  |  EOS: {eos_names[e]}\n"
    t += f"T = {T:.2f} K  |  P = {P:.3f} bar  |  Tr = {T / g[TC]:.3f}  |  Pr = {P / g[PC]:.3f}\n"
    t += f"Roots calculated: {len(Z)}\n\n"
    if len(Z) == 3:
        t += state_block("Vapor Root (Z_max)", evaluate(e, g, T, P, Z[2]))
        t += state_block("Liquid Root (Z_min)", evaluate(e, g, T, P, Z[0]))
        t += f">> Stable Phase: {phase} (Lowest Fugacity Criterion)\n\n"
    else:
        t += f">> Phase: {phase}\n" + state_block("Properties", st)
    H = ideal_h(g, T) + st["HR"]
    S = ideal_s(g, T, P) + st["SR"]
    t += "--- Total Properties (Ref: 298.15 K, 1 bar Ideal Gas) ---\n"
    t += f"  Ideal Enthalpy H^ig       = {ideal_h(g, T):.2f} J/mol\n"
    t += f"  Real Enthalpy  H          = {H:.2f} J/mol\n"
    t += f"  Ideal Entropy  S^ig       = {ideal_s(g, T, P):.4f} J/(mol K)\n"
    t += f"  Real Entropy   S          = {S:.4f} J/(mol K)\n"
    return t


def change_text(g, e, T1, P1, T2, P2):
    s1, ph1 = stable_state(e, g, T1, P1)
    s2, ph2 = stable_state(e, g, T2, P2)
    dHig = ideal_h(g, T2) - ideal_h(g, T1)
    dSig = ideal_s(g, T2, P2) - ideal_s(g, T1, P1)
    dH = dHig + s2["HR"] - s1["HR"]
    dS = dSig + s2["SR"] - s1["SR"]
    t = f"State Transformation for {g[NAME]} ({eos_names[e]})\n"
    t += "=" * 54 + "\n"
    t += f"Parameter            State 1            State 2\n"
    t += "-" * 54 + "\n"
    t += f"T (K)             {T1:14.2f}     {T2:14.2f}\n"
    t += f"P (bar)           {P1:14.3f}     {P2:14.3f}\n"
    t += f"Phase             {ph1:>14}     {ph2:>14}\n"
    t += f"Z                 {s1['Z']:14.5f}     {s2['Z']:14.5f}\n"
    t += f"H^R (J/mol)       {s1['HR']:14.2f}     {s2['HR']:14.2f}\n"
    t += f"S^R (J/mol K)     {s1['SR']:14.4f}     {s2['SR']:14.4f}\n"
    t += f"f (bar)           {s1['f']:14.4f}     {s2['f']:14.4f}\n"
    t += "-" * 54 + "\n"
    t += f"dH (Real) = {dH:.2f} J/mol   (dH^ig = {dHig:.2f} J/mol)\n"
    t += f"dS (Real) = {dS:.4f} J/(mol K) (dS^ig = {dSig:.4f} J/(mol K))\n"
    return t, s1, s2


def table_rows(g, e, T, P1, P2, N):
    rows = []
    for i in range(N):
        P = P1 + (P2 - P1) * i / (N - 1)
        s, ph = stable_state(e, g, T, P)
        rows.append((P, ph, s["Z"], s["V"], s["HR"], s["SR"], s["phi"], s["f"]))
    return rows


def table_text(g, e, T, rows):
    t = f"Isotherm Scan: {g[NAME]} | {eos_names[e]} | T = {T:.2f} K\n"
    t += "-" * 84 + "\n"
    t += f"{'P(bar)':>8} {'Phase':>7} {'Z':>8} {'V(cm3/mol)':>12} {'HR(J/mol)':>11} {'SR(J/mol K)':>12} {'phi':>8} {'f(bar)':>10}\n"
    t += "-" * 84 + "\n"
    for r in rows:
        t += f"{r[0]:8.2f} {r[1]:>7} {r[2]:8.4f} {r[3]:12.2f} {r[4]:11.1f} {r[5]:12.3f} {r[6]:8.4f} {r[7]:10.2f}\n"
    return t


def compare_text(g, T, P):
    t = f"EOS Benchmark Comparison: {g[NAME]} at T = {T:.2f} K, P = {P:.2f} bar\n"
    t += "-" * 76 + "\n"
    t += f"{'EOS Model':<22} {'Phase':>8} {'Z':>9} {'V(cm3/mol)':>13} {'HR(J/mol)':>11} {'phi':>8}\n"
    t += "-" * 76 + "\n"
    res = []
    for k in range(4):
        s, ph = stable_state(k, g, T, P)
        t += f"{eos_names[k]:<22} {ph:>8} {s['Z']:9.4f} {s['V']:13.2f} {s['HR']:11.1f} {s['phi']:8.4f}\n"
        res.append((eos_names[k], s['Z'], s['V']))
    return t, res


def saturation_text(g, e, T):
    found, Ps, v, l = saturation(e, g, T)
    if not found:
        return f"No VLE equilibrium detected for {g[NAME]} at T = {T:.2f} K.\n(T must be below Tc = {g[TC]:.2f} K).", None
    t = f"VLE Saturation: {g[NAME]} ({eos_names[e]})\n"
    t += f"Temperature T = {T:.2f} K (Tc = {g[TC]:.2f} K)\n"
    t += f"Psat = {Ps:.4f} bar\n"
    t += f"Heat of Vaporization dH_vap = {v['HR'] - l['HR']:.2f} J/mol\n\n"
    t += state_block("Saturated Vapor Phase", v)
    t += state_block("Saturated Liquid Phase", l)
    return t, Ps


# ----------------------------------------------------------------------------
#  PART 3 : DESKTOP GUI & EMBEDDED GRAPH RENDERING
# ----------------------------------------------------------------------------
class ThermoDesktopApp:
    def __init__(self, root):
        self.root = root
        self.root.title("ThermoSim - Cubic Equations of State Suite")
        self.root.geometry("1180x740")
        self.root.minsize(980, 620)
        self.root.configure(bg="#0f172a")

        self.setup_styles()
        self.build_ui()
        self.root.after(150, self.draw_placeholder)

    def setup_styles(self):
        style = ttk.Style()
        style.theme_use("clam")

        style.configure(".", background="#0f172a", foreground="#f8fafc", font=("Segoe UI", 9))
        style.configure("TNotebook", background="#0f172a", borderwidth=0)
        style.configure("TNotebook.Tab", background="#1e293b", foreground="#94a3b8",
                        padding=[12, 6], font=("Segoe UI", 9, "bold"))
        style.map("TNotebook.Tab",
                  background=[("selected", "#0ea5e9")],
                  foreground=[("selected", "#ffffff")])

        style.configure("TFrame", background="#0f172a")
        style.configure("Card.TFrame", background="#1e293b", relief="flat")
        style.configure("TLabel", background="#1e293b", foreground="#cbd5e1", font=("Segoe UI", 9))
        style.configure("Header.TLabel", background="#0f172a", foreground="#f8fafc", font=("Segoe UI", 16, "bold"))
        style.configure("SubHeader.TLabel", background="#0f172a", foreground="#94a3b8", font=("Segoe UI", 9))

        style.configure("Action.TButton", font=("Segoe UI", 9, "bold"),
                        background="#0ea5e9", foreground="#ffffff", borderwidth=0, padding=6)
        style.map("Action.TButton", background=[("active", "#0284c7")])

        style.configure("Secondary.TButton", font=("Segoe UI", 9),
                        background="#334155", foreground="#ffffff", borderwidth=0, padding=6)
        style.map("Secondary.TButton", background=[("active", "#475569")])

    def build_ui(self):
        # 1. Header
        header = tk.Frame(self.root, bg="#0f172a")
        header.pack(fill="x", padx=16, pady=(10, 6))

        ttk.Label(header, text="Thermodynamics Simulator", style="Header.TLabel").pack(anchor="w")
        ttk.Label(header, text="Pure Component Thermodynamics, Phase Equilibria & Residual Properties", style="SubHeader.TLabel").pack(anchor="w")

        # 2. Control Bar
        ctrl_card = ttk.Frame(self.root, style="Card.TFrame", padding=(12, 8))
        ctrl_card.pack(fill="x", padx=14, pady=4)

        ttk.Label(ctrl_card, text="Select Gas:").pack(side="left", padx=(4, 6))
        self.gas_box = ttk.Combobox(ctrl_card, values=[g[NAME] for g in gases], state="readonly", width=16)
        self.gas_box.current(2)  # Propane
        self.gas_box.pack(side="left", padx=(0, 20))

        ttk.Label(ctrl_card, text="Equation of State:").pack(side="left", padx=(4, 6))
        self.eos_box = ttk.Combobox(ctrl_card, values=eos_names, state="readonly", width=22)
        self.eos_box.current(3)  # Peng-Robinson
        self.eos_box.pack(side="left", padx=(0, 10))

        # 3. Main Workspace
        main_paned = ttk.PanedWindow(self.root, orient="horizontal")
        main_paned.pack(fill="both", expand=True, padx=14, pady=8)

        # Left Column: Tabs & Text Output
        left_box = ttk.Frame(main_paned, style="TFrame")
        main_paned.add(left_box, weight=1)

        self.book = ttk.Notebook(left_box)
        self.book.pack(fill="x", pady=(0, 8))
        self.build_tabs()

        out_card = ttk.Frame(left_box, style="Card.TFrame", padding=6)
        out_card.pack(fill="both", expand=True)

        ttk.Label(out_card, text="Calculated Output:", font=("Segoe UI", 9, "bold")).pack(anchor="w", pady=(0, 4))
        self.out = tk.Text(out_card, font=("Consolas", 9), wrap="none", bg="#0a0f1d", fg="#38bdf8",
                           insertbackground="white", bd=0, relief="flat", height=14)
        sy = ttk.Scrollbar(out_card, orient="vertical", command=self.out.yview)
        sx = ttk.Scrollbar(out_card, orient="horizontal", command=self.out.xview)
        self.out.configure(yscrollcommand=sy.set, xscrollcommand=sx.set)
        sy.pack(side="right", fill="y")
        sx.pack(side="bottom", fill="x")
        self.out.pack(fill="both", expand=True)

        # Right Column: Canvas Plots
        right_box = ttk.Frame(main_paned, style="Card.TFrame", padding=8)
        main_paned.add(right_box, weight=1)

        ttk.Label(right_box, text="Visual Thermodynamic Graph:", font=("Segoe UI", 9, "bold")).pack(anchor="w", pady=(0, 4))
        self.canvas = tk.Canvas(right_box, bg="#0a0f1d", highlightthickness=1, highlightbackground="#334155")
        self.canvas.pack(fill="both", expand=True)

        self.write_output("Choose gas & EOS, fill input parameters above and click calculate.\n"
                          "Thermodynamic results and visual curves will update automatically.")

    def build_tabs(self):
        f1, self.v1 = self.create_tab_inputs("State (T, P)", [("T (K)", "300"), ("P (bar)", "5")])
        ttk.Button(f1, text="Calculate State", style="Action.TButton",
                   command=lambda: self.safe(self.run_state)).grid(row=1, column=2, padx=10, pady=4)

        f2, self.v2 = self.create_tab_inputs("dH & dS", [("T1 (K)", "300"), ("P1 (bar)", "1"),
                                                         ("T2 (K)", "400"), ("P2 (bar)", "20")])
        ttk.Button(f2, text="Compute Delta", style="Action.TButton",
                   command=lambda: self.safe(self.run_change)).grid(row=1, column=4, padx=10, pady=4)

        f3, self.v3 = self.create_tab_inputs("Isotherm (P Scan)", [("T (K)", "350"), ("P start (bar)", "1"),
                                                                   ("P end (bar)", "50"), ("Points", "20")])
        ttk.Button(f3, text="Scan & Plot", style="Action.TButton",
                   command=lambda: self.safe(self.run_table)).grid(row=1, column=4, padx=6, pady=4)
        ttk.Button(f3, text="Export CSV", style="Secondary.TButton",
                   command=lambda: self.safe(self.export_csv)).grid(row=1, column=5, padx=6, pady=4)

        f4, self.v4 = self.create_tab_inputs("Compare EOS", [("T (K)", "350"), ("P (bar)", "30")])
        ttk.Button(f4, text="Compare All", style="Action.TButton",
                   command=lambda: self.safe(self.run_compare)).grid(row=1, column=2, padx=10, pady=4)

        f5, self.v5 = self.create_tab_inputs("VLE Saturation", [("T (K)", "300")])
        ttk.Button(f5, text="Find Psat & Plot Curve", style="Action.TButton",
                   command=lambda: self.safe(self.run_sat)).grid(row=1, column=1, padx=10, pady=4)

    def create_tab_inputs(self, title, fields):
        frame = ttk.Frame(self.book, style="Card.TFrame", padding=10)
        self.book.add(frame, text=f" {title} ")
        vars_list = []
        for col, (label_txt, default_val) in enumerate(fields):
            ttk.Label(frame, text=label_txt).grid(row=0, column=col, sticky="w", padx=4, pady=2)
            var = tk.StringVar(value=default_val)
            ent = tk.Entry(frame, textvariable=var, width=11, bg="#0a0f1d", fg="#f8fafc",
                           insertbackground="white", bd=1, relief="solid")
            ent.grid(row=1, column=col, padx=4, pady=2)
            vars_list.append(var)
        return frame, vars_list

    def write_output(self, text):
        self.out.config(state="normal")
        self.out.delete("1.0", "end")
        self.out.insert("end", text)
        self.out.config(state="disabled")

    def parse_float(self, var, lo, hi, name):
        try:
            val = float(var.get())
        except ValueError:
            raise ValueError(f"'{name}' must be numeric.")
        if val < lo or val > hi:
            raise ValueError(f"'{name}' must be in range [{lo}, {hi}].")
        return val

    def safe(self, action):
        try:
            action()
        except ValueError as err:
            messagebox.showerror("Validation Error", str(err))

    def cur_gas(self):
        return gases[self.gas_box.current()]

    def cur_eos(self):
        return self.eos_box.current()

    # ---- GRAPH RENDERING (CANVAS) ----
    def draw_placeholder(self):
        self.canvas.delete("all")
        w = self.canvas.winfo_width() or 480
        h = self.canvas.winfo_height() or 400
        self.canvas.create_text(w / 2, h / 2, text="Run any calculation to display\nreal-time graph here.",
                                fill="#64748b", font=("Segoe UI", 11), justify="center")

    def draw_line_plot(self, x_pts, y_pts, title, xlabel, ylabel, highlight_pt=None):
        self.canvas.delete("all")
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        if w < 100 or h < 100:
            return

        pad_left, pad_right, pad_top, pad_bot = 65, 30, 45, 50
        gw = w - pad_left - pad_right
        gh = h - pad_top - pad_bot

        self.canvas.create_text(w / 2, 20, text=title, fill="#f8fafc", font=("Segoe UI", 11, "bold"))

        xmin, xmax = min(x_pts), max(x_pts)
        ymin, ymax = min(y_pts), max(y_pts)
        if xmin == xmax:
            xmax += 1.0
        if ymin == ymax:
            ymax += 1.0

        # Axis Lines
        self.canvas.create_line(pad_left, pad_top + gh, pad_left + gw, pad_top + gh, fill="#475569", width=2)
        self.canvas.create_line(pad_left, pad_top, pad_left, pad_top + gh, fill="#475569", width=2)

        # Grids & Ticks
        for i in range(5):
            gx = pad_left + gw * (i / 4.0)
            gy = pad_top + gh * (1.0 - i / 4.0)
            self.canvas.create_line(gx, pad_top, gx, pad_top + gh, fill="#1e293b", dash=(2, 4))
            self.canvas.create_line(pad_left, gy, pad_left + gw, gy, fill="#1e293b", dash=(2, 4))
            xv = xmin + (xmax - xmin) * (i / 4.0)
            yv = ymin + (ymax - ymin) * (i / 4.0)
            self.canvas.create_text(gx, pad_top + gh + 15, text=f"{xv:.2f}", fill="#94a3b8", font=("Segoe UI", 8))
            self.canvas.create_text(pad_left - 8, gy, text=f"{yv:.2f}", fill="#94a3b8", font=("Segoe UI", 8), anchor="e")

        self.canvas.create_text(w / 2, h - 12, text=xlabel, fill="#94a3b8", font=("Segoe UI", 9))
        self.canvas.create_text(18, h / 2, text=ylabel, fill="#94a3b8", font=("Segoe UI", 9), angle=90)

        # Draw Points & Curve
        pixel_points = []
        for x, y in zip(x_pts, y_pts):
            px = pad_left + (x - xmin) / (xmax - xmin) * gw
            py = pad_top + gh - (y - ymin) / (ymax - ymin) * gh
            pixel_points.append((px, py))

        for i in range(len(pixel_points) - 1):
            p1, p2 = pixel_points[i], pixel_points[i + 1]
            self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill="#38bdf8", width=3, smooth=True)

        for px, py in pixel_points:
            self.canvas.create_oval(px - 3, py - 3, px + 3, py + 3, fill="#0ea5e9", outline="#f8fafc")

        if highlight_pt:
            hx, hy = highlight_pt
            px = pad_left + (hx - xmin) / (xmax - xmin) * gw
            py = pad_top + gh - (hy - ymin) / (ymax - ymin) * gh
            self.canvas.create_oval(px - 6, py - 6, px + 6, py + 6, fill="#ef4444", outline="#ffffff", width=2)
            self.canvas.create_text(px + 10, py - 10, text=f"State: ({hx:.2f}, {hy:.3f})", fill="#ef4444",
                                    font=("Segoe UI", 8, "bold"), anchor="w")

    def draw_bar_plot(self, categories, values, title, ylabel):
        self.canvas.delete("all")
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()
        if w < 100 or h < 100:
            return

        pad_left, pad_right, pad_top, pad_bot = 65, 30, 45, 50
        gw = w - pad_left - pad_right
        gh = h - pad_top - pad_bot

        self.canvas.create_text(w / 2, 20, text=title, fill="#f8fafc", font=("Segoe UI", 11, "bold"))
        self.canvas.create_line(pad_left, pad_top + gh, pad_left + gw, pad_top + gh, fill="#475569", width=2)

        min_v = min(0.0, min(values))
        max_v = max(values) * 1.15
        if max_v == min_v:
            max_v += 1.0

        n = len(categories)
        bar_w = (gw / n) * 0.52
        spacing = gw / n
        colors = ["#38bdf8", "#818cf8", "#c084fc", "#f43f5e"]

        for i, (cat, val) in enumerate(zip(categories, values)):
            bx = pad_left + i * spacing + spacing / 2
            bar_h = (val - min_v) / (max_v - min_v) * gh
            by = pad_top + gh - bar_h

            c = colors[i % len(colors)]
            self.canvas.create_rectangle(bx - bar_w / 2, by, bx + bar_w / 2, pad_top + gh, fill=c, outline="")
            self.canvas.create_text(bx, by - 10, text=f"{val:.3f}", fill="#f8fafc", font=("Segoe UI", 8, "bold"))
            self.canvas.create_text(bx, pad_top + gh + 18, text=cat, fill="#94a3b8", font=("Segoe UI", 8))

        self.canvas.create_text(18, h / 2, text=ylabel, fill="#94a3b8", font=("Segoe UI", 9), angle=90)

    # ---- HANDLERS ----
    def run_state(self):
        T = self.parse_float(self.v1[0], 1, 3000, "T")
        P = self.parse_float(self.v1[1], 1e-5, 5000, "P")
        g = self.cur_gas()
        e = self.cur_eos()

        self.write_output(state_text(g, e, T, P))

        p_vals = [max(0.1, P * f) for f in [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5]]
        z_vals = [stable_state(e, g, T, p)[0]["Z"] for p in p_vals]
        cur_z = stable_state(e, g, T, P)[0]["Z"]

        self.draw_line_plot(p_vals, z_vals, f"Z vs Pressure (T = {T:.1f} K)", "Pressure P (bar)", "Compressibility Z", (P, cur_z))

    def run_change(self):
        T1 = self.parse_float(self.v2[0], 1, 3000, "T1")
        P1 = self.parse_float(self.v2[1], 1e-5, 5000, "P1")
        T2 = self.parse_float(self.v2[2], 1, 3000, "T2")
        P2 = self.parse_float(self.v2[3], 1e-5, 5000, "P2")
        g = self.cur_gas()
        e = self.cur_eos()

        txt, s1, s2 = change_text(g, e, T1, P1, T2, P2)
        self.write_output(txt)
        self.draw_bar_plot(["State 1", "State 2"], [s1["HR"], s2["HR"]], "Residual Enthalpy H^R (J/mol)", "H^R (J/mol)")

    def run_table(self):
        T = self.parse_float(self.v3[0], 1, 3000, "T")
        P1 = self.parse_float(self.v3[1], 1e-5, 5000, "P start")
        P2 = self.parse_float(self.v3[2], P1, 5000, "P end")
        N = int(self.parse_float(self.v3[3], 3, 100, "Points"))
        g = self.cur_gas()
        e = self.cur_eos()

        rows = table_rows(g, e, T, P1, P2, N)
        self.write_output(table_text(g, e, T, rows))

        pressures = [r[0] for r in rows]
        z_vals = [r[2] for r in rows]
        self.draw_line_plot(pressures, z_vals, f"{g[NAME]} Isotherm: Z vs P (T={T:.1f} K)", "Pressure P (bar)", "Compressibility Z")

    def export_csv(self):
        T = self.parse_float(self.v3[0], 1, 3000, "T")
        P1 = self.parse_float(self.v3[1], 1e-5, 5000, "P start")
        P2 = self.parse_float(self.v3[2], P1, 5000, "P end")
        N = int(self.parse_float(self.v3[3], 3, 100, "Points"))
        g = self.cur_gas()
        e = self.cur_eos()

        rows = table_rows(g, e, T, P1, P2, N)
        path = filedialog.asksaveasfilename(defaultextension=".csv", initialfile="thermo_data.csv",
                                            filetypes=[("CSV file", "*.csv")])
        if not path:
            return
        with open(path, "w") as f:
            f.write("P_bar,Phase,Z,V_cm3_mol,HR_J_mol,SR_J_molK,phi,f_bar\n")
            for r in rows:
                f.write(f"{r[0]:.6g},{r[1]},{r[2]:.6g},{r[3]:.6g},{r[4]:.6g},{r[5]:.6g},{r[6]:.6g},{r[7]:.6g}\n")
        messagebox.showinfo("Export Successful", f"Saved {len(rows)} data rows to {path}")

    def run_compare(self):
        T = self.parse_float(self.v4[0], 1, 3000, "T")
        P = self.parse_float(self.v4[1], 1e-5, 5000, "P")
        g = self.cur_gas()

        txt, res = compare_text(g, T, P)
        self.write_output(txt)

        names = ["vdW", "RK", "SRK", "PR"]
        z_vals = [r[1] for r in res]
        self.draw_bar_plot(names, z_vals, f"Compressibility Z Comparison for {g[NAME]}", "Z")

    def run_sat(self):
        T = self.parse_float(self.v5[0], 1, 3000, "T")
        g = self.cur_gas()
        e = self.cur_eos()

        txt, psat = saturation_text(g, e, T)
        self.write_output(txt)
        if psat is None:
            return

        t_min = max(g[TC] * 0.45, 50.0)
        t_max = g[TC] * 0.98
        t_range = [t_min + (t_max - t_min) * i / 14.0 for i in range(15)]
        p_sats, valid_t = [], []
        for temp in t_range:
            ok, p_val, _, _ = saturation(e, g, temp)
            if ok:
                p_sats.append(p_val)
                valid_t.append(temp)

        if valid_t:
            self.draw_line_plot(valid_t, p_sats, f"Vapor Pressure Curve: {g[NAME]}", "Temperature (K)", "Psat (bar)", (T, psat))


# ----------------------------------------------------------------------------
#  LAUNCHER
# ----------------------------------------------------------------------------
def main():
    root = tk.Tk()
    app = ThermoDesktopApp(root)

    # Force window to top and take focus on startup
    root.lift()
    root.focus_force()
    root.attributes("-topmost", True)
    root.after_idle(root.attributes, "-topmost", False)

    root.mainloop()


if __name__ == "__main__":
    main()