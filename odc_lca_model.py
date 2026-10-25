#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
odc_lca_model.py
Modelo parametrico do no orbital de 1 MW (SSO 700 km) e calculos da linha de
base terrestre do artigo "Data Centers Orbitais de IA: Estrutura de LCA...".

Reproduz: Eq. 2-10, 13-15 e Tabelas 2-5 (massa, confiabilidade, agua, teto de
emissoes incorporadas), OAT e Monte Carlo (N=3000, semente 1).
Todos os parametros de entrada sao HIPOTESES (Tabela A1 do artigo).
Uso:  python odc_lca_model.py
"""
import math
import random
import json

# ----------------------------- constantes ---------------------------------
MU = 398600.4418          # km^3/s^2
RE = 6378.137             # km
J2 = 1.08263e-3
G0_SUN = 1361.0           # W/m^2
ECC_EARTH = 0.0167
OBLIQ = math.radians(23.44)
SIGMA = 5.670374e-8
G_STD = 9.80665
YEAR_S = 365.2422 * 86400.0

# ----------------------------- premissas (Tabela A1) ----------------------
DEFAULTS = dict(
    Pload=1.05e6,    # W      P1  (1 MW de TI + 5%)
    eta_c=0.97,      #        P2
    eta_rt=0.92,     #        P2
    DoD=0.35,        #        P3
    wb=165.0,        # Wh/kg  P3
    eta_cell=0.295,  #        P4 (BOL)
    k_arr=0.80,      #        P4
    rho_pv=3.2,      # kg/m2  P5
    q=360.0,         # W/m2   P6
    rho_rad=6.0,     # kg/m2  P6
    kIT=6.0,         # kg/kW  P7
    kst=0.15,        #        P8
    mga=0.20,        #        P9
    Wwing=30.0,      # m      P11
    hub=5.0,         # m      P11
    hspec=17.5,      # N.m.s/kg  P10
    redund=2.0,      #        P10
    Mhw=300.0,       # kg     P12
    rho_air=1e-13,   # kg/m3  P13
    Cd=2.2,          #        P13
    Isp=1500.0,      # s      P13
    prop_res=0.10,   #        P13
)

# faixas para OAT e Monte Carlo (P18)
RANGES = {
    "eta_cell": (0.27, 0.32), "k_arr": (0.70, 0.90), "rho_pv": (2.5, 4.5),
    "DoD": (0.25, 0.50), "wb": (120.0, 200.0), "q": (300.0, 420.0),
    "rho_rad": (4.0, 9.0), "kIT": (6.0, 12.0), "kst": (0.10, 0.25),
    "mga": (0.20, 0.30), "redund": (1.0, 3.0), "rho_air": (3e-14, 5e-13),
    "Pload": (1.05e6, 1.10e6),
}
MC_MODES = {"kIT": 8.0, "redund": 2.0}   # demais: modo = caso-base


# ----------------------------- geometria orbital --------------------------
def orbit_geometry(h_km=700.0):
    """Eq. 2: periodo e inclinacao heliossincrona."""
    r = RE + h_km
    T = 2.0 * math.pi * math.sqrt(r ** 3 / MU)
    n = 2.0 * math.pi / T
    omega_sso = 2.0 * math.pi / YEAR_S
    cosi = -omega_sso / (1.5 * n * J2 * (RE / r) ** 2)
    return dict(h=h_km, r=r, T_min=T / 60.0, n=n,
                inc_deg=math.degrees(math.acos(cosi)),
                beta_star_deg=math.degrees(math.asin(RE / r)))


def eclipse_min(beta_deg, orb):
    """Eq. 3: duracao do eclipse (sombra cilindrica), em minutos."""
    if abs(beta_deg) >= orb["beta_star_deg"]:
        return 0.0
    h, r = orb["h"], orb["r"]
    x = math.sqrt(h * h + 2.0 * RE * h) / (r * math.cos(math.radians(beta_deg)))
    return orb["T_min"] / math.pi * math.acos(min(1.0, x))


def beta_series_dawn_dusk(orb):
    """beta(t) ao longo do ano para LTAN 06:00 (Sol em movimento medio)."""
    sini = math.sin(math.radians(orb["inc_deg"]))
    cosi = math.cos(math.radians(orb["inc_deg"]))
    out = []
    for d in range(365):
        lam = 2.0 * math.pi * (d - 79.5) / 365.25
        sd = math.sin(OBLIQ) * math.sin(lam)
        cd = math.sqrt(1.0 - sd * sd)
        sb = cd * sini + sd * cosi
        out.append(math.degrees(math.asin(max(-1.0, min(1.0, sb)))))
    return out


def eclipse_stats(orb):
    betas = beta_series_dawn_dusk(orb)
    ecl = [eclipse_min(b, orb) for b in betas]
    ndays = sum(1 for e in ecl if e > 0.0)
    orb_per_day = 86400.0 / (orb["T_min"] * 60.0)
    absb = [abs(b) for b in betas]
    return dict(tecl_max_A=max(ecl), tecl_B=eclipse_min(0.0, orb),
                ndays_eclipse_A=ndays, beta_min=min(absb), beta_max=max(absb),
                orbits_10y=10.0 * 365.25 * orb_per_day,
                orbits_eclipse_A_10y=10.0 * ndays * orb_per_day)


# ----------------------------- dimensionamento ----------------------------
def size_node(p, t_ecl_min, orb):
    """Eq. 4-9: potencia, bateria, arranjo, radiador, ADCS, massa, arrasto."""
    G = G0_SUN / (1.0 + ECC_EARTH) ** 2          # irradiancia de afelio
    T_min, r, n = orb["T_min"], orb["r"], orb["n"]
    t_sun = T_min - t_ecl_min
    E_ecl = p["Pload"] * t_ecl_min / 60.0         # Wh
    E_bat = E_ecl / p["DoD"]
    M_bat = E_bat / p["wb"]
    P_gen = p["Pload"] / p["eta_c"] + \
        (E_ecl / (p["eta_c"] * p["eta_rt"])) / (t_sun / 60.0)
    A_pv = P_gen / (p["eta_cell"] * G * p["k_arr"])
    M_pv = A_pv * p["rho_pv"]
    Q_avg = P_gen * t_sun / T_min
    A_rad = Q_avg / (2.0 * p["q"])
    M_rad = A_rad * p["rho_rad"]
    M_it = p["kIT"] * 1000.0
    S = M_pv + M_rad + M_bat + M_it
    M_st = p["kst"] * S

    # Eq. 8: gradiente de gravidade e momento ciclico (duas asas simetricas)
    m_wing = M_pv / 2.0
    L = (A_pv / 2.0) / p["Wwing"]
    yc = p["hub"] + L / 2.0
    I_xx = 2.0 * m_wing * (L * L / 12.0 + yc * yc)
    T_gg = 1.5 * MU * 1e9 / ((r * 1e3) ** 3) * I_xx
    h_cyc = T_gg / (2.0 * n)
    M_adcs = p["redund"] * h_cyc / p["hspec"]

    # Eq. 9: massa seca e de lancamento (iteracao de ponto fixo p/ arrasto)
    M_dry = (S + M_st + M_adcs + p["Mhw"]) * (1.0 + p["mga"])
    v = math.sqrt(MU / r) * 1000.0
    F_drag = 0.5 * p["rho_air"] * v * v * p["Cd"] * A_pv
    life = 10.0 * 365.25 * 86400.0
    M_w, dv, m_p = M_dry, 0.0, 0.0
    for _ in range(60):
        dv = F_drag / M_w * life
        m_p = M_w * (1.0 - math.exp(-dv / (G_STD * p["Isp"]))) * (1.0 + p["prop_res"])
        M_w = M_dry + m_p
    return dict(t_ecl=t_ecl_min, t_sun=t_sun, E_ecl_kWh=E_ecl / 1e3,
                E_bat_kWh=E_bat / 1e3, M_bat=M_bat, P_gen_MW=P_gen / 1e6,
                A_pv=A_pv, M_pv=M_pv, A_rad=A_rad, M_rad=M_rad, M_it=M_it,
                M_st=M_st, I_xx=I_xx, T_gg=T_gg, h_cyc=h_cyc, M_adcs=M_adcs,
                M_hw=p["Mhw"], M_dry=M_dry, F_drag=F_drag, dv_drag=dv,
                m_prop=m_p, M_launch=M_w)


def radiator_check(orb, p):
    """Secao 3.3: fluxo liquido por face com IR planetario e albedo."""
    rho = math.asin(RE / orb["r"])
    FE = (2.0 / math.pi) * (rho / 2.0 - math.sin(2.0 * rho) / 4.0)
    q_ideal = 0.88 * 0.94 * SIGMA * 328.0 ** 4
    abs_ir = 0.88 * FE * 237.0
    abs_alb = 0.2 * FE * 0.30 * 1361.0 * 0.5
    q_env = q_ideal - abs_ir - abs_alb
    return dict(q_ideal=q_ideal, F_earth=FE, q_env=q_env, margin=p["q"] / q_env)


def perturbations(res):
    """Ordens de grandeza (hipoteses: bracos e dipolo residual da Sec. 3.4)."""
    return dict(SRP_N=4.56e-6 * 1.3 * res["A_pv"],
                aero_torque_Nm=res["F_drag"] * 5.0,
                mag_torque_Nm=1000.0 * 3.5e-5,
                gravity_gradient_Nm=res["T_gg"])


# ----------------------------- sensibilidade ------------------------------
def oat_sensitivity(orb, stats, base=None):
    base = base or DEFAULTS
    out = {}
    for k, (lo, hi) in RANGES.items():
        row = []
        for val in (lo, hi):
            pp = dict(base)
            pp[k] = val
            a = size_node(pp, stats["tecl_max_A"], orb)["M_launch"] / 1000.0
            b = size_node(pp, stats["tecl_B"], orb)["M_launch"] / 1000.0
            row.append((round(a, 1), round(b, 1)))
        out[k] = row
    return out


def monte_carlo(orb, stats, base=None, n=3000, seed=1):
    base = base or DEFAULTS
    random.seed(seed)
    A, B = [], []
    for _ in range(n):
        pp = dict(base)
        for k, (lo, hi) in RANGES.items():
            mode = MC_MODES.get(k, base[k])
            pp[k] = random.triangular(lo, hi, mode)
        A.append(size_node(pp, stats["tecl_max_A"], orb)["M_launch"] / 1000.0)
        B.append(size_node(pp, stats["tecl_B"], orb)["M_launch"] / 1000.0)
    A.sort()
    B.sort()

    def pct(a, q):
        return round(a[int(q * (len(a) - 1))], 1)
    return dict(A=[pct(A, 0.05), pct(A, 0.5), pct(A, 0.95)],
                B=[pct(B, 0.05), pct(B, 0.5), pct(B, 0.95)])


# ----------------------------- lancamento ---------------------------------
def launch_co2_floor():
    """Eq. 13: piso de CO2 de combustao por kg de carga (Tsiolkovsky)."""
    c_rp1 = 0.8593 * 44.009 / 12.011      # kg CO2 / kg RP-1 (86% C)
    c_ch4 = 44.009 / 16.043               # kg CO2 / kg CH4

    def floor(dv, isp, xf, c):
        return (math.exp(dv * 1000.0 / (G_STD * isp)) - 1.0) * xf * c
    out = {}
    for dv in (9.0, 9.5):
        for isp in (330.0, 365.0):
            out["RP1_dv%.1f_isp%d" % (dv, isp)] = round(floor(dv, isp, 1.0 / 3.6, c_rp1), 2)
        for isp in (350.0, 380.0):
            out["CH4_dv%.1f_isp%d" % (dv, isp)] = round(floor(dv, isp, 1.0 / 4.6, c_ch4), 2)
    return out


# ----------------------------- confiabilidade -----------------------------
def failure_rate_llama3(n_unexpected=419, n_gpu=16384, days=54):
    """lambda por GPU (h^-1) a partir de Dubey et al. (2024) [36]."""
    return n_unexpected / (n_gpu * days * 24.0)


def binom_tail(n, p, k):
    s = 0.0
    lp, lq = math.log(p), math.log(1.0 - p)
    for j in range(k, n + 1):
        s += math.exp(math.lgamma(n + 1) - math.lgamma(j + 1)
                      - math.lgamma(n - j + 1) + j * lp + (n - j) * lq)
    return s


def n_min(n_req, p_s, conf=0.99):
    lo, hi = n_req, n_req * 4
    while lo < hi:
        mid = (lo + hi) // 2
        if binom_tail(mid, p_s, n_req) >= conf:
            hi = mid
        else:
            lo = mid + 1
    return lo


def reliability_table(n_req=1000, hours=87600.0):
    """Eq. 10 / Tabela 3."""
    lam = failure_rate_llama3()
    rows = []
    for phi in (0.05, 0.10, 0.20, 0.30):
        ps = math.exp(-phi * lam * hours)
        nm = n_min(n_req, ps)
        rows.append(dict(phi=phi, p_s=round(ps, 3), n_min=nm,
                         oversize_pct=round(100.0 * (nm / n_req - 1.0))))
    return dict(lambda_per_h=lam, rows=rows)


# ----------------------------- linha de base terrestre --------------------
def water_tables():
    """Eq. 14-15 / Tabela 4 e agua indireta da geracao [8]."""
    e_it_gwh = 1.0 * 87600.0 / 1000.0
    e_fac = (e_it_gwh * 1.2, e_it_gwh * 1.5)
    direct = {}
    for name, wue in (("base_0.32-0.67", (0.32, 0.67)),
                      ("autores_1-2", (1.0, 2.0)),
                      ("extremo_5", (5.0, 5.0))):
        direct[name] = (round(e_it_gwh * wue[0], 3), round(e_it_gwh * wue[1], 3))
    indirect = {}
    for name, f in (("eolica_0.28", 0.28), ("carvao_2.33", 2.33),
                    ("hidreletrica_70", 70.0)):
        indirect[name] = (round(e_fac[0] * f, 1), round(e_fac[1] * f, 1))
    return dict(E_IT_GWh=e_it_gwh, E_fac_GWh=e_fac,
                direct_Mlitres=direct, indirect_Mlitres=indirect)


def carbon_ceiling(m_launch_A, m_launch_B):
    """Tabela 5: teto de G_inc (G_ODC,yr = 0), em t CO2-eq e kg/kg."""
    e_fac = (87.6 * 1.2, 87.6 * 1.5)       # GWh
    rows = []
    for I in (450, 350, 220, 50):          # g CO2/kWh  (GWh*g/kWh = t)
        tot = (e_fac[0] * I, e_fac[1] * I)
        rows.append(dict(
            I_gCO2_kWh=I,
            annual_t=(round(tot[0] / 10.0), round(tot[1] / 10.0)),
            total_10y_t=(round(tot[0]), round(tot[1])),
            kg_per_kg_A=(round(tot[0] * 1000.0 / m_launch_A, 1),
                         round(tot[1] * 1000.0 / m_launch_A, 1)),
            kg_per_kg_B=(round(tot[0] * 1000.0 / m_launch_B, 1),
                         round(tot[1] * 1000.0 / m_launch_B, 1))))
    return rows


def payback_constant(G_inc, I0_kg_kWh, PUE, G_odc_yr=0.0):
    """Eq. 18 (intensidade constante). G em kg CO2-eq; retorna anos."""
    E_it_yr = 8.76e6                                   # kWh/ano
    den = E_it_yr * PUE * I0_kg_kWh - G_odc_yr
    return G_inc / den if den > 0 else float("inf")


def payback_exp_decay(G_inc, A, gamma):
    """Eq. 20 (G_ODC,yr = 0; A em kg/ano; gamma em 1/ano)."""
    x = gamma * G_inc / A
    return -math.log(1.0 - x) / gamma if x < 1.0 else float("inf")


# ----------------------------- fechamento discreto ------------------------
def soc_simulation(res, p, orb, dt=10.0, n_orb=5, soc0=0.65):
    """Verificacao em tempo discreto do SOC (Secao 3.2)."""
    E_bat = res["E_bat_kWh"] * 1e3
    soc = soc0 * E_bat
    t_sun_s = res["t_sun"] * 60.0
    n_steps = int(orb["T_min"] * 60.0 / dt)
    P_gen = res["P_gen_MW"] * 1e6
    out = []
    for _ in range(n_orb):
        lo = hi = soc
        for k in range(n_steps):
            if k * dt < t_sun_s:
                ch = (P_gen - p["Pload"] / p["eta_c"]) * p["eta_c"] * p["eta_rt"]
                soc = min(E_bat, soc + ch * dt / 3600.0)
            else:
                soc = soc - p["Pload"] * dt / 3600.0
            lo, hi = min(lo, soc), max(hi, soc)
        out.append((round(lo / E_bat, 3), round(hi / E_bat, 3), round(soc / E_bat, 3)))
    return out


# ----------------------------- orquestracao e testes ----------------------
def run_all(p=None, mc_n=3000):
    p = p or DEFAULTS
    orb = orbit_geometry(700.0)
    st = eclipse_stats(orb)
    A = size_node(p, st["tecl_max_A"], orb)
    B = size_node(p, st["tecl_B"], orb)
    return dict(orbit=orb, eclipse=st, caseA=A, caseB=B,
                radiator=radiator_check(orb, p), perturb_B=perturbations(B),
                oat=oat_sensitivity(orb, st, p),
                monte_carlo=monte_carlo(orb, st, p, n=mc_n),
                launch_floor=launch_co2_floor(),
                reliability=reliability_table(), water=water_tables(),
                ceiling=carbon_ceiling(A["M_launch"], B["M_launch"]),
                soc_A=soc_simulation(A, p, orb), soc_B=soc_simulation(B, p, orb))


def self_check(R):
    """Compara com os valores publicados no artigo."""
    exp = [
        ("T_orb_min", 98.77, R["orbit"]["T_min"], 0.01),
        ("inc_deg", 98.19, R["orbit"]["inc_deg"], 0.01),
        ("beta_star_deg", 64.3, R["orbit"]["beta_star_deg"], 0.05),
        ("Tecl_A_min", 18.8, R["eclipse"]["tecl_max_A"], 0.05),
        ("Tecl_B_min", 35.3, R["eclipse"]["tecl_B"], 0.05),
        ("dias_eclipse_A", 83, R["eclipse"]["ndays_eclipse_A"], 0.5),
        ("Mdry_A_kg", 51752, R["caseA"]["M_dry"], 2.0),
        ("Mdry_B_kg", 67504, R["caseB"]["M_dry"], 2.0),
        ("Mlaunch_A_kg", 52388, R["caseA"]["M_launch"], 2.0),
        ("Mlaunch_B_kg", 68316, R["caseB"]["M_launch"], 2.0),
        ("MC_A_P50_t", 62.1, R["monte_carlo"]["A"][1], 0.05),
        ("MC_B_P50_t", 79.8, R["monte_carlo"]["B"][1], 0.05),
        ("n_min_phi5", 1114, R["reliability"]["rows"][0]["n_min"], 0.5),
        ("n_min_phi30", 1760, R["reliability"]["rows"][3]["n_min"], 0.5),
        ("floor_RP1_max", 15.6, R["launch_floor"]["RP1_dv9.5_isp330"], 0.05),
        ("floor_CH4_min", 6.08, R["launch_floor"]["CH4_dv9.0_isp380"], 0.05),
    ]
    return [dict(item=n, expected=e, got=round(g, 3), ok=abs(e - g) <= t)
            for (n, e, g, t) in exp]


def main():
    R = run_all()
    print("=== Geometria ===")
    print(json.dumps(R["orbit"], indent=1))
    print(json.dumps(R["eclipse"], indent=1))
    print("\n=== Tabela 2: orcamento (A = LTAN 06:00; B = LTAN 12:00) ===")
    for k in R["caseA"]:
        print("%-12s A=%14.3f   B=%14.3f" % (k, R["caseA"][k], R["caseB"][k]))
    print("\n=== Radiador ===")
    print(json.dumps(R["radiator"], indent=1))
    print("\n=== Perturbacoes (caso B) ===")
    print(json.dumps(R["perturb_B"], indent=1))
    print("\n=== Tabela 3: confiabilidade ===")
    print("lambda = %.3e 1/h" % R["reliability"]["lambda_per_h"])
    for row in R["reliability"]["rows"]:
        print(row)
    print("\n=== Tabela 4: agua ===")
    print(json.dumps(R["water"], indent=1))
    print("\n=== Tabela 5: teto de G_inc ===")
    for row in R["ceiling"]:
        print(row)
    print("\n=== Piso de CO2 de lancamento (kg CO2/kg) ===")
    print(json.dumps(R["launch_floor"], indent=1))
    print("\n=== OAT (massa de lancamento, t) [(A,B) em lo, hi] ===")
    for k, v in R["oat"].items():
        print("%-9s %s" % (k, v))
    print("\n=== Monte Carlo P5/P50/P95 (t) ===")
    print(R["monte_carlo"])
    print("\n=== SOC discreto (min, max, final) na 5a orbita ===")
    print("A:", R["soc_A"][-1], " B:", R["soc_B"][-1])
    print("\n=== Conferencia com o artigo ===")
    chk = self_check(R)
    for c in chk:
        print("%-16s esperado=%-9s obtido=%-10s %s" %
              (c["item"], c["expected"], c["got"], "OK" if c["ok"] else "FALHA"))
    print("TODOS OK:", all(c["ok"] for c in chk))


if __name__ == "__main__":
    main()
