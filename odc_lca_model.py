#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
odc_lca_model.py
Modelo parametrico do no orbital de 1 MW (SSO 700 km) e calculos da linha de
base terrestre do artigo "Data Centers Orbitais de IA: Estrutura de LCA...".

Reproduz as Tabelas de massa, confiabilidade, agua, teto de emissoes
incorporadas, sensibilidade OAT e Monte Carlo (N=3000, semente 1) do artigo.
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

# constantes do radiador e das perturbacoes (P20 da Tabela A1)
EPS_IR = 0.88             # emissividade
ETA_FIN = 0.94            # eficiencia de aleta
T_RAD_K = 328.0           # temperatura do painel
IR_PLANETA = 237.0        # W/m^2, IR planetario medio
ABS_SOLAR = 0.2           # absortividade solar do radiador
ALBEDO = 0.30
FATOR_ILUM = 0.5
CP_SRP = 4.56e-6          # N/m^2
CR_SRP = 1.3
BRACO_AERO = 5.0          # m
DIPOLO_RES = 1000.0       # A.m^2
B_TERRA = 3.5e-5          # T

# intensidades de lancamento + reentrada, Ohs et al. [18], v3 (kg CO2e/kg)
I_LR = {"Falcon9": 51.7, "Starship": 34.4}

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

# faixas para OAT e Monte Carlo (P18); rho_rad ampliada para 4-22 kg/m2 para
# cobrir o equivalente de Ohs et al. [19] (~21,8 kg/m2)
RANGES = {
    "eta_cell": (0.27, 0.32), "k_arr": (0.70, 0.90), "rho_pv": (2.5, 4.5),
    "DoD": (0.25, 0.50), "wb": (120.0, 200.0), "q": (300.0, 420.0),
    "rho_rad": (4.0, 22.0), "kIT": (6.0, 12.0), "kst": (0.10, 0.25),
    "mga": (0.20, 0.30), "redund": (1.0, 3.0), "rho_air": (3e-14, 5e-13),
    "Pload": (1.05e6, 1.10e6),
}
MC_MODES = {"kIT": 8.0, "redund": 2.0}   # demais: modo = caso-base


# ----------------------------- geometria orbital --------------------------
def orbit_geometry(h_km=700.0):
    """Periodo e inclinacao heliossincrona (eq:orb)."""
    r = RE + h_km
    T = 2.0 * math.pi * math.sqrt(r ** 3 / MU)
    n = 2.0 * math.pi / T
    omega_sso = 2.0 * math.pi / YEAR_S
    cosi = -omega_sso / (1.5 * n * J2 * (RE / r) ** 2)
    return dict(h=h_km, r=r, T_min=T / 60.0, n=n,
                inc_deg=math.degrees(math.acos(cosi)),
                beta_star_deg=math.degrees(math.asin(RE / r)))


def eclipse_min(beta_deg, orb):
    """Duracao do eclipse, sombra cilindrica, em minutos (eq:ecl)."""
    if abs(beta_deg) >= orb["beta_star_deg"]:
        return 0.0
    h, r = orb["h"], orb["r"]
    x = math.sqrt(h * h + 2.0 * RE * h) / (r * math.cos(math.radians(beta_deg)))
    return orb["T_min"] / math.pi * math.acos(min(1.0, x))


def beta_series_dawn_dusk(orb, ltan_h=6.0):
    """beta(t) anual COM SINAL (graus), Sol em movimento medio.
    LTAN 06:00 => sin(Omega-alpha) = -1 ; LTAN 18:00 => +1."""
    s_node = math.sin(math.radians((ltan_h - 12.0) * 15.0))
    sini = math.sin(math.radians(orb["inc_deg"]))
    cosi = math.cos(math.radians(orb["inc_deg"]))
    out = []
    for d in range(365):
        lam = 2.0 * math.pi * (d - 79.5) / 365.25
        sd = math.sin(OBLIQ) * math.sin(lam)
        cd = math.sqrt(1.0 - sd * sd)
        sb = sini * cd * s_node + cosi * sd
        out.append(math.degrees(math.asin(max(-1.0, min(1.0, sb)))))
    return out


def eclipse_runs(orb, ltan_h=6.0):
    """Janelas (dia inicial, dia final) com eclipse; dia 0 = 1 de janeiro.
    Um dia final > 364 indica que a temporada atravessa a virada do ano."""
    betas = beta_series_dawn_dusk(orb, ltan_h)
    dias = [d for d, b in enumerate(betas) if eclipse_min(b, orb) > 0.0]
    s = set(dias)
    runs = []
    for a in dias:
        if (a - 1) % 365 not in s:
            b = a
            while (b + 1) % 365 in s:
                b += 1
            runs.append((a, b))
    return runs


def eclipse_stats(orb, ltan_h=6.0):
    betas = beta_series_dawn_dusk(orb, ltan_h)
    ecl = [eclipse_min(b, orb) for b in betas]
    ndays = sum(1 for e in ecl if e > 0.0)
    orb_per_day = 86400.0 / (orb["T_min"] * 60.0)
    absb = [abs(b) for b in betas]
    return dict(tecl_max_A=max(ecl), tecl_B=eclipse_min(0.0, orb),
                ndays_eclipse_A=ndays, beta_min=min(absb), beta_max=max(absb),
                orbits_10y=10.0 * 365.25 * orb_per_day,
                orbits_eclipse_A_10y=10.0 * ndays * orb_per_day)


# ----------------------------- dimensionamento ----------------------------
def size_node(p, t_ecl_min, orb, scale_it=False, rad_kg_per_kw=None):
    """Potencia, bateria, arranjo, radiador, ADCS, massa e arrasto
    (eq:bat a eq:mass). rad_kg_per_kw substitui rho_rad (calibracao [19]);
    scale_it escala M_IT com Pload."""
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
    if rad_kg_per_kw is None:
        M_rad = A_rad * p["rho_rad"]
    else:
        M_rad = rad_kg_per_kw * Q_avg / 1000.0
    M_it = p["kIT"] * 1000.0 * ((p["Pload"] / 1.05e6) if scale_it else 1.0)
    S = M_pv + M_rad + M_bat + M_it
    M_st = p["kst"] * S

    # gradiente de gravidade e momento ciclico (duas asas simetricas)
    m_wing = M_pv / 2.0
    L = (A_pv / 2.0) / p["Wwing"]
    yc = p["hub"] + L / 2.0
    I_xx = 2.0 * m_wing * (L * L / 12.0 + yc * yc)
    T_gg = 1.5 * MU * 1e9 / ((r * 1e3) ** 3) * I_xx
    h_cyc = T_gg / (2.0 * n)
    M_adcs = p["redund"] * h_cyc / p["hspec"]

    # massa seca e de lancamento (ponto fixo para o propelente de arrasto)
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
                A_pv=A_pv, M_pv=M_pv, A_rad=A_rad, M_rad=M_rad, Q_avg=Q_avg,
                M_it=M_it, M_st=M_st, I_xx=I_xx, T_gg=T_gg, h_cyc=h_cyc,
                M_adcs=M_adcs, M_hw=p["Mhw"], M_dry=M_dry, F_drag=F_drag,
                dv_drag=dv, m_prop=m_p, M_launch=M_w)


def radiator_check(orb, p):
    """Fluxo liquido por face com IR planetario e albedo (Secao termica)."""
    rho = math.asin(RE / orb["r"])
    FE = (2.0 / math.pi) * (rho / 2.0 - math.sin(2.0 * rho) / 4.0)
    q_ideal = EPS_IR * ETA_FIN * SIGMA * T_RAD_K ** 4
    abs_ir = EPS_IR * FE * IR_PLANETA
    abs_alb = ABS_SOLAR * FE * ALBEDO * G0_SUN * FATOR_ILUM
    q_env = q_ideal - abs_ir - abs_alb
    return dict(q_ideal=q_ideal, F_earth=FE, q_env=q_env, margin=p["q"] / q_env)


def perturbations(res):
    """Ordens de grandeza das perturbacoes (hipoteses da Secao de atitude)."""
    return dict(SRP_N=CP_SRP * CR_SRP * res["A_pv"],
                aero_torque_Nm=res["F_drag"] * BRACO_AERO,
                mag_torque_Nm=DIPOLO_RES * B_TERRA,
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
    rng = random.Random(seed)             # gerador isolado (nao altera o global)
    A, B = [], []
    for _ in range(n):
        pp = dict(base)
        for k, (lo, hi) in RANGES.items():
            mode = MC_MODES.get(k, base[k])
            pp[k] = rng.triangular(lo, hi, mode)
        A.append(size_node(pp, stats["tecl_max_A"], orb)["M_launch"] / 1000.0)
        B.append(size_node(pp, stats["tecl_B"], orb)["M_launch"] / 1000.0)
    A.sort()
    B.sort()

    def pct(a, q):
        return round(a[int(q * (len(a) - 1))], 1)
    return dict(A=[pct(A, 0.05), pct(A, 0.5), pct(A, 0.95)],
                B=[pct(B, 0.05), pct(B, 0.5), pct(B, 0.95)])


def monte_carlo_rho_rad(orb, stats, base=None, n=3000, seed=1):
    """Sensibilidade do Monte Carlo a distribuicao de rho_rad (tabela mcvar)."""
    base = base or DEFAULTS
    variantes = {
        "tri_4_9_moda_6": ("tri", 4.0, 9.0, 6.0),
        "tri_4_22_moda_6": ("tri", 4.0, 22.0, 6.0),
        "tri_4_22_moda_12": ("tri", 4.0, 22.0, 12.0),
        "tri_4_22_moda_21.8": ("tri", 4.0, 22.0, 21.8),
        "uni_4_22": ("uni", 4.0, 22.0, None),
    }
    out = {}
    for nome, (tipo, lo_r, hi_r, modo_r) in variantes.items():
        rng = random.Random(seed)
        A, B = [], []
        for _ in range(n):
            pp = dict(base)
            for k, (lo, hi) in RANGES.items():
                if k == "rho_rad":
                    if tipo == "uni":
                        pp[k] = rng.uniform(lo_r, hi_r)
                    else:
                        pp[k] = rng.triangular(lo_r, hi_r, modo_r)
                else:
                    pp[k] = rng.triangular(lo, hi, MC_MODES.get(k, base[k]))
            A.append(size_node(pp, stats["tecl_max_A"], orb)["M_launch"] / 1000.0)
            B.append(size_node(pp, stats["tecl_B"], orb)["M_launch"] / 1000.0)
        A.sort()
        B.sort()

        def pct(a, q):
            return round(a[int(q * (len(a) - 1))], 1)
        out[nome] = dict(A=[pct(A, 0.05), pct(A, 0.5), pct(A, 0.95)],
                         B=[pct(B, 0.05), pct(B, 0.5), pct(B, 0.95)])
    return out


# ----------------------------- lancamento ---------------------------------
def launch_co2_floor():
    """Piso de CO2 de combustao por kg de carga (eq:floor, Tsiolkovsky)."""
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
    if p >= 1.0:
        return 1.0
    if p <= 0.0:
        return 0.0
    s = 0.0
    lp, lq = math.log(p), math.log(1.0 - p)
    for j in range(k, n + 1):
        s += math.exp(math.lgamma(n + 1) - math.lgamma(j + 1)
                      - math.lgamma(n - j + 1) + j * lp + (n - j) * lq)
    return s


def n_min(n_req, p_s, conf=0.99):
    lo, hi = n_req, n_req * 4
    if binom_tail(hi, p_s, n_req) < conf:
        raise ValueError("n_min: confianca nao atingida com ate 4x unidades "
                         "(p_s=%.3f)" % p_s)
    while lo < hi:
        mid = (lo + hi) // 2
        if binom_tail(mid, p_s, n_req) >= conf:
            hi = mid
        else:
            lo = mid + 1
    return lo


def reliability_table(n_req=1000, hours=87600.0):
    """eq:rel / Tabela de sobredimensionamento."""
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
    """Agua direta e agua indireta da geracao [8] (milhoes de litros)."""
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
    """Teto de G_inc (G_ODC,yr = 0), em t CO2-eq e kg/kg."""
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
    """eq:payback (intensidade constante). G em kg CO2-eq; retorna anos."""
    E_it_yr = 8.76e6                                   # kWh/ano
    den = E_it_yr * PUE * I0_kg_kWh - G_odc_yr
    return G_inc / den if den > 0 else float("inf")


def payback_exp_decay(G_inc, A, gamma):
    """eq:paybackclosed (G_ODC,yr = 0; A em kg/ano; gamma em 1/ano)."""
    x = gamma * G_inc / A
    return -math.log(1.0 - x) / gamma if x < 1.0 else float("inf")


# ----------------------------- fechamento discreto ------------------------
def soc_simulation(res, p, orb, dt=10.0, n_orb=5, soc0=0.65):
    """Verificacao em tempo discreto do SOC."""
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


# ----------------------------- comparacao com [18] e [19] -----------------
def calibracao_ref19(orb, stats, base=None):
    """Massa de lancamento ao trocar parametros pelos de Ohs et al. [19]:
    manta 3,41 kg/m2; bateria 3,75 kg/kWh; radiador 3,03 m2/kW x 10 kg/m2."""
    base = base or DEFAULTS
    casos = (("base", base, None),
             ("manta_19", dict(base, rho_pv=3.41), None),
             ("bateria_19", dict(base, wb=1000.0 / 3.75), None),
             ("radiador_19", base, 30.3),
             ("todos_19", dict(base, rho_pv=3.41, wb=1000.0 / 3.75), 30.3))
    out = {}
    for nome, par, rk in casos:
        a = size_node(par, stats["tecl_max_A"], orb, rad_kg_per_kw=rk)
        b = size_node(par, stats["tecl_B"], orb, rad_kg_per_kw=rk)
        out[nome] = dict(A_t=round(a["M_launch"] / 1000.0, 1),
                         B_t=round(b["M_launch"] / 1000.0, 1),
                         M_launch_A_kg=a["M_launch"], M_launch_B_kg=b["M_launch"])
    return out


def lr_vs_teto(mA_kg, mB_kg):
    """Lancamento+reentrada de [18] como fracao do teto de G_inc."""
    e_fac = (87.6 * 1.2, 87.6 * 1.5)       # GWh
    out = dict(emissao_t={}, fracao_do_teto_pct={})
    for nome, i in I_LR.items():
        out["emissao_t"][nome] = (round(mA_kg * i / 1000.0), round(mB_kg * i / 1000.0))
    for I in (450, 50):                    # g CO2/kWh
        fr = [100.0 * m * i / 1000.0 / (e * I)
              for m in (mA_kg, mB_kg) for i in I_LR.values() for e in e_fac]
        out["fracao_do_teto_pct"][I] = (round(min(fr), 1), round(max(fr), 1))
    return out


# ----------------------------- orquestracao e testes ----------------------
def run_all(p=None, mc_n=3000, variantes=True):
    p = p or DEFAULTS
    orb = orbit_geometry(700.0)
    st = eclipse_stats(orb)
    A = size_node(p, st["tecl_max_A"], orb)
    B = size_node(p, st["tecl_B"], orb)
    R = dict(orbit=orb, eclipse=st, eclipse_runs=eclipse_runs(orb),
             caseA=A, caseB=B, radiator=radiator_check(orb, p),
             perturb_B=perturbations(B), oat=oat_sensitivity(orb, st, p),
             monte_carlo=monte_carlo(orb, st, p, n=mc_n),
             launch_floor=launch_co2_floor(), reliability=reliability_table(),
             water=water_tables(),
             ceiling=carbon_ceiling(A["M_launch"], B["M_launch"]),
             soc_A=soc_simulation(A, p, orb), soc_B=soc_simulation(B, p, orb),
             cal19=calibracao_ref19(orb, st, p),
             lr=lr_vs_teto(A["M_launch"], B["M_launch"]))
    if variantes:
        R["mc_rho_rad"] = monte_carlo_rho_rad(orb, st, p, n=mc_n)
    return R


def self_check(R):
    """Compara com os valores publicados no artigo."""
    mc = R["monte_carlo"]
    exp = [
        ("T_orb_min", 98.77, R["orbit"]["T_min"], 0.01),
        ("inc_deg", 98.19, R["orbit"]["inc_deg"], 0.01),
        ("beta_star_deg", 64.3, R["orbit"]["beta_star_deg"], 0.05),
        ("Tecl_A_min", 18.8, R["eclipse"]["tecl_max_A"], 0.05),
        ("Tecl_B_min", 35.3, R["eclipse"]["tecl_B"], 0.05),
        ("dias_eclipse_A", 83, R["eclipse"]["ndays_eclipse_A"], 0.5),
        ("beta_min_A", 58.4, R["eclipse"]["beta_min"], 0.05),
        ("beta_max_A", 89.9, R["eclipse"]["beta_max"], 0.05),
        ("inicio_eclipse_A_dia", 312, R["eclipse_runs"][0][0], 0.5),
        ("Mdry_A_kg", 51752, R["caseA"]["M_dry"], 2.0),
        ("Mdry_B_kg", 67504, R["caseB"]["M_dry"], 2.0),
        ("Mlaunch_A_kg", 52388, R["caseA"]["M_launch"], 2.0),
        ("Mlaunch_B_kg", 68316, R["caseB"]["M_launch"], 2.0),
        ("MC_A_P5_t", 57.2, mc["A"][0], 0.05),
        ("MC_A_P50_t", 71.0, mc["A"][1], 0.05),
        ("MC_A_P95_t", 91.3, mc["A"][2], 0.05),
        ("MC_B_P5_t", 73.0, mc["B"][0], 0.05),
        ("MC_B_P50_t", 89.3, mc["B"][1], 0.05),
        ("MC_B_P95_t", 111.1, mc["B"][2], 0.05),
        ("n_min_phi5", 1114, R["reliability"]["rows"][0]["n_min"], 0.5),
        ("n_min_phi30", 1760, R["reliability"]["rows"][3]["n_min"], 0.5),
        ("floor_RP1_max", 15.6, R["launch_floor"]["RP1_dv9.5_isp330"], 0.05),
        ("floor_CH4_min", 6.08, R["launch_floor"]["CH4_dv9.0_isp380"], 0.05),
        ("rad_q_env", 484.8, R["radiator"]["q_env"], 0.05),
        ("cal19_rad_A", 85.7, R["cal19"]["radiador_19"]["A_t"], 0.05),
        ("cal19_rad_B", 102.2, R["cal19"]["radiador_19"]["B_t"], 0.05),
        ("cal19_todos_A", 84.2, R["cal19"]["todos_19"]["A_t"], 0.05),
        ("cal19_todos_B", 98.6, R["cal19"]["todos_19"]["B_t"], 0.05),
        ("LR_F9_A_t", 2708, R["lr"]["emissao_t"]["Falcon9"][0], 1.0),
        ("LR_St_B_t", 2350, R["lr"]["emissao_t"]["Starship"][1], 1.0),
        ("LR_frac450_min", 3.0, R["lr"]["fracao_do_teto_pct"][450][0], 0.05),
        ("LR_frac50_max", 67.2, R["lr"]["fracao_do_teto_pct"][50][1], 0.05),
    ]
    if "mc_rho_rad" in R:
        v = R["mc_rho_rad"]
        exp += [
            ("MCvar_4_9_A_P50", 62.1, v["tri_4_9_moda_6"]["A"][1], 0.05),
            ("MCvar_unif_A_P50", 77.4, v["uni_4_22"]["A"][1], 0.05),
            ("MCvar_unif_B_P95", 118.8, v["uni_4_22"]["B"][2], 0.05),
            ("MCvar_m218_B_P50", 102.8, v["tri_4_22_moda_21.8"]["B"][1], 0.05),
        ]
    return [dict(item=n, expected=e, got=round(g, 3), ok=abs(e - g) <= t)
            for (n, e, g, t) in exp]


def main():
    R = run_all()
    print("=== Geometria ===")
    print(json.dumps(R["orbit"], indent=1))
    print(json.dumps(R["eclipse"], indent=1))
    print("janelas de eclipse (caso A, dia 0 = 1 jan):", R["eclipse_runs"])
    print("\n=== Orcamento (A = LTAN 06:00; B = LTAN 12:00) ===")
    for k in R["caseA"]:
        print("%-12s A=%14.3f   B=%14.3f" % (k, R["caseA"][k], R["caseB"][k]))
    print("\n=== Radiador ===")
    print(json.dumps(R["radiator"], indent=1))
    print("\n=== Perturbacoes (caso B) ===")
    print(json.dumps(R["perturb_B"], indent=1))
    print("\n=== Confiabilidade ===")
    print("lambda = %.3e 1/h" % R["reliability"]["lambda_per_h"])
    for row in R["reliability"]["rows"]:
        print(row)
    print("\n=== Agua ===")
    print(json.dumps(R["water"], indent=1))
    print("\n=== Teto de G_inc ===")
    for row in R["ceiling"]:
        print(row)
    print("\n=== Piso de CO2 de lancamento (kg CO2/kg) ===")
    print(json.dumps(R["launch_floor"], indent=1))
    print("\n=== OAT (massa de lancamento, t): [(A,B) em lo, hi] ===")
    for k, v in R["oat"].items():
        print("%-9s %s" % (k, v))
    print("\n=== Monte Carlo P5/P50/P95 (t) ===")
    print(R["monte_carlo"])
    print("\n=== Monte Carlo: variantes de rho_rad ===")
    for k, v in R["mc_rho_rad"].items():
        print("%-20s %s" % (k, v))
    print("\n=== SOC discreto (min, max, final) na 5a orbita ===")
    print("A:", R["soc_A"][-1], " B:", R["soc_B"][-1])
    print("\n=== Comparacao com [19] (massa de lancamento, t) ===")
    for k, v in R["cal19"].items():
        print("%-12s A=%5.1f  B=%5.1f" % (k, v["A_t"], v["B_t"]))
    print("\n=== Lancamento+reentrada [18] vs teto ===")
    print(json.dumps(R["lr"], indent=1))
    print("\n=== Conferencia com o artigo ===")
    chk = self_check(R)
    for c in chk:
        print("%-22s esperado=%-9s obtido=%-10s %s" %
              (c["item"], c["expected"], c["got"], "OK" if c["ok"] else "FALHA"))
    print("TODOS OK:", all(c["ok"] for c in chk))


if __name__ == "__main__":
    main()
