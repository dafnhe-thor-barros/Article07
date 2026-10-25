#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
odc_lca_figures.py
Figuras do artigo "Data Centers Orbitais de IA: Estrutura de LCA...":
  1) fig_beta_eclipse.pdf/png    - |beta|(t) anual e eclipse por orbita (caso A)
  2) fig_eclipse_vs_beta.pdf/png - eclipse em funcao de |beta| (h = 700 km)
  3) fig_tornado.pdf/png         - tornado OAT da massa de lancamento (A e B)

Requer: matplotlib, numpy e odc_lca_model.py na mesma pasta.
Uso:
    python odc_lca_figures.py              # grava na pasta atual
    python odc_lca_figures.py figuras      # grava na pasta "figuras"
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from odc_lca_model import (DEFAULTS, RANGES, orbit_geometry, eclipse_min,
                           beta_series_dawn_dusk, eclipse_stats, size_node)

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.axisbelow": True})

MONTH_START = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
MONTH_NAME = ["Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
              "Jul", "Ago", "Set", "Out", "Nov", "Dez"]

# rotulo (mathtext), fator de escala para exibicao, unidade
LABELS = {
    "rho_pv":   (r"$\rho_{PV}$", 1.0, "kg/m$^2$"),
    "rho_rad":  (r"$\rho_{rad}$", 1.0, "kg/m$^2$"),
    "kIT":      (r"$M_{IT}$ espec.", 1.0, "kg/kW"),
    "k_arr":    (r"$k_{arr}$", 1.0, ""),
    "kst":      (r"$k_{st}$ (estrutura)", 1.0, ""),
    "DoD":      (r"DoD", 1.0, ""),
    "eta_cell": (r"$\eta_{cell}$", 1.0, ""),
    "wb":       (r"$e_b$", 1.0, "Wh/kg"),
    "q":        (r"$q$ (radiador)", 1.0, "W/m$^2$"),
    "mga":      (r"margem de massa", 1.0, ""),
    "redund":   (r"$N_{red}$ (CMG)", 1.0, ""),
    "rho_air":  (r"$\rho_{atm}$ (700 km)", 1.0, "kg/m$^3$"),
    "Pload":    (r"$P_{load}$", 1e-6, "MW"),
}


def fmt(v):
    """Numero com virgula decimal."""
    return ("%g" % v).replace(".", ",")


def month_axis(ax):
    ax.set_xticks(MONTH_START)
    ax.set_xticklabels(MONTH_NAME)
    ax.set_xlim(0, 364)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------- Figura 1: |beta|(t) e eclipse (caso A) ------------------
def fig_beta_eclipse(orb, st):
    betas = np.abs(np.array(beta_series_dawn_dusk(orb)))   # |beta|
    days = np.arange(len(betas))
    ecl = np.array([eclipse_min(b, orb) for b in betas])
    bstar = orb["beta_star_deg"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1.2, 1]))
    ax1.plot(days, betas, lw=1.6, color="tab:blue",
             label=r"$|\beta|(t)$, caso A (LTAN 06:00)")
    ax1.axhline(bstar, ls="--", lw=1.1, color="tab:red",
                label=r"$\beta^*=%s^\circ$" % fmt(round(bstar, 1)))
    ax1.fill_between(days, betas, bstar, where=betas < bstar,
                     color="tab:orange", alpha=0.3, label="dias com eclipse")
    imin = int(np.argmin(betas))
    ax1.annotate("mín. %s°" % fmt(round(float(betas[imin]), 1)),
                 xy=(imin, betas[imin]), xytext=(imin - 90, betas[imin] + 6),
                 arrowprops=dict(arrowstyle="->", lw=0.8))
    ax1.set_ylabel(r"$|\beta|$ (graus)")
    ax1.set_ylim(55, 92)
    ax1.legend(loc="lower center", fontsize=8)

    ax2.plot(days, ecl, lw=1.6, color="tab:green",
             label="caso A: eclipse por órbita")
    ax2.axhline(st["tecl_B"], ls=":", lw=1.4, color="tab:purple",
                label="caso B (LTAN 12:00): %s min" % fmt(round(st["tecl_B"], 1)))
    ax2.set_ylabel("Eclipse máximo (min)")
    ax2.set_ylim(0, 40)
    ax2.legend(loc="center", fontsize=8)
    month_axis(ax2)
    ax2.set_xlabel("Dia do ano")
    fig.suptitle("SSO a 700 km: ângulo beta anual e duração do eclipse",
                 fontsize=10)
    save(fig, "fig_beta_eclipse")


# ---------------- Figura 2: eclipse em funcao de beta ---------------------
def fig_eclipse_vs_beta(orb, st):
    bstar = orb["beta_star_deg"]
    bgrid = np.linspace(0.0, 75.0, 751)
    tgrid = np.array([eclipse_min(b, orb) for b in bgrid])
    bmin = min(abs(b) for b in beta_series_dawn_dusk(orb))

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(bgrid, tgrid, lw=1.8, color="tab:blue")
    ax.axvline(bstar, ls="--", lw=1.1, color="tab:red",
               label=r"$\beta^*=%s^\circ$" % fmt(round(bstar, 1)))
    ax.plot([bmin], [st["tecl_max_A"]], "o", color="tab:green",
            label="caso A (máx. anual): %s min" % fmt(round(st["tecl_max_A"], 1)))
    ax.plot([0.0], [st["tecl_B"]], "s", color="tab:purple",
            label=r"caso B ($\beta\approx0$): %s min" % fmt(round(st["tecl_B"], 1)))
    ax.set_xlabel(r"$|\beta|$ (graus)")
    ax.set_ylabel("Duração do eclipse (min)")
    ax.set_xlim(0, 75)
    ax.set_ylim(0, 40)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(r"Eclipse em função de $\beta$ (h = 700 km, sombra cilíndrica)",
                 fontsize=10)
    save(fig, "fig_eclipse_vs_beta")


# ---------------- Figura 3: tornado OAT -----------------------------------
def tornado_rows(orb, t_ecl, base):
    m0 = size_node(base, t_ecl, orb)["M_launch"] / 1000.0
    rows = []
    for k, (lo, hi) in RANGES.items():
        pl = dict(base)
        pl[k] = lo
        ph = dict(base)
        ph[k] = hi
        ml = size_node(pl, t_ecl, orb)["M_launch"] / 1000.0
        mh = size_node(ph, t_ecl, orb)["M_launch"] / 1000.0
        rows.append((k, lo, hi, ml, mh, abs(mh - ml)))
    rows.sort(key=lambda r: r[5])          # maior efeito fica no topo do barh
    return m0, rows


def draw_tornado(ax, m0, rows, titulo, top_n=10):
    rows = rows[-top_n:]
    ylabels = []
    for i, (k, lo, hi, ml, mh, _) in enumerate(rows):
        ax.barh(i, ml - m0, left=m0, color="tab:blue", height=0.6)
        ax.barh(i, mh - m0, left=m0, color="tab:red", height=0.6)
        lab, sc, un = LABELS.get(k, (k, 1.0, ""))
        ylabels.append("%s\n[%s–%s %s]" % (lab, fmt(lo * sc), fmt(hi * sc), un))
    ax.axvline(m0, color="k", lw=1.0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xlabel("Massa de lançamento (t)")
    ax.set_title("%s (base = %s t)" % (titulo, fmt(round(m0, 1))), fontsize=9)


def fig_tornado(orb, st, base=None):
    base = base or DEFAULTS
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.2))
    m0a, ra = tornado_rows(orb, st["tecl_max_A"], base)
    m0b, rb = tornado_rows(orb, st["tecl_B"], base)
    draw_tornado(axes[0], m0a, ra, "Caso A (LTAN 06:00)")
    draw_tornado(axes[1], m0b, rb, r"Caso B (LTAN $\approx$ 12:00)")
    handles = [plt.Rectangle((0, 0), 1, 1, color="tab:blue"),
               plt.Rectangle((0, 0), 1, 1, color="tab:red")]
    fig.legend(handles, ["valor mínimo da faixa", "valor máximo da faixa"],
               loc="lower center", ncol=2, fontsize=8, frameon=False)
    fig.suptitle("Sensibilidade um-por-vez (OAT) da massa de lançamento",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    save(fig, "fig_tornado")


def main():
    orb = orbit_geometry(700.0)
    st = eclipse_stats(orb)
    fig_beta_eclipse(orb, st)
    fig_eclipse_vs_beta(orb, st)
    fig_tornado(orb, st)
    print("Geradas em '%s': fig_beta_eclipse, fig_eclipse_vs_beta, fig_tornado "
          "(.pdf/.png)" % OUT)
    print("Eclipse A max = %.2f min | B = %.2f min | dias com eclipse (A) = %d"
          % (st["tecl_max_A"], st["tecl_B"], st["ndays_eclipse_A"]))


if __name__ == "__main__":
    main()
