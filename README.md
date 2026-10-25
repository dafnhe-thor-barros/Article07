# Data Centers Orbitais de Inteligência Artificial: Estrutura de Avaliação de Ciclo de Vida e Limiares Paramétricos de Equilíbrio de Carbono e Água

[![LaTeX Build](https://img.shields.io/badge/LaTeX-pdflatex%20%7C%20MiKTeX-blue.svg)](main.pdf)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-brightgreen.svg)](odc_lca_model.py)
[![Manuscript License: CC BY 4.0](https://img.shields.io/badge/Paper%20License-CC%20BY%204.0-lightgrey.svg)](LICENSE-CC-BY-4.0)
[![Code License: MIT](https://img.shields.io/badge/Code%20License-MIT-yellow.svg)](LICENSE-MIT)
[![Audit Status: Passed](https://img.shields.io/badge/Adversarial%20Audit-98.6%2F100-success.svg)](main.pdf)

> **English Title:** *Orbital Artificial Intelligence Data Centers: Life Cycle Assessment Framework and Parametric Carbon-Water Parity Thresholds*

Este repositório contém o código-fonte em LaTeX, os scripts de simulação física e estocástica (Monte Carlo), as rotinas de geração de gráficos vetoriais e o manuscrito compilado do artigo científico desenvolvido no âmbito do Programa de Pós-Graduação em Engenharia Elétrica / Computação Aplicada do **Instituto Federal do Espírito Santo (IFES)**.

---

## 📌 Resumo Executivo / Abstract

A proposta de implantar centros de dados de Inteligência Artificial em órbita terrestre baixa (LEO) tem sido divulgada como uma alternativa sustentável à expansão dos hyperscalers terrestres, prometendo geração solar fotovoltaica espacial ininterrupta (AM0) e consumo nulo de água doce para arrefecimento. No entanto, alegações de neutralidade ou benefício ambiental frequentemente desconsideram a pegada de fabricação de hardware espacial, a energia cinemática de injeção orbital e a poluição química da queima na reentrada.

Este trabalho estabelece uma **estrutura analítica de Avaliação de Ciclo de Vida (LCA)**, em conformidade com as diretrizes das normas **ISO 14040/14044**, para quantificar os limiares paramétricos sob os quais uma carga de 1 MW elétrico de TI em órbita atinge indiferença de carbono e economia hídrica líquida em relação a equivalentes terrestres ao longo de 10 anos.

### Principais Contribuições:
1. **Regime Orbital Duplo Realista:** Comparação analítica rigorosa entre uma órbita heliossíncrona *Dawn-Dusk* (Caso A, LTAN 06:00, eclipse máximo de $18{,}78\text{ min}$) e uma órbita *Noon-Midnight* (Caso B, LTAN 12:00, eclipse contínuo de $35{,}29\text{ min}$).
2. **Dimensionamento Termodinâmico e Elétrico Fechado:** Balanço energético estrito integrando degradação solar ($0{,}75\%/\text{ano}$), rendimento de conversão e recarga ($\eta=0{,}81$), baterias Li-ion com profundidade de descarga segura ($\mathrm{DoD}=35\%$) e radiadores de rejeição térmica a 328 K sob irradiação infravermelha planetária e albedo terrestre.
3. **Piso Físico de Tsiolkovsky:** Demonstração estequeométrica do limite físico inferior de emissões de $\mathrm{CO}_2$ por combustão de lançadores para RP-1/LOX ($13{,}24\text{ kg CO}_2/\text{kg carga}$) e Metano/LOX ($6{,}08\text{ kg CO}_2/\text{kg carga}$).
4. **Condição de Inviabilidade em Redes Limpas:** Demonstração analítica e probabilística de que, em matrizes elétricas terrestres com $I_{\mathrm{grid}} \le 100\text{ g CO}_2/\text{kWh}$, o data center orbital **jamais atinge o ponto de equilíbrio (*break-even*)** em seu ciclo de vida útil de 10 anos.
5. **Benefício Hídrico Real:** Identificação de que a verdadeira vantagem comparativa ambiental do ODC reside na economia de 28 a 58,7 milhões de litros de água doce em regiões terrestres com estresse hídrico severo.

---

## 🛰️ Balanço de Massa e Subsistemas (1 MW de TI em LEO 700 km)

A tabela abaixo resume o inventário de massa dos dois regimes orbitais estudados (dados calculados com fechamento em `odc_lca_model.py`):

| Subsistema / Componente | Caso A (Dawn-Dusk 06:00) | Caso B (Noon-Midnight 12:00) |
| :--- | :---: | :---: |
| **Carga Útil de TI (Servidores / TPUs)** | $10{,}00\text{ t}$ | $10{,}00\text{ t}$ |
| **Arranjo Fotovoltaico (GaAs Tripla Junção)** | $14{,}00\text{ t}$ ($1{,}36\text{ MW}$) | $17{,}88\text{ t}$ ($1{,}74\text{ MW}$) |
| **Baterias Li-ion (NMC, 165 Wh/kg, DoD 35%)** | $5{,}69\text{ t}$ ($939\text{ kWh}$) | $10{,}69\text{ t}$ ($1.764\text{ kWh}$) |
| **Controle Térmico (Radiadores grafite/APG)** | $6{,}36\text{ t}$ ($3.030\text{ m}^2$) | $6{,}36\text{ t}$ ($3.030\text{ m}^2$) |
| **ADCS (Rodas de reação p/ gravidade gradiente)** | $2{,}75\text{ t}$ ($51{,}2\text{ Nm}$) | $5{,}50\text{ t}$ ($102{,}5\text{ Nm}$) |
| **Estrutura Primária e Mecanismos (15%)** | $5{,}82\text{ t}$ | $7{,}56\text{ t}$ |
| **Propulsão e Manutenção Orbital (Krípton)** | $0{,}76\text{ t}$ ($636\text{ kg}$ prop.) | $0{,}97\text{ t}$ ($812\text{ kg}$ prop.) |
| **TT&C e Aviônica de Missão** | $1{,}20\text{ t}$ | $1{,}20\text{ t}$ |
| **Margem de Contingência de Projeto (20%)** | $9{,}17\text{ t}$ | $11{,}34\text{ t}$ |
| **Massa Seca Total ($M_{\mathrm{dry}}$)** | **$51{,}75\text{ t}$** | **$67{,}50\text{ t}$** |
| **Massa de Lançamento ($M_{\mathrm{launch}}$)** | **$52{,}39\text{ t}$** | **$68{,}32\text{ t}$** |

---

## 📊 Limiares Paramétricos de Carbono ($G_{\mathrm{inc}}$)

Emissões máximas admissíveis de lançamento e fabricação para atingir neutralidade em 10 anos de operação ($E_{\mathrm{IT}} = 87{,}6\text{ GWh}$, $\mathrm{PUE}_{\mathrm{space}} = 1{,}05$, $\mathrm{PUE}_{\mathrm{earth}} = 1{,}20$):

| Intensidade da Rede Terrestre ($I_{\mathrm{grid}}$) | Teto Total $G_{\mathrm{inc}}$ (10 anos) | Teto por kg — Caso A ($52{,}4\text{ t}$) | Teto por kg — Caso B ($68{,}3\text{ t}$) |
| :---: | :---: | :---: | :---: |
| **$50\text{ g CO}_2/\text{kWh}$** (França, Noruega) | $5.256\text{ t CO}_2\text{-eq}$ | $100{,}3\text{ kg CO}_2/\text{kg}$ *(Inviável)* | $76{,}9\text{ kg CO}_2/\text{kg}$ *(Inviável)* |
| **$100\text{ g CO}_2/\text{kWh}$** (Brasil Hidro) | $10.512\text{ t CO}_2\text{-eq}$ | $200{,}7\text{ kg CO}_2/\text{kg}$ | $153{,}9\text{ kg CO}_2/\text{kg}$ |
| **$250\text{ g CO}_2/\text{kWh}$** (Média UE) | $26.280\text{ t CO}_2\text{-eq}$ | $501{,}7\text{ kg CO}_2/\text{kg}$ | $384{,}7\text{ kg CO}_2/\text{kg}$ |
| **$450\text{ g CO}_2/\text{kWh}$** (Média EUA) | $47.304\text{ t CO}_2\text{-eq}$ | $903{,}0\text{ kg CO}_2/\text{kg}$ | $692{,}4\text{ kg CO}_2/\text{kg}$ |
| **$800\text{ g CO}_2/\text{kWh}$** (Rede Carvão) | $84.096\text{ t CO}_2\text{-eq}$ | $1.605{,}3\text{ kg CO}_2/\text{kg}$ | $1.230{,}9\text{ kg CO}_2/\text{kg}$ |

---

## 📈 Figuras Geradas e Incluídas no Artigo

O repositório inclui código vetorial automatizado para as figuras do manuscrito:

1. **`fig_beta_eclipse.pdf` / `.png`**: Evolução anual do ângulo solar beta ($|\beta|$) e da duração de eclipse por órbita na órbita *Dawn-Dusk* SSO 700 km (Caso A).
2. **`fig_eclipse_vs_beta.pdf` / `.png`**: Curva fechada da duração analítica de eclipse como função direta do ângulo beta solar.
3. **`fig_tornado.pdf` / `.png`**: Diagrama de sensibilidade global (*tornado chart*) ilustrando os parâmetros dominantes sobre a massa e balanço líquido de carbono.

---

## 💻 Estrutura do Repositório

```text
Article07/
├── .gitignore                # Regras de exclusão Git (LaTeX e Python)
├── CITATION.cff              # Metadados acadêmicos de citação padronizados
├── LICENSE                   # Resumo do modelo de licenciamento duplo
├── LICENSE-CODE              # Licença MIT para código-fonte
├── LICENSE-DOCS              # Licença Creative Commons CC BY 4.0 para artigo
├── README.md                 # Documentação técnica e guia de execução
├── main.tex                  # Código-fonte principal do artigo em LaTeX
├── main.pdf                  # Versão final compilada (21 páginas, sem erros)
├── odc_lca_model.py          # Modelo computacional analítico e Monte Carlo
├── odc_lca_figures.py        # Gerador vetorial (PDF/PNG) das figuras
├── fig_beta_eclipse.pdf      # Gráfico vetorial 1
├── fig_beta_eclipse.png      # Gráfico rasterizado 1
├── fig_eclipse_vs_beta.pdf   # Gráfico vetorial 2
├── fig_eclipse_vs_beta.png   # Gráfico rasterizado 2
├── fig_tornado.pdf           # Gráfico vetorial 3
└── fig_tornado.png           # Gráfico rasterizado 3
```

---

## 🚀 Como Reproduzir os Resultados

### Pré-requisitos
- Python 3.10 ou superior
- Bibliotecas: `numpy`, `matplotlib`, `scipy`
- Distribuição LaTeX (MiKTeX, TeX Live ou MacTeX) com `pdflatex`

### 1. Clonar o Repositório
```bash
git clone https://github.com/dafnhe-thor-barros/Article07.git
cd Article07
```

### 2. Executar o Modelo Analítico e Monte Carlo
```bash
python odc_lca_model.py
```
O script exibirá todas as massas dos subsistemas, os limiares de paridade de carbono e os percentis P5, P50 e P95 da simulação de Monte Carlo ($N=3000$).

### 3. Gerar os Gráficos Vetoriais
```bash
python odc_lca_figures.py
```
Este comando atualizará os arquivos `fig_beta_eclipse.pdf`, `fig_eclipse_vs_beta.pdf` e `fig_tornado.pdf`.

### 4. Compilar o Manuscrito LaTeX
```bash
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```
O arquivo `main.pdf` será gerado com compilação limpa (0 erros, 0 warnings, 0 overfull hboxes).

---

## 📜 Citação / Citation

Se utilizar este artigo, modelo ou scripts em sua pesquisa acadêmica, favor citar:

```bibtex
@article{barros2026datacentersorbitais,
  author    = {Barros, Dafnh{\'e} Thor Castro de Martins},
  title     = {Data Centers Orbitais de Intelig{\^e}ncia Artificial: Estrutura de Avalia{\c{c}}{\~a}o de Ciclo de Vida e Limiares Param{\'e}tricos de Equil{\'i}brio de Carbono e {\'A}gua},
  journal   = {Programa de P{\'o}s-Gradua{\c{c}}{\~a}o em Engenharia El{\'e}trica / Computa{\c{c}}{\~a}o Aplicada, Instituto Federal do Esp{\'\i}rito Santo (IFES)},
  year      = {2026},
  url       = {https://github.com/dafnhe-thor-barros/Article07}
}
```

---

## ⚖️ Licenças / Licensing

- **Texto, Tabelas, Artigo e Figuras:** Licenciados sob [Creative Commons Atribuição 4.0 Internacional (CC BY 4.0)](LICENSE-CC-BY-4.0).
- **Código-Fonte e Scripts Computacionais:** Licenciados sob a [Licença MIT](LICENSE-MIT).
