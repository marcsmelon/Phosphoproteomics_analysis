#!/usr/bin/env python3
"""
Compare kinase activity: Human Heart vs iPSC (combined + individual)
Stouffer's method to combine iPSC z-scores across 4 genotypes.
Three outputs:
  1. Correlation scatter: Human vs iPSC combined
  2. Volcano: iPSC combined activity vs combined p-value
  3. 2x2 grid: Individual genotype vs Human
"""

import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import norm, spearmanr
from adjustText import adjust_text

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.size'] = 6

# ══════════════════════════════════════════════════════════════════════════════
# SETTINGS — change TIMEPOINT to "D30" for the second analysis
# ══════════════════════════════════════════════════════════════════════════════
TIMEPOINT = "D15"   # ← change to "D30" for second run

BASE_DIR = "/Users/marcmelon/Downloads/Python/phosphonetworks"

HUMAN_FILE = os.path.join(BASE_DIR, "outputs/Human_Heart_C2/kinase_activity_Human_Heart_C2.csv")

IPSC_FILES = {
    "TPM1":  os.path.join(BASE_DIR, f"outputs/TPM1_D15/kinase_activity_TPM1_D15.csv"),
    "TTN":   os.path.join(BASE_DIR, f"outputs/TTN_D15/kinase_activity_TTN_D15.csv"),
    "TNNT2": os.path.join(BASE_DIR, f"outputs/TNNT2_D15/kinase_activity_TNNT2_D15.csv"),
    "DSP":   os.path.join(BASE_DIR, f"outputs/DSP_D15/kinase_activity_DSP_D15.csv"),
}

OUT_DIR = os.path.join(BASE_DIR, f"outputs/Human_vs_iPSC_{TIMEPOINT}")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Colours ──
DCM_COL     = "#F5CD6A"
HEALTHY_COL = "#3A2044"
NS_COL      = "#CCCCCC"
DISCORDANT  = "#A8A8D0"

# ── Thresholds ──
SIG_Z = 1.0    # colour if |activity| > 1 on both human and iPSC sides (no p-value filter)
SIG_P = 0.05   # used only for volcano label selection

# ── Figure dimensions ──
FIG_W = 189.34 / 72
FIG_H = 169.38 / 72

# ══════════════════════════════════════════════════════════════════════════════
# HELPER
# ══════════════════════════════════════════════════════════════════════════════
def save_fig(fig, name):
    for ext, dpi in [('png', 300), ('pdf', None), ('jpeg', 300)]:
        path = os.path.join(OUT_DIR, f"{name}.{ext}")
        if dpi: fig.savefig(path, dpi=dpi, bbox_inches='tight')
        else:   fig.savefig(path, bbox_inches='tight')

def assign_colour(row, human_col, ipsc_col, human_p=None, ipsc_p=None):
    sig_h = abs(row[human_col]) > SIG_Z
    sig_i = abs(row[ipsc_col])  > SIG_Z
    if sig_h and sig_i:  # activity > 1 on both sides, no p-value filter
        dir_h = "dcm" if row[human_col] > 0 else "healthy"
        dir_i = "dcm" if row[ipsc_col]  > 0 else "healthy"
        if dir_h == dir_i:
            return DCM_COL if dir_h == "dcm" else HEALTHY_COL
        return DISCORDANT
    return NS_COL

# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA — filter to combined resource
# ══════════════════════════════════════════════════════════════════════════════
print(f"Loading data for {TIMEPOINT}...")

human = pd.read_csv(HUMAN_FILE)
if "kinase" in human.columns and "UniProt" not in human.columns:
    human = human.rename(columns={"kinase": "UniProt"})
if "pvalue" not in human.columns:
    human["pvalue"] = 2 * norm.sf(np.abs(human["activity"]))
# Keep best activity per kinase across ALL resources
human["abs_activity"] = human["activity"].abs()
human = human.sort_values("abs_activity", ascending=False).drop_duplicates(subset="UniProt")
human = human[["UniProt", "gene", "activity", "pvalue"]].copy()
human = human.rename(columns={"activity": "activity_human", "pvalue": "pvalue_human"})

ipsc_dfs = {}
for geno, fpath in IPSC_FILES.items():
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        if "kinase" in df.columns and "UniProt" not in df.columns:
            df = df.rename(columns={"kinase": "UniProt"})
        if "pvalue" not in df.columns:
            df["pvalue"] = 2 * norm.sf(np.abs(df["activity"]))
        # Keep best activity per kinase across ALL resources
        df["abs_activity"] = df["activity"].abs()
        df = df.sort_values("abs_activity", ascending=False).drop_duplicates(subset="UniProt")
        df = df[["UniProt", "activity", "pvalue"]].copy()
        df = df.rename(columns={"activity": f"activity_{geno}", "pvalue": f"pvalue_{geno}"})
        ipsc_dfs[geno] = df
        print(f"  Loaded {geno}: {len(df)} kinases (best across all resources)")
    else:
        print(f"  WARNING: {fpath} not found")

# ══════════════════════════════════════════════════════════════════════════════
# MERGE ALL
# ══════════════════════════════════════════════════════════════════════════════
merged = human.copy()
for geno, df in ipsc_dfs.items():
    merged = merged.merge(df, on="UniProt", how="inner")

geno_names = list(ipsc_dfs.keys())
n_geno = len(geno_names)
activity_cols = [f"activity_{g}" for g in geno_names]
print(f"  Kinases in common: {len(merged)}")

# ══════════════════════════════════════════════════════════════════════════════
# STOUFFER'S COMBINED Z-SCORE
# ══════════════════════════════════════════════════════════════════════════════
print("Computing Stouffer's combined z-score...")

merged["activity_iPSC_mean"]     = merged[activity_cols].mean(axis=1)
merged["activity_iPSC_stouffer"] = merged["activity_iPSC_mean"] * np.sqrt(n_geno)
merged["pvalue_iPSC"]            = 2 * norm.sf(np.abs(merged["activity_iPSC_stouffer"]))
merged["neg_log10_p_iPSC"]      = -np.log10(merged["pvalue_iPSC"].clip(lower=1e-300))

# Consistency across genotypes
merged["n_positive"]  = (merged[activity_cols] > 0).sum(axis=1)
merged["n_negative"]  = (merged[activity_cols] < 0).sum(axis=1)
merged["consistency"] = merged[["n_positive", "n_negative"]].max(axis=1)

# Concordance colour
merged["colour"] = merged.apply(
    lambda r: assign_colour(r, "activity_human", "activity_iPSC_stouffer",
                            "pvalue_human", "pvalue_iPSC"), axis=1)

# Save table
save_cols = ["gene", "UniProt", "activity_human", "pvalue_human",
             "activity_iPSC_stouffer", "pvalue_iPSC", "consistency"] + activity_cols
merged[save_cols].sort_values("activity_iPSC_stouffer", ascending=False).to_csv(
    os.path.join(OUT_DIR, f"kinase_comparison_human_vs_iPSC_{TIMEPOINT}.csv"), index=False)

conc_dcm = (merged["colour"] == DCM_COL).sum()
conc_hlt = (merged["colour"] == HEALTHY_COL).sum()
disc     = (merged["colour"] == DISCORDANT).sum()
print(f"  Concordant DCM: {conc_dcm} | Concordant Healthy: {conc_hlt} | Discordant: {disc}")

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1: CORRELATION SCATTER — Human vs iPSC combined
# ══════════════════════════════════════════════════════════════════════════════
print("\nPlot 1: Correlation scatter...")

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

for col, zorder in [(NS_COL, 1), (DISCORDANT, 2), (HEALTHY_COL, 3), (DCM_COL, 3)]:
    sub = merged[merged["colour"] == col]
    if len(sub) == 0: continue
    ax.scatter(sub["activity_human"], sub["activity_iPSC_stouffer"],
               color=col, alpha=0.7 if col != NS_COL else 0.3,
               s=4 if col != NS_COL else 2,
               edgecolors='none', zorder=zorder)

# Label top concordant hits
top_conc = merged[merged["colour"].isin([DCM_COL, HEALTHY_COL])].copy()
top_conc["abs_sum"] = top_conc["activity_human"].abs() + top_conc["activity_iPSC_stouffer"].abs()
top_label = top_conc.nlargest(15, "abs_sum")
texts = []
for _, row in top_label.iterrows():
    texts.append(ax.text(row["activity_human"], row["activity_iPSC_stouffer"],
                         f"$\\it{{{row['gene']}}}$", fontsize=6, fontfamily='Arial'))
if texts:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.3, color='grey'))

ax.axhline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
ax.axvline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])]
ax.plot(lims, lims, '--', color='grey', linewidth=0.5, alpha=0.5)

corr, pval = spearmanr(merged["activity_human"], merged["activity_iPSC_stouffer"])
ax.text(0.03, 0.97, f'ρ = {corr:.3f}\np = {pval:.1e}',
        transform=ax.transAxes, fontsize=6, va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='grey', linewidth=0.3, alpha=0.8))

ax.set_xlabel('Human Heart Activity (z-score)', fontsize=6)
ax.set_ylabel(f'iPSC Combined Activity — {TIMEPOINT} (Stouffer z)', fontsize=6)
ax.tick_params(labelsize=6, width=0.4, length=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for s in ['left', 'bottom']: ax.spines[s].set_linewidth(0.5)

legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor=DCM_COL,     markersize=4, label='Concordant DCM'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=HEALTHY_COL, markersize=4, label='Concordant Donor'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=DISCORDANT,  markersize=4, label='Discordant'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor=NS_COL,      markersize=4, label='NS'),
]
ax.legend(handles=legend_elements, fontsize=6, loc='lower right',
          frameon=False, handletextpad=0.3)

plt.tight_layout()
save_fig(fig, f"scatter_human_vs_iPSC_{TIMEPOINT}")
plt.close()

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2: VOLCANO — iPSC combined kinase activity
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 2: Volcano — iPSC combined...")

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

for col, zorder in [(NS_COL, 1), (DISCORDANT, 2), (HEALTHY_COL, 3), (DCM_COL, 3)]:
    sub = merged[merged["colour"] == col]
    if len(sub) == 0: continue
    ax.scatter(sub["activity_iPSC_stouffer"], sub["neg_log10_p_iPSC"],
               color=col, alpha=0.7 if col != NS_COL else 0.3,
               s=4 if col != NS_COL else 2,
               edgecolors='none', zorder=zorder)

sig_iPSC = merged[(merged["pvalue_iPSC"] < SIG_P) & (merged["activity_iPSC_stouffer"].abs() > SIG_Z)]
top_sig = sig_iPSC.nlargest(10, "neg_log10_p_iPSC")
texts = []
for _, row in top_sig.iterrows():
    texts.append(ax.text(row["activity_iPSC_stouffer"], row["neg_log10_p_iPSC"],
                         f"$\\it{{{row['gene']}}}$", fontsize=6, fontfamily='Arial'))
if texts:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.3, color='grey'))

ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=0.5, alpha=0.7)
ax.axvline(0, color='black', linestyle='-', linewidth=0.5, alpha=0.7)

ax.set_xlabel(f'iPSC Combined Kinase Activity — {TIMEPOINT} (Stouffer z)', fontsize=6)
ax.set_ylabel('$-$log$_{10}$(p-value)', fontsize=6)
ax.tick_params(labelsize=6, width=0.4, length=2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
for s in ['left', 'bottom']: ax.spines[s].set_linewidth(0.5)
ax.legend(handles=legend_elements, fontsize=6, loc='upper right',
          frameon=False, handletextpad=0.3)

plt.tight_layout()
save_fig(fig, f"volcano_iPSC_combined_{TIMEPOINT}")
plt.close()
# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3: 2x2 GRID — Individual genotype vs Human (independent merges)
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 3: Individual genotype scatter grid (2x2)...")

fig, axes = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
fig.subplots_adjust(hspace=0.5, wspace=0.4,
                    left=0.12, right=0.97, top=0.93, bottom=0.08)

pair_cache = {}  # cache pairs for reuse in Plot 4

for idx, geno in enumerate(geno_names):
    ax = axes[idx // 2, idx % 2]

    # Fresh load for this genotype only — all resources, best per kinase
    ipsc_single = pd.read_csv(IPSC_FILES[geno])
    if "kinase" in ipsc_single.columns and "UniProt" not in ipsc_single.columns:
        ipsc_single = ipsc_single.rename(columns={"kinase": "UniProt"})
    if "pvalue" not in ipsc_single.columns:
        ipsc_single["pvalue"] = 2 * norm.sf(np.abs(ipsc_single["activity"]))
    ipsc_single["abs_activity"] = ipsc_single["activity"].abs()
    ipsc_single = ipsc_single.sort_values("abs_activity", ascending=False).drop_duplicates(subset="UniProt")
    ipsc_single = ipsc_single[["UniProt", "activity", "pvalue"]].copy()
    ipsc_single = ipsc_single.rename(columns={"activity": "activity_iPSC", "pvalue": "pvalue_iPSC"})

    # Reload human fresh
    human_fresh = pd.read_csv(HUMAN_FILE)
    if "kinase" in human_fresh.columns and "UniProt" not in human_fresh.columns:
        human_fresh = human_fresh.rename(columns={"kinase": "UniProt"})
    if "pvalue" not in human_fresh.columns:
        human_fresh["pvalue"] = 2 * norm.sf(np.abs(human_fresh["activity"]))
    human_fresh["abs_activity"] = human_fresh["activity"].abs()
    human_fresh = human_fresh.sort_values("abs_activity", ascending=False).drop_duplicates(subset="UniProt")
    human_fresh = human_fresh[["UniProt", "gene", "activity", "pvalue"]].copy()
    human_fresh = human_fresh.rename(columns={"activity": "activity_human", "pvalue": "pvalue_human"})

    # Merge just these two
    pair = human_fresh.merge(ipsc_single, on="UniProt", how="inner")

    # Assign colours
    pair["colour"] = pair.apply(
        lambda r: assign_colour(r, "activity_human", "activity_iPSC",
                                "pvalue_human", "pvalue_iPSC"), axis=1)
    pair_cache[geno] = pair  # save for Plot 4

    for col, zorder in [(NS_COL, 1), (DISCORDANT, 2), (HEALTHY_COL, 3), (DCM_COL, 3)]:
        sub = pair[pair["colour"] == col]
        if len(sub) == 0: continue
        ax.scatter(sub["activity_human"], sub["activity_iPSC"],
                   color=col, alpha=0.7 if col != NS_COL else 0.5,
                   s=2 if col != NS_COL else 2,
                   edgecolors='none', zorder=zorder)

    # Top 3 from DCM + top 3 from Donor labelled separately
    texts = []
    for grp_col in [DCM_COL, HEALTHY_COL]:
        top_g = pair[pair["colour"] == grp_col].copy()
        if len(top_g) == 0:
            continue
        top_g["abs_sum"] = top_g["activity_human"].abs() + top_g["activity_iPSC"].abs()
        top_label = top_g.nlargest(3, "abs_sum")
        for _, row in top_label.iterrows():
            texts.append(ax.text(row["activity_human"], row["activity_iPSC"],
                                 f"$\\it{{{row['gene']}}}$", fontsize=6, fontfamily='Arial'))
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.2, color='grey'))

    ax.axhline(0, color='black', linestyle='-', linewidth=0.3, alpha=0.5)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.3, alpha=0.5)

    corr, pval = spearmanr(pair["activity_human"], pair["activity_iPSC"])
    ax.text(0.03, 0.97, f'ρ={corr:.2f}\nn={len(pair)}',
            transform=ax.transAxes, fontsize=6, va='top', ha='left')

    ax.set_title(geno, fontsize=6, fontweight='bold', pad=2)
    ax.tick_params(labelsize=6, width=0.3, length=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for s in ['left', 'bottom']: ax.spines[s].set_linewidth(0.3)

fig.text(0.52, 0.01, 'Human Heart Activity (z-score)', ha='center', fontsize=6)
fig.text(0.01, 0.5, f'iPSC Activity — {TIMEPOINT} (z-score)', va='center',
         rotation='vertical', fontsize=6)
fig.legend(handles=legend_elements, fontsize=6, loc='lower center',
           frameon=False, bbox_to_anchor=(0.52, -0.03), ncol=4)

plt.tight_layout()
save_fig(fig, f"scatter_individual_genotypes_2x2_{TIMEPOINT}")
plt.close()

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 4: 1x4 VERTICAL — Individual genotype vs Human
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 4: Individual genotype scatter grid (1x4 vertical)...")

fig, axes = plt.subplots(4, 1, figsize=(FIG_W / 2, FIG_H * 2))
fig.subplots_adjust(hspace=0.6,
                    left=0.15, right=0.97, top=0.96, bottom=0.06)

for idx, geno in enumerate(geno_names):
    ax = axes[idx]
    pair = pair_cache[geno]  # reuse from Plot 3

    for col, zorder in [(NS_COL, 1), (DISCORDANT, 2), (HEALTHY_COL, 3), (DCM_COL, 3)]:
        sub = pair[pair["colour"] == col]
        if len(sub) == 0: continue
        ax.scatter(sub["activity_human"], sub["activity_iPSC"],
                   color=col, alpha=0.7 if col != NS_COL else 0.5,
                   s=2 if col != NS_COL else 2,
                   edgecolors='none', zorder=zorder)

    # Top 3 from DCM + top 3 from Donor labelled separately
    texts = []
    for grp_col in [DCM_COL, HEALTHY_COL]:
        top_g = pair[pair["colour"] == grp_col].copy()
        if len(top_g) == 0:
            continue
        top_g["abs_sum"] = top_g["activity_human"].abs() + top_g["activity_iPSC"].abs()
        top_label = top_g.nlargest(3, "abs_sum")
        for _, row in top_label.iterrows():
            texts.append(ax.text(row["activity_human"], row["activity_iPSC"],
                                 f"$\\it{{{row['gene']}}}$", fontsize=6, fontfamily='Arial'))
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.2, color='grey'))

    ax.axhline(0, color='black', linestyle='-', linewidth=0.3, alpha=0.5)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.3, alpha=0.5)

    corr, pval = spearmanr(pair["activity_human"], pair["activity_iPSC"])
    ax.text(0.03, 0.97, f'ρ={corr:.2f}\nn={len(pair)}',
            transform=ax.transAxes, fontsize=6, va='top', ha='left')

    ax.set_title(geno, fontsize=6, fontweight='bold', pad=2)
    ax.tick_params(labelsize=6, width=0.3, length=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for s in ['left', 'bottom']: ax.spines[s].set_linewidth(0.3)

fig.text(0.52, 0.01, 'Human Heart Activity (z-score)', ha='center', fontsize=6)
fig.text(0.01, 0.5, f'iPSC Activity — {TIMEPOINT} (z-score)', va='center',
         rotation='vertical', fontsize=6)
fig.legend(handles=legend_elements, fontsize=6, loc='lower center',
           frameon=False, bbox_to_anchor=(0.52, -0.02), ncol=4)

plt.tight_layout()
save_fig(fig, f"scatter_individual_genotypes_1x4_{TIMEPOINT}")
plt.close()

print(f"\n✅ Done! Outputs in: {OUT_DIR}/")
