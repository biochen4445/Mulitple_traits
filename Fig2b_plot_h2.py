import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from adjustText import adjust_text

F = 'Mutli_Traits_TableS_2026_09_28.xlsx'
OUT = 'Fig2b_h2SNP'

st1 = pd.read_excel(F, sheet_name='ST1', header=1).iloc[:, :4]
st4 = pd.read_excel(F, sheet_name='ST4', header=1).iloc[:, :2]
st1.columns = ['Trait', 'Full', 'Category', 'N']
st4.columns = ['Trait', 'h2se']

key = lambda s: re.sub(r'[\s_]', '', str(s)).lower()
st1['k'] = st1.Trait.map(key)
st4['k'] = st4.Trait.map(key)
st1['Trait'] = st1.Trait.str.strip()
d = st1.merge(st4[['k', 'h2se']], on='k', how='inner')
assert len(d) == 95, len(d)
d = d[~d.Trait.isin(['BT', 'Anticcp_Ab'])].copy()
d['h2'] = d.h2se.str.extract(r'^(-?[\d.]+)').astype(float)
d['se'] = d.h2se.str.extract(r'\(([\d.]+)\)').astype(float)
d['N'] = d.N.astype(str).str.extract(r'^(\d+)').astype(int)

order = ['Anthropometric', 'Hematology', 'Metabolism', 'Tumor marker', 'Vital sign',
         'Echocardiography', 'Cardiac marker', 'Kidney', 'Electrolyte', 'Coagulation',
         'Liver', 'Inflammatory', 'Hormone', 'Protein', 'Ophthalmology']
assert set(order) == set(d.Category)
colors = {
    'Anthropometric': '#7B8FD6', 'Hematology': '#B07CC6', 'Metabolism': '#5FAE88',
    'Tumor marker': '#A07A9C', 'Vital sign': '#EE9A76', 'Echocardiography': '#D96C7B',
    'Cardiac marker': '#C4564F', 'Kidney': '#C9B45A', 'Electrolyte': '#4FA7B8',
    'Coagulation': '#8C7BE0', 'Liver': '#9CC66E', 'Inflammatory': '#8F8F8F',
    'Hormone': '#D9A441', 'Protein': '#6FA3E0', 'Ophthalmology': '#3F8F8F'}

def size(n):
    return 360 if n < 1e5 else 720 if n < 2e5 else 1140 if n < 3e5 else 1680

fig, ax = plt.subplots(figsize=(44, 20))
x0, gap = 0.0, 1.2
texts, xs, ys, ticks = [], [], [], []
for c in order:
    g = d[d.Category == c].sort_values('h2', ascending=False)
    n = len(g)
    x = x0 + np.arange(n)
    ax.add_patch(Rectangle((x0 - 0.6, -1), n + 0.2, 2, color=colors[c], alpha=0.12, lw=0, zorder=0))
    ax.scatter(x, g.h2, s=[size(v) for v in g.N], color=colors[c], alpha=0.55,
               edgecolor=colors[c], linewidth=2, zorder=3)
    for xi, (_, r) in zip(x, g.iterrows()):
        texts.append(ax.text(xi, r.h2, r.Trait, fontsize=36, zorder=4))
        xs.append(xi); ys.append(r.h2)
    ticks.append((x0 + (n - 1) / 2, c))
    x0 += n + gap



ax.set_xlim(-1.2, x0 - gap + 0.6)
ax.set_ylim(-0.01, 0.43)
ax.set_xticks([t[0] for t in ticks])
ax.set_xticklabels([t[1] for t in ticks], rotation=65, ha='right', rotation_mode='anchor', fontsize=51, fontweight='bold')
ax.set_ylabel('h² SNP', fontsize=54, fontweight='bold')
ax.tick_params(axis='y', labelsize=39, width=2, length=10)
ax.tick_params(axis='x', width=2, length=10)
for t in ax.get_yticklabels():
    t.set_fontweight('bold')
for sp in ax.spines.values():
    sp.set_linewidth(2)

adjust_text(texts, x=xs, y=ys, ax=ax, expand=(1.15, 1.45),
            force_text=(0.4, 0.8), force_static=(0.3, 0.6),
            arrowprops=dict(arrowstyle='-', color='grey', lw=1.5), time_lim=60)

for lab, n in [('<100k', 5e4), ('>100k', 1.5e5), ('>200k', 2.5e5), ('>300k', 3.5e5)]:
    ax.scatter([], [], s=size(n), color='grey', alpha=0.55, label=lab)
leg = ax.legend(title='Sample Size', loc='upper left', ncol=4, frameon=True, framealpha=0.9, edgecolor='none',
          fontsize=36, title_fontsize=39, columnspacing=1.5, handletextpad=0.4, borderpad=0.6)
leg.get_title().set_fontweight('bold')
leg._legend_box.align = 'left'

plt.tight_layout()
fig.savefig(OUT + '.png', dpi=300)
fig.savefig(OUT + '.pdf')
d.drop(columns='k').to_csv(OUT + '_data.csv', index=False)
print(d.groupby('Category').size().reindex(order))
print(d.sort_values('h2').head(4)[['Trait', 'h2', 'se', 'N']])
