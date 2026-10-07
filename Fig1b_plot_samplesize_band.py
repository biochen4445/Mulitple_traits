import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator

F = 'Mutli_Traits_TableS_2026_09_28.xlsx'
OUT = 'Fig1b_SampleSize_15categories_band'

import logging; logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
plt.rcParams['font.family'] = ['Arial', 'Liberation Sans', 'DejaVu Sans']

d = pd.read_excel(F, sheet_name='ST1', header=1).iloc[:, :4]
d.columns = ['Trait', 'Full', 'Category', 'N']
d = d.dropna(subset=['Trait'])
d['Trait'] = d.Trait.str.strip()
d['N'] = d.N.astype(str).str.extract(r'^(\d+)')[0].astype(int)
assert len(d) == 95

# Same category colours as the h2 SNP figure (plot_h2.py)
colors = {
    'Anthropometric': '#7B8FD6', 'Hematology': '#B07CC6', 'Metabolism': '#5FAE88',
    'Tumor marker': '#A07A9C', 'Vital sign': '#EE9A76', 'Echocardiography': '#D96C7B',
    'Cardiac marker': '#C4564F', 'Kidney': '#C9B45A', 'Electrolyte': '#4FA7B8',
    'Coagulation': '#8C7BE0', 'Liver': '#9CC66E', 'Inflammatory': '#8F8F8F',
    'Hormone': '#D9A441', 'Protein': '#6FA3E0', 'Ophthalmology': '#3F8F8F'}

# Row order and display labels as in the original sample-size figure (top to bottom)
rows = [('Hematology', 'Hematology'), ('Metabolism', 'Metabolism'),
        ('Echocardiography', 'Echocardiography'), ('Inflammatory', 'Inflammation'),
        ('Kidney', 'Kidney'), ('Liver', 'Liver'), ('Coagulation', 'Coagulation'),
        ('Vital sign', 'Vital signs'), ('Electrolyte', 'Electrolytes'),
        ('Hormone', 'Hormones'), ('Anthropometric', 'Anthropometric'),
        ('Protein', 'Protein'), ('Ophthalmology', 'Ophthalmology'),
        ('Tumor marker', 'Tumor markers'), ('Cardiac marker', 'Cardiac marker')]
assert set(c for c, _ in rows) == set(d.Category)

rng = np.random.default_rng(1)
fig, ax = plt.subplots(figsize=(10, 4.7))
nrow = len(rows)
ypos = {}
for i, (c, lab) in enumerate(rows):
    y = nrow - 1 - i
    g = d[d.Category == c]
    jit = rng.uniform(-0.28, 0.28, len(g)) if len(g) > 1 else np.zeros(1)
    ax.axhspan(y - 0.42, y + 0.42, color=colors[c], alpha=0.12, lw=0, zorder=0)
    ax.scatter(g.N, y + jit, s=55, color=colors[c], edgecolor='white', linewidth=0.8, zorder=3)
    for t, yy in zip(g.Trait, y + jit):
        ypos[t] = yy
    ax.text(1.005, y, str(len(g)), transform=ax.get_yaxis_transform(), va='center', ha='left',
            fontsize=14, color='#555555')
ax.text(1.005, nrow - 0.35, 'Traits', transform=ax.get_yaxis_transform(), va='bottom', ha='left',
        fontsize=14, color='#555555')

ax.set_xscale('log')
ax.set_xlim(1000, 400000)
ax.set_ylim(-0.6, nrow - 0.4)
ax.set_yticks(range(nrow))
ax.set_yticklabels([lab for _, lab in rows][::-1], fontsize=15)
ax.tick_params(axis='y', length=0, pad=6)
ax.xaxis.set_major_locator(FixedLocator([1e3, 1e4, 1e5]))
ax.set_xticklabels(['1,000', '10,000', '100,000'], fontsize=14)
ax.tick_params(axis='x', which='major', length=6)
ax.tick_params(axis='x', which='minor', length=3)
ax.set_xlabel('Individuals with the measurement', fontsize=15)
for v in [1e4, 1e5]:
    ax.axvline(v, color='#D0D0D0', lw=1, zorder=1)
for sp in ['top', 'right', 'left']:
    ax.spines[sp].set_visible(False)
ax.spines['bottom'].set_linewidth(1.2)

def note(trait, text, xytext, color, ha='left'):
    r = d[d.Trait == trait].iloc[0]
    ax.annotate(text, xy=(r.N, ypos[trait]), xytext=xytext, textcoords='data',
                fontsize=13, color=color, ha=ha, va='center',
                arrowprops=dict(arrowstyle='-', color=color, lw=0.9, shrinkA=2, shrinkB=4), zorder=2)

# Annotations removed (2026-10-05)

plt.tight_layout()
fig.savefig(OUT + '.png', dpi=300)
fig.savefig(OUT + '.pdf')
print(d.groupby('Category').size())
