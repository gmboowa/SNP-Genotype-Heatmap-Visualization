
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# (1) Recreate the genotype matrix from the template:
genotypes = np.array([
    [0,   np.nan, 2, 1, 0, 0, 0, 2, 1, 2],
    [0,   np.nan, 2, 0, 0, 0, 0, 1, 1, 0],
    [1,   0,      0, 0, 1, 2, 0, 1, 1, 0],
    [1,   0,      0, np.nan, np.nan, 2, 0, 0, 1, 1],
    [0,   1,      0, np.nan, 0,      1, 0, 1, 1, 0],
    [np.nan, 2,  np.nan, 2, 1, np.nan, 0, 0, 0, 0]
])
samples = [f"Sample_{i}" for i in range(1, 7)]
snps    = [f"SNP_{i}"    for i in range(1, 11)]

# (2) Build a 3‐color colormap (0→blue, 1→orange, 2→red) and set NA→light gray:
cmap = ListedColormap(['#1f77b4', '#ff7f0e', '#d62728'])
cmap.set_bad(color='lightgray')

# (3) Plot the heatmap:
fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(genotypes, cmap=cmap, vmin=0, vmax=2)

# (4) Annotate each cell with its genotype (skip text for NA):
for i in range(genotypes.shape[0]):
    for j in range(genotypes.shape[1]):
        if not np.isnan(genotypes[i, j]):
            val = int(genotypes[i, j])
            # Use white text on red (val=2), black text otherwise
            txt_color = 'white' if val == 2 else 'black'
            ax.text(j, i, str(val), ha='center', va='center', color=txt_color)

# (5) Axis labels, ticks, and title:
ax.set_xticks(np.arange(len(snps)))
ax.set_xticklabels(snps, rotation=45, ha='right')
ax.set_yticks(np.arange(len(samples)))
ax.set_yticklabels(samples)
ax.set_title("Heatmap of SNP Genotypes (0=Ref, 1=Het, 2=Alt, NA=Gray)")

# (6) Legend showing color→genotype mapping:
legend_elements = [
    Patch(facecolor='#1f77b4', edgecolor='k', label='0 (Ref)'),
    Patch(facecolor='#ff7f0e', edgecolor='k', label='1 (Het)'),
    Patch(facecolor='#d62728', edgecolor='k', label='2 (Alt)'),
    Patch(facecolor='lightgray', edgecolor='k', label='NA')
]
ax.legend(handles=legend_elements, title="Genotype",
          bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
# Save the figure (you can download it from '/mnt/data/categorical_heatmap.png')
plt.savefig('./categorical_heatmap.png')
plt.show()
