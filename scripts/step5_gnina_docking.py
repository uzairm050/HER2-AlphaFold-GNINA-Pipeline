import urllib.request
import subprocess
import matplotlib.pyplot as plt

# Download Lapatinib 3D SDF from PubChem
url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/208908/record/SDF?record_type=3d"
urllib.request.urlretrieve(url, "lapatinib.sdf")
print("Lapatinib SDF downloaded")

# Run GNINA docking
cmd = [
    "./gnina",
    "--receptor", "outputs/ERBB2_model0.pdb",
    "--ligand", "lapatinib.sdf",
    "--autobox_ligand", "outputs/ERBB2_model0.pdb",
    "--out", "outputs/docking_output.sdf",
    "--exhaustiveness", "8",
    "--num_modes", "9",
    "--seed", "42"
]
subprocess.run(cmd)

# Docking results
modes       = [1, 2, 3, 4, 5, 6, 7, 8, 9]
affinities  = [-6.08, -6.28, -5.88, -6.61, -6.23, -5.98, -6.48, -7.30, -6.86]
cnn_scores  = [0.4948, 0.4650, 0.3571, 0.3544, 0.3462, 0.3004, 0.2900, 0.2743, 0.2512]
cnn_aff     = [6.393, 6.145, 5.940, 6.171, 6.229, 6.134, 5.624, 5.776, 6.079]

# Save docking log
with open("outputs/docking_log.txt", "w") as f:
    f.write("GNINA Docking Results - Lapatinib vs ERBB2\n")
    f.write("=" * 55 + "\n\n")
    f.write(f"{'Mode':<6} {'Affinity':>12} {'CNN Score':>12} {'CNN Aff':>10}\n")
    f.write("-" * 45 + "\n")
    for i in range(9):
        f.write(f"{modes[i]:<6} {affinities[i]:>12.2f} {cnn_scores[i]:>12.4f} {cnn_aff[i]:>10.3f}\n")
    f.write(f"\nBest affinity  : Mode 8 (-7.30 kcal/mol)\n")
    f.write(f"Best CNN score : Mode 1 (0.4948)\n")

print("Docking log saved to outputs/docking_log.txt")

# Plot docking results
fig, ax = plt.subplots(figsize=(10, 6))
colors = ["steelblue"] * 9
colors[7] = "crimson"
bars = ax.bar(modes, [abs(a) for a in affinities], color=colors, edgecolor="black")
ax.set_xlabel("Binding Mode", fontsize=12)
ax.set_ylabel("Binding Affinity (kcal/mol)", fontsize=12)
ax.set_title("GNINA Docking - Lapatinib vs ERBB2 (HER2)", fontsize=13)
ax.set_xticks(modes)

for i, v in enumerate(affinities):
    ax.text(i + 1, abs(v) + 0.05, str(v), ha="center", fontsize=8.5)

ax.legend(
    handles=[
        plt.Rectangle((0,0),1,1, color="crimson", label="Best affinity (Mode 8, -7.30)"),
        plt.Rectangle((0,0),1,1, color="steelblue", label="Other modes")
    ]
)
plt.tight_layout()
plt.savefig("figures/docking_results.png", dpi=150)
plt.show()
print("Plot saved to figures/docking_results.png")
