from Bio.PDB import MMCIFParser, PDBIO
import json
import os

# Convert AlphaFold CIF to PDB
cif_file = "fold_2026_05_31_15_02_model_0.cif"
pdb_file = "outputs/ERBB2_model0.pdb"

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("ERBB2", cif_file)

io = PDBIO()
io.set_structure(structure)
io.save(pdb_file)

print(f"Converted {cif_file} to {pdb_file}")

# Read and report confidence scores
conf_file = "fold_2026_05_31_15_02_summary_confidences_0.json"
with open(conf_file) as f:
    conf = json.load(f)

report = f"""AlphaFold2 Structure Prediction Report
=======================================

Gene       : ERBB2
Protein    : HER2 (receptor tyrosine-protein kinase erbB-2)
Tool       : AlphaFold Server (alphafoldserver.com)
Model      : model_0 (best ranked)
Date       : 2026-05-31

Confidence Metrics
------------------
pTM              : {conf.get('ptm', 'N/A')}
Chain pTM        : {conf.get('chain_ptm', ['N/A'])[0]}
Ranking Score    : {conf.get('ranking_score', 'N/A')}
Has Clash        : {conf.get('has_clash', 'N/A')}
Recycles         : {conf.get('num_recycles', 'N/A')}
Fraction Disordered: {conf.get('fraction_disordered', 'N/A')}

Interpretation
--------------
pTM of 0.55 indicates acceptable overall fold confidence for a large
1255 aa protein spanning multiple structural domains including the
extracellular domain and the intracellular kinase domain.
"""

with open("outputs/AlphaFold2_report.txt", "w") as f:
    f.write(report)

print(report)
print("Report saved to outputs/AlphaFold2_report.txt")
