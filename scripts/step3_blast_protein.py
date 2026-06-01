from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = "your@email.com"

# Fetch protein sequence
handle = Entrez.efetch(db="protein", id="NP_004439", rettype="fasta", retmode="text")
protein_record = SeqIO.read(handle, "fasta")

with open("outputs/ERBB2_protein.fasta", "w") as f:
    SeqIO.write(protein_record, f, "fasta")

print("Protein name:", protein_record.description)
print("Length:", len(protein_record.seq), "amino acids")
print("Saved to outputs/ERBB2_protein.fasta")

# Run BLAST
with open("outputs/ERBB2_DNA.fasta") as f:
    lines = f.readlines()
dna_seq = "".join([l.strip() for l in lines if not l.startswith(">")])

print("\nRunning BLAST (this takes 1-2 minutes)...")
result_handle = NCBIWWW.qblast("blastn", "nr", dna_seq[:1000])
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)

with open("outputs/ERBB2_BLAST_results.txt", "w") as out:
    out.write("BLAST Results - ERBB2 DNA Query\n")
    out.write("=" * 50 + "\n\n")
    for i, alignment in enumerate(blast_record.alignments[:5]):
        hsp = alignment.hsps[0]
        out.write(f"Hit {i+1}: {alignment.title[:80]}\n")
        out.write(f"Score: {hsp.score} | E-value: {hsp.expect} | Identity: {hsp.identities}\n\n")

print("BLAST results saved to outputs/ERBB2_BLAST_results.txt")
