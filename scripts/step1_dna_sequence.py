from Bio import Entrez, SeqIO

Entrez.email = "your@email.com"

handle = Entrez.efetch(
    db="nucleotide",
    id="NC_000017.11",
    rettype="fasta",
    retmode="text",
    seq_start=39688094,
    seq_stop=39728658
)

record = SeqIO.read(handle, "fasta")

with open("outputs/ERBB2_DNA.fasta", "w") as f:
    SeqIO.write(record, f, "fasta")

print("ID:", record.id)
print("Length:", len(record.seq), "bp")
print("Saved to outputs/ERBB2_DNA.fasta")
