input_file  = "outputs/ERBB2_DNA.fasta"
output_file = "outputs/ERBB2_mRNA.fasta"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            outfile.write(line)
        else:
            rna_seq = line.strip().upper().replace("T", "U")
            outfile.write(rna_seq + "\n")

print(f"mRNA FASTA saved as: {output_file}")

# Preview conversion
with open(input_file) as f:
    lines = f.readlines()
dna_preview = [l.strip() for l in lines if not l.startswith(">")][0][:60]
mrna_preview = dna_preview.replace("T", "U")
print("DNA  :", dna_preview)
print("mRNA :", mrna_preview)
