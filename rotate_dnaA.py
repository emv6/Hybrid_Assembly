from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Input files
fasta_file = "23EV599.fasta"
gff_file = "23EV599.gff"
output_file = "23EV599_rotated.fasta"

# Find dnaA in GFF
dnaA_found = False
with open(gff_file) as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        if "dnaa" in line.lower():
            fields = line.strip().split("\t")
            contig = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            dnaA_found = True
            break

if not dnaA_found:
    raise ValueError("dnaA not found in GFF.")

# Determine rotation point
rotation_point = end if strand == "-" else start
print(f"Rotating contig '{contig}' at position {rotation_point} (strand: {strand})")

# Rotate and optionally reverse-complement
rotated_records = []
for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id == contig:
        idx = rotation_point - 1  # 0-based
        seq = record.seq[idx:] + record.seq[:idx]
        if strand == "-":
            seq = seq.reverse_complement()
        rotated = SeqRecord(seq, id=record.id, description="rotated to dnaA")
        rotated_records.append(rotated)
    else:
        rotated_records.append(record)

SeqIO.write(rotated_records, output_file, "fasta")
print(f"Rotated genome saved to: {output_file}")
