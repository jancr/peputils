

def fasta_to_protein_hash(fasta_file):
    proteins = {}
    if isinstance(fasta_file, str):
        fasta_file = open(fasta_file)

    acc = None
    for line in fasta_file:
        if line.startswith(">"):
            acc = line.rstrip("\r\n>").lstrip('>').split(' ')[0]
            proteins[acc] = ""
        else:
            proteins[acc] += line.rstrip()

    if isinstance(fasta_file, str):
        fasta_file.close()
    return proteins
