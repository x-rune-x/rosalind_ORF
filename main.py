# Create a class that takes a DNA sequence and ID. The DnaSeq objects will have all the open reading frames (ORF)
# stored as attributes.
class DnaSeq:
    def __init__(self, fasta_id, seq):
        self.fasta_id = fasta_id
        self.seq = seq
        self.length = len(seq)
        self.compliment = self.get_compliment()
        self.orfs = [
            self.seq,
            self.seq[1:],
            self.seq[2:],
            self.compliment,
            self.compliment[1:],
            self.compliment[2:]
        ]

    # object method to create a complementary strand of DNA.
    def get_compliment(self):
        nuc_pairs = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }

        # DNA strands are antiparallel so 5' to 3' runs in opposite directions for the new strand. This loop iterates
        # through the initial strand from back to front and replaces each base with its complement.
        compliment = ""
        for nucleotide in self.seq:
            compliment += nuc_pairs[nucleotide]

        compliment = compliment[::-1]

        return compliment


# Divide DNA sequences into groups of 3 to form codons.
def find_codons(input_seq):
    codon_list = []
    current_codon = ""
    for pos in range(len(input_seq)):
        current_codon += input_seq[pos]
        if len(current_codon) == 3:
            codon_list.append(current_codon)
            current_codon = ""

    return codon_list


# Read a list of codons and form a protein sequence.
def translate(codon_list):

    dna_codon = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "stop",
        "TAG": "stop",
        "TGT": "C",
        "TGC": "C",
        "TGA": "stop",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G"
    }

    prot_seq = []
    for codon in codon_list:
        prot_seq.append(dna_codon[codon])

    return prot_seq


# Read through a protein sequence to find all substrings that begin with a methionine and end with a stop codon. These
# are protein coding segments.
def protein_strings(prot_seq):
    prot_strings = []
    # Go through each amino acid residue in the sequence and check if it is an M.
    for amino_acid in range(len(prot_seq)):
        # When an M is found, read through the sequence from the position of the M and adding the following residues to
        # a string.
        if prot_seq[amino_acid] == "M":
            prot_string = ''
            # If a stop is encountered, add the string to a list of protein coding strings, exit the loop and resume
            # searching the sequence for methionines.
            for aa in range(amino_acid, len(prot_seq)):
                if prot_seq[aa] != "stop":
                    prot_string += prot_seq[aa]
                else:
                    prot_strings.append(prot_string)
                    break

    # Problem calls for distinct protein strings
    prot_strings = list(dict.fromkeys(prot_strings))
    return prot_strings


# Rosalind supplies the sequences in FASTA format as a text file with the ID on the first line and the DNA sequence on
# subsequent lines. This function takes those attributes and uses them as parameters to create a DnaSeq object.
def get_sequence(input_file):
    file = open(input_file)
    fasta_id = file.readline().strip()
    lines = file.readlines()
    fasta_seq = ''
    for line in lines:
        fasta_seq += line.strip()
    seq_obj = DnaSeq(fasta_id, fasta_seq)
    file.close()

    # Now the DnaSeq object has been created, the function reads through each ORF, finds the codons, translates them to
    # a protein sequence and searches that sequences for protein coding segments.
    viable_strings = []
    for orf in seq_obj.orfs:
        codons = find_codons(orf)
        protein_seq = translate(codons)
        strings = protein_strings(protein_seq)
        for string in strings:
            # Separate layer of duplicate removal.
            if string not in viable_strings:
                viable_strings.append(string)

    # Convert the list of protein coding sequences from list form to a string to suit Rosalind.
    string_output = ''
    for vs in viable_strings:
        string_output += vs + '\n'

    return string_output


solution = get_sequence("rosalind_orf(3).txt")
solution_file = open("solution_file.txt", "w")
solution_file.write(solution)
print(solution)
