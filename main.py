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
            self.compliment[2:],
        ]

    def get_compliment(self):
        nuc_pairs = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }

        compliment = ""
        for nucleotide in self.seq:
            compliment += nuc_pairs[nucleotide]

        return compliment


def find_codons(input_seq):
    codon_list = []
    current_codon = ""
    for pos in range(len(input_seq)):
        current_codon += input_seq[pos]
        if len(current_codon) == 3:
            codon_list.append(current_codon)
            # print(current_codon)
            current_codon = ""

    return codon_list


def translate(self):

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
        "TAA": "Stop",
        "TAG": "Stop",
        "TGT": "C",
        "TGC": "C",
        "TGA": "Stop",
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

    prot_seq = ""
    for codon in codon_list:
        prot_seq += dna_codon[codon]

    return prot_seq

# prot_strings = []
# for amino_acid in range(len(prot_seq)):
#    if prot_seq[amino_acid] == "M":
#        prot_string = prot_seq[amino_acid, prot_seq.find("Stop")]
#        prot_strings.append(prot_string)
#
# return prot_strings


sequence = DnaSeq("1a", "TGACTGACTGAC")
print(sequence.compliment)
print(find_codons(sequence.compliment))
