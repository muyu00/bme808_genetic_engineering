from collections import Counter
from Bio.Seq import Seq
from Bio.Data import CodonTable

# find the most common amino acid in the translated protein line
# take in the protien and nucleotide sequence (changing U to T since RNA)
def most_common(p_seq, n_seq):
    # codon table from library for RNA (since U has replaced all the T bases)
    table = CodonTable.unambiguous_rna_by_name["Standard"]

    # count the amino acids
    aa_counts = Counter(p_seq) # dictionary of all the amino acids and its occurance
    most_common_aa, aa_count = aa_counts.most_common(1)[0] # () will return the most frequent element and its count

    # Count codons that encode this amino acid
    codon_counter = Counter() # dictonary to store the counting of amino acids

    # loop through the length of the nucleotide sequence, and count the occurances of codons which attribute to protein sequence
    for i in range(0, len(n_seq)-2, 3):
        codon = n_seq[i:i+3]
        # if the amino acid that the codon sequence creates is the most frequent, add to dictionary
        if table.forward_table.get(codon) == most_common_aa:
            codon_counter[codon] += 1
    
    # Print results
    print("\nMost common amino acid:", most_common_aa, "with count:", aa_count)
    print("Codons contributing to", most_common_aa, ":")
    for codon, count in codon_counter.items():
        print(f"{codon}: {count}")
    
# read and count the number of occurances of A, T, C, G
def read_and_count_bases(f):
    # read the file; first line and spaces should be removed
    next(f) # skips the first line with the info ">"

    # print out the file's content
    mRNA = f.read()
    mRNA = mRNA.replace("\n", "")

    lengthA = mRNA.count('A')
    lengthT = mRNA.count('T')
    lengthC = mRNA.count('C')
    lengthG = mRNA.count('G')

    print("Total Length of the Nucleotide: ", len(mRNA))
    print("Number of A: ", lengthA)
    print("Number of T: ", lengthT)
    print("Number of C: ", lengthC)
    print("Number of G: ", lengthG)

    return lengthA, lengthT, len(mRNA) # for question 3

# calculate the percentage of A and T
def proportion(A, T, tot):
    # calculate the percentage of the nucleotide sequence which is A and T
    prop = (A+T)/tot
    prop = prop*100
    
    print("the proprotion of the nucleotide that is A/T is: ", prop)

# count the number of times CAT, CGT, CCT, CTT, and CxT occurs in the nucleotide sequence
# count in frame 1, by threes
def find_CxT(f):
    next(f) # skip the title line ">"

    seq = f.read() # read the remainder of the file
    seq = seq.replace("\n", "") # replace the line breaks to make one long sequence

    # count by threes (frame 1)
    tot_count = 0
    CAT_count = 0
    CCT_count = 0
    CGT_count = 0
    CTT_count = 0
    
    # for x starts at 0, is still in range of the length of the string w/ increments of 3
    for x in range(0, len(seq), 3):
    # if it does, count, if not, continue the for
        if seq[x:x+3] == "CAT":
            tot_count += 1 
            CAT_count += 1
        if seq[x:x+3] == "CCT":
            tot_count += 1
            CCT_count += 1
        if seq[x:x+3] == "CGT":
            tot_count += 1
            CGT_count += 1
        if seq[x:x+3] == "CTT":
            tot_count += 1
            CTT_count += 1
    
    print("Number of CAT in the sequence: ", CAT_count)
    print("Number of CCT in the sequence: ", CCT_count)
    print("Number of CGT in the sequence: ", CGT_count)
    print("Number of CTT in the sequence: ", CTT_count)
    print("Total CxT in the sequence: ", tot_count)

# convert the nucleotide line (replacing T with U) into protien sequence (AA)
def translateSQ(f):
    next(f)  # skip FASTA header

    # read entire nucleotide sequence
    seq = f.read().replace("\n", "")

    # translate T to U
    seq = seq.replace("T", "U")

    # error checker -- ensure length is divisible by 3 (frame 1)
    seq = seq[:len(seq) - (len(seq) % 3)]

    # library function -  translate the protein sequence based on the three pairs
    protein_seq = Seq(seq).translate()

    # display first 3–4 lines (assume standard where its 60 aa per line - not four lines of FASTA files)
    # this information was googled: https://thebumblingbiochemist.com/365-days-of-science/fasta/
    print("First 3–4 lines of translated amino acid sequence:")
    for i in range(0, 240, 60):
        print(protein_seq[i:i+60])

    return protein_seq, seq

# main code to call functions for each question
if __name__ == "__main__":
    # Question 1
    # open file and count the number of bases, total length of the sequence
    # return A, T and total length for question 3
    txt_file = open('Plasmodium Falciparum Sodium Channel.txt', 'r')   # or put in the path to the file you downloaded
    print("Information on Plasmodium Falciparum Sodium Channel")
    lenA_PFSC, lenT_PFSC, tot_PFSC = read_and_count_bases(txt_file)
    print("\n")

    # Question 2
    # open the homo sapein nucleotide, count the total length of the sequence and bases
    # return A, T and total length for question 3
    txt_file2 = open('sequence.fasta','r')
    print("Information on Homo_Sapien nucleotides")
    lenA_HS, lenT_HS, tot_HS = read_and_count_bases(txt_file2)

    # Question 3
    # determine the proportion of the nucleotide that is A/T
    print("Question 3:")
    proportion(lenA_PFSC, lenT_PFSC, tot_PFSC)
    proportion(lenA_HS, lenT_HS, tot_HS)
    print("\n")

    # Question 4
    # count the number of CAT, CCT, CGT, CTT and CxT (where x is any nucleotide)
    # assume reading from frame 1
    print("Question 4:")
    txt_file2 = open('sequence.fasta','r')
    find_CxT(txt_file2)

    # Question 5
    # translation of the Homo-Sapien Sequence
    txt_file2 = open('sequence.fasta','r')
    protein_seq_HS, seq_tot = translateSQ(txt_file2)

    # Question 6
    # determine the most common amino acid, and the nucleotides which build the amino acid
    most_common(protein_seq_HS, seq_tot)
