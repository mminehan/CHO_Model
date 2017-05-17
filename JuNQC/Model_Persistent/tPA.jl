include("Include.jl")

#from http://www.ebi.ac.uk/ena/data/view/AAA60111
mRNA_sequence = "TCGCCCAGCCAGGAAATCCATGCCCGATTCAGAAGAGGAGCCAGATCTTACCAAGTGATC
TGCAGAGATGAAAAAACGCAGATGATATACCAGCAACATCAGTCATGGCTGCGCCCTGTG
CTCAGAAGCAACCGGGTGGAATATTGCTGGTGCAACAGTGGCAGGGCACAGTGCCACTCA
GTGCCTGTCAAAAGTTGCAGCGAGCCAAGGTGTTTCAACGGGGGCACCTGCCAGCAGGCC
CTGTACTTCTCAGATTTCGTGTGCCAGTGCCCCGAAGGATTTGCTGGGAAGTGCTGTGAA
ATAGATACCAGGGCCACGTGCTACGAGGACCAGGGCATCAGCTACAGGGGCACGTGGAGC
ACAGCGGAGAGTGGCGCCGAGTGCACCAACTGGAACAGCAGCGCGTTGGCCCAGAAGCCC
TACAGCGGGCGGAGGCCAGACGCCATCAGGCTGGGCCTGGGGAACCACAACTACTGCAGA
AACCCAGATCGAGACTCAAAGCCCTGGTGCTACGTCTTTAAGGCGGGGAAGTACAGCTCA
GAGTTCTGCAGCACCCCTGCCTGCTCTGAGGGAAACAGTGACTGCTACTTTGGGAATGGG
TCAGCCTACCGTGGCACGCACAGCCTCACCGAGTCGGGTGCCTCCTGCCTCCCGTGGAAT
TCCATGATCCTGATAGGCAAGGTTTACACAGCACAGAACCCCAGTGCCCAGGCACTGGGC
CTGGGCAAACATAATTACTGCCGGAATCCTGATGGGGATGCCAAGCCCTGGTGCCACGTG
CTGAAGAACCGCAGGCTGACGTGGGAGTACTGTGATGTGCCCTCCTGCTCCACCTGCGGC
CTGAGACAGTACAGCCAGCCTCAGTTTCGCATCAAAGGAGGGCTCTTCGCCGACATCGCC
TCCCACCCCTGGCAGGCTGCCATCTTTGCCAAGCACAGGAGGTCGCCCGGAGAGCGGTTC
CTGTGCGGGGGCATACTCATCAGCTCCTGCTGGATTCTCTCTGCCGCCCACTGCTTCCAG
GAGAGGTTTCCGCCCCACCACCTGACGGTGATCTTGGGCAGAACATACCGGGTGGTCCCT
GGCGAGGAGGAGCAGAAATTTGAAGTCGAAAAATACATTGTCCATAAGGAATTCGATGAT
GACACTTACGACAATGACATTGCGCTGCTGCAGCTGAAATCGGATTCGTCCCGCTGTGCC
CAGGAGAGCAGCGTGGTCCGCACTGTGTGCCTTCCCCCGGCGGACCTGCAGCTGCCGGAC
TGGACGGAGTGTGAGCTCTCCGGCTACGGCAAGCATGAGGCCTTGTCTCCTTTCTATTCG
GAGCGGCTGAAGGAGGCTCATGTCAGACTGTACCCATCCAGCCGCTGCACATCACAACAT
TTACTTAACAGAACAGTCACCGACAACATGCTGTGTGCTGGAGACACTCGGAGCGGCGGG
CCCCAGGCAAACTTGCACGACGCCTGCCAGGGCGATTCGGGAGGCCCCCTGGTGTGTCTG
AACGATGGCCGCATGACTTTGGTGGGCATCATCAGCTGGGGCCTGGGCTGTGGACAGAAG
GATGTCCCGGGTGTGTACACCAAGGTTACCAACTACCTAGACTGGATTCGTGACAACATG
CGACCGTGA"

#From http://www.uniprot.org/uniprot/P00750.fasta
protein_sequence = "MDAMKRGLCCVLLLCGAVFVSPSQEIHARFRRGARSYQVICRDEKTQMIYQQHQSWLRPV
LRSNRVEYCWCNSGRAQCHSVPVKSCSEPRCFNGGTCQQALYFSDFVCQCPEGFAGKCCE
IDTRATCYEDQGISYRGTWSTAESGAECTNWNSSALAQKPYSGRRPDAIRLGLGNHNYCR
NPDRDSKPWCYVFKAGKYSSEFCSTPACSEGNSDCYFGNGSAYRGTHSLTESGASCLPWN
SMILIGKVYTAQNPSAQALGLGKHNYCRNPDGDAKPWCHVLKNRRLTWEYCDVPSCSTCG
LRQYSQPQFRIKGGLFADIASHPWQAAIFAKHRRSPGERFLCGGILISSCWILSAAHCFQ
ERFPPHHLTVILGRTYRVVPGEEEQKFEVEKYIVHKEFDDDTYDNDIALLQLKSDSSRCA
QESSVVRTVCLPPADLQLPDWTECELSGYGKHEALSPFYSERLKEAHVRLYPSSRCTSQH
LLNRTVTDNMLCAGDTRSGGPQANLHDACQGDSGGPLVCLNDGRMTLVGIISWGLGCGQK
DVPGVYTKVTNYLDWIRDNMRP"

#initialize nucleotide count
nucleotide_count = zeros(1:4)

#count number of each amino acid
nucleotide_count[1] = length(searchall(mRNA_sequence,"A")) #Adenine
nucleotide_count[2] = length(searchall(mRNA_sequence,"C")) #Cytosine
nucleotide_count[3] = length(searchall(mRNA_sequence,"G")) #Guanine
nucleotide_count[4] = length(searchall(mRNA_sequence,"T")) #Uracil

mRNA_length = sum(nucleotide_count)
#initialize amino acid count array
AA_count = zeros(1:20)

#count number of each amino acid
AA_count[1] = length(searchall(protein_sequence,"G")) #Glycine
AA_count[2] = length(searchall(protein_sequence,"A")) #Alanine
AA_count[3] = length(searchall(protein_sequence,"R")) #Arginine
AA_count[4] = length(searchall(protein_sequence,"N")) #Asparagine
AA_count[5] = length(searchall(protein_sequence,"D")) #Aspartic Acid
AA_count[6] = length(searchall(protein_sequence,"C")) #Cysteine
AA_count[7] = length(searchall(protein_sequence,"E")) #Glutamic Acid
AA_count[8] = length(searchall(protein_sequence,"H")) #Histidine
AA_count[9] = length(searchall(protein_sequence,"I")) #Isoleucine
AA_count[10] = length(searchall(protein_sequence,"L")) #Leucine
AA_count[11] = length(searchall(protein_sequence,"K")) #Lysine
AA_count[12] = length(searchall(protein_sequence,"M")) #Methionine
AA_count[13] = length(searchall(protein_sequence,"F")) #Phenylalanine
AA_count[14] = length(searchall(protein_sequence,"P")) #Proline
AA_count[15] = length(searchall(protein_sequence,"S")) #Serine
AA_count[16] = length(searchall(protein_sequence,"T")) #Threonine
AA_count[17] = length(searchall(protein_sequence,"W")) #Tryptophan
AA_count[18] = length(searchall(protein_sequence,"Y")) #Tyrosine
AA_count[19] = length(searchall(protein_sequence,"V")) #Valine
AA_count[20] = length(searchall(protein_sequence,"Q")) #Glutamine

protein_length = sum(AA_count)

#Get species ids
a_id = "s_0773_c_1"
c_id = "s_0059_c_1"
g_id = "s_0226_c_1"
u_id = "s_0431_c_1"

AA_ind = Array{String}(20)

AA_ind[1] = "s_0316_c_1"
AA_ind[2] = "s_0349_c_1"
AA_ind[3] = "s_0004_c_1"
AA_ind[4] = "s_0612_c_1"
AA_ind[5] = "s_0746_c_1"
AA_ind[6] = "s_0754_c_1"
AA_ind[7] = "s_0679_c_1"
AA_ind[8] = "s_0201_c_1"
AA_ind[9] = "s_0366_c_1"
AA_ind[10] = "s_0487_c_1"
AA_ind[11] = "s_0750_c_1"
AA_ind[12] = "s_0391_c_1"
AA_ind[13] = "s_0013_c_1"
AA_ind[14] = "s_0629_c_1"
AA_ind[15] = "s_0562_c_1"
AA_ind[16] = "s_0381_c_1"
AA_ind[17] = "s_0946_c_1"
AA_ind[18] = "s_0037_c_1"
AA_ind[19] = "s_0066_c_1"
AA_ind[20] = "s_0425_c_1"


# Reactions
TX = string(nucleotide_count[1])*"*"*a_id*"+"*string(nucleotide_count[2])*"*"*c_id*"+"*string(nucleotide_count[3])*"*"*g_id*"+"*string(nucleotide_count[4])*"*"*u_id*", product_mRNA,0,inf"

TL = string(AA_count[1])*"*"*AA_ind[1]*"+"*string(AA_count[2])*"*"*AA_ind[2]*"+"*string(AA_count[3])*"*"*AA_ind[3]*"+"*string(AA_count[4])*"*"*AA_ind[4]*"+"*string(AA_count[5])*"*"*AA_ind[5]*"+"*string(AA_count[6])*"*"*AA_ind[6]*"+"*string(AA_count[7])*"*"*AA_ind[7]*"+"*string(AA_count[8])*"*"*AA_ind[8]*"+"*string(AA_count[9])*"*"*AA_ind[9]*"+"*string(AA_count[10])*"*"*AA_ind[10]*"+"*string(AA_count[11])*"*"*AA_ind[11]*"+"*string(AA_count[12])*"*"*AA_ind[12]*"+"*string(AA_count[13])*"*"*AA_ind[13]*"+"*string(AA_count[14])*"*"*AA_ind[14]*"+"*string(AA_count[15])*"*"*AA_ind[15]*"+"*string(AA_count[16])*"*"*AA_ind[16]*"+"*string(AA_count[17])*"*"*AA_ind[17]*"+"*string(AA_count[18])*"*"*AA_ind[18]*"+"*string(AA_count[19])*"*"*AA_ind[19]*"+"*string(AA_count[20])*"*"*AA_ind[20]*",product,0,inf"
