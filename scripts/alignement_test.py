from Bio import SeqIO
import pysam
from colorama import Fore, Style

bam_input_file = '/mnt/e/2018_Xuyu_project/BAMS/BRCA_RNAseq_Realign/1f578ca2-3646-44f7-a1d3-76491344ac27_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'

fasta_input_file = '/mnt/e/2018_Xuyu_project/fasta/Homo_sapiens_L1.L1HS.fa'

fasta_sequences = SeqIO.parse(open(fasta_input_file),'fasta')

l1_seq = None

for fasta in fasta_sequences:
    l1_seq = str(fasta.seq)

samfile = pysam.AlignmentFile(bam_input_file, "rb")

bam_reads = []

for read in samfile.fetch():
    start_offset = 100
    if start_offset < read.reference_start <= start_offset + 100:
        read_str = ''
        cigar = read.cigartuples
        start = read.reference_start - start_offset
        seq = read.seq
        if cigar[0][0] == 4:
            start -= cigar[0][1]

        for type, length in cigar:
            if type == 0:
                read_str += Fore.GREEN
            elif type == 1:
                read_str += Fore.BLUE
            elif type == 2:
                read_str += Fore.MAGENTA
                seq = '-'*length + seq
                print(cigar)
            elif type == 3:
                read_str += Fore.YELLOW
            elif type == 4:
                read_str += Fore.RED
            elif type == 5:
                read_str += Fore.WHITE
            elif type == 6:
                read_str += Fore.WHITE
            elif type == 7:
                read_str += Fore.WHITE
            elif type == 8:
                read_str += Fore.WHITE
            read_str += seq[0:length]
            seq = seq[length:]

        read_str += Fore.WHITE
        bam_reads.append((start, read_str))

samfile.close()
bam_reads.sort(key=lambda x: x[0])
print('\n\n\n\n\n')
for read in bam_reads:
    print('.'*(25+read[0]) + '{}'.format(read[1]))
print(Style.RESET_ALL + '-'*120)
print('-'*25 + '{}'.format(l1_seq[start_offset:start_offset + 100]))


# M	BAM_CMATCH	0
# I	BAM_CINS	1
# D	BAM_CDEL	2
# N	BAM_CREF_SKIP	3
# S	BAM_CSOFT_CLIP	4
# H	BAM_CHARD_CLIP	5
# P	BAM_CPAD	6
# =	BAM_CEQUAL	7
# X	BAM_CDIFF	8
# B	BAM_CBACK	9
