from Bio import SeqIO
from Bio import Seq
import pysam
from colorama import Fore, Style

bam_input_file = '../test_data/3e25dd86-256f-4b4a-bd54-8d8e83d47e37_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'

fasta_input_file = '../fasta/Homo_sapiens_L1.L1HS.fa'

fasta_sequences = SeqIO.parse(open(fasta_input_file),'fasta')

l1_seq = None

for fasta in fasta_sequences:
    l1_seq = str(fasta.seq)

samfile = pysam.AlignmentFile(bam_input_file, "rb")

# ┌┐└┘├┤┬┴┼─│

def print_grid(grid, reference):
    totals = [x['A'] + x['C'] + x['G'] + x['T'] for x in grid]
    grid = [
            {
                'A': x['A']/totals[i],
                'C': x['C']/totals[i],
                'G': x['G']/totals[i],
                'T': x['T']/totals[i]
                } for i, x in enumerate(grid)]
            # print('┌───────────┬' + '─────┬'*len(reference[:-1]) + '─────┐')
    # print('│ Reference │' + '│'.join(['  {}  '.format(x) for x in reference]) + '│')
    # print('├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤')
    # print('│     A     │' + '│'.join(['{:>4.0%} '.format(x['A']) for x in grid]) + '│')
    # print('├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤')
    # print('│     C     │' + '│'.join(['{:>4.0%} '.format(x['C']) for x in grid]) + '│')
    # print('├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤')
    # print('│     G     │' + '│'.join(['{:>4.0%} '.format(x['G']) for x in grid]) + '│')
    # print('├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤')
    # print('│     T     │' + '│'.join(['{:>4.0%} '.format(x['T']) for x in grid]) + '│')
    # print('└───────────┴' + '─────┴'*len(reference[:-1]) + '─────┘')
    tmp_string = ''
    tmp_string += '┌───────────┬' + '─────┬'*len(reference[:-1]) + '─────┐\n'
    tmp_string += '│ Reference │' + '│'.join(['  {}  '.format(x) for x in reference]) + '│\n'
    tmp_string += '├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤\n'
    tmp_string += '│     A     │' + '│'.join(['{:>4.0%} '.format(x['A']) for x in grid]) + '│\n'
    tmp_string += '├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤\n'
    tmp_string += '│     C     │' + '│'.join(['{:>4.0%} '.format(x['C']) for x in grid]) + '│\n'
    tmp_string += '├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤\n'
    tmp_string += '│     G     │' + '│'.join(['{:>4.0%} '.format(x['G']) for x in grid]) + '│\n'
    tmp_string += '├───────────┼' + '─────┼'*len(reference[:-1]) + '─────┤\n'
    tmp_string += '│     T     │' + '│'.join(['{:>4.0%} '.format(x['T']) for x in grid]) + '│\n'
    tmp_string += '└───────────┴' + '─────┴'*len(reference[:-1]) + '─────┘\n'
    return tmp_string

def compare_reads(start_pos, seqlength):
    start_pos -= 1
    bam_reads = []
    samfile.reset()
    grid = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(seqlength)]
    print('grid len {}'.format(len(grid)))
    print('seq  len {}'.format(seqlength))
    matched_reads = 0
    other_reads = 0
    counter = 0
    mismatch_seq = []
    for read in samfile.fetch():

        start_offset = start_pos - 1
        ref_s = read.reference_start - read.query_alignment_start
        ref_e = ref_s + read.query_length

        if start_offset - 50 < ref_s <= start_offset + seqlength:
            counter += 1
            read_str = ''
            cigar = read.cigartuples
            start = read.reference_start - start_offset
            seq = read.query_sequence
            if cigar[0][0] == 4:
                start -= cigar[0][1]
            index = 0
            len_so_far = 0

            ref_seq = l1_seq[start_pos:start_pos + seqlength].upper()
            r_match = ref_seq[:max(start_pos, ref_s) - start_pos] + seq[max(start_pos, ref_s) - ref_s:start_pos + seqlength - ref_s] + ref_seq[len(ref_seq) - (start_pos + seqlength - min(start_pos + seqlength, ref_e)):]
            # print(
                    # 'seq :',
                    # '.'*(max(start_pos, ref_s) - start_pos) + \
                            # seq[max(start_pos, ref_s) - ref_s:start_pos + seqlength - ref_s] + \
                            # '.'*(start_pos + seqlength - min(start_pos + seqlength, ref_e))
                            # ,
                    # max(start_pos, ref_s) - start_pos,
                    # start_pos + seqlength - min(start_pos + seqlength, ref_e)
                # )
            if r_match != ref_seq:
                mismatch_seq.append(Seq.translate(r_match))
                other_reads += 1
            else:
                matched_reads += 1
    print(mismatch_seq)
    # for seq in mismatch_seq:
        # print(Seq.translate(seq))
    print('Total: {}; Matched: {}; Mismatched: {}'.format(counter, matched_reads, other_reads))
    print('-'*100)



def print_reads(start_pos, seqlength):
    bam_reads = []
    samfile.reset()
    grid = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(seqlength)]
    print('grid len {}'.format(len(grid)))
    print('seq  len {}'.format(seqlength))
    for read in samfile.fetch():
        start_offset = start_pos - 1

        if start_offset - 50 < read.reference_start <= start_offset + seqlength:
            read_str = ''
            cigar = read.cigartuples
            start = read.reference_start - start_offset
            seq = read.seq
            if cigar[0][0] == 4:
                start -= cigar[0][1]
            index = 0
            len_so_far = 0
            for type, length in cigar:
                if type == 0:
                    for c_idx in range(length):
                        if start_offset <= len_so_far + c_idx + read.reference_start < start_offset + seqlength:
                            try:
                                grid[len_so_far + c_idx + read.reference_start - start_offset][seq[c_idx].upper()] += 1
                            except:
                                print('error {}'.format(len_so_far + c_idx + read.reference_start - start_offset))
                                print('grid_len {}'.format(len(grid)))

                    # grid[index] += 1
                    # index += 1
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
                len_so_far += length

            read_str += Fore.WHITE
            bam_reads.append((start, read_str))
    # return print_grid(grid, l1_seq[start_offset:start_offset + seqlength].upper())
# samfile.close()
    # bam_reads.sort(key=lambda x: x[0])
    print('\n\n\n\n\n')
    for read in bam_reads:
        print('.'*(50+read[0]) + '{}'.format(read[1]))
    print(Style.RESET_ALL + '-'*120)
    print('-'*0 + '{}'.format(Fore.WHITE + l1_seq[start_offset - 50: start_offset] + Fore.GREEN + l1_seq[start_offset:start_offset + seqlength].upper() + Fore.WHITE + l1_seq[start_offset+seqlength:start_offset+seqlength + 20]))


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


if __name__ in "__main__":
    print_reads(2621, 20)
