from Bio.Blast.Applications import NcbitblastnCommandline
import pandas as pd
import xml.etree.ElementTree as ET
from alignment_test import print_reads, compare_reads

fasta_input_file = '../fasta/Homo_sapiens_L1.L1HS.fa'
query = '../test_data/query.fa'
#
blastx_cline = NcbitblastnCommandline(query=query, subject=fasta_input_file, word_size=2,
                                      outfmt=5)
stdout, stderr = blastx_cline()
#
# print(stdout)
# print(stderr)

with open('../test_data/example.xml', 'w') as fh:
    fh.write(stdout)

root = ET.fromstring(stdout)
i = 0
with open('../test_data/test_results.txt', 'w') as fh:
    for iterations in root.findall('./BlastOutput_iterations/Iteration'):
        # input()
        print('\t{}'.format(iterations.find('Iteration_query-def').text))
        # print('\t{}'.format(iterations.find('Iteration_message').text))
        for hit in iterations.findall('./Iteration_hits/Hit'):
            # fh.write('\t\t{}\n'.format(hit.find('./Hit_hsps/Hsp/Hsp_qseq').text))
            # fh.write('\t\t{}\n'.format(hit.find('./Hit_hsps/Hsp/Hsp_midline').text))
            # fh.write('\t\t{}\n'.format(hit.find('./Hit_hsps/Hsp/Hsp_hseq').text))
            print('\t\tqseq: {}'.format(hit.find('./Hit_hsps/Hsp/Hsp_qseq').text))
            print('\t\t      {}'.format(hit.find('./Hit_hsps/Hsp/Hsp_midline').text))
            print('\t\thseq: {}'.format(hit.find('./Hit_hsps/Hsp/Hsp_hseq').text))
            start = int(hit.find('./Hit_hsps/Hsp/Hsp_hit-from').text)
            end = int(hit.find('./Hit_hsps/Hsp/Hsp_hit-to').text) + 1
            print('\t\t{}'.format(start))
            print('\t\t{}'.format(end))
            # TODO: Make sure offsets are correct
            tmp_string = compare_reads(start, end-start)
            # print(print_reads(start, end-start))
            # fh.write(tmp_string)

        # if i == 1:
            # break
        # i += 1

#


# BUILDS THE fasta file from the tsv file
# csv_file = '/mnt/d/projects/2018_xuya/test_data/ORF2.BRCA_OV.CPTAC.ensembl.stringent.txt'
# df = pd.read_csv(csv_file, sep='\t')
# with open(query, 'w') as fh:
#     for i, v in enumerate(df['peptide'].values):
#         fh.write('>SEQ{:0>2}\n{}\n'.format(i, v))
