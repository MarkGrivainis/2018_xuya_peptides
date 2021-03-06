import pandas as pd
import sys
from Bio import Seq, SeqIO
import pprint

from my_packages.bam.view import View
import my_packages.blast.blast_text as pblast


AA_dict = {
    'Alanine': ('Ala', 'A'),
    'Arginine': ('Arg', 'R'),
    'Asparagine': ('Asn', 'N'),
    'Aspartic_acid': ('Asp', 'D'),
    'Cysteine': ('Cys', 'C'),
    'Glutamic_acid': ('Glu', 'E'),
    'Glutamine': ('Gln', 'Q'),
    'Glycine': ('Gly', 'G'),
    'Histidine': ('His', 'H'),
    'Hydroxyproline': ('Hyp', 'O'),
    'Isoleucine': ('Ile', 'I'),
    'Leucine': ('Leu', 'L'),
    'Lysine': ('Lys', 'K'),
    'Methionine': ('Met', 'M'),
    'Phenylalanine': ('Phe', 'F'),
    'Proline': ('Pro', 'P'),
    'Pyroglutamatic': ('Glp', 'U'),
    'Serine': ('Ser', 'S'),
    'Threonine': ('Thr', 'T'),
    'Tryptophan': ('Trp', 'W'),
    'Tyrosine': ('Tyr', 'Y'),
    'Valine': ('Val', 'V'),
    'stop_codon': ('*', '*'),
}

AA_lookup = {x[1]: idx+1 for idx, x in enumerate(AA_dict.values())}
#
# ORF2_START = View.ORF2_START
# fasta_input_file = './fasta/Homo_sapiens_L1.L1HS.fa'
#
# l1_seq = None
#
# fasta_sequences = SeqIO.parse(open(fasta_input_file), 'fasta')
#
#
# for fasta in fasta_sequences:
#     l1_seq = Seq.translate(str(fasta.seq[View.ORF2_START:5815]))
#     tst_seq = str(fasta.seq[ORF2_START:5815])
#
#
#
# column = 1
# if len(sys.argv) > 1:
#     column = int(sys.argv[1]) if sys.argv[1].isdigit() else column
#
# bam_input_file = (
#     './test_data/3e25dd86-256f-4b4a-bd54-8d8e83d47e37_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
#     # './test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_'
#     # 'gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
# )
# reader = View(bam_input_file, 0, 200)
#
# lookup_table = (
#     './test_data/tables/'
#     'TCGA-BRCA.UUID_Barcode.final_tumor.txt')
#
# TCGA_PEP_TABLE = (
#     './test_data/tables/'
#     'all_ORF2_wide.txt'
# )
# # bam_id = bam_input_file.split('/')[-1].split('_')[0]
# #
# # bam_line = filter(
# #     lambda x: x[0].split('_')[0] == bam_id,
# #     map(
# #         lambda x: x.strip().split('\t'),
# #         open(lookup_table, 'r').readlines()
# #     )
# # )
# # tcga_ids = []
# # for id in bam_line:
# #     tcga_ids.append(id[2])
# # id = '-'.join(tcga_ids[0].split('-')[1:3])
# # tcga_pep_positive = pd.read_table(TCGA_PEP_TABLE)
# # peptides = tcga_pep_positive['breast' + '.' + id].dropna().index
# # pep_dict = pblast.build_peptide_dict(peptides)
# # blast_xml_tree = pblast.blast_fasta('\n'.join(['>{}\n{}'.format(k, v) for k,v in pep_dict.items()]), fasta_input_file)
# # blast_results = pblast.process_blast_output(blast_xml_tree)
# # blast_full_length = {k: [z for z in v if len(z[1]['query']) == len(pep_dict[k])] for k, v in blast_results.items()}
#
# bam_directories = {'breast': '/run/media/mark/Scripts/Data/2018_Xuyu_project/BAMS/BRCA_RNAseq_Realign/',
#                    'ovarian': '/run/media/mark/Scripts/Data/2018_Xuyu_project/BAMS/OV_RNAseq_Realign/'}
# # bam_input_file = (
# #         '../test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_'
# #         'gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
# #         )
# # bam_id = bam_input_file.split('/')[-1].split('_')[0]
# # lookup_table = (
# #         '../test_data/tables/'
# #         'TCGA-OV.final_UUID_barcode.cleaned.uniqueBam.txt')
# lookup_table = (
#     './test_data/tables/'
#     'TCGA-BRCA.UUID_Barcode.final_tumor.txt')
#
# TCGA_PEP_TABLE = (
#     './test_data/tables/'
#     'all_ORF2_wide.txt'
# )
# tcga_pep_positive = pd.read_table(TCGA_PEP_TABLE)
# samples = list(tcga_pep_positive)
# print(samples)
# sample_type, selected_sample = samples[column].split('.')
# pep_dict = pblast.build_peptide_dict(list(tcga_pep_positive['.'.join((sample_type, selected_sample))].dropna().index))
# bam_line = filter(
#     lambda x: selected_sample in x[2],
#     map(
#         lambda x: x.strip().split('\t'),
#         open(lookup_table, 'r').readlines()
#     )
# )
# matches = []
# for x in bam_line:
#     x[0] = x[0].split('.')[0] + '.Aligned.sortedByCoord.out.bam'
#     matches.append(x)
# print(matches)
# bam_files = []
# for id in matches:
#     bam_files.append(bam_directories[sample_type] + id[0])
# print(bam_files)
# blast_xml_tree = pblast.blast_fasta('\n'.join(['>{}\n{}'.format(k, v) for k,v in pep_dict.items()]))
# blast_results = pblast.process_blast_output(blast_xml_tree)
# print(blast_results)
# blast_full_length = {k: [z for z in v if len(z[1]['query']) == len(pep_dict[k])] for k, v in blast_results.items()}
# print(blast_full_length)
#
# k = 0
# offset = 0
#
# height, width = stdscr.getmaxyx()
#
# reader.set_region('L1HS', 0, width-22)
# AA_from_RNA = reader.update()
# AA_from_RNA_seq = reader.change_read(0)
def get_paths():
    # bam_directories = {'breast': '/run/media/mark/Scripts/Data/2018_Xuyu_project/BAMS/BRCA_RNAseq_Realign/',
    #                    'ovarian': '/run/media/mark/Scripts/Data/2018_Xuyu_project/BAMS/OV_RNAseq_Realign/'}
    bam_directories = {'breast': '/mnt/d/Data/2018_Xuyu_project/BAMS/BRCA_RNAseq_Realign/',
                       'ovarian': '/mnt/d/Data/2018_Xuyu_project/BAMS/OV_RNAseq_Realign/'}
    # bam_input_file = (
    #     './test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
    # )
    lookup_table = (
        './test_data/tables/'
        'TCGA-BRCA.UUID_Barcode.final_tumor.txt')

    TCGA_PEP_TABLE = (
        './test_data/tables/'
        'all_ORF2_wide.txt'
    )
    return bam_directories, lookup_table, TCGA_PEP_TABLE


def main():
    bam_directories, lookup_table, TCGA_PEP_TABLE = get_paths()
    fasta_input_file = './fasta/Homo_sapiens_L1.L1HS.fa'

    l1_seq = None

    fasta_sequences = SeqIO.parse(open(fasta_input_file), 'fasta')

    for fasta in fasta_sequences:
        l1_seq = Seq.translate(str(fasta.seq[View.ORF2_START:View.ORF2_END]))
        tst_seq = str(fasta.seq[View.ORF2_START:View.ORF2_END])

    column = 0
    if len(sys.argv) > 1:
        assert sys.argv[1].isdigit()
        column = int(sys.argv[1])


    l1 = l1_seq[:]
    statistics = {
        'filename': None,
        'length': len(l1),
        'rna_con_match': 0,
        'rna_con_mismatch': 0,
        'con_no_rna': 0,
        'total_aa': 0,
        'aa_with_rna': 0,
        'aa_without_rna': 0
    }
    df = pd.DataFrame(columns=statistics.keys())
    tcga_pep_positive = pd.read_table(TCGA_PEP_TABLE)
    samples = list(tcga_pep_positive)
    for column, val in enumerate(samples):
        sample_type, selected_sample = samples[column].split('.')
        pep_dict = pblast.build_peptide_dict(list(tcga_pep_positive['.'.join((sample_type, selected_sample))].dropna().index))
        bam_line = filter(
            lambda x: selected_sample in x[2],
            map(
                lambda x: x.strip().split('\t'),
                open(lookup_table, 'r').readlines()
            )
        )
        matches = []
        for x in bam_line:
            x[0] = x[0].split('.')[0] + '.Aligned.sortedByCoord.out.bam'
            matches.append(x)

        for bam_input_file, id, tcga_id in matches:
            reader = View(bam_directories[sample_type] + bam_input_file, 0, len(l1_seq))
            AA_from_RNA = reader.update()



            statistics = {
                'filename': bam_input_file,
                'length': len(l1),
                'rna_con_match': 0,
                'rna_con_mismatch': 0,
                'con_no_rna': 0,
                'total_aa': 0,
                'aa_with_rna': 0,
                'aa_without_rna': 0
            }
            blast_xml_tree = pblast.blast_fasta('\n'.join(['>{}\n{}'.format(k, v) for k,v in pep_dict.items()]))
            blast_results = pblast.process_blast_output(blast_xml_tree)
            blast_full_length = {k: [z for z in v if len(z[1]['query']) == len(pep_dict[k])] for k, v in blast_results.items()}
            for seq, matches in blast_full_length.items():
                for match in matches:
                    translated_start_pos = (match[0] - View.ORF2_START)//3
                    for i, AA in enumerate(AA_from_RNA[translated_start_pos:translated_start_pos + len(match[1]['query'])]):
                        statistics['total_aa'] += 1
                        if match[1]['query'][i] in AA:
                            statistics['aa_with_rna'] += 1
                        else:
                            statistics['aa_without_rna'] += 1
            seq = ''
            for i, AA_set in enumerate(AA_from_RNA):
                if l1[i] in AA_set:
                    statistics['rna_con_match'] += 1
                    AA_set.remove(l1[i])
                    seq += l1[i]
                else:
                    statistics['con_no_rna'] += 1
                    seq += ' '
                if AA_set:
                    statistics['rna_con_mismatch'] += 1
            df = df.append(statistics, ignore_index=True)
    print(df)
    df.to_csv('./test_data/rna_and_peptide_counts.csv', index=False)

if __name__ == '__main__':
    main()
