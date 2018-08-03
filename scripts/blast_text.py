from Bio.Blast.Applications import NcbitblastnCommandline
import xml.etree.ElementTree as ET
import pandas as pd


def build_fasta(sequences):
    """TODO: Docstring for build_fasta.

    :sequences: TODO
    :returns: TODO

    """
    pass


def blast_fasta(fasta, reference):
    """TODO: Docstring for blast_fasta.

    :fasta: Fasta file containing peptide sequences
    :reference: reference sequence to blast against (subject)
    :returns: XML tree of blast output

    """
    fasta_input_file = '../fasta/Homo_sapiens_L1.L1HS.fa'
    query = '../test_data/query_71c5ab4f.fa'
    #
    blastx_cline = NcbitblastnCommandline(
            query=query,
            subject=fasta_input_file,
            word_size=2,
            outfmt=5
            )
    stdout, stderr = blastx_cline()
    return ET.fromstring(stdout)


def process_blast_output(xml_tree):
    """TODO: Docstring for process_blast_output.

    :xml_tree: TODO
    :returns: TODO

    """
    for iterations in xml_tree.findall('./BlastOutput_iterations/Iteration'):
        print('\t{}'.format(iterations.find('Iteration_query-def').text))
        for hit in iterations.findall('./Iteration_hits/Hit'):
            start = int(hit.find('./Hit_hsps/Hsp/Hsp_hit-from').text)
            # end = int(hit.find('./Hit_hsps/Hsp/Hsp_hit-to').text) + 1
            print(start, hit.find('./Hit_hsps/Hsp/Hsp_qseq').text)


# BUILDS THE fasta file from the tsv file
# csv_file = ('/mnt/d/projects/2018_xuya/test_data/'
# 'ORF2.BRCA_OV.CPTAC.ensembl.stringent.txt')
# df = pd.read_csv(csv_file, sep='\t')
# with open(query, 'w') as fh:
#     for i, v in enumerate(df['peptide'].values):
#         fh.write('>SEQ{:0>2}\n{}\n'.format(i, v))

if __name__ == "__main__":
    bam_input_file = (
            '../test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_'
            'gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
            )
    bam_id = bam_input_file.split('/')[-1].split('_')[0]
    # lookup_table = (
    #         '../test_data/tables/'
    #         'TCGA-OV.final_UUID_barcode.cleaned.uniqueBam.txt')
    lookup_table = (
            '../test_data/tables/'
            'TCGA-BRCA.UUID_Barcode.final_tumor.txt')

    TCGA_PEP_TABLE = (
            '../test_data/tables/'
            'all_ORF2_wide.txt'
            )

    bam_line = filter(
            lambda x: x[0].split('_')[0] == bam_id,
            map(
                lambda x: x.strip().split('\t'),
                open(lookup_table, 'r').readlines()
                )
            )
    tcga_ids = []
    for id in bam_line:
        tcga_ids.append(id[2])
    print(tcga_ids)
    id = '-'.join(tcga_ids[0].split('-')[1:3])
    tcga_pep_positive = pd.read_table(TCGA_PEP_TABLE)
    # print(tcga_pep_positive)
    peptides = tcga_pep_positive['breast' + '.' + id].dropna().index
    print(peptides)
    for pep in peptides:
        print(pep)
    # main()
