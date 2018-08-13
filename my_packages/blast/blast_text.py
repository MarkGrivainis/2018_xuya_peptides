from Bio.Blast.Applications import NcbitblastnCommandline
from typing import List, Dict
import xml.etree.ElementTree as ET
import pandas as pd


def build_peptide_dict(sequences: List[str]) -> Dict[str, str]:
    """TODO: Docstring for build_fasta.

    :sequences: TODO
    :returns: TODO

    """
    return {'SEQ{:0>2}'.format(i): v for i, v in enumerate(sequences)}


def blast_fasta(query: str, reference: str =None):
    """TODO: Docstring for blast_fasta.

    :fasta: Fasta file containing peptide sequences
    :reference: reference sequence to blast against (subject)
    :returns: XML tree of blast output

    """
    if not reference:
        reference = '../fasta/Homo_sapiens_L1.L1HS.fa'
    #
    blastx_cline = NcbitblastnCommandline(
            # query=query,
            subject=reference,
            word_size=2,
            outfmt=5
            )
    # stdout, stderr = blastx_cline()
    stdout, stderr = blastx_cline(stdin=query)
    return ET.fromstring(stdout)


def process_blast_output(xml_tree):
    """TODO: Docstring for process_blast_output.

    :xml_tree: TODO
    :returns: TODO

    """
    results = {}
    for iterations in xml_tree.findall('./BlastOutput_iterations/Iteration'):
        seq = iterations.find('Iteration_query-def').text
        if seq not in results:
            results[seq] = []
        for hit in iterations.findall('./Iteration_hits/Hit'):
            for h in hit.findall('./Hit_hsps/Hsp'):
                start = int(h.find('./Hsp_hit-from').text)
                qseq = h.find('./Hsp_qseq').text
                mlin = h.find('./Hsp_midline').text
                hseq = h.find('./Hsp_hseq').text
                results[seq].append(
                    (start,
                     {
                         'query': qseq,
                         'mid': mlin,
                         'hit': hseq
                     }
                     ))
    return results


def main():
    bam_directories = {'breast': '/run/media/mark/Scripts/Data/2018_Xuyu_project/BAMS/BRCA_RNAseq_Realign/',
                       'ovarian': '/run/media/mark/Scripts/Data/2018_Xuyu_project/BAMS/OV_RNAseq_Realign/'}
    # bam_input_file = (
    #         '../test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_'
    #         'gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
    #         )
    # bam_id = bam_input_file.split('/')[-1].split('_')[0]
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

    tcga_pep_positive = pd.read_table(TCGA_PEP_TABLE)
    samples = list(tcga_pep_positive)
    print(samples)
    sample_type, selected_sample = samples[1].split('.')
    pep_dict = build_peptide_dict(list(tcga_pep_positive['.'.join((sample_type, selected_sample))].dropna().index))
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
    print(matches)
    bam_files = []
    for id in matches:
        bam_files.append(bam_directories[sample_type] + id[0])
    print(bam_files)
    blast_xml_tree = blast_fasta('\n'.join(['>{}\n{}'.format(k, v) for k,v in pep_dict.items()]))
    blast_results = process_blast_output(blast_xml_tree)
    print(blast_results)
    blast_full_length = {k: [z for z in v if len(z[1]['query']) == len(pep_dict[k])] for k, v in blast_results.items()}
    print(blast_full_length)


if __name__ == "__main__":
    main()

