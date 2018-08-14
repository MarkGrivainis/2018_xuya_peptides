import curses

import pandas as pd
import sys
from Bio import Seq, SeqIO

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
ORF2_START = View.ORF2_START
fasta_input_file = './fasta/Homo_sapiens_L1.L1HS.fa'

l1_seq = None

fasta_sequences = SeqIO.parse(open(fasta_input_file), 'fasta')


for fasta in fasta_sequences:
    l1_seq = Seq.translate(str(fasta.seq[View.ORF2_START:5815]))
    tst_seq = str(fasta.seq[ORF2_START:5815])


def draw_menu(stdscr):
    column = 1
    if len(sys.argv) > 1:
        column = int(sys.argv[1]) if sys.argv[1].isdigit() else column

    bam_input_file = (
        './test_data/3e25dd86-256f-4b4a-bd54-8d8e83d47e37_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
        # './test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_'
        # 'gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
    )
    reader = View(bam_input_file, 0, 200)

    lookup_table = (
        './test_data/tables/'
        'TCGA-BRCA.UUID_Barcode.final_tumor.txt')

    TCGA_PEP_TABLE = (
        './test_data/tables/'
        'all_ORF2_wide.txt'
    )
    # bam_id = bam_input_file.split('/')[-1].split('_')[0]
    #
    # bam_line = filter(
    #     lambda x: x[0].split('_')[0] == bam_id,
    #     map(
    #         lambda x: x.strip().split('\t'),
    #         open(lookup_table, 'r').readlines()
    #     )
    # )
    # tcga_ids = []
    # for id in bam_line:
    #     tcga_ids.append(id[2])
    # id = '-'.join(tcga_ids[0].split('-')[1:3])
    # tcga_pep_positive = pd.read_table(TCGA_PEP_TABLE)
    # peptides = tcga_pep_positive['breast' + '.' + id].dropna().index
    # pep_dict = pblast.build_peptide_dict(peptides)
    # blast_xml_tree = pblast.blast_fasta('\n'.join(['>{}\n{}'.format(k, v) for k,v in pep_dict.items()]), fasta_input_file)
    # blast_results = pblast.process_blast_output(blast_xml_tree)
    # blast_full_length = {k: [z for z in v if len(z[1]['query']) == len(pep_dict[k])] for k, v in blast_results.items()}

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
        './test_data/tables/'
        'TCGA-BRCA.UUID_Barcode.final_tumor.txt')

    TCGA_PEP_TABLE = (
        './test_data/tables/'
        'all_ORF2_wide.txt'
    )
    tcga_pep_positive = pd.read_table(TCGA_PEP_TABLE)
    samples = list(tcga_pep_positive)
    print(samples)
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
    print(matches)
    bam_files = []
    for id in matches:
        bam_files.append(bam_directories[sample_type] + id[0])
    print(bam_files)
    blast_xml_tree = pblast.blast_fasta('\n'.join(['>{}\n{}'.format(k, v) for k,v in pep_dict.items()]))
    blast_results = pblast.process_blast_output(blast_xml_tree)
    print(blast_results)
    blast_full_length = {k: [z for z in v if len(z[1]['query']) == len(pep_dict[k])] for k, v in blast_results.items()}
    print(blast_full_length)

    k = 0
    offset = 0

    height, width = stdscr.getmaxyx()

    reader.set_region('L1HS', 0, width-22)
    AA_from_RNA = reader.update()
    AA_from_RNA_seq = reader.change_read(0)

    # Clear and refresh the screen for a blank canvas
    stdscr.clear()
    stdscr.refresh()

    # Start colors in curses
    curses.start_color()
    curses.init_pair(1, curses.COLOR_CYAN, curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_BLACK, curses.COLOR_WHITE)
    curses.init_pair(4, curses.COLOR_BLACK, curses.COLOR_BLACK)
    curses.init_pair(5, curses.COLOR_BLACK, curses.COLOR_WHITE)
    curses.init_pair(6, curses.COLOR_MAGENTA, curses.COLOR_BLACK)
    curses.init_pair(7, curses.COLOR_MAGENTA, curses.COLOR_WHITE)
    curses.init_pair(8, curses.COLOR_GREEN, curses.COLOR_BLACK)
    curses.init_pair(9, curses.COLOR_GREEN, curses.COLOR_WHITE)
    curses.init_pair(10, curses.COLOR_BLUE, curses.COLOR_BLACK)
    curses.init_pair(11, curses.COLOR_BLUE, curses.COLOR_WHITE)
    print(curses.COLORS)

    # Loop where k is the last character pressed
    while (k != ord('q')):
        # Initialization
        stdscr.clear()

        height, width = stdscr.getmaxyx()

        if k == curses.KEY_RIGHT:
            offset += 50
            offset = max(0, offset)
            offset = min(offset, len(l1_seq) - width + 22)
            reader.set_offset(offset)
            AA_from_RNA = reader.update()
            AA_from_RNA_seq = reader.change_read(0)

        elif k == curses.KEY_LEFT:
            offset -= 50
            offset = max(0, offset)
            offset = min(offset, len(l1_seq) - width + 22)
            reader.set_offset(offset)
            AA_from_RNA = reader.update()
            AA_from_RNA_seq = reader.change_read(0)

        elif k == curses.KEY_DOWN:
            AA_from_RNA_seq = reader.change_read(-1)

        elif k == curses.KEY_UP:
            AA_from_RNA_seq = reader.change_read(1)

        # AA_from_RNA, AA_from_RNA_seq = update_reads(offset, width)

        # stdscr.addstr(50, 50,
        #         str((reads_list[read_idx].reference_start -
        #             reads_list[read_idx].query_alignment_start -
        #             (ORF2_START + offset*3) - translate_offset)/3 + 1))
        # stdscr.addstr(51, 50,
        #         str(reads_list[read_idx].reference_start -
        #             reads_list[read_idx].query_alignment_start))
        # stdscr.addstr(52, 50, str(ORF2_START + offset*3))
        # stdscr.addstr(53, 50, str(translate_offset))
        # stdscr.addstr(54, 50, str(reads_list[read_idx].reference_length))
        # stdscr.addstr(55, 50, str(read_idx))
        # stdscr.addstr(56, 50, str(len(reads_list)))
        # stdscr.addstr(56, 50, reads_list[read_idx].seq)
        # stdscr.addstr(57, 50, reads_list[read_idx].cigarstring)

        statusbarstr = "Press 'q' to exit | STATUS BAR | Pos: {}, {}".format(
                0, 0)

        # Rendering the header
        whstr = "Width: {}, Height: {}".format(width, height)
        stdscr.addstr(0, 0, whstr, curses.color_pair(1))
        offsetstr = "Offset: {}".format(offset)
        stdscr.addstr(1, 0, offsetstr, curses.color_pair(1))

        # Render status bar
        stdscr.attron(curses.color_pair(3))
        stdscr.addstr(height-1, 0, statusbarstr)
        stdscr.addstr(height-1,
                      len(statusbarstr), " " * (width - len(statusbarstr) - 1))
        stdscr.attroff(curses.color_pair(3))

        # Turning on attributes for title
        stdscr.attron(curses.color_pair(2))
        stdscr.attron(curses.A_BOLD)

        # Rendering title
        for idx, AA in enumerate(AA_dict.values()):
            stdscr.addch(idx+7, 1, AA[1])

        box1 = stdscr.subwin(len(AA_lookup) + 4, int(width-20), 4, 5)
        box1.box()
        for idx in range(1, width-21):
            box1.addch(1, idx, l1_seq[idx-1+offset])

        box1.addstr(2, 1, ''.join(AA_from_RNA_seq))

        for idx, AA in enumerate(l1_seq[offset: offset + width - 22]):
            idx_color = {}
            idx_color[AA] = idx_color.get(AA, 4) + 1
            # TODO make sure that there are no duplicates at each position
            for match_list in blast_full_length.values():
                for pep1 in match_list:
                    if (pep1[0]-ORF2_START)//3 <= \
                        idx+offset <= \
                            (pep1[0]-ORF2_START)//3 + len(pep1[1]['query']) - 1:
                        pep = pep1[1]['query'][idx+offset-(pep1[0]-ORF2_START)//3]
                        idx_color[pep] = idx_color.get(pep, 4) + 4
            for con_AA in AA_from_RNA[idx]:
                idx_color[con_AA] = idx_color.get(con_AA, 4) + 2
            for k, v in idx_color.items():
                box1.attron(curses.color_pair(v))
                if k != 'X':
                    box1.addch(AA_lookup[k]+2, idx+1, '■')
                box1.attroff(curses.color_pair(v))

        # consensus
        stdscr.attron(curses.color_pair(1))
        stdscr.addstr(len(AA_lookup) + 8, 8, 'Consensus:')
        stdscr.attroff(curses.color_pair(1))
        stdscr.attron(curses.color_pair(5))
        stdscr.addch(len(AA_lookup) + 8, 19, '■')
        stdscr.attroff(curses.color_pair(5))
        # RNA
        stdscr.attron(curses.color_pair(1))
        stdscr.addstr(len(AA_lookup) + 8, 21, '| RNA:')
        stdscr.attroff(curses.color_pair(1))
        stdscr.attron(curses.color_pair(6))
        stdscr.addch(len(AA_lookup) + 8, 28, '■')
        stdscr.attroff(curses.color_pair(6))
        # RNA + Consensus

        # stdscr.attron(curses.color_pair(7))
        # stdscr.addch(len(AA_lookup) + 8, 53, '■')
        # stdscr.attroff(curses.color_pair(7))

        # # RNA + Peptide + Consensus
        # stdscr.attron(curses.color_pair(9))
        # stdscr.addch(len(AA_lookup) + 8, 28, '■')
        # stdscr.attroff(curses.color_pair(9))
        # # Peptide
        stdscr.attron(curses.color_pair(1))
        stdscr.addstr(len(AA_lookup) + 8, 30, '| Peptide:')
        stdscr.attroff(curses.color_pair(1))
        stdscr.attron(curses.color_pair(8))
        stdscr.addch(len(AA_lookup) + 8, 41, '■')
        stdscr.attroff(curses.color_pair(8))
        # # Peptide + Consensus
        # stdscr.attron(curses.color_pair(11))
        # stdscr.addch(len(AA_lookup) + 8, 32, '■')
        # stdscr.attroff(curses.color_pair(11))
        # # RNA + Peptide
        stdscr.attron(curses.color_pair(1))
        stdscr.addstr(len(AA_lookup) + 8, 43, '| RNA + Peptide:')
        stdscr.attroff(curses.color_pair(1))
        stdscr.attron(curses.color_pair(10))
        stdscr.addch(len(AA_lookup) + 8, 60, '■')
        stdscr.attroff(curses.color_pair(10))

        # draw_sequences
        start_y = len(AA_lookup) + 10
        for match_list in blast_full_length.values():
            for pep1 in match_list:
                stdscr.addstr(start_y, 20, pep1[1]['query'])
                start_y += 1
                stdscr.addstr(start_y, 20, pep1[1]['mid'])
                start_y += 1
                stdscr.addstr(start_y, 20, pep1[1]['hit'])
                start_y += 2
        # Turning off attributes for title
        stdscr.attroff(curses.color_pair(2))
        stdscr.attroff(curses.A_BOLD)

        #  Refresh the screen
        stdscr.refresh()

        # Wait for next input
        k = stdscr.getch()


def main():

    curses.wrapper(draw_menu)


if __name__ == "__main__":
    main()
