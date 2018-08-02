import sys,os
import curses
from Bio import Seq, SeqIO

from view import View

    # if read.reference_length <= 50:
        # print(read.reference_start - reads_list[read_idx].query_alignment_start, read.query_sequence, read.reference_length, read.cigarstring)
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
    'undetermined': ('X', 'X')
}


AA_lookup = {x[1]: idx+1 for idx, x in enumerate(AA_dict.values())}
ORF2_START = View.ORF2_START
fasta_input_file = '../fasta/Homo_sapiens_L1.L1HS.fa'

l1_seq = None

fasta_sequences = SeqIO.parse(open(fasta_input_file),'fasta')
pep1 = (3392, 'KSPGPDGFIAEFYK')
for fasta in fasta_sequences:
    l1_seq = Seq.translate(str(fasta.seq[View.ORF2_START:5815]))
    tst_seq = str(fasta.seq[ORF2_START:5815])

bam_input_file = '../test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
reader = View(bam_input_file, 0, 200)

def draw_menu(stdscr):
    k = 0
    read_idx = 0
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

        # stdscr.addstr(50, 50, str((reads_list[read_idx].reference_start - reads_list[read_idx].query_alignment_start - (ORF2_START + offset*3) - translate_offset)/3 + 1))
        # stdscr.addstr(51, 50, str(reads_list[read_idx].reference_start - reads_list[read_idx].query_alignment_start))
        # stdscr.addstr(52, 50, str(ORF2_START + offset*3))
        # stdscr.addstr(53, 50, str(translate_offset))
        # stdscr.addstr(54, 50, str(reads_list[read_idx].reference_length))
        # stdscr.addstr(55, 50, str(read_idx))
        # stdscr.addstr(56, 50, str(len(reads_list)))
        # stdscr.addstr(56, 50, reads_list[read_idx].seq)
        # stdscr.addstr(57, 50, reads_list[read_idx].cigarstring)

        statusbarstr = "Press 'q' to exit | STATUS BAR | Pos: {}, {}".format(0, 0)

        # Rendering the header
        whstr = "Width: {}, Height: {}".format(width, height)
        stdscr.addstr(0, 0, whstr, curses.color_pair(1))
        offsetstr = "Offset: {}".format(offset)
        stdscr.addstr(1,0, offsetstr, curses.color_pair(1))

        # Render status bar
        stdscr.attron(curses.color_pair(3))
        stdscr.addstr(height-1, 0, statusbarstr)
        stdscr.addstr(height-1, len(statusbarstr), " " * (width - len(statusbarstr) - 1))
        stdscr.attroff(curses.color_pair(3))

        # Turning on attributes for title
        stdscr.attron(curses.color_pair(2))
        stdscr.attron(curses.A_BOLD)

        # Rendering title
        for idx, AA in enumerate(AA_dict.values()):
            stdscr.addch(idx+7, 1, AA[1])

        box1 = stdscr.subwin(len(AA_lookup) + 4, int(width-20), 4, 5)
        box1.box()
        for idx in range(1,width-21):
            box1.addch(1, idx, l1_seq[idx-1+offset])
            # box1.addch(AA_lookup[l1_seq[idx-1+offset]], idx, l1_seq[idx-1+offset])

        box1.addstr(2, 1, ''.join(AA_from_RNA_seq))

        for idx, AA in enumerate(l1_seq[offset: offset + width - 22]):
            color = 4
            pep = ''
            rna = False
            con = False
            idx_color = {}
            idx_color[AA] = idx_color.get(AA, 4) + 1
            # TODO make sure that there are no duplicates at each position
            if (pep1[0]-ORF2_START)//3 <= idx+offset <= (pep1[0]-ORF2_START)//3 + len(pep1[1]) - 1 :
                pep = pep1[1][idx+offset-(pep1[0]-ORF2_START)//3]
                idx_color[pep] = idx_color.get(pep, 4) + 4
            for con_AA in AA_from_RNA[idx]:
                idx_color[con_AA] = idx_color.get(con_AA, 4) + 2
            for k,v in idx_color.items():
                box1.attron(curses.color_pair(v))
                box1.addch(AA_lookup[k]+2, idx+1, '■')
                box1.attroff(curses.color_pair(v))

        # for idx, AAs in enumerate(AA_from_RNA):
            # plotted_ref = False
            # for AA in AAs:
                # if l1_seq[idx+offset] == AA:
                    # color = 6
                    # plotted_ref = True
                # elif (pep1[0]-ORF2_START)//3 <= idx+offset <= (pep1[0]-ORF2_START)//3 + len(pep1[1]) :
                    # color=7
                # else:
                    # color = 5
                # box1.attron(curses.color_pair(color))
                # box1.addch(AA_lookup[AA]+2, idx+1, '■')
                # box1.attroff(curses.color_pair(color))
            # if not plotted_ref:
                # box1.attron(curses.color_pair(4))
                # box1.addch(AA_lookup[l1_seq[idx+offset]]+2, idx+1, ' ')
                # box1.attroff(curses.color_pair(4))

        # Turning off attributes for title
        stdscr.attroff(curses.color_pair(2))
        stdscr.attroff(curses.A_BOLD)


        # Refresh the screen
        stdscr.refresh()

        # Wait for next input
        k = stdscr.getch()

def main():
    # pass
    curses.wrapper(draw_menu)

if __name__ == "__main__":
    main()


