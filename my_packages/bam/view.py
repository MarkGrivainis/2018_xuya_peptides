from Bio import SeqIO
from Bio import Seq
from .bamreader import ReadAligner

class View(ReadAligner):
    """Update the characters that will be shown on the screen"""


    __slots__ = ('current_read', 'loaded_reads')

    def __init__(self, filename, offset, width):
        """TODO: to be defined. """
        super().__init__(filename, offset, width)
        self.current_read = 0
        self.loaded_reads = []

    def update(self):
        AA_from_RNA = [set() for x in range(self.width)]
        self.loaded_reads = []

        for read in self.reads:
            if read.reference_length <= 50:
                self.loaded_reads.append(read)
                translate_offset = (read.reference_start - read.query_alignment_start - self.ORF2_START)%3
                read_start = (read.reference_start - read.query_alignment_start - (self.ORF2_START + self.offset*3) - translate_offset)//3
                read_translated_seq = Seq.translate('N'*translate_offset + read.seq)
                for idx, char in enumerate(read_translated_seq):
                    if read_start + idx >= 0 and read_start + idx < len(AA_from_RNA):
                        AA_from_RNA[read_start + idx].add(char)

        self.loaded_reads.sort(key=lambda x: x.reference_start - x.query_alignment_start)

        return AA_from_RNA

    def change_read(self, incriment):
        self.current_read += incriment
        self.current_read = max(0, self.current_read)
        self.current_read = min(self.current_read, len(self.loaded_reads) - 1)

        AA_from_RNA_seq = [' ' for x in range(self.width)]

        translate_offset = (self.loaded_reads[self.current_read].reference_start - self.loaded_reads[self.current_read].query_alignment_start - self.ORF2_START)%3
        read_start = (self.loaded_reads[self.current_read].reference_start - self.loaded_reads[self.current_read].query_alignment_start - (self.ORF2_START + self.offset*3) - translate_offset)//3
        read_translated_seq = Seq.translate('N'*translate_offset + self.loaded_reads[self.current_read].seq)
        for idx, char in enumerate(read_translated_seq):
            if read_start + idx >= 0 and read_start + idx < len(AA_from_RNA_seq):
                AA_from_RNA_seq[read_start + idx] = char
        return AA_from_RNA_seq

    def set_offset(self, offset):
        """Function to update the offset and load the relavant reads

        :offset: TODO

        """
        self.offset = offset
        self.reads = self.samfile.fetch(self.ref, self.ORF2_START + self.offset, self.ORF2_START + (self.width + self.offset)*3)
        self.current_read = 0

if __name__ == "__main__":
    bam_input_file = '../test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
    v = View(bam_input_file, 0, 100)
    v.update()
