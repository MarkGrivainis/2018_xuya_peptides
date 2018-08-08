import pysam


class ReadAligner(object):

    """Object which aligns reads using genomic offsets for
    visualization purposes"""

    ORF2_START = 1987
    __slots__ = ('offset', 'width', 'samfile', 'reads', 'ref')

    def __init__(self, filename, offset=0, width=100):
        self.samfile = pysam.AlignmentFile(filename, "rb")
        self.offset = offset
        self.width = width
        self.set_region(
                'L1HS',
                self.ORF2_START + offset,
                (width + offset)*3 + self.ORF2_START
                )

    def set_region(self, ref, offset, width):
        """Function to set which region of the consensus is currently being viewed

        :ref: TODO
        :start: TODO
        :end: TODO
        :returns: TODO

        """
        self.ref = ref
        self.offset = offset
        self.width = width
        self.reads = self.samfile.fetch(
                self.ref,
                self.ORF2_START + self.offset,
                self.ORF2_START + (self.width + self.offset)*3
                )

    def set_offset(self, offset):
        """Function to update the offset and load the relavant reads

        :offset: TODO

        """
        self.offset = offset
        self.reads = self.samfile.fetch(
                self.ref,
                self.ORF2_START + self.offset,
                self.ORF2_START + (self.width + self.offset)*3
                )

    def __iter__(self):
        return self.reads


if __name__ == "__main__":
    bam_input_file = (
            '../test_data/71c5ab4f-ce13-432d-9a90-807ec33cf891'
            '_gdc_realn_rehead.Aligned.sortedByCoord.out.bam'
            )
    reader = ReadAligner(bam_input_file)
    for read in reader:
        print(read)
    reader.set_region('L1HS', 0, 50)
    print('-'*50)
    for read in reader:
        print(read)
