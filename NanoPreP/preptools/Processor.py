from NanoPreP.seqtools.SeqFastq import SeqFastq

class Processor:
    def trimmer(
        read: SeqFastq,
        trim_poly: bool,
        trim_adapt: bool
    ) -> None:
        a = 0
        b = None
        # get 5' cutpoint
        if read.annot.poly5 > 0 and trim_poly:
            a = read.annot.ploc5 + read.annot.poly5
            read.annot.poly5 = read.annot.ploc5 = 0 # 0: trimmed
        elif read.annot.ploc5 > 0 and trim_adapt:
            a = read.annot.ploc5
            read.annot.ploc5 = 0  # 0: trimmed
        else:
            a = 0

        # get 3' cutpoint
        if read.annot.poly3 > 0 and trim_poly:
            b = read.annot.ploc3 - read.annot.poly3
            read.annot.poly3 = read.annot.ploc3 = 0  # 0: trimmed
        elif read.annot.ploc3 > 0 and trim_adapt:
            b = read.annot.ploc3
            read.annot.ploc3 = 0  # 0: trimmed
        else:
            b = None

        # trim read
        read.seq = read.seq[a:b]
        read.qual = read.qual[a:b]

        return

    def orientor(read: SeqFastq, to: int) -> None:
        # strand unknown
        if read.annot.strand == 0:
            return

        # re-orient
        if read.annot.strand * read.annot.orientation * to < 0:
            read.reverse_complement()

        return
