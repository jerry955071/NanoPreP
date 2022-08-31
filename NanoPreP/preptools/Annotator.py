from seqtools.SeqFastq import SeqFastq, SeqAnnot
from aligntools.edlibAligner import edlibAligner as aligner
from preptools.polyFinder import polyFinder
import re

class Annotator(object):
    def __init__(
        self,
        p5_sense: str,
        p3_sense: str,
        isl5: tuple,
        isl3: tuple,
        pid_isl: float,
        pid_body: float,
        w: int,
        k: int
    ) -> None:
        ### parse primer sequences
        prog5 = re.compile("(?P<p>[A-Z]+)(\{(?P<max_n>[0-9]*)\}(?P<n>[A-Z]))*")
        prog3 = re.compile("((?P<n>[A-Z])\{(?P<max_n>[0-9]*)\})*(?P<p>[A-Z]+)")
        polymers = {}
        ## p5_sense
        p5_sense, n, max_n = \
            prog5.match(p5_sense).group("p", "n", "max_n")
        if n and max_n:
            polymers[(5, 1)] = (n, int(max_n))

        ## p3_sense
        p3_sense, n, max_n = \
            prog3.match(p3_sense).group("p", "n", "max_n")
        if n and max_n:
            polymers[(3, 1)] = (n, int(max_n))

        ## p5_antisense
        p5_anti = SeqFastq.reverse_complement_static(p3_sense)
        if (3, 1) in polymers.keys():
            n, max_n = polymers[(3, 1)]
            polymers[(5, -1)] = (SeqFastq.reverse_complement_static(n), max_n)

        ## p3_antisense
        p3_anti = SeqFastq.reverse_complement_static(p5_sense)
        if (5, 1) in polymers.keys():
            n, max_n = polymers[(5, 1)]
            polymers[(3, -1)] = (SeqFastq.reverse_complement_static(n), max_n)

        # assign attrs
        self.p5_sense = p5_sense
        self.p3_sense = p3_sense
        self.p5_anti = p5_anti
        self.p3_anti = p3_anti
        self.isl5 = isl5
        self.isl3 = isl3
        self.pid_isl = pid_isl
        self.pid_body = pid_body
        self.w = w  # window size for polymer trimming
        self.k = k  # required number of polymers in the window
        self.poly = polymers

        pass

    def annotate(self, read: SeqFastq):
        # reset annot
        read.annot = SeqAnnot()

        # detect fusion
        name, res = aligner.bestAlign(
            {
                "p5_sense": self.p5_sense,
                "p3_sense": self.p3_sense,
                "p5_anti": self.p5_anti,
                "p3_anti": self.p3_anti
            },
            read.seq[self.isl5[1]:self.isl3[0]],
            mode="HW",
            task="locations",
            pid=self.pid_body
        )
        if res["pid"] >= self.pid_body:
            read.annot.fusion = 1
            return

        # align 5' primers to 5' isl
        strand, res = aligner.bestAlign(
            {1: self.p5_sense, -1: self.p5_anti},
            read.seq[self.isl5[0]:self.isl5[1]],
            mode="HW",
            task="locations",
            pid=self.pid_isl
        )
        if res["pid"] >= self.pid_isl:
            read.annot.ploc5 = res["locations"][-1][-1] + self.isl5[0]
            read.annot.strand += round(strand * res["pid"] * .5, 2)

        # align 3' primers to 3' isl
        strand, res = aligner.bestAlign(
            {1: self.p3_sense, -1: self.p3_anti},
            read.seq[self.isl3[0]:self.isl3[1]],
            mode="HW",
            task="locations",
            pid=self.pid_isl
        )
        if res["pid"] >= self.pid_isl:
            read.annot.ploc3 = res["locations"][0][0] + \
                self.isl3[0] + len(read.seq)
            read.annot.strand += round(strand * res["pid"] * .5, 2)

        # identify full-length reads
        if read.annot.ploc5 > 0 and read.annot.ploc3 > 0:
            read.annot.full_length = 1

        # identify homopolymers next to the primers
        if self.poly:
            for (end, strand), (n, max_n) in self.poly.items():
                # strand == read.annot.strand
                if strand * read.annot.strand > 0:
                    # get finder
                    finder = {
                        5: polyFinder.lfind,
                        3: polyFinder.rfind
                    }[end]
                    finder(read, n, max_n, k=self.k, w=self.w)

        return
