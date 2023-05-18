from NanoPreP.seqtools.SeqFastq import SeqFastq, SeqAnnot
from NanoPreP.aligntools.edlibAligner import edlibAligner as aligner
from NanoPreP.preptools.polyFinder import polyFinder
from collections import namedtuple
from typing import Tuple, Dict
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
        # parser of primer sequences
        prog5 = re.compile("(?P<p>[A-Z]+)((?P<n>[A-Z])\{(?P<max_n>[0-9]*)\})*")
        prog3 = re.compile("((?P<n>[A-Z])\{(?P<max_n>[0-9]*)\})*(?P<p>[A-Z]+)")
        polymers = {}

        # p5_sense -> p5_sense, n, max_n
        p5_sense, n, max_n = \
            prog5.match(p5_sense).group("p", "n", "max_n")
        if n and max_n:
            polymers[(5, 1)] = (n, int(max_n))
        
        # p3_sense -> p3_sense, n, max_n
        p3_sense, n, max_n = \
            prog3.match(p3_sense).group("p", "n", "max_n")
        if n and max_n:
            polymers[(3, 1)] = (n, int(max_n))

        # reverse_complement(p3_sense) -> p5_antisense
        p5_anti = SeqFastq.reverse_complement_static(p3_sense)
        if (3, 1) in polymers.keys():
            n, max_n = polymers[(3, 1)]
            polymers[(5, -1)] = (SeqFastq.reverse_complement_static(n), max_n)

        # reverse_complement(p5_sense) -> p3_antisense
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
        return

    # annotate features on SeqFastq
    def annotate(self, read: SeqFastq, always: bool = False) -> None:
        # reset annot
        read.annot = SeqAnnot()

        # detect fusion
        if len(read) > \
                (self.isl5[1] - self.isl5[0] + self.isl3[1] - self.isl3[0]):
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
            if res["pid"] > 0:
                read.annot.fusion = 1
                if not always:
                    return

        # align 5' primers to 5' isl
        strand, res = aligner.bestAlign(
            {1: self.p5_sense, -1: self.p5_anti},
            read.seq[self.isl5[0]:self.isl5[1]],
            mode="HW",
            task="locations",
            pid=self.pid_isl,
            tie_breaking="right"
        )
        if res["pid"] > 0:
            read.annot.ploc5 = res["location"][-1] + self.isl5[0]
            read.annot.strand += round(strand * res["pid"] * .5, 2)
            strand5 = strand

        # align 3' primers to 3' isl
        strand, res = aligner.bestAlign(
            {1: self.p3_sense, -1: self.p3_anti},
            read.seq[self.isl3[0]:self.isl3[1]],
            mode="HW",
            task="locations",
            pid=self.pid_isl,
            tie_breaking="left"
        )
        if res["pid"] > 0:
            read.annot.ploc3 = res["location"][0] + \
                self.isl3[0] + len(read.seq)
            read.annot.strand += round(strand * res["pid"] * .5, 2)
            strand3 = strand

        # identify full-length reads meeting the following criteria:
        # 1. Both 5' and 3' primers are identified
        # 2. Both 5' and 3' primers are from the same strand
        if read.annot.ploc5 > 0 and \
            read.annot.ploc3 > 0 and \
                strand5 * strand3 > 0:
            read.annot.full_length = 1

        # identify homopolymers with `polyFinder`
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

    
class Optimizer:
    def __init__(
            self,
            p5_sense: str,
            p3_sense: str,
            p5_anti: str = None,
            p3_anti: str = None,
            isl5: tuple = None,
            isl3: tuple = None
        ) -> None:
        self.p5_sense = p5_sense
        self.p3_sense = p3_sense
        self.p5_anti = p5_anti if p5_anti else SeqFastq.reverse_complement_static(p3_sense)
        self.p3_anti = p3_anti if p3_anti else SeqFastq.reverse_complement_static(p5_sense)
        self.isl5 = isl5
        self.isl3 = isl3
        return 
    
    def alignPNbests(self, read: SeqFastq, plen: int, nbests: int) -> Dict[str, namedtuple]:
        """_summary_

        Args:
            read (SeqFastq): _description_
            plen (int): _description_
            nbests (int): _description_

        Returns:
            dict: {"primer_name": SimpleResult(pid, location)}
        """
        SimpleResult = namedtuple("SimpleResult", "pid, location")
        out = {}
        out["p5_sense"] = \
            [SimpleResult(i["pid"], i["location"][1]) for i in 
            aligner.ntopAligns(
                self.p5_sense[-plen:],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        out["p5_anti"] = \
            [SimpleResult(i["pid"], i["location"][1]) for i in 
            aligner.ntopAligns(
                self.p5_anti[-plen:],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        out["p3_sense"] =  \
            [SimpleResult(i["pid"], len(read.seq) - i["location"][0]) for i in 
            aligner.ntopAligns(
                self.p3_sense[:plen],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        out["p3_anti"] = \
            [SimpleResult(i["pid"], len(read.seq) - i["location"][0]) for i in 
            aligner.ntopAligns(
                self.p3_anti[:plen],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        return out
    
    
    def alignP2P(self, plen: int) -> dict:
        """Align primer to primers (5' to 5'; 3' to 3')

        Args:
            plen (int): the length of the primer substring to be aligned

        Returns:
            dict: {"p5": pid, "p3": pid}
        """
        out = {}
        out["p5"] = aligner.singleAlign(
                self.p5_sense[-plen:],
                self.p5_anti[-plen:],
                "HW",
                "locations",
                -1
            )["pid"]
        out["p3"] = aligner.singleAlign(
                self.p3_sense[:plen],
                self.p3_anti[:plen],
                "HW",
                "locations",
                -1
            )["pid"]
    
        return out
    
    
    def alignInB(self, read: SeqFastq) -> Dict[str, Tuple[float, float]]:
        """Align primer to ISL and read body

        Args:
            read (SeqFastq): _description_

        Returns:
            Dict[Tuple[float, float]]: 
            {
                "p5": (pid_isl, pid_body),
                "p3": (pid_isl, pid_body)
            }
        """
        
        out = {}
        # align 5' primers to 5' isl
        strand, res = aligner.bestAlign(
            {1: self.p5_sense, -1: self.p5_anti},
            read.seq[self.isl5[0]:self.isl5[1]],
            mode="HW",
            task="locations",
            pid=-1
        )
        pid_body = aligner.singleAlign(
            {1: self.p5_sense, -1: self.p5_anti}[strand],
            read.seq[self.isl5[1]:self.isl3[0]],
            mode="HW",
            task="locations",
            pid=-1
        )["pid"]
        out["p5"] = (res["pid"], pid_body)

        # align 3' primers to 3' isl
        strand, res = aligner.bestAlign(
            {1: self.p3_sense, -1: self.p3_anti},
            read.seq[self.isl3[0]:self.isl3[1]],
            mode="HW",
            task="locations",
            pid=-1,
        )
        pid_body = aligner.singleAlign(
            {1: self.p3_sense, -1: self.p3_anti}[strand],
            read.seq[self.isl5[1]:self.isl3[0]],
            mode="HW",
            task="locations",
            pid=-1
        )["pid"]
        out["p3"] = (res["pid"], pid_body)
        
        return out