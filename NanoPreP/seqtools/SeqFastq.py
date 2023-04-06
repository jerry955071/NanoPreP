"""Object to represent NanoPreP-styled FASTQ records"""
from typing import Tuple
import numpy as np
import re

class SeqAnnot(object):
    """Object to represent features (annotation) on NanoPreP-styled records"""
    def __init__(
        self,
        strand: float = 0,
        orientation: int = 1,
        full_length: int = 0,
        fusion: int = 0,
        ploc5: int = -1,
        ploc3: int = -1,
        poly5: int = -1,
        poly3: int = -1
    ) -> None:
        self.strand = strand  # 0: unknown, > 0: sense, < 0: antisense
        self.orientation = orientation  # 1: generic, -1: filpped
        self.full_length = full_length  # 1: full-length, 0: non-full-length
        self.fusion = fusion  # 1: fusion, 0: non-fusion
        self.ploc5 = ploc5  # 0: trimmed, -1: not found, int > 0: primer location on `seq`
        self.ploc3 = ploc3  # 0: trimmed, -1: not found, int > 0: primer location on `seq`
        self.poly5 = poly5  # 0: trimmed, -1: not found, int > 0: polymer length
        self.poly3 = poly3  # 0: trimmed, -1: not found, int > 0: polymer length
        pass

    def __str__(self) -> str:
        return (
            "strand=%.2f "
            "full_length=%s "
            "fusion=%s "
            "ploc5=%s "
            "ploc3=%s "
            "poly5=%s "
            "poly3=%s"
        ) % (
            self.strand * self.orientation,
            self.full_length,
            self.fusion,
            self.ploc5,
            self.ploc3,
            self.poly5,
            self.poly3
        )

    def __repr__(self) -> str:
        return (
            "SeqAnnot("
            "strand=%r, "
            "orientation=%r, "
            "full_length=%r, "
            "fusion=%r, "
            "ploc5=%r, "
            "ploc3=%r, "
            "poly5=%r, "
            "poly3=%r)"
        ) % (
            self.strand,
            self.orientation,
            self.full_length,
            self.fusion,
            self.ploc5,
            self.ploc3,
            self.poly5,
            self.poly3
        )

    @staticmethod
    def from_id(x: str) -> Tuple[str, "SeqAnnot"]:
        program = re.compile(
            r"(?P<prefix>.*)"
            r"strand=(?P<strand>-?\d+\.\d*) "
            r"full_length=(?P<full_length>[01]) "
            r"fusion=(?P<fusion>[01]) "
            r"ploc5=(?P<ploc5>-?\d+) "
            r"ploc3=(?P<ploc3>-?\d+) "
            r"poly5=(?P<poly5>-?\d+) "
            r"poly3=(?P<poly3>-?\d+)"
            r"(?P<suffix>.*)"
        )

        # match `x` with `program`
        res = program.match(x)
        if res:
            return res.group("prefix") + res.group("suffix"), \
                SeqAnnot(
                    strand=float(res.group("strand")),
                    orientation=1,
                    full_length=int(res.group("full_length")),
                    fusion=int(res.group("fusion")),
                    ploc5=int(res.group("ploc5")),
                    ploc3=int(res.group("ploc3")),
                    poly5=int(res.group("poly5")),
                    poly3=int(res.group("poly3"))
            )
        else:
            return x, SeqAnnot()


class SeqFastq(object):
    """Object to represent NanoPreP-styled FASTQ record"""
    def __init__(
        self,
        id: str = "",
        seq: str = "",
        id2: str = "",
        qual: str = "",
        annot: SeqAnnot = None
    ) -> None:
        self.id = id
        self.seq = seq
        self.id2 = id2
        self.qual = qual
        self.annot = annot if annot else SeqAnnot()
        pass

    def __str__(self) -> str:
        space = " " if self.id else ""
        return "@%s%s\n%s\n+%s\n%s\n" % (
            self.id + space,
            self.annot,
            self.seq,
            self.id2,
            self.qual
        )

    def __repr__(self) -> str:
        return str(self)

    def __len__(self) -> int:
        return len(self.seq)

    @staticmethod
    def reverse_complement_static(sequence:str) -> str:
        base_complment = {
            "A": "T", "T": "A",
            "C": "G", "G": "C",
            "N": "N"
        }
        return "".join(
            [base_complment[i] for i in sequence[::-1]]
        )

    def reverse_complement(self) -> None:
        # seq
        self.seq = SeqFastq.reverse_complement_static(self.seq)

        # qual
        self.qual = self.qual[::-1]

        # annot
        self.annot.orientation = {
            1: -1,
            -1: 1
        }[self.annot.orientation]

        # record new annotation
        new = {}
        if self.annot.ploc5 > 0:
            new["ploc3"] = len(self.seq) - self.annot.ploc5
        else:
            new["ploc3"] = self.annot.ploc5
        if self.annot.ploc3 > 0:
            new["ploc5"] = len(self.seq) - self.annot.ploc3
        else:
            new["ploc5"] = self.annot.ploc3
        if self.annot.poly5 > 0:
            new["poly3"] = self.annot.poly5
        else:
            new["poly3"] = self.annot.poly5
        if self.annot.poly3 > 0:
            new["poly5"] = self.annot.poly3
        else:
            new["poly5"] = self.annot.poly3

        # update new annotation
        self.annot.ploc5 = new["ploc5"]
        self.annot.ploc3 = new["ploc3"]
        self.annot.poly5 = new["poly5"]
        self.annot.poly3 = new["poly3"]
        return

    @staticmethod
    def parse(id, seq, id2, qual) -> object:
        # id
        id = id.lstrip("@").rstrip("\n")
        id, annot = SeqAnnot.from_id(id)
        # seq
        seq = seq.strip()
        # qual
        qual = qual.strip()
        # id2
        id2 = id2.lstrip("+").rstrip("\n")
        return SeqFastq(id, seq, id2, qual, annot)

    @staticmethod
    def meanq(read) -> float:
        if len(read) == 0:
            return 0
        # ascii -> dec
        quals = [ord(i) - 33 for i in read.qual]
        quals = np.fromiter(quals, dtype=int, count=len(quals))
        # q -> p -> mean p
        mean_p = np.mean(np.power(10, quals/-10))
        # mean p -> mean q
        mean_q = -10 * np.log10(mean_p)
        return mean_q
