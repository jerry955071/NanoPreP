from NanoPreP.seqtools.SeqFastq import SeqFastq

class polyFinder:
    # find the length of homopolymers `n` next to ploc5
    def lfind(
        read: SeqFastq,
        n: str,
        max_n: int,
        k: int,
        w: int
    ) -> None:
        # start from `ploc5`
        i = read.annot.ploc5
        # if `ploc5` unavailable
        if i == -1:
            return

        # window-sliding algorithm
        while read.seq[i:i+w].count(n) >= k:
            # exceed `max_n`
            if i + w - read.annot.ploc5 > max_n:
                break
            # slide 1 position
            i += 1

        # report poly
        read.annot.poly5 = i + w - 1 - read.annot.ploc5
        return

    # find the length of homopolymers `n` next to ploc3
    def rfind(
        read: SeqFastq,
        n: str,
        max_n: int,
        k: int,
        w: int
    ) -> None:
        # start from `ploc3`
        i = read.annot.ploc3
        # if `ploc3` unavailable
        if i == -1:
            return

        # window-sliding algorithm
        while read.seq[i-w:i].count(n) >= k:
            # exceed `max_n`
            if read.annot.ploc3 - i + w > max_n:
                break
            # slide 1 position
            i -= 1

        # report poly
        read.annot.poly3 = read.annot.ploc3 - i + w - 1
        return
