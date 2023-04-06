from NanoPreP.seqtools.SeqFastq import SeqFastq
from io import TextIOWrapper

class FastqIO:
    # FASTQ generator
    def read(handle: TextIOWrapper) -> SeqFastq:
        # reading from the handle
        for line in handle:
            yield SeqFastq.parse(
                line,
                next(handle),
                next(handle),
                next(handle)
            )

    # write SeqFastq to files in FASTQ format
    def write(handle: TextIOWrapper, record: SeqFastq) -> None:
        handle.write(str(record))
        return
