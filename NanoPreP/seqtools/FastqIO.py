from NanoPreP.seqtools.SeqFastq import SeqFastq
from io import TextIOWrapper
from pathlib import Path
from datetime import datetime
from typing import List
import gzip, random, sys


class FastqIO:
    # Batch read FASTQ records
    def batch_read(file: str, from_nth: int, batch_size: int) -> List[SeqFastq]:
        with FastqIO.openg(file, "r") as handle:
            line_generator = (i for i in handle)
            line_counter = -1
            res = []
            while True:
                new_line = next(line_generator, None)
                if new_line is None:
                    break
                line_counter += 1
                if line_counter // 4 >= from_nth:
                    res.append(SeqFastq.parse(
                        new_line,
                        next(line_generator),
                        next(line_generator),
                        next(line_generator)
                    ))
                    line_counter += 3
                    if len(res) >= batch_size:
                        break
        return res
    
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
    
    # sample FASTQ file
    def sample(file: str, n: int, seed: int = 42) -> list:
        # get number of lines
        with FastqIO.openg(file, "r") as handle:
            for nlines, line in enumerate(handle):
                pass
               
        # sample lines
        num_reads = (nlines + 1) // 4 # integer division to avoid float
        random.seed(seed)
        randidx = random.sample(range(num_reads), min(n, num_reads))
        randidx = set(randidx)
        
        # get FASTQ records
        res = []
        fq_lines = []
        keep_lines = 0
        with FastqIO.openg(file, "r") as handle_in:
            for line_num, line in enumerate(handle_in):
                if line_num % 4 == 0:
                    idx = line_num / 4
                    if idx in randidx:
                        keep_lines = 4
                if keep_lines > 0:
                    fq_lines.append(line)
                    keep_lines -= 1
                    if keep_lines == 0:
                        res.append(SeqFastq.parse(*fq_lines))
                        fq_lines = []
        return res, num_reads
    
    @staticmethod
    def openg(p:str, mode:str):
        p = Path(p) if not isinstance(p, Path) else p
        if p.suffix == ".gz":
            return gzip.open(p, mode + "t")
        else:
            return open(p, mode)
        

class FastqIndexIO:
    def __init__(self, file: str) -> None:
        self.file = file
        self.names, self.offsets = FastqIndexIO.fqidx(file)
        return
        
    @staticmethod
    def openg(p:str, mode:str):
        p = Path(p) if not isinstance(p, Path) else p
        if p.suffix == ".gz":
            return gzip.open(p, mode + "t")
        else:
            return open(p, mode)
        
    # index FASTQ file
    @staticmethod
    def fqidx(file: str) -> dict:
        # index the handle
        ordered_keys = []
        offsets = {}
        with FastqIndexIO.openg(file, "r") as handle:
            while True:
                offset = handle.tell()
                identifier = handle.readline().strip()
                if len(identifier) == 0:
                    break
                ordered_keys.append(identifier[1:])
                offsets[identifier[1:]] = offset
                handle.readline()
                handle.readline()
                handle.readline()
        return ordered_keys, offsets
    
    # get SeqFastq by index
    def get(self, name: str) -> SeqFastq:
        with FastqIndexIO.openg(self.file, "r") as handle:
            handle.seek(self.offsets[name])
            return SeqFastq.parse(
                handle.readline(),
                handle.readline(),
                handle.readline(),
                handle.readline()
            )
    
    # iter
    def __iter__(self):
        for name in self.names:
            yield self.get(name)
        return