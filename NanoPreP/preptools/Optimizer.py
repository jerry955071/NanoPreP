from NanoPreP.seqtools.SeqFastq import SeqFastq
from NanoPreP.aligntools.edlibAligner import edlibAligner as aligner
from typing import Tuple, Dict, List
from multiprocessing import Pool
import numpy as np
import pandas as pd
import re

class Optimizer:
    def __init__(
            self,
            p5_sense: str,
            p3_sense: str,
            p5_anti: str = None,
            p3_anti: str = None
        ) -> None:
        # parser of primer sequences
        prog5 = re.compile("(?P<p>[A-Z]+)((?P<n>[A-Z])\{(?P<max_n>[0-9]*)\})*")
        prog3 = re.compile("((?P<n>[A-Z])\{(?P<max_n>[0-9]*)\})*(?P<p>[A-Z]+)")
        
        # p5_sense -> p5_sense, n, max_n
        # p3_sense -> p3_sense, n, max_n
        p5_sense, poly_n5, max_n5 = \
            prog5.match(p5_sense).group("p", "n", "max_n")
        p3_sense, poly_n3, max_n3 = \
            prog3.match(p3_sense).group("p", "n", "max_n")
        
        # set values
        self.p5_sense = p5_sense
        self.p3_sense = p3_sense
        self.p5_anti = p5_anti if p5_anti else SeqFastq.reverse_complement_static(p3_sense)
        self.p3_anti = p3_anti if p3_anti else SeqFastq.reverse_complement_static(p5_sense)
        self.poly_n5 = poly_n5
        self.poly_n3 = poly_n3
        self.max_n5 = max_n5 if max_n5 else 0
        self.max_n3 = max_n3 if max_n3 else 0
        return
    
    def optimize(
            self,
            fq_iter:iter,
            plens: List[float],
            n_iqr:List[int],
            processes: int,
            target: str,
            beta: float
        ) -> Dict[str, object]:
        out = {
            "left": self.optimizeLIP(fq_iter, plens, "left", n_iqr, processes, target, beta),
            "right": self.optimizeLIP(fq_iter, plens, "right", n_iqr, processes, target, beta)
        }
        if self.poly_n5:
            out["left"]["seq"] = \
                self.p5_sense[-int(round(len(self.p5_sense) * out["left"]["plen"])):] \
                + self.poly_n5 \
                + "{" + str(self.max_n5) + "}"
        else:
            out["left"]["seq"] = \
                self.p5_sense[-int(round(len(self.p5_sense) * out["left"]["plen"])):]
                
        if self.poly_n3:
            out["right"]["seq"] = \
                self.poly_n3 \
                + "{" + str(self.max_n3) + "}" \
                + self.p3_sense[:int(round(len(self.p3_sense) * out["right"]["plen"]))]
        else:
            out["right"]["seq"] = \
                self.p3_sense[:int(round(len(self.p3_sense) * out["right"]["plen"]))]
            
        return out
    
    def optimizeLIP(
            self,
            fq_iter:iter,
            plens: List[float],
            site: str,
            n_iqr: List[int],
            processes: int,
            target: str,
            beta: float
        ) -> Dict[str, object]:
        """Find the optimal primer length (int), pid cutoff (float), and ideal searching location Tuple[int, int]

        Args:
            fq_iter (iter): an iterator to fastq file

        Returns:
            Dict[str, object]:  optimzed parameters
        """
        out = {target: -1}
        # iter over primer lengths (ascending)
        plens.sort()
        task = []
        for plen in plens:
            task.append((fq_iter, plen, site, n_iqr, target, beta))
        
        with Pool(processes) as pool:
            results = pool.starmap(self.optimizePI, task)
        
        for res in results:
            if res[target] > out[target]:
                out = res
        
        return out
    
    
    def optimizePI(
            self,
            fq_iter: iter,
            plen: float, 
            site: str,
            n_iqr: List[int],
            target: str,
            beta: float
        ) -> Dict[str, object]:
        """Find the optimal pid cutoff (float) and ideal searching location (int) at certain primer length

        Args:
            fq_iter (iter): an iterator to fastq file
            plen (float): length of the primer to use
            site (str): "left" or "right"
            
        Returns:
            Dict[str, object]: optimzed parameters
        """
        out = {target: -1}
        # get pid, location, and class dataframe
        data = self.getPLC(fq_iter, plen, site)
        
        # get pid targets
        pids = data["pid"].unique()
        pids.sort()
        pids = pids[::-1]
        
        # get loc targets
        med_loc = np.median(data[data["cls"] == "positive"]["loc"])
        iqr_loc = np.quantile(
            data[data["cls"] == "positive"]["loc"], .75, method="median_unbiased"
            ) - np.quantile(
            data[data["cls"] == "positive"]["loc"], .25, method="median_unbiased"
        )
        isls = [med_loc + n * iqr_loc for n in n_iqr]
        
        # iter over pid and loc targets
        for pid in pids:
            for location in isls:
                res = self.calculateAP(data, pid, location, beta)
                if res[target] > out[target]:
                    out = res
        out["plen"] = plen
        return out
    
    
    @staticmethod
    def calculateAP(
            data: pd.DataFrame,
            pid_cutoff: int,
            isl_cutoff: int,
            beta: float
        ) -> Dict[str, int or float]:
        tp = sum((data["cls"] == "positive") & (data["pid"] >= pid_cutoff) & (data["loc"] <= isl_cutoff))
        fp = sum((data["cls"] == "negative") & (data["pid"] >= pid_cutoff) & (data["loc"] <= isl_cutoff))
        tn = sum(data["cls"] == "negative") - fp
        fn = sum(data["cls"] == "positive") - tp
        accr = (tp + tn) / (tp + tn + fp + fn)
        prec = tp / (tp + fp) if (tp + fp) > 0 else -1
        recall = tp / (tp + fn) if (tp + fn) > 0 else -1
        fscore = (1 + beta ** 2) * prec * recall / (((beta ** 2) * prec) + recall) if (prec > 0) and (recall > 0) else -1
        return {
            "pid": pid_cutoff,
            "loc": isl_cutoff,
            "accr": accr,
            "prec": prec,
            "recall": recall,
            "fscore": fscore
        }
    
    # get PID, Location, and Class dataframe of 5'/3' primers
    def getPLC(
            self,
            fq_iter: iter,
            plen: float,
            site: str
        ) -> pd.DataFrame:
        # configure nbest_aligner
        nbest_aligner = {
            "left": self.nbestsl,
            "right": self.nbestsr
        }[site]
        
        # initialize data
        data = {
            "pid": [],
            "loc": [],
            "cls": []
        }
        
        # iter over reads in fq_iter
        for read in fq_iter:
            # get [pid, location, class] of top 2 best aligned primers
            res = nbest_aligner(read, plen, 2)
            data["pid"] += [i[0] for i in res]
            data["loc"] += [i[1] for i in res]
            data["cls"] += ["positive", "negative"]
             
        # return a pd.dataframe 
        return pd.DataFrame(data)
    
    
    # n best 5' primers
    def nbestsl(
            self,
            read: SeqFastq,
            plen: float,
            nbests: int
        ) -> List[Tuple[int, int]]:
        """Get the best n alignments of the primers to a read

        Args:
            read (SeqFastq): the read to be aligned
            plen (float): length of the primer to use
            nbests (int): number of best alignments to return

        Returns:
            List[Tuple[int, int]]
        """

        l = int(round(len(self.p5_sense) * plen))
        sense = \
            [(i["pid"], i["location"][1]) for i in 
            aligner.ntopAligns(
                self.p5_sense[-l:],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        l = int(round(len(self.p5_anti) * plen))
        anti = \
            [(i["pid"], i["location"][1]) for i in 
            aligner.ntopAligns(
                self.p5_anti[-l:],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        if sense[0][0] > anti[0][0]:
            return sense
        else:
            return anti
        
        
    # n best 3' primers
    def nbestsr(
            self,
            read: SeqFastq,
            plen: int,
            nbests: int
        ) -> List[Tuple[int, int]]:
        """Get the best n alignments of 3' primers to a read

        Args:
            read (SeqFastq): the read to be aligned
            plen (int): length of the primer to use
            nbests (int): number of best alignments to return

        Returns:
            List[SimpleResult(pid, location)]
        """
        l = int(round(len(self.p3_sense) * plen))
        sense =  \
            [(i["pid"], len(read.seq) - i["location"][0]) for i in 
            aligner.ntopAligns(
                self.p3_sense[:l],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        anti = \
            [(i["pid"], len(read.seq) - i["location"][0]) for i in 
            aligner.ntopAligns(
                self.p3_anti[:l],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        if sense[0][0] > anti[0][0]:
            return sense
        else:
            return anti
    