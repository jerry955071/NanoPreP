from NanoPreP.seqtools.SeqFastq import SeqFastq
from NanoPreP.aligntools.edlibAligner import edlibAligner as aligner
from collections import namedtuple, Counter
from typing import Tuple, Dict, List
from multiprocessing import Pool
import numpy as np
import pandas as pd
import sys
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
        p5_sense, n, max_n = \
            prog5.match(p5_sense).group("p", "n", "max_n")
        
        # p3_sense -> p3_sense, n, max_n
        p3_sense, n, max_n = \
            prog3.match(p3_sense).group("p", "n", "max_n")

        self.p5_sense = p5_sense
        self.p3_sense = p3_sense
        self.p5_anti = p5_anti if p5_anti else SeqFastq.reverse_complement_static(p3_sense)
        self.p3_anti = p3_anti if p3_anti else SeqFastq.reverse_complement_static(p5_sense)
        return
    
    
    def optimze(
            self,
            fq_iter:iter,
            plens: List[int],
            n_iqr:List[int],
            processes: int,
            target: str,
            beta: float
        ) -> Dict[str, object]:
        out = {
            "left": self.optimizeLIP(fq_iter, plens, "left", n_iqr, processes, target, beta),
            "right": self.optimizeLIP(fq_iter, plens, "right", n_iqr, processes, target, beta)
        }
        return out
    
    def optimizeLIP(
            self,
            fq_iter:iter,
            plens: List[int],
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
        for plen in plens:
            res = self.optimizePI(fq_iter, plen, site, n_iqr, processes, target, beta)
            if res[target] > out[target]:
                out = res
                out["plen"] = plen
        return out
    
    
    def optimizePI(self, fq_iter: iter, plen: int, site: str, n_iqr: List[int], processes: int, target: str, beta: float) -> Dict[str, object]:
        """Find the optimal pid cutoff (float) and ideal searching location (int) at certain primer length

        Args:
            fq_iter (iter): an iterator to fastq file
            plen (int): length of the primer to use
            site (str): "left" or "right"
            
        Returns:
            Dict[str, object]: optimzed parameters
        """
        print(f"plen = {plen}: aligning {site} primers to the read", file=sys.stderr)
        out = {target: -1}
        # get pid, location, and class dataframe
        data = self.getPLC(fq_iter, plen, site, processes)
        
        # get pid targets
        pids = data["pid"].unique()
        pids.sort()
        pids = pids[::-1]
        
        # get loc targets
        med_loc = np.median(data[data["cls"] == "positive"]["loc"])
        iqr_loc = np.quantile(data[data["cls"] == "positive"]["loc"], .75) - np.quantile(data[data["cls"] == "positive"]["loc"], .25)
        isls = [med_loc + n * iqr_loc for n in n_iqr]
        
        # iter over pid and loc targets
        print_target = target if not target == "fscore" else f"{target} (beta = {beta})"
        print(
            f"plen = {plen}: optimizing {print_target}. PID cutoffs: {pids[0]}-{pids[-1]}; ISL cutoffs: {med_loc} + n * {iqr_loc}",
            file=sys.stderr
        )
        tasks = []
        for pid in pids:
            for location in isls:
                tasks.append((data, pid, location, beta))
                # tp = sum((data["cls"] == "positive") & (data["pid"] >= pid) & (data["loc"] <= location))
                # fp = sum((data["cls"] == "negative") & (data["pid"] >= pid) & (data["loc"] <= location))
                # tn = sum(data["cls"] == "negative") - fp
                # fn = sum(data["cls"] == "positive") - tp
                # accr = (tp + tn) / (tp + tn + fp + fn)
                # prec = tp / (tp + fp) if (tp + fp) > 0 else -1
                # if accr > out["accr"]:
                #     out["pid"] = pid
                #     out["loc"] = location
                #     out["accr"] = accr
                #     out["prec"] = prec
                # elif (accr == out["accr"]) and (prec > out["prec"]):
                #     out["pid"] = pid
                #     out["loc"] = location
                #     out["accr"] = accr
                #     out["prec"] = prec
                # # else:
                # #     print(f"pid = {pid}: stop trying at location = {location}, accuracy = {accr}, precision = {prec}", file=sys.stderr)
                # #     break
        with Pool(processes) as pool:
            result = pool.starmap_async(self.calculateAP, tasks)
            r = result.get()
        
        for res in r:
            if res[target] > out[target]:
                out = res
            # elif (res["accr"] == out["accr"]) and (res["prec"] > out["prec"]):
            #     out = res
        return out
    
    
    @staticmethod
    def calculateAP(data: pd.DataFrame, pid_cutoff: int, isl_cutoff: int, beta: float) -> Dict[str, int or float]:
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
    def getPLC(self, fq_iter: iter, plen: int, site: str, processes: int) -> pd.DataFrame:
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
        with Pool(processes) as pool:
            # for read in fq_iter:
            # get [pid, location, class] of top 2 best aligned primers
            # res = nbest_aligner(read, plen, 2)
            # data["pid"] += [i.pid for i in res]
            # data["loc"] += [i.location for i in res]
            # data["cls"] += ["positive", "negative"]
            result = [pool.apply_async(nbest_aligner, (read, plen, 2)) for read in fq_iter]
            data = []
            for i in result:
                data += i.get()
            data = pd.DataFrame(data, columns=["pid", "loc"])
            data["cls"] = ["positive", "negative"] * (len(data) // 2) # `//` to avoid float
        # return as a pd.dataframe 
        return pd.DataFrame(data)
    
    
    # n best 5' primers
    def nbestsl(self, read: SeqFastq, plen: int, nbests: int) -> List[Tuple[int, int]]:
        """Get the best n alignments of the primers to a read

        Args:
            read (SeqFastq): the read to be aligned
            plen (int): length of the primer to use
            nbests (int): number of best alignments to return

        Returns:
            List[Tuple[int, int]]
        """
        sense = \
            [(i["pid"], i["location"][1]) for i in 
            aligner.ntopAligns(
                self.p5_sense[-plen:],
                read.seq,
                "HW",
                "locations",
                -1,
                nbests
            )
        ]
        anti = \
            [(i["pid"], i["location"][1]) for i in 
            aligner.ntopAligns(
                self.p5_anti[-plen:],
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
    def nbestsr(self, read: SeqFastq, plen: int, nbests: int) -> List[Tuple[int, int]]:
        """Get the best n alignments of 3' primers to a read

        Args:
            read (SeqFastq): the read to be aligned
            plen (int): length of the primer to use
            nbests (int): number of best alignments to return

        Returns:
            List[SimpleResult(pid, location)]
        """
        sense =  \
            [(i["pid"], len(read.seq) - i["location"][0]) for i in 
            aligner.ntopAligns(
                self.p3_sense[:plen],
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
                self.p3_anti[:plen],
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
    
    # def alignP2P(self, plen: int) -> dict:
    #     """Align primer to primers (5' to 5'; 3' to 3')

    #     Args:
    #         plen (int): the length of the primer substring to be aligned

    #     Returns:
    #         dict: {"p5": pid, "p3": pid}
    #     """
    #     out = {}
    #     out["p5"] = aligner.singleAlign(
    #             self.p5_sense[-plen:],
    #             self.p5_anti[-plen:],
    #             "HW",
    #             "locations",
    #             -1
    #         )["pid"]
    #     out["p3"] = aligner.singleAlign(
    #             self.p3_sense[:plen],
    #             self.p3_anti[:plen],
    #             "HW",
    #             "locations",
    #             -1
    #         )["pid"]
    #     return out
    
    
    # def alignInB(self, read: SeqFastq) -> Dict[str, Tuple[float, float]]:
    #     """Align primer to ISL and read body

    #     Args:
    #         read (SeqFastq): _description_

    #     Returns:
    #         Dict[Tuple[float, float]]: 
    #         {
    #             "p5": (pid_isl, pid_body),
    #             "p3": (pid_isl, pid_body)
    #         }
    #     """
        
    #     out = {}
    #     # align 5' primers to 5' isl
    #     strand, res = aligner.bestAlign(
    #         {1: self.p5_sense, -1: self.p5_anti},
    #         read.seq[self.isl5[0]:self.isl5[1]],
    #         mode="HW",
    #         task="locations",
    #         pid=-1
    #     )
    #     pid_body = aligner.singleAlign(
    #         {1: self.p5_sense, -1: self.p5_anti}[strand],
    #         read.seq[self.isl5[1]:self.isl3[0]],
    #         mode="HW",
    #         task="locations",
    #         pid=-1
    #     )["pid"]
    #     out["p5"] = (res["pid"], pid_body)

    #     # align 3' primers to 3' isl
    #     strand, res = aligner.bestAlign(
    #         {1: self.p3_sense, -1: self.p3_anti},
    #         read.seq[self.isl3[0]:self.isl3[1]],
    #         mode="HW",
    #         task="locations",
    #         pid=-1,
    #     )
    #     pid_body = aligner.singleAlign(
    #         {1: self.p3_sense, -1: self.p3_anti}[strand],
    #         read.seq[self.isl5[1]:self.isl3[0]],
    #         mode="HW",
    #         task="locations",
    #         pid=-1
    #     )["pid"]
    #     out["p3"] = (res["pid"], pid_body)
        
    #     return out