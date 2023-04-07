"""Caller of edlib.align()"""
from typing import Tuple
import edlib

class edlibAligner:
    def singleAlign(
        query: str,
        target: str,
        mode: str,
        task: str,
        pid: float,
        tie_breaking: str = "middle"
    ) -> dict:
        # call edlib for the alignment
        res = edlib.align(
            query,
            target,
            mode,
            task,
            round((1 - pid) * len(query))
        )

        # add 2 fields: `pid` and `location` to `res`
        if res["editDistance"] >= 0:
            # pid
            res["pid"] = 1 - res["editDistance"] / len(query)
            
            # tie-breaking if edlib.align report more than one location
            idx = {
                "left": 0,
                "middle": len(res["locations"]) // 2,
                "right": -1
            }
            res["location"] = res.pop("locations")[idx[tie_breaking]] 
        else:
            res["pid"] = -1
            res["location"] = (-1, -1)
        
        return res

    def bestAlign(
        querys: dict,
        target: str,
        mode: str,
        task: str,
        pid: float,
        tie_breaking: str = "middle"
    ) -> Tuple[str, dict]:
        res = name = None
        # iterate over querys to find the best-aligned query
        for qname, query in querys.items():
            new = edlibAligner.singleAlign(
                query,
                target,
                mode,
                task,
                pid,
                tie_breaking
            )
            
            # the first alingment result
            if not res:
                res = new
                name = qname
                continue

            # if the new alignment result is better
            if new["pid"] > res["pid"]:
                res = new
                name = qname

        return name, res
