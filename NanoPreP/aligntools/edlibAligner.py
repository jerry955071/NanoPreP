"""Caller of edlib.align()"""
import edlib


class edlibAligner:
    def singleAlign(
        query: str,
        target: str,
        mode: str,
        task: str,
        pid: float
    ):
        # call edlib for the alignment
        res = edlib.align(
            query,
            target,
            mode,
            task,
            round((1 - pid) * len(query))
        )

        # add pid to result
        if res["editDistance"] >= 0:
            res["pid"] = 1 - res["editDistance"] / len(query)
        else:
            res["pid"] = 0

        return res

    def bestAlign(
        querys: dict,
        target: str,
        mode: str,
        task: str,
        pid: float
    ):
        res = name = None
        for qname, query in querys.items():
            new = edlibAligner.singleAlign(query, target, mode, task, pid)
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
