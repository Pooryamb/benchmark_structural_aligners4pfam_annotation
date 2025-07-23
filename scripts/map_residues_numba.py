import os
import sys
import glob
import pandas as pd
import numpy as np
import numba as nb

@nb.njit
def map_residues_numba(qseq, tseq, qstart=1, tstart=1):
    n = len(qseq)
    mapping = np.zeros((n, 2), dtype=np.int16)
    qpos, tpos = qstart, tstart
    count = 0
    for i in range(n):
        qc = qseq[i]
        tc = tseq[i]
        if qc != ord('-') and tc != ord('-'):
            mapping[count, 0] = qpos
            mapping[count, 1] = tpos
            count += 1
            qpos += 1
            tpos += 1
        elif qc != ord('-'):
            qpos += 1
        elif tc != ord('-'):
            tpos += 1
    return mapping[:count]
