#!/usr/bin/env python3
"""
Time one saved least-squares system with SPQR, the way LSsurf's iterate_fit
solves it, at several thread counts -- the same script on the ADE and on a
DPS worker, so the two are compared on identical work (docs/plan_dps_speed.sh
D2).  No ATL11 read, no S3 traffic during the timing.

The system is iteration 0 of E1340_N-2420's quarterly prelim fit (1.81M x
477k, nnz 7.7M), saved from a local run on 2026-09-24: A0.npz is Ip_r TCinv
G (CSC), b0.npy the weighted right-hand side, x0.npy SPQR's solution, which
every timed solve is checked against.

Also times two single-threaded kernels that need no SuiteSparse:
  dgemm1   a 2000x2000 matrix product on one BLAS thread (best of 5):
           per-core floating-point speed
  spmv     200 products A @ x on the saved A: memory bandwidth

What the answers mean (H1/H2 in the plan):
  every number slower than the ADE by one factor  -> slower cores (H1)
  1-thread QR close to the ADE, 4-thread QR not   -> threads do not get
                                                     cores (H2); compare the
                                                     cpu lines run.sh prints
Usage:
  bench_solve.py <source> [--threads 1,2,4]
    source  a directory, local or s3://, holding A0.npz, b0.npy and x0.npy
Prints one line per measurement and a final BENCH: summary line.
"""
import argparse
import os
import sys
import tempfile
import time

import numpy as np
import scipy.sparse as sp
from threadpoolctl import threadpool_limits

FILES = ('A0.npz', 'b0.npy', 'x0.npy')


def fetch(source, dest):
    """Local paths to the three files, copying them from s3:// if needed."""
    if not source.startswith('s3://'):
        return [os.path.join(source, f) for f in FILES]
    import s3fs
    fs = s3fs.S3FileSystem()
    paths = []
    for f in FILES:
        path = os.path.join(dest, f)
        fs.get(f'{source.rstrip("/")}/{f}', path)
        paths.append(path)
    return paths


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('source')
    parser.add_argument('--threads', default='1,2,4')
    args = parser.parse_args(argv)
    threads = [int(t) for t in args.threads.split(',')]

    import sparseqr
    with tempfile.TemporaryDirectory() as tmp:
        t = time.time()
        a_path, b_path, x_path = fetch(args.source, tmp)
        A = sp.load_npz(a_path).tocsc()
        b = np.load(b_path).ravel()
        x_ref = np.load(x_path).ravel()
        print(f'bench: loaded {args.source} in {time.time() - t:.1f} s: '
              f'A {A.shape} nnz {A.nnz}', flush=True)

    summary = {}
    for n in threads:
        with threadpool_limits(limits={'openmp': n, 'blas': n}):
            t = time.time()
            x = np.asarray(sparseqr.solve(A.tocoo(), b, ordering=6)).ravel()
            dt = time.time() - t
        err = np.max(np.abs(x - x_ref))
        print(f'bench: SPQR solve, {n} thread(s): {dt:.1f} s  max|x - x0| = {err:.1e}',
              flush=True)
        # a wrong answer is not a timing: say so and fail
        if not err < 1e-6:
            print(f'bench: ERROR: the solve does not reproduce x0 (max diff {err:.1e})',
                  file=sys.stderr)
            return 1
        summary[f'qr_t{n}'] = f'{dt:.1f}'

    rng = np.random.default_rng(0)
    M = rng.standard_normal((2000, 2000))
    best = np.inf
    with threadpool_limits(limits={'blas': 1}):
        for _ in range(5):
            t = time.time()
            M @ M
            best = min(best, time.time() - t)
    print(f'bench: dgemm 2000x2000, 1 thread, best of 5: {best:.3f} s', flush=True)
    summary['dgemm1'] = f'{best:.3f}'

    x = np.ones(A.shape[1])
    t = time.time()
    for _ in range(200):
        A @ x
    dt = time.time() - t
    print(f'bench: 200 x (A @ x): {dt:.1f} s', flush=True)
    summary['spmv'] = f'{dt:.1f}'

    print('BENCH: ' + ' '.join(f'{k}={v}' for k, v in summary.items()), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
