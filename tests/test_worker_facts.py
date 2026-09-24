"""
docs/plan_dps_speed.sh D1-D2: the job reports its worker (worker_facts.py),
the CPU each step got (run_with_rusage.py), and the bench step's timings
(bench_solve.py); ogc_jobs parses all three for collect_jobs.py.
"""
import os
import subprocess
import sys

import numpy as np
import pytest
import scipy.sparse as sp

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'maap'))
from ogc_jobs import parse_bench, parse_cpu, parse_worker  # noqa: E402
import bench_solve  # noqa: E402

WORKER_KEYS = {'instance_type', 'instance_id', 'az', 'cpu', 'vcpus', 'affinity',
               'threads_per_core', 'cpu_quota', 'mem_gib', 'mem_limit_gib',
               'blas', 'load1'}


def test_worker_facts_prints_one_parseable_line():
    out = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'worker_facts.py')],
                         capture_output=True, text=True, timeout=60)
    assert out.returncode == 0
    lines = [ln for ln in out.stdout.splitlines() if ln.startswith('WORKER: ')]
    assert len(lines) == 1
    facts = parse_worker(out.stdout)
    assert set(facts) == WORKER_KEYS
    # every value is one token, so the line stays key=value parseable
    assert all(v and ' ' not in v for v in facts.values())
    assert int(facts['vcpus']) >= 1


@pytest.mark.parametrize('code', [0, 3])
def test_rusage_reports_cpu_and_keeps_the_exit_status(code):
    burn = 'import time\nt=time.process_time()\nwhile time.process_time()-t<0.5: pass\n'
    out = subprocess.run([sys.executable, os.path.join(REPO, 'scripts', 'run_with_rusage.py'),
                          'burn', sys.executable, '-c', burn + f'raise SystemExit({code})'],
                         capture_output=True, text=True, timeout=60)
    assert out.returncode == code
    assert '=== rusage [burn]: elapsed' in out.stdout
    cpu = parse_cpu(out.stdout)['burn']
    # a busy loop gets about one core; allow a loaded test machine
    assert float(cpu['cpu_s']) >= 0.4
    assert 0.2 < float(cpu['cores']) <= 1.5


def test_parse_cpu_reads_na_fields():
    text = ('=== cpu [error]: cpu_time 12.0 s = na cores avg, '
            'load1 min/med/max na, steal na, throttled na\n')
    assert parse_cpu(text) == {'error': {'cpu_s': '12.0', 'cores': 'na', 'load': 'na',
                                         'steal_s': 'na', 'throttled_s': 'na',
                                         'n_throttled': '-'}}


def test_jobs_before_d1_parse_empty():
    old = '=== rusage [fit]: elapsed 1213.5 s, peak RSS 5.15 GiB (before this step: 0.00 GiB), exit 0\n'
    assert parse_worker(old) == {} and parse_cpu(old) == {} and parse_bench(old) == {}


def small_system(tmp_path, seed=0):
    rng = np.random.default_rng(seed)
    A = (sp.random(400, 120, density=0.05, random_state=rng) + sp.eye(400, 120)).tocsc()
    b = rng.standard_normal(400)
    x = np.linalg.lstsq(A.toarray(), b, rcond=None)[0]
    sp.save_npz(tmp_path / 'A0.npz', A)
    np.save(tmp_path / 'b0.npy', b)
    np.save(tmp_path / 'x0.npy', x)
    return x


def test_bench_times_and_checks_the_solve(tmp_path, capsys):
    small_system(tmp_path)
    assert bench_solve.main([str(tmp_path), '--threads', '1,2']) == 0
    bench = parse_bench(capsys.readouterr().out)
    assert set(bench) == {'qr_t1', 'qr_t2', 'dgemm1', 'spmv'}


def test_bench_fails_on_a_wrong_answer(tmp_path):
    x = small_system(tmp_path)
    np.save(tmp_path / 'x0.npy', x + 1)
    assert bench_solve.main([str(tmp_path), '--threads', '1']) == 1
