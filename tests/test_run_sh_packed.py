"""
run.sh --tiles: several prelim/matched tiles in one job (plan_pack_tiles.sh K4).

run.sh itself is run, with `conda` stubbed on PATH: the stub runs the command
directly, ATL11_to_ATL15.py is a fake that writes (or does not write) the tile
named by --xy0, and s3_tiles.py's get/put copy within a local "bucket"
directory.  So the argument handling, lanes, per-tile isolation and
TILE_STATUS lines are tested; nothing is solved and nothing leaves the machine.
"""
import os
import re
import stat
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

CONDA = r'''#!/usr/bin/env bash
# conda run --no-capture-output -n <env> <cmd...>
shift 4
if [ "$1" = python ] && [[ "${2:-}" == */s3_tiles.py ]]; then
    shift 2; exec "$FAKE_BIN/s3_tiles" "$@"
fi
if [ "$1" = python ] && [[ "${2:-}" == */worker_facts.py ]]; then
    echo "WORKER: fake"; exit 0
fi
exec "$@"
'''

S3_TILES = r'''#!/usr/bin/env bash
# put <file> <prefix> <step>   |   get <prefix> prelim <x0> <y0> <spacing> <dir>
set -e
root=$FAKE_BUCKET/${2#s3://}
case "$1" in
    put) root=$FAKE_BUCKET/${3#s3://}; mkdir -p "$root/$4"; cp "$2" "$root/$4/" ;;
    get) awk -v x="$4" -v y="$5" -v d="$6" 'BEGIN{for(i=-1;i<=1;i++)for(j=-1;j<=1;j++)
             printf "E%d_N%d.h5\n", int((x+i*d)/1000), int((y+j*d)/1000)}' |
         while read -r f; do
             if [ -f "$root/$3/$f" ]; then cp "$root/$3/$f" "$7/"; else echo "missing $f"; fi
         done ;;
esac
'''

# --xy0 X Y (not matched), --base_directory (last wins), --out_name, --data_file, --prelim,
# --matched, --calc_error_for_xy.  FAKE_NODATA / FAKE_FAIL: "x,y x,y ..."
SOLVER = r'''#!/usr/bin/env python3
import os, sys, time
a = sys.argv[1:]
def last(opt):
    v = [a[i + 1] for i, s in enumerate(a) if s == opt]
    return v[-1] if v else None
if '--xy0' in a:
    i = a.index('--xy0'); x, y = a[i + 1], a[i + 2]
    xy = '%s,%s' % (x, y)
    name = 'E%d_N%d.h5' % (int(float(x) / 1000), int(float(y) / 1000))
else:   # matched: the tile is named by --out_name
    name = os.path.basename(last('--out_name')); xy = None
threads = [s for s in a if s.startswith('--THREADS=')][0]
with open(os.environ['FAKE_CALLS'], 'a') as fh:
    fh.write('%s %s %s %s\n' % (os.getpid(), name, threads,
             'error' if '--calc_error_for_xy' in a else ('matched' if '--matched' in a else 'fit')))
time.sleep(float(os.environ.get('FAKE_SLEEP', '0')))
if xy in os.environ.get('FAKE_FAIL', '').split():
    print('fake solver: failing %s' % name); sys.exit(1)
if '--matched' in a:
    assert os.path.isfile(last('--data_file')), last('--data_file')
    open(last('--out_name'), 'w').write('matched')
elif '--calc_error_for_xy' not in a:
    if xy in os.environ.get('FAKE_NODATA', '').split():
        print('fake solver: no data for %s' % name); sys.exit(0)
    d = os.path.join(last('--base_directory'), 'prelim')
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, name), 'w').write('prelim')
'''


@pytest.fixture
def job(tmp_path):
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    for name, text in [('conda', CONDA), ('s3_tiles', S3_TILES), ('ATL11_to_ATL15.py', SOLVER)]:
        f = bin_dir / name
        f.write_text(text)
        f.chmod(f.stat().st_mode | stat.S_IXUSR)
    args = tmp_path / 'input_args_XX.txt'
    args.write_text('-W=40000\n--tile_spacing=40000\n-b=/nowhere\n')
    work = tmp_path / 'job'
    work.mkdir()
    bucket = tmp_path / 'bucket'
    env = {k: v for k, v in os.environ.items() if not k.startswith(('MAAP_', 'ATL1415_', 'FAKE_'))}
    env.update(PATH='%s:%s' % (bin_dir, env['PATH']), FAKE_BIN=str(bin_dir),
               FAKE_BUCKET=str(bucket), FAKE_CALLS=str(tmp_path / 'calls'),
               ATL1415_THREADS='')

    def run(*argv, **fake):
        e = dict(env, **{k: str(v) for k, v in fake.items()})
        if not e['ATL1415_THREADS']:
            del e['ATL1415_THREADS']
        p = subprocess.run([str(REPO / 'run.sh'), '--args_file', str(args), *argv],
                           cwd=work, env=e, capture_output=True, text=True)
        calls = (tmp_path / 'calls').read_text().split('\n')[:-1] if (tmp_path / 'calls').exists() else []
        return p, calls

    run.work, run.bucket = work, bucket
    return run


def statuses(out):
    """the summary's TILE_STATUS lines (after the 'packed job summary' banner)"""
    summary = out.split('packed job summary', 1)[1]
    return dict(re.findall(r'^TILE_STATUS (\S+) (\S+)$', summary, re.M))


TILES = '0,0;40000,0|0,40000;-40000,-40000'


def test_packed_prelim_ok_nodata_failed(job):
    p, calls = job('--x0', '0', '--y0', '0', '--step', 'prelim', '--tiles', TILES,
                   '--tile_prefix', 's3://b/run',
                   FAKE_NODATA='40000,0', FAKE_FAIL='0,40000')
    assert p.returncode == 1, p.stdout + p.stderr
    assert statuses(p.stdout) == {'E0_N0': 'ok', 'E40_N0': 'nodata',
                                  'E0_N40': 'failed', 'E-40_N-40': 'ok'}
    assert 'TILES: 2 ok, 1 nodata, 1 failed, of 4' in p.stdout
    # the ok tiles were uploaded, the others were not
    assert sorted(f.name for f in (job.bucket / 'b/run/prelim').iterdir()) == \
        ['E-40_N-40.h5', 'E0_N0.h5']
    # a failed tile does not stop the tile after it in its lane, and only the
    # tiles with a fit get an error step
    kinds = sorted((c.split()[1], c.split()[3]) for c in calls)
    assert kinds == sorted([('E0_N0.h5', 'fit'), ('E0_N0.h5', 'error'),
                            ('E40_N0.h5', 'fit'), ('E0_N40.h5', 'fit'),
                            ('E-40_N-40.h5', 'fit'), ('E-40_N-40.h5', 'error')])
    for name in ['E0_N0', 'E40_N0', 'E0_N40', 'E-40_N-40']:
        assert (job.work / 'output/tile_logs' / (name + '.log')).is_file()


def test_packed_all_ok_exits_0_and_lanes_overlap(job):
    p, calls = job('--x0', '0', '--y0', '0', '--step', 'prelim', '--tiles', TILES, FAKE_SLEEP=1)
    assert p.returncode == 0, p.stdout + p.stderr
    assert set(statuses(p.stdout).values()) == {'ok'}
    # 2 lanes: the first tile of each lane starts before either lane's first
    # tile is done; threads split between them (ATL1415_THREADS unset)
    starts = re.findall(r'^TILE_START (\S+) lane (\d)', p.stdout, re.M)
    assert {s[1] for s in starts} == {'1', '2'}
    first_end = p.stdout.index('TILE_END')
    assert p.stdout.index('TILE_START E0_N40') < first_end
    assert p.stdout.index('TILE_START E0_N0') < first_end


def test_packed_threads_override(job):
    p, calls = job('--x0', '0', '--y0', '0', '--step', 'prelim', '--tiles', TILES,
                   ATL1415_THREADS=3)
    assert p.returncode == 0, p.stdout + p.stderr
    assert {c.split()[2] for c in calls} == {'--THREADS=3'}


def test_packed_matched_reads_neighbours_per_tile(job):
    prelim = job.bucket / 'b/run/prelim'
    prelim.mkdir(parents=True)
    for name in ['E0_N0', 'E40_N0', 'E0_N40']:
        (prelim / (name + '.h5')).write_text('prelim')
    p, calls = job('--x0', '0', '--y0', '0', '--step', 'matched',
                   '--tiles', '0,0;40000,0|80000,80000', '--tile_prefix', 's3://b/run')
    # 80000,80000 has no prelim tile of its own: that tile fails, the rest do not
    assert p.returncode == 1, p.stdout + p.stderr
    assert statuses(p.stdout) == {'E0_N0': 'ok', 'E40_N0': 'ok', 'E80_N80': 'failed'}
    assert sorted(f.name for f in (job.bucket / 'b/run/matched').iterdir()) == ['E0_N0.h5', 'E40_N0.h5']
    # each tile fetched its own 3x3 into its own work directory
    assert sorted(f.name for f in (job.work / 'tile_work/E40_N0/input/prelim').iterdir()) == \
        ['E0_N0.h5', 'E0_N40.h5', 'E40_N0.h5']


def test_single_tile_unchanged(job):
    p, calls = job('--x0', '0', '--y0', '0', '--step', 'prelim', '--tiles', '-')
    assert p.returncode == 0, p.stdout + p.stderr
    assert 'TILE_STATUS' not in p.stdout and 'xy0         : 0 0' in p.stdout
    assert (job.work / 'output/prelim/E0_N0.h5').is_file()
    assert [c.split()[3] for c in calls] == ['fit', 'error']


@pytest.mark.parametrize('argv, message', [
    (['--x0', '1', '--y0', '0', '--tiles', '0,0;40000,0'], 'must be the first tile'),
    (['--x0', '0', '--y0', '0', '--tiles', '0,0;0,0'], 'appears twice'),
    (['--x0', '0', '--y0', '0', '--tiles', '0,0||40000,0'], 'empty lane'),
    (['--x0', '0', '--y0', '0', '--tiles', '0,0;40000'], 'is not <x>,<y>'),
    (['--x0', '0', '--y0', '0', '--step', 'mosaic', '--tiles', '0,0'], 'for steps prelim and matched'),
], ids=['x0_not_first', 'duplicate', 'empty_lane', 'bad_pair', 'mosaic_step'])
def test_bad_tiles_stop_before_any_solve(job, argv, message):
    if '--step' not in argv:
        argv = argv + ['--step', 'prelim']
    p, calls = job(*argv)
    assert p.returncode == 2 and message in p.stderr, p.stdout + p.stderr
    assert calls == []
