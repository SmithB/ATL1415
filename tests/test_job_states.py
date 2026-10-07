"""
scripts/maap/ogc_jobs.JobStates: job statuses from list calls, with
get_job_status only for jobs that left the active lists.  A fake MAAP client;
no network.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts', 'maap'))
from ogc_jobs import JobStates  # noqa: E402


class Resp:
    def __init__(self, body, code=200):
        self.body, self.status_code = body, code

    def json(self):
        return self.body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f'HTTP {self.status_code}')


class FakeMaap:
    """jobs: {id: status}; pages capped at `cap` like the real API (250)"""
    def __init__(self, jobs, cap=3):
        self.jobs, self.cap, self.status_calls, self.list_calls, self.down = jobs, cap, 0, 0, False

    def list_jobs(self, status=None, page_size=10, offset=0, get_job_details=True):
        self.list_calls += 1
        if self.down:
            return Resp({}, 502)
        ids = sorted(j for j, s in self.jobs.items() if s == status)
        page = ids[offset:offset + min(page_size, self.cap)]
        return Resp({'jobs': [{'jobID': j, 'status': status} for j in page]})

    def get_job_status(self, j):
        self.status_calls += 1
        return Resp({'status': self.jobs[j]})


def test_active_jobs_cost_list_calls_only():
    jobs = {f'j{i}': 'running' if i % 2 else 'accepted' for i in range(10)}
    m = FakeMaap(jobs)
    st = JobStates(m).poll(list(jobs))
    assert st == jobs and m.status_calls == 0
    assert m.list_calls == 6            # per status: 5 jobs at cap 3 = pages of 3, 2, then an empty one


def test_finished_jobs_are_asked_once():
    jobs = {f'j{i}': 'running' for i in range(6)}
    m = FakeMaap(jobs)
    js = JobStates(m)
    js.poll(list(jobs))
    jobs['j0'], jobs['j1'] = 'successful', 'failed'
    st = js.poll(list(jobs))
    assert st['j0'] == 'successful' and st['j1'] == 'failed' and m.status_calls == 2
    js.poll(list(jobs))
    assert m.status_calls == 2          # never asked again


def test_budget_leaves_the_rest_unknown_and_prefers_jobs_that_just_left():
    jobs = {f'j{i}': 'successful' for i in range(5)}
    jobs['run'] = 'running'
    m = FakeMaap(jobs)
    js = JobStates(m)
    st = js.poll(list(jobs), max_status_calls=2)
    assert sum(v == 'unknown' for v in st.values()) == 3 and m.status_calls == 2
    jobs['run'] = 'failed'
    st = js.poll(list(jobs), max_status_calls=1)
    assert st['run'] == 'failed'        # the one that left the lists goes first


def test_listing_outage_keeps_last_states():
    jobs = {'a': 'running', 'b': 'accepted'}
    m = FakeMaap(jobs)
    js = JobStates(m)
    js.poll(list(jobs))
    m.down = True
    jobs['a'] = 'successful'
    st = js.poll(list(jobs))
    assert st == {'a': 'running', 'b': 'accepted'} and js.error and m.status_calls == 0
