"""
time_coverage_* root attributes (ATL1415_attrs_meta.set_time_range).

The duration must agree with the start and end attributes written beside it:
start + duration == end, in whole seconds.  The bug this guards against
(docs/plan_monthly_on_maap.sh M0b) computed
int((datetime_start-datetime_end).seconds), which reversed the operands and
then kept only the within-day part of the negative timedelta, so a 7.5-year
product reported 54385 s -- about 15 hours.  No files here: set_time_range
only writes root_info and two METADATA/Extent attributes.
"""
from datetime import datetime
from types import SimpleNamespace

import pytest

from ATL1415.ATL1415_attrs_meta import set_time_range

FMT = '%Y-%m-%dT%H:%M:%S.%fZ'


class FakeExtent:
    def __init__(self):
        self.attrs = {}

    def setncattr(self, name, value):
        self.attrs[name] = value


class FakeDst(dict):
    """Only /METADATA/Extent is touched by set_time_range."""
    def __init__(self):
        super().__init__({'/METADATA/Extent': FakeExtent()})


def run(t_crop, region='IS'):
    dst = FakeDst()
    root_info = {}
    set_time_range(dst, root_info,
                   SimpleNamespace(t_crop=t_crop, region=region, verbose=False))
    return dst, root_info


@pytest.mark.parametrize('region', ['IS', 'GL', 'AA'])
def test_duration_matches_the_start_and_end_written(region):
    # the invariant that actually matters: the three attributes agree
    _, root_info = run((2019.0, 2026.5), region=region)
    start = datetime.strptime(root_info['time_coverage_start'], FMT)
    end = datetime.strptime(root_info['time_coverage_end'], FMT)
    assert root_info['time_coverage_duration'] == int((end - start).total_seconds())


def test_known_span_of_the_0332_quarterly():
    # --t_crop=2019,2026.5 with region IS, the five 0332 files re-written in M0b
    _, root_info = run((2019.0, 2026.5))
    assert root_info['time_coverage_duration'] == 236681615
    assert root_info['time_coverage_start'] == '2019-01-01T00:06:25.000000Z'
    assert root_info['time_coverage_end'] == '2026-07-02T09:00:00.000000Z'


def test_duration_is_positive_and_not_truncated_to_within_a_day():
    # the shape of the old bug: positive, and larger than one day
    _, root_info = run((2019.0, 2026.5))
    assert root_info['time_coverage_duration'] > 86400
    assert root_info['time_coverage_duration'] != 54385


def test_one_year_span():
    # a whole year is 365.25 days by the function's own epoch arithmetic,
    # less the region offset that set_time_range adds to the start
    _, root_info = run((2020.0, 2021.0))
    offset = ord('I') * 3 + ord('S') * 2
    assert root_info['time_coverage_duration'] == int(365.25 * 86400) - offset


def test_duration_is_a_plain_int_of_seconds():
    # rel005 carries 2.17e8 and main computes seconds too: numeric, not ISO 8601
    _, root_info = run((2019.0, 2026.5))
    assert isinstance(root_info['time_coverage_duration'], int)


def test_extent_range_attributes_mirror_the_root_ones():
    dst, root_info = run((2019.0, 2026.5))
    extent = dst['/METADATA/Extent'].attrs
    assert extent['rangeBeginningDateTime'] == root_info['time_coverage_start']
    assert extent['rangeEndingDateTime'] == root_info['time_coverage_end']
