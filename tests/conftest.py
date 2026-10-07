import pytest


@pytest.fixture(autouse=True)
def _fresh_tile_records():
    """ATL1415.tile_meta caches each tiles_dir's records for the life of the
    process (one writer run); tests rewrite tiles in place, so start clean."""
    from ATL1415 import tile_meta
    tile_meta._CACHE.clear()
    yield
    tile_meta._CACHE.clear()
