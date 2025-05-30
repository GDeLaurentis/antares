from pathlib import Path

cache_dir = Path.home() / ".cache" / "antares"
cache_dir.mkdir(parents=True, exist_ok=True)

CACHE_PATH = cache_dir
DISKCACHE_SIZE_LIMIT_IN_GB = 32

from .version import __version__  # noqa
from .terms.terms import Terms    # noqa
