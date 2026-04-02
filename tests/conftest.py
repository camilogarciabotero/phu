"""pytest configuration and fixtures."""

import sys
from pathlib import Path

# Ensure src/ is in sys.path for imports
src_path = Path(__file__).parent.parent / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))
