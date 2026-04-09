from __future__ import annotations

from collections.abc import Callable
from typing import Any

import click


def run_click_task(label: str, func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """Run a blocking task with a minimal Click progress indicator."""
    with click.progressbar(
        length=1,
        label=label,
        show_eta=False,
        show_percent=False,
        show_pos=False,
    ) as bar:
        result = func(*args, **kwargs)
        bar.update(1)
    return result