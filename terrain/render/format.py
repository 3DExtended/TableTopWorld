"""Serialize render specs to stable text for tests and inspection."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from terrain.render.spec import FlowerRenderSpec


def format_render_spec(spec: FlowerRenderSpec) -> str:
    """Deterministic, human-readable text (JSON with sorted keys)."""
    return json.dumps(spec.to_dict(), indent=2, sort_keys=True) + "\n"
