#!/usr/bin/env python3
"""Write a shields-style coverage SVG and append a GitHub Actions step summary."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def _text_width(text: str) -> int:
    return 7 * len(text) + 10


def _color(percent: float) -> str:
    if percent < 70:
        return "#e05d44"
    if percent < 90:
        return "#dfb317"
    return "#4c1"


def _badge_svg(percent: float) -> str:
    label = "coverage"
    message = f"{percent:.1f}%"
    left = _text_width(label)
    right = _text_width(message)
    width = left + right
    color = _color(percent)
    aria = f"{label}: {message}"
    font = "DejaVu Sans,Verdana,Geneva,sans-serif"
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="20"
     role="img" aria-label="{aria}">
  <title>{aria}</title>
  <linearGradient id="s" x2="0" y2="100%">
    <stop offset="0" stop-color="#bbb" stop-opacity=".1"/>
    <stop offset="1" stop-opacity=".1"/>
  </linearGradient>
  <mask id="m"><rect width="{width}" height="20" rx="3" fill="#fff"/></mask>
  <g mask="url(#m)">
    <rect width="{left}" height="20" fill="#555"/>
    <rect x="{left}" width="{right}" height="20" fill="{color}"/>
    <rect width="{width}" height="20" fill="url(#s)"/>
  </g>
  <g fill="#fff" text-anchor="middle" font-family="{font}" font-size="11">
    <text x="{left / 2}" y="15" fill="#010101" fill-opacity=".3">{label}</text>
    <text x="{left / 2}" y="14">{label}</text>
    <text x="{left + right / 2}" y="15" fill="#010101" fill-opacity=".3">{message}</text>
    <text x="{left + right / 2}" y="14">{message}</text>
  </g>
</svg>
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_json", type=Path)
    parser.add_argument("badge_svg", type=Path)
    args = parser.parse_args()

    data = json.loads(args.summary_json.read_text(encoding="utf-8"))
    percent = float(data["line_percent"])
    message = f"{percent:.1f}%"

    args.badge_svg.parent.mkdir(parents=True, exist_ok=True)
    args.badge_svg.write_text(_badge_svg(percent), encoding="utf-8")

    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        with open(summary_path, "a", encoding="utf-8") as fh:
            fh.write(f"## Coverage\n\n**Lines:** {message}\n\n")
            fh.write("HTML report is attached as the `coverage-html` artifact.\n")

    print(f"line_coverage={percent:.1f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
