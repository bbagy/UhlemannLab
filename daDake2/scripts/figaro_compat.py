#!/usr/bin/env python3
"""Run legacy FIGARO under modern Python versions."""

import collections
import collections.abc
import os
import runpy
import sys


for name in (
    "Callable",
    "Iterable",
    "Mapping",
    "MutableMapping",
    "Sequence",
    "MutableSequence",
    "Set",
    "MutableSet",
):
    if not hasattr(collections, name):
        setattr(collections, name, getattr(collections.abc, name))

if len(sys.argv) < 2:
    raise SystemExit("Usage: figaro_compat.py /path/to/figaro.py [FIGARO arguments]")

target = sys.argv[1]
sys.argv = sys.argv[1:]
sys.path.insert(0, os.path.dirname(os.path.abspath(target)))
runpy.run_path(target, run_name="__main__")
