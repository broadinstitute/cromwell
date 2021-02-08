#!/usr/bin/env python3

from pathlib import Path

"""
Can be used to read files from the resources directory, like:
 (RESOURCES / filename).read_text()
"""
RESOURCES = Path(__file__).parent.parent / "resources"
