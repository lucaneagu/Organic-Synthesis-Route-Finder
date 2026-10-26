"""Entry point for ``python -m chemsynth``."""

import os
import sys

from .cli import main

try:
    code = main()
except BrokenPipeError:
    # Someone piped us into `head` or `less` and closed the pipe early.
    os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    code = 0
except KeyboardInterrupt:
    code = 130

raise SystemExit(code)
