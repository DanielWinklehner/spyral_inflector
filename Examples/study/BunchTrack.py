"""Deck wrapper (no space charge): python -m spyral_inflector.tracking.bunch (original in legacy/)."""
import sys

from spyral_inflector.tracking.bunch import main

if __name__ == "__main__":
    main([a for a in sys.argv[1:] if a not in ("--sc", "--no-sc")])
