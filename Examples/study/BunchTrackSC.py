"""Deck wrapper (space charge on unless --no-sc): python -m spyral_inflector.tracking.bunch --sc
(original in legacy/)."""
import sys

from spyral_inflector.tracking.bunch import main

if __name__ == "__main__":
    argv = [a for a in sys.argv[1:] if a != "--sc"]
    if "--no-sc" in argv:
        argv.remove("--no-sc")
    else:
        argv.append("--sc")
    main(argv)
