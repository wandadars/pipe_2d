#!/usr/bin/env python3
import sys

from pipe_grid import ChannelGrid


def print_usage():
    print("Usage: python3 channel_grid.py nI nJ nK xMin xMax yMin yMax zMin zMax")


def main():
    if len(sys.argv) < 2 or sys.argv[1].lower() in ("-h", "--help", "help"):
        print_usage()
        sys.exit(0)

    try:
        grid = ChannelGrid()
        grid.write_grid()
        print("\nGrid written to file.grd")
    except ValueError as error:
        print(error)
        sys.exit(1)


if __name__ == "__main__":
    main()
