#!/usr/bin/env python3
# A pure Python implementation of rpm2cpio
import sys
import struct
import os

def main(rpm_path):
    with open(rpm_path, "rb") as f:
        data = f.read()

    lead = data[0:96]
    if lead[0:4] != b"\xed\xab\xee\xdb":
        sys.stderr.write("Not an RPM file\n")
        sys.exit(1)

    # Find start of cpio archive
    cpio_start = data.find(b"070701")
    if cpio_start == -1:
        sys.stderr.write("CPIO archive not found\n")
        sys.exit(1)

    archive = data[cpio_start:]
    sys.stdout.buffer.write(archive)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: rpm2cpio.py file.rpm")
        sys.exit(1)
    main(sys.argv[1])
