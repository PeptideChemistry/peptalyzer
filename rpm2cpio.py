#!/usr/bin/env python3
import sys
import os
from cpioextract import extract_cpio_archive

def main(rpm_path):
    with open(rpm_path, "rb") as f:
        data = f.read()

    cpio_start = data.find(b"070701")
    if cpio_start == -1:
        sys.stderr.write("CPIO archive not found\n")
        sys.exit(1)

    archive = data[cpio_start:]

    out_dir = os.path.join(os.getcwd(), "extracted")
    os.makedirs(out_dir, exist_ok=True)

    with open("payload.cpio", "wb") as out:
        out.write(archive)

    extract_cpio_archive("payload.cpio", out_dir)
    os.remove("payload.cpio")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: rpm2cpio.py file.rpm")
        sys.exit(1)
    main(sys.argv[1])
