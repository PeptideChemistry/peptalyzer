#!/usr/bin/env python3
import sys
import os

def extract_cpio(archive_data, out_dir):
    pos = 0
    while pos + 110 <= len(archive_data):
        header = archive_data[pos:pos+110]
        if header[:6] != b'070701':
            break  # Invalid header, stop

        name_size = int(header[94:102], 16)
        file_size = int(header[54:62], 16)

        pos += 110
        name = archive_data[pos:pos+name_size].rstrip(b'\x00').decode()
        pos += ((name_size + 3) // 4) * 4  # Align

        if name == "TRAILER!!!":
            break

        file_path = os.path.join(out_dir, name)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        with open(file_path, "wb") as f:
            f.write(archive_data[pos:pos+file_size])

        pos += ((file_size + 3) // 4) * 4  # Align

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
    extract_cpio(archive, out_dir)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: rpm2cpio.py file.rpm")
        sys.exit(1)
    main(sys.argv[1])
