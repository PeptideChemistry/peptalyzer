#!/bin/bash
set -e

echo "✅ Starting wkhtmltopdf setup"

mkdir -p tmp bin
cd tmp

# Download CentOS-compatible RPM
curl -LO https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6-1/wkhtmltox-0.12.6-1.centos7.x86_64.rpm

# Extract RPM using Python script
python3 ../rpm2cpio.py wkhtmltox-0.12.6-1.centos7.x86_64.rpm

# Copy the binary
cp extracted/usr/local/bin/wkhtmltopdf ../bin/wkhtmltopdf
chmod +x ../bin/wkhtmltopdf

# Cleanup
cd ..
rm -rf tmp extracted