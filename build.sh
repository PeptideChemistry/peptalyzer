#!/bin/bash

set -e  # Exit if any command fails

echo "✅ Starting build.sh for wkhtmltopdf setup"

mkdir -p bin
cd bin

# Download and extract the CentOS version of wkhtmltopdf
curl -LO https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6-1/wkhtmltox-0.12.6-1.centos7.x86_64.rpm

# Install rpm2cpio if not present
command -v rpm2cpio >/dev/null 2>&1 || {
  echo >&2 "❌ rpm2cpio is not available"; exit 1;
}

# Extract binary
rpm2cpio wkhtmltox-0.12.6-1.centos7.x86_64.rpm | cpio -idmv

# Copy binary to project bin dir
cp usr/local/bin/wkhtmltopdf ../bin/wkhtmltopdf
chmod +x ../bin/wkhtmltopdf

# Cleanup
cd ..
rm -rf bin/usr
rm -f bin/wkhtmltox-0.12.6-1.centos7.x86_64.rpm

echo "✅ wkhtmltopdf installed at: /app/bin/wkhtmltopdf"