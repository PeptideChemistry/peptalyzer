#!/bin/bash
set -e

echo "✅ Starting wkhtmltopdf setup"

# Create temp extraction dir
mkdir -p tmp bin
cd tmp

# Download CentOS-compatible RPM
curl -LO https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6-1/wkhtmltox-0.12.6-1.centos7.x86_64.rpm

# Check if rpm2cpio is available
if ! command -v rpm2cpio &> /dev/null; then
  echo "❌ rpm2cpio not found — PDF export will fail!"
  exit 1
fi

# Extract RPM contents
rpm2cpio wkhtmltox-0.12.6-1.centos7.x86_64.rpm | cpio -idmv

# Copy binary to final bin dir
cp usr/local/bin/wkhtmltopdf ../b