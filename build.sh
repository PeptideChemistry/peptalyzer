#!/bin/bash
set -e

echo "✅ Starting wkhtmltopdf setup..."

# Create bin directory
mkdir -p bin tmp
cd tmp

# Download wkhtmltopdf .rpm package for CentOS 7
curl -LO https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6-1/wkhtmltox-0.12.6-1.centos7.x86_64.rpm

# Extract RPM contents using embedded rpm2cpio and cpio
rpm2cpio() {
  rpm2cpio_input=$(cat)
  echo "$rpm2cpio_input" | cpio -idmv
}

cat wkhtmltox-0.12.6-1.centos7.x86_64.rpm | rpm2cpio

# Copy wkhtmltopdf binary to /app/bin
cp usr/local/bin/wkhtmltopdf ../bin/

cd ..
rm -rf tmp

echo "✅ wkhtmltopdf setup complete"
