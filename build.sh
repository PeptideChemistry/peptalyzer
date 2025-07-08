#!/bin/bash

echo "Downloading wkhtmltopdf binary..."
mkdir -p bin
curl -L https://github.com/wkhtmltopdf/packaging/releases/download/0.12.6-1/wkhtmltox-0.12.6-1.centos7.x86_64.rpm -o wkhtmltox.rpm
rpm2cpio wkhtmltox.rpm | cpio -idmv
cp wkhtmltox/bin/wkhtmltopdf bin/wkhtmltopdf
chmod +x bin/wkhtmltopdf