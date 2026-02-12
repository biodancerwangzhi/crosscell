#!/bin/bash
# 预初始化 basilisk 环境
# 修复 GnuTLS 问题并初始化 zellkonverter 的 Python 环境

set -e

echo "=== Fixing git config for GnuTLS ==="
git config --global http.postBuffer 524288000
git config --global http.version HTTP/1.1

echo "=== Initializing basilisk environment for zellkonverter ==="
R -e '
library(zellkonverter)
library(basilisk)
env <- zellkonverter:::zellkonverterAnnDataEnv()
proc <- basiliskStart(env)
cat("basilisk environment initialized OK\n")
basiliskStop(proc)
'

echo "=== Done ==="
