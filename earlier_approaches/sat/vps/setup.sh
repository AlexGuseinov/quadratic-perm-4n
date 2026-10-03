#!/bin/bash
# One-time setup on a fresh Ubuntu VPS (tested plan; run as a normal user with sudo).
set -e
sudo apt-get update -y
sudo apt-get install -y build-essential git python3 python3-pip
pip3 install --break-system-packages python-sat || pip3 install python-sat
mkdir -p ~/tools && cd ~/tools
[ -d kissat ] || git clone --depth 1 https://github.com/arminbiere/kissat.git
(cd kissat && ./configure && make -j"$(nproc)")
[ -d cadical ] || git clone --depth 1 https://github.com/arminbiere/cadical.git
(cd cadical && ./configure && make -j"$(nproc)")
[ -d lrat-trim ] || git clone --depth 1 https://github.com/arminbiere/lrat-trim.git
(cd lrat-trim && ./configure && make) || echo "lrat-trim build failed (optional)"
echo "tools ready in ~/tools"
