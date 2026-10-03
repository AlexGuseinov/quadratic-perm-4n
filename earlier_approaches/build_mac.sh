#!/bin/bash
# Builds drat-trim and constrank8, then verifies the n = 8 Desarguesian DRAT proof.
# On macOS it picks an SDK that the installed linker can use (newer SDKs can be
# unreadable by an older linker: "tapi error: malformed file").
# Usage:  bash earlier_approaches/build_mac.sh
set -e
cd "$(dirname "$0")"

TMP=$(mktemp -d)
echo 'int main(void){return 0;}' > "$TMP/t.c"

SDK=""
if cc "$TMP/t.c" -o "$TMP/t" 2>/dev/null; then
  echo "Default toolchain works."
else
  for s in $(ls -rd /Library/Developer/CommandLineTools/SDKs/MacOSX*.sdk 2>/dev/null); do
    if SDKROOT="$s" cc "$TMP/t.c" -o "$TMP/t" 2>/dev/null; then SDK="$s"; break; fi
  done
  if [ -z "$SDK" ]; then
    echo "No working SDK found. Update Command Line Tools (softwareupdate --list)."
    exit 1
  fi
  echo "Using SDK: $SDK"
fi
rm -rf "$TMP"

build() {  # build <args...>
  if [ -n "$SDK" ]; then SDKROOT="$SDK" cc "$@"; else cc "$@"; fi
}

[ -d drat-trim ] || git clone --depth 1 https://github.com/marijnheule/drat-trim.git
build -std=c99 -O2 -o drat-trim/drat-trim drat-trim/drat-trim.c
echo "built drat-trim/drat-trim"

build -O3 -o structured/constrank8 structured/constrank8.c
echo "built structured/constrank8"

echo
echo "== DRAT check of the n = 8 Desarguesian case (expect: s VERIFIED)"
./drat-trim/drat-trim sat/cnf/d8_desarg.cnf sat/proofs/d8_desarg.drat | tr "\r" "\n" | grep -E "^s |verification time"
