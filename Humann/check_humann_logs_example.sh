#!/usr/bin/env bash
set -euo pipefail

cd /media/uhlemann/core5/01_MG/20260409_DEAPIM30

tail -n 100 humann3_out/6_logs/DPM10014_S213.humann.log || true
echo
tail -n 100 humann3_out/1_humann3_out/DPM10014_S213_humann_temp/DPM10014_S213.log || true
