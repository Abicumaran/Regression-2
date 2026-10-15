#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
python3 -m pip install -r requirements.txt
exec python3 -m streamlit run app.py
