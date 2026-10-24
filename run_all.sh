#!/usr/bin/env bash
# Run the full pipeline from the project root.
#   PY=path/to/python RS=Rscript ./run_all.sh
# Step 0 needs internet (downloads about 5 MB from figshare). Skip it with SKIP_DOWNLOAD=1.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-python3}"
RS="${RS:-Rscript}"
if [ "${SKIP_DOWNLOAD:-0}" != "1" ]; then
  "$PY" -I python/00_fetch_soilhealthdb.py data_raw
fi
for s in 01_load_clean 02_effect_sizes 03_random_effects 04_moderators 05_bias_sensitivity 07_great_plains_subset 06_figures 99_session_info; do
  echo "---- $s"
  "$RS" "R/$s.R"
done
echo "Done. Tables in outputs/tables, figures in outputs/figures, logs in outputs/logs."
