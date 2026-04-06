#!/bin/bash

INPUT_DIR="/eos/user/a/achoubey/htobbgg/WH_Signal_File"


for MASS in 15 20 25 30 35 40 45 50 55 60
do
    INPUT_FILE="${INPUT_DIR}/WH_2024_M${MASS}.root"
    OUTPUT_FILE="curated_signal_m${MASS}.root"

    echo "Running mass point M=${MASS}"
    python3 curate_signal_rdf_v2.py "${INPUT_FILE}" "${OUTPUT_FILE}" "${MASS}"
done
