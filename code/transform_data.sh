#!/bin/bash

JSON_FILE="data/output/transformed_data.json"

if [ -f "$JSON_FILE" ]; then
    rm "$JSON_FILE"
fi

python3 code/1_csv_to_json.py
python3 code/2_xlsx_to_json.py
python3 code/3_fix_problems.py
python3 code/4_rs_to_json.py