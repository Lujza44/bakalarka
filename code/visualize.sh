#!/bin/bash

RAW_TABLE1="data/output/raw_table1.csv"
RAW_TABLE2="data/output/raw_table2.csv"
TABLES="data/output/Tables.xlsx"

remove_file_if_exists() {
    if [ -f "$1" ]; then
        rm "$1"
    fi
}

remove_file_if_exists "$RAW_TABLE1"
remove_file_if_exists "$RAW_TABLE2"
remove_file_if_exists "$TABLES"

python3 code/5_table1_vis_prep.py
python3 code/6_table2_vis_prep.py
python3 code/7_visualize.py