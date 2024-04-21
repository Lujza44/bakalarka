import pandas as pd
import xlsxwriter

data = pd.read_csv('data/vis/raw_vis.csv')

workbook = xlsxwriter.Workbook('data/vis/Table_S1.xlsx', {'nan_inf_to_errors': True})
worksheet = workbook.add_worksheet()

red_format = workbook.add_format({'color': 'red', 'underline': 1, 'font_name': 'Courier New'})
default_format = workbook.add_format({'color': 'black', 'font_name': 'Courier New'})

def color_red(sequence, indices):
    parts = []
    last_index = 0
    for i in range(len(sequence)):
        if i in indices:
            if last_index != i:
                parts.append(default_format)
                parts.append(sequence[last_index:i])
            parts.append(red_format)
            parts.append(sequence[i])
            last_index = i + 1
    if last_index < len(sequence):
        parts.append(default_format)
        parts.append(sequence[last_index:])
    return parts

max_widths = {}

def write_and_track_width(row, col, value, cell_format=default_format):
    if cell_format:
        worksheet.write(row, col, value, cell_format)
    else:
        worksheet.write(row, col, value)
    max_widths[col] = max(max_widths.get(col, 0), len(str(value)))

skipped_columns = [7, 9]

write_col_index = 0
for col_num, value in enumerate(data.columns):
    if col_num in skipped_columns:
        continue
    write_and_track_width(0, write_col_index, value)
    write_col_index += 1

for row_index, row in data.iterrows():
    write_col_index = 0
    for col_index, value in enumerate(row):
        if pd.isna(value):
            value = ''
        if col_index in skipped_columns:
            continue
        elif col_index == 6 and pd.notna(row['5\' SNP indexes']):
            snps_5_indexes = [int(num.strip()) for num in row['5\' SNP indexes'].split(",")]
            worksheet.write_rich_string(row_index + 1, write_col_index, *color_red(value, snps_5_indexes))
            max_widths[write_col_index] = max(max_widths.get(write_col_index, 0), len(str(value)))
        elif col_index == 8 and pd.notna(row['3\' SNP indexes']):
            snps_3_indexes = [int(num.strip()) for num in row['3\' SNP indexes'].split(",")]
            worksheet.write_rich_string(row_index + 1, write_col_index, *color_red(value, snps_3_indexes))
            max_widths[write_col_index] = max(max_widths.get(write_col_index, 0), len(str(value)))
        else:
            write_and_track_width(row_index + 1, write_col_index, value)
        write_col_index += 1

for col_index, width in max_widths.items():
    worksheet.set_column(col_index, col_index, width + 1)
worksheet.set_column('G:G', 105)  
worksheet.set_column('H:H', 105)

workbook.close()
