import pandas as pd
import xlsxwriter
import csv


workbook = xlsxwriter.Workbook('data/vis/Tables.xlsx', {'nan_inf_to_errors': True})


# TABULKA S1
worksheet1 = workbook.add_worksheet("Table S1")

data = pd.read_csv('data/vis/raw_vis.csv')

# formatovanie prveho sheetu
red_format = workbook.add_format({'color': 'red', 'underline': 1, 'font_name': 'Courier New'})
default_format = workbook.add_format({'color': 'black', 'font_name': 'Courier New'})

# funkcia, ktora zafarbi vsetky najdene SNP na cerveno
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

# funkcia, ktora pomaha nastavit dostatocnu sirku stlpcov
def write_and_track_width(row, col, value, cell_format=default_format):
    if cell_format:
        worksheet1.write(row, col, value, cell_format)
    else:
        worksheet1.write(row, col, value)
    max_widths[col] = max(max_widths.get(col, 0), len(str(value)))


skipped_columns = [7, 9]

# zapisanie heads okrem dvoch stlpcov, ktore boli len pomocne a vynechame ich
write_col_index = 0
for col_num, value in enumerate(data.columns):
    if col_num in skipped_columns:
        continue
    write_and_track_width(0, write_col_index, value)
    write_col_index += 1

# prepisanie naformatovanych dat, pomocne stlpce s indexami SNPov sa neprepisuju
for row_index, row in data.iterrows():
    write_col_index = 0
    for col_index, value in enumerate(row):
        if pd.isna(value):
            value = ''
        if col_index in skipped_columns:
            continue
        elif col_index == 6 and pd.notna(row['5\' SNP indexes']):
            snps_5_indexes = [int(num.strip()) for num in row['5\' SNP indexes'].split(",")]
            worksheet1.write_rich_string(row_index + 1, write_col_index, *color_red(value, snps_5_indexes))
            max_widths[write_col_index] = max(max_widths.get(write_col_index, 0), len(str(value)))
        elif col_index == 8 and pd.notna(row['3\' SNP indexes']):
            snps_3_indexes = [int(num.strip()) for num in row['3\' SNP indexes'].split(",")]
            worksheet1.write_rich_string(row_index + 1, write_col_index, *color_red(value, snps_3_indexes))
            max_widths[write_col_index] = max(max_widths.get(write_col_index, 0), len(str(value)))
        else:
            write_and_track_width(row_index + 1, write_col_index, value)
        write_col_index += 1

# nastavenie potrebnej sirky stlpcov
for col_index, width in max_widths.items():
    worksheet1.set_column(col_index, col_index, width + 1)
worksheet1.set_column('G:G', 105)  
worksheet1.set_column('H:H', 105)


# TABULKA S6
worksheet2 = workbook.add_worksheet('Table S6')

# formatovanie
format_A = workbook.add_format({'bg_color': '#92d050', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
format_G = workbook.add_format({'bg_color': '#ffff00', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
format_C = workbook.add_format({'bg_color': '#adb9ca', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
format_T = workbook.add_format({'bg_color': '#ffc7ce', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
vertical_format = workbook.add_format({'rotation': 90, 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
center_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'}) 
default_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'color': 'black', 'font_name': 'Courier New'})

# prepisanie sformatovanych dat z pripravneho suboru do tabulky
with open('data/vis/raw_vis6.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    column_widths = {}  
    for row_idx, row in enumerate(reader):
        if row[1] == 'Allele' or row[2] == 'GRCh38 coordinates':
            worksheet2.set_row(row_idx, 100)
        if row[2] == 'Distance from repeat region':
             worksheet2.set_row(row_idx, 35)
        for col_idx, value in enumerate(row):
            try: 
                float_value = float(value)
                int_value = int(float_value)
                if float_value == int_value:
                    value = str(int_value)
                else:
                    value = str(float_value)
            except ValueError:
                pass
            if col_idx >= 3:
                if row[2] == 'GRCh38 coordinates' or row[2] == 'Distance from repeat region' or (row[1] == 'Allele' and not value.isdigit()):
                    cell_format = vertical_format
                elif value == 'A':
                    cell_format = format_A
                elif value == 'G':
                    cell_format = format_G
                elif value == 'C':
                    cell_format = format_C
                elif value == 'T':
                    cell_format = format_T
                else:
                    cell_format = center_format
            else:
                cell_format = default_format  
            worksheet2.write(row_idx, col_idx, value, cell_format)

            if col_idx >= 3:
                column_widths[col_idx] = 2.5 # nastavenie sirky stlpcov na velmi uzku, pre prehladnost sekvencii
            else:
                column_widths[col_idx] = max(column_widths.get(col_idx, 0), len(value))

# nastavenie sirky stlpcov
for col_idx, width in column_widths.items():
    worksheet2.set_column(col_idx, col_idx, width)

workbook.close()