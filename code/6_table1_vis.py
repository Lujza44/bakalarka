import pandas as pd
import xlsxwriter

data = pd.read_csv('data/raw_vis.csv')

workbook = xlsxwriter.Workbook('data/Table_S1.xlsx', {'nan_inf_to_errors': True})
worksheet = workbook.add_worksheet()

red_format = workbook.add_format({'color': 'red', 'underline': 1})
default_format = workbook.add_format({'color': 'black'})

# zapisanie heads
write_col_index = 0
for col_num, value in enumerate(data.columns):
    if col_num in [7, 9]:
        continue  # preskocanie 7. a 9. stlpca, ktore boli len pomocne
    worksheet.write(0, write_col_index, value)
    write_col_index += 1

# funkcia, ktora farbi SNPs v sekvenciach
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

# citanie dat z CSV a zapisovanie sformatovanych dat do vyslednej XLSX tabulky
for row_index, row in data.iterrows():
    write_col_index = 0
    for col_index, value in enumerate(row):
        if pd.isna(value):
            value = ''
        if col_index == 6 and pd.notna(row['5\' SNP indexes']):
            snps_5_indexes = [int(num.strip()) for num in row['5\' SNP indexes'].split(",")]
            worksheet.write_rich_string(row_index + 1, write_col_index, *color_red(value, snps_5_indexes))
        elif col_index == 8 and pd.notna(row['3\' SNP indexes']):
            snps_3_indexes = [int(num.strip()) for num in row['3\' SNP indexes'].split(",")]
            worksheet.write_rich_string(row_index + 1, write_col_index, *color_red(value, snps_3_indexes))
        elif col_index in [7, 9]:
            continue
        else:
            worksheet.write(row_index + 1, write_col_index, value)
        if col_index not in [7, 9]:
            write_col_index += 1

workbook.close()