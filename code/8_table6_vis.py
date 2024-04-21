import csv
import xlsxwriter

workbook = xlsxwriter.Workbook('data/vis/Table_S6.xlsx')
worksheet = workbook.add_worksheet()

format_A = workbook.add_format({'bg_color': '#92d050', 'align': 'center', 'font_name': 'Courier New'})
format_G = workbook.add_format({'bg_color': '#ffff00', 'align': 'center', 'font_name': 'Courier New'})
format_C = workbook.add_format({'bg_color': '#adb9ca', 'align': 'center', 'font_name': 'Courier New'})
format_T = workbook.add_format({'bg_color': '#ffc7ce', 'align': 'center', 'font_name': 'Courier New'})
vertical_format = workbook.add_format({'rotation': 90, 'align': 'center', 'font_name': 'Courier New'})
center_format = workbook.add_format({'align': 'center', 'font_name': 'Courier New'}) 
default_format = workbook.add_format({'color': 'black', 'font_name': 'Courier New'})

with open('data/vis/raw_vis6.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    column_widths = {}  
    for row_idx, row in enumerate(reader):
        is_special_row = (len(row) > 1 and row[1] == 'Allele') or (len(row) > 2 and row[2] == 'GRCh38 coordinates')
        is_allele_row = (len(row) > 1 and row[1] == 'Allele')
        vertical_grch38 = (len(row) > 2 and row[2] == 'GRCh38 coordinates')

        if is_special_row:
            worksheet.set_row(row_idx, 100)  

        for col_idx, value in enumerate(row):
            if col_idx >= 3:
                if vertical_grch38 and col_idx > 2:
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
                if value == 'A':
                    cell_format = format_A
                elif value == 'G':
                    cell_format = format_G
                elif value == 'C':
                    cell_format = format_C
                elif value == 'T':
                    cell_format = format_T
                else:
                    cell_format = default_format  
            worksheet.write(row_idx, col_idx, value, cell_format)

            if col_idx >= 3:
                column_widths[col_idx] = 3 
            else:
                column_widths[col_idx] = max(column_widths.get(col_idx, 0), len(value))

for col_idx, width in column_widths.items():
    worksheet.set_column(col_idx, col_idx, width)

workbook.close()