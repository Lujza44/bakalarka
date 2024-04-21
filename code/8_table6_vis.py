import csv
import xlsxwriter

workbook = xlsxwriter.Workbook('data/vis/Table_S6.xlsx')
worksheet = workbook.add_worksheet()

format_A = workbook.add_format({'bg_color': '#92d050', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
format_G = workbook.add_format({'bg_color': '#ffff00', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
format_C = workbook.add_format({'bg_color': '#adb9ca', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
format_T = workbook.add_format({'bg_color': '#ffc7ce', 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
vertical_format = workbook.add_format({'rotation': 90, 'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'})
center_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'font_name': 'Courier New'}) 
default_format = workbook.add_format({'align': 'center', 'valign': 'vcenter', 'color': 'black', 'font_name': 'Courier New'})

with open('data/vis/raw_vis6.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    column_widths = {}  
    for row_idx, row in enumerate(reader):
        if row[1] == 'Allele' or row[2] == 'GRCh38 coordinates':
            worksheet.set_row(row_idx, 100)
        if row[2] == 'Distance from repeat region':
             worksheet.set_row(row_idx, 35)
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
            worksheet.write(row_idx, col_idx, value, cell_format)

            if col_idx >= 3:
                column_widths[col_idx] = 2.5
            else:
                column_widths[col_idx] = max(column_widths.get(col_idx, 0), len(value))

for col_idx, width in column_widths.items():
    worksheet.set_column(col_idx, col_idx, width)

workbook.close()