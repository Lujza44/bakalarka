import csv
import xlsxwriter

# Create a new Excel workbook and add a worksheet.
workbook = xlsxwriter.Workbook('data/vis/Table_S6.xlsx')
worksheet = workbook.add_worksheet()

# Define cell formats for different nucleotides and vertical text.
format_A = workbook.add_format({'bg_color': '#92d050', 'align': 'center'})
format_G = workbook.add_format({'bg_color': '#ffff00', 'align': 'center'})
format_C = workbook.add_format({'bg_color': '#adb9ca', 'align': 'center'})
format_T = workbook.add_format({'bg_color': '#ffc7ce', 'align': 'center'})
vertical_format = workbook.add_format({'rotation': 90, 'align': 'center'})
center_format = workbook.add_format({'align': 'center'})  # Generic center format for other cells

# Open the CSV file for reading.
with open('data/vis/raw_vis6.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    column_widths = {}  # Dictionary to store the maximum width for each column.

    # Iterate through each row in the CSV file.
    for row_idx, row in enumerate(reader):
        # Check for 'Allele' at index 1 or 'GRCh38 coordinates' at index 2
        is_special_row = (len(row) > 1 and row[1] == 'Allele') or (len(row) > 2 and row[2] == 'GRCh38 coordinates')
        is_allele_row = (len(row) > 1 and row[1] == 'Allele')
        vertical_grch38 = (len(row) > 2 and row[2] == 'GRCh38 coordinates')

        # Set row height for special rows
        if is_special_row:
            worksheet.set_row(row_idx, 100)  # Set height to 150 pixels

        # Iterate through each cell in the row.
        for col_idx, value in enumerate(row):
            # Apply vertical or center format based on conditions.
            if col_idx >= 3:
                if vertical_grch38 and col_idx > 2:
                    cell_format = vertical_format  # Use vertical format if it meets the condition
                elif value == 'A':
                    cell_format = format_A
                elif value == 'G':
                    cell_format = format_G
                elif value == 'C':
                    cell_format = format_C
                elif value == 'T':
                    cell_format = format_T
                else:
                    cell_format = center_format  # Use generic center format for non-nucleotide values
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
                    cell_format = None  # Default format for other cells

            # Write the cell to the worksheet, applying the correct format.
            worksheet.write(row_idx, col_idx, value, cell_format)

            # Set fixed width for columns from index 3 onwards or calculate width for other columns.
            if col_idx >= 3:
                column_widths[col_idx] = 3  # Set fixed width for columns from index 3 onwards.
            else:
                column_widths[col_idx] = max(column_widths.get(col_idx, 0), len(value))

# Set the column widths based on the maximum observed width for each column or a fixed width of 40.
for col_idx, width in column_widths.items():
    worksheet.set_column(col_idx, col_idx, width)

# Close the workbook.
workbook.close()
