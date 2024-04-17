import xlsxwriter

# DNA sequence
sequence = "AAAAGGCTACGTT"

# Create a workbook and add a worksheet
workbook = xlsxwriter.Workbook('data/Table_S6.xlsx')
worksheet = workbook.add_worksheet()

# Define formats for each nucleotide
format_A = workbook.add_format({'bg_color': '#a8fc83'})
format_G = workbook.add_format({'bg_color': '#f6fc83'})
format_C = workbook.add_format({'bg_color': '#93aac7'})
format_T = workbook.add_format({'bg_color': '#fc83f6'})

# Map nucleotides to formats
formats = {'A': format_A, 'G': format_G, 'C': format_C, 'T': format_T}

column_width = 3
worksheet.set_column(0, len(sequence) - 1, column_width)

# Write sequence to worksheet, letter by letter
for index, letter in enumerate(sequence):
    worksheet.write(0, index, letter, formats[letter]) # vzdy zapise do 0. riadku, index = column

# Close the workbook to create the Excel file
workbook.close()