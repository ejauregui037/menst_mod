import pandas as pd
import openpyxl

# convert each sheet in the excel file to a csv file
# note: "Cover" csv will be nonsensical
for sheet_name, df in pd.read_excel(r"250407 YWS & CLUE Data For Lauren (1).xlsx", sheet_name=None).items():
    df.to_csv(f'output_{sheet_name}.csv', index=False, encoding='utf-8')
