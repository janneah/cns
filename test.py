import pandas as pd

file = pd.read_table("/home/janneae/cns/all_chr.ldscore")
print(file.head())
print(file.dtypes)
# file["CHR"] = pd.to_numeric(file["CHR"])

# file = file.sort_values('CHR')

# file.to_csv("/home/janneae/cns/all_chr.ldscore")