import pandas as pd
import peppy
import argparse
import yaml

p = peppy.Project(".tests/config/pep.yaml")
sample_metadata = p.sample_table.rename(columns = {"sample_name": "alias"})
print(sample_metadata)

with open(".tests/config/config.yaml") as h:
    conf = yaml.full_load(h)
print(conf)

file_name_template = "results/data/{run}_R{pair}_clean.fastq.gz"
print(file_name_template)

with pd.ExcelWriter('test.xlsx') as writer:
    sample_metadata.to_excel(writer, sheet_name='sample_metadata')
    ##df2.to_excel(writer, sheet_name='Sheet2')
