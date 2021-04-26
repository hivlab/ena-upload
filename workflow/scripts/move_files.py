import pandas as pd
import peppy
import os
import shutil

pep = peppy.Project(".tests/config/pep.yaml")
sample_metadata = pep.sample_table.rename(columns={"sample_name": "alias"})

grouped = sample_metadata.groupby("experiment")

for name, group in grouped:
    dirname = f".tests/results/{name}/data"
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    for file in group.file:
        shutil.move(
            os.path.join(
                ".tests/results/data", file.replace("001.fastq.gz", "clean.fastq.gz")
            ),
            dirname,
        )
