import xml.etree.ElementTree as ET
from Bio import SeqIO
import re
import gzip
import os

def get_alias(identifier):
    parts = identifier.split("/")
    assert len(parts) == 5 and parts[0] == "SARS-CoV-2" and parts[1] == "human"
    return(parts[3])

# Define variables
results = "/Users/taavi/Projects/sars-cov2-et-results/results/"
fasta_path = os.path.join(results, "consensus_masked.fa")
xml_path = "results/PRJEB42961_samples.xml"

ASSEMBLY_TYPE = "COVID-19 outbreak"
PROGRAM = "doi:10.5281/zenodo.4542000"
PLATFORM = "Illumina"

# Parse ENA xml
tree = ET.parse(xml_path)
root = tree.getroot()

study = {}
for x in root.iter('SAMPLE'):
    SUBMITTER_ID = x.find('IDENTIFIERS/SUBMITTER_ID').text
    ENA_STUDY = x.findall('SAMPLE_LINKS/SAMPLE_LINK/XREF_LINK/[DB="ENA-STUDY"]')[0][1].text
    ENA_RUN = x.findall('SAMPLE_LINKS/SAMPLE_LINK/XREF_LINK/[DB="ENA-RUN"]')[0][1].text
    PRIMARY_ID = x.find('IDENTIFIERS/PRIMARY_ID').text
    with open(os.path.join(results, SUBMITTER_ID, "covstats.txt"), "r") as covstats:
        COVERAGE = round(float(covstats.readlines()[1].split()[1]))
    study.update({SUBMITTER_ID : {"STUDY": ENA_STUDY, "SAMPLE": PRIMARY_ID, "RUN_REF": ENA_RUN,"ASSEMBLY_TYPE": ASSEMBLY_TYPE, "PROGRAM": PROGRAM, "PLATFORM": PLATFORM, "COVERAGE": COVERAGE}})

# Parse sequences
seq_dict = SeqIO.index(fasta_path, "fasta", key_function=get_alias)
fasta = {}
for sample in list(study.keys()):
        record = seq_dict[sample]
        splits = record.id.split("/")
        record.id = f"{splits[3]}_{splits[2]}_{splits[4]}"
        record.description = ""
        record.seq = record.seq.strip("N")
        FASTA = f"output/{sample}.fasta.gz"
        CHR_LIST = f"output/{sample}.chrom.gz"
        with gzip.open(FASTA, "wt") as h:
            h.write(record.format("fasta"))
        with gzip.open(CHR_LIST, "wt") as h:
            h.write(f"{record.id}   1   Chromosome")
        fasta.update({sample: {"FASTA": os.path.basename(FASTA), "NAME": record.id, "CHROMOSOME_LIST": os.path.basename(CHR_LIST)}})

manifest = {}
for key in set().union(study, fasta):
    if key in study: manifest.setdefault(key, {}).update(study[key])
    if key in fasta: manifest.setdefault(key, {}).update(fasta[key])

for s, m in manifest.items():
    with open(f"output/{s}_manifest.txt", "w") as h:
        for k, v in m.items():
            h.write(f"{k}   {v}\n")
 