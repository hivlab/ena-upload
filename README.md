[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4574997.svg)](https://doi.org/10.5281/zenodo.4574997)

# ena-upload

Snakemake workflow to trim and upload Sars-Cov-2 reads to ENA. This workflow uses [ena-upload-cli](https://github.com/usegalaxy-eu/ena-upload-cli) for upload and updated *process_xlsx.py* script for metadata table preparation from [Galaxys ena_upload tool](https://github.com/galaxyproject/tools-iuc/tree/master/tools/ena_upload). Raw reads are adapter-trimmed using [bbduk.sh](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) tool and mapped to host/human reference genome using [bbmap.sh](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/). 

## Install

- Download and install miniconda3: <https://docs.conda.io/en/latest/miniconda.html>.
- Clone this repo and create conda environment

```bash
git clone https://github.com/avilab/ena-upload.git
cd ena-upload
conda env create -f workflow/envs/environment.yml
```

- Install R **pepr** separately, as *r-pepr* is not available on anaconda.org.

```bash
conda activate ena_upload
Rscript -e "install.packages('pepr', repos='https://cloud.r-project.org/')"
```

## Running

Update files in `config` directory:

- `metadata.csv` sample metadata (sample_name, collection_date, experiment id, runs/raw read file names)

- `pep.yaml` **peppy** configuration file to provide full paths to raw reads.

- `config.yaml` study wide configuration.

- Check upload rule flags in `Snakefile`, by default, TEST submission (--dev flag in extra parameter) with 'add' action is performed.

- Create `.secret.yml` file with [Webin](https://www.ebi.ac.uk/ena/submit/sra/#home) credentials.

- Test run:

```bash
snakemake -n --use-conda
```

- Run:

```bash
snakemake -j --use-conda
```

- Hold your breath.

## Workflow graph

```bash
snakemake --dag -d .tests/ | dot -Tsvg > resources/rulegraph.svg
```

![rulegraph](resources/rulegraph.svg)
