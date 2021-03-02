# ena-upload

Snakemake workflow to trim and upload Sars-Cov-2 reads to ENA.

## Install

- Download and install miniconda3: <https://docs.conda.io/en/latest/miniconda.html>.
- Clone this repo and create conda environment

```bash
git clone https://github.com/avilab/ena-upload.git
cd ena-upload
conda env create -f workflow/envs/environment.yml
```

- Install R **pepr**

```bash
conda activate ena_upload
Rscript -e 'Sys.unsetenv("GITHUB_PAT");remotes::install_github("pepkit/pepr", dependencies = FALSE)'
```

## Running

Update files in `config` directory:

- `metadata.csv` sample metadata (sample_name, collection_date, experiment id)

- `runs.csv` sample_name mapped to run name(s) and to file names of raw reads.

- `pep.yaml` **peppy** configuration file to provide full paths to raw reads.

- `config.yaml` study wide configuration.

- Check upload rule flags in `Snakefile`, by default, TEST submission with 'add' action is performed.

- Hold your breath.

## Workflow graph

```bash
snakemake --dag -d .tests/ | dot -Tsvg > resources/rulegraph.svg
```

![rulegraph](resources/rulegraph.svg)
