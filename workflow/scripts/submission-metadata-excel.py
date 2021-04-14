import pandas as pd
import peppy
import argparse
import yaml
import re
import os


def mergeDict(dict1, dict2):
    """ Merge dictionaries and keep values of common keys in list"""
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
            dict3[key] = [value, dict1[key]]
    return dict3


class ChainedAssignment:
    def __init__(self, chained=None):
        acceptable = [None, "warn", "raise"]
        assert chained in acceptable, "chained must be in " + str(acceptable)
        self.swcw = chained

    def __enter__(self):
        self.saved_swcw = pd.options.mode.chained_assignment
        pd.options.mode.chained_assignment = self.swcw
        return self

    def __exit__(self, *args):
        pd.options.mode.chained_assignment = self.saved_swcw


pep = peppy.Project(snakemake.input.pep)
sample_metadata = pep.sample_table.rename(columns={"sample_name": "alias"})

with open(snakemake.input.config) as h:
    conf = yaml.full_load(h)

file_name_template = snakemake.params[0]

# ENA_study
ena_study_cols = ["alias", "title", "study_type", "study_abstract"]
ena_study_comments = [
    "Unique identificator for a study. This is used to link experiments to the study.",
    "Title of the study as would be used in a publication.",
    "The STUDY_TYPE presents a controlled vocabulary for expressing the overall purpose of the study.",
    "Briefly describes the goals, purpose, and scope of the Study.  This need not be listed if it can be inherited from a referenced publication.",
]

study_head = dict(zip(ena_study_cols, ena_study_comments))
ena_study = conf["ena_study"]
study = pd.DataFrame(mergeDict(ena_study, study_head))[ena_study_cols]

# ENA_sample
ena_sample_cols = [
    "alias",
    "title",
    "scientific_name",
    "sample_description",
    "collection date",
    "geographic location (country and/or sea)",
    "host common name",
    "host health state",
    "host sex",
    "host scientific name",
    "collector name",
    "collecting institution",
    "isolate",
]
ena_sample_comments = [
    "Unique identificator for each sample. This is used to link experiments to the samples.",
    "Short text that can be used to call out sample records in search results or in displays.",
    "Scientific name of sample that distinguishes its taxonomy. Please use a name or synonym that is tracked in the INSDC Taxonomy database. Also, this field can be used to confirm the TAXON_ID setting.",
    "Free-form text describing the sample, its origin, and its method of isolation.",
    "",
    "The geographical origin of the sample as defined by the country or sea. Country or sea names should be chosen from the INSDC country list (http://insdc.org/country.html).",
    "common name of the host, e.g. human",
    "Health status of the host at the time of sample collection.",
    "Gender or sex of the host.",
    "Scientific name of the natural (as opposed to laboratory) host to the organism from which sample was obtained.",
    "Name of the person who collected the specimen. Example: John Smith",
    "Name of the institution to which the person collecting the specimen belongs. Format: Institute Name, Institute Address",
    "individual isolate from which the sample was obtained",
]
sample_head = dict(zip(ena_sample_cols, [[i] for i in ena_sample_comments]))
sample_conf = conf["ena_sample"]
ena_sample = sample_metadata[["alias", "collection_date"]]
ena_sample = ena_sample.drop_duplicates()

with ChainedAssignment():
    for k, v in sample_conf.items():
        ena_sample.loc[:, k] = v

isolate = ena_sample.apply(
    lambda row: f"SARS-CoV-2/human/{row['geographic location (country and/or sea)']}/{row.alias}/{row.collection_date[0:4]}",
    axis=1,
)
ena_sample = ena_sample.assign(isolate=isolate.values)
ena_sample.rename(columns={"collection_date": "collection date"}, inplace=True)
samples = pd.concat([pd.DataFrame(sample_head), ena_sample], sort=False)[
    ena_sample_cols
]

# ENA_experiment
ena_experiment_cols = [
    "alias",
    "title",
    "study_alias",
    "sample_alias",
    "design_description",
    "library_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "insert_size",
    "library_construction_protocol",
    "platform",
    "instrument_model",
]
ena_experiment_comments = [
    "Unique identificator for each experiment. This is used to link runs to experiments.",
    "Short text that can be used to call out experiment records in searches or in displays. This element is technically optional but should be used for all new records.",
    "from study_metadata",
    "from sample_metadata",
    "Goal and setup of the individual library including library was constructed.",
    "The submitter's name for this library.",
    "Sequencing technique intended for this library.",
    "The LIBRARY_SOURCE specifies the type of source material that is being sequenced.",
    "Method used to enrich the target in the sequence library preparation",
    "LIBRARY_LAYOUT specifies whether to expect single, paired, or other configuration of reads. In the case of paired reads, information about the relative distance and orientation is specified.",
    "Insert size for paired reads",
    "Free form text describing the protocol by which the sequencing library was constructed.",
    "The PLATFORM record selects which sequencing platform and platform-specific runtime parameters. This will be determined by the Center. Nor required if 'instrument_model' is provided.",
    "Model of the sequencing instrument.",
]

experiment_head = dict(zip(ena_experiment_cols, [[i] for i in ena_experiment_comments]))
experiment_conf = conf["ena_experiment"]
ena_experiment = sample_metadata[["alias", "experiment"]]
ena_experiment = ena_experiment.drop_duplicates()
ena_experiment.rename(columns={"alias": "sample_alias"}, inplace=True)
alias = ena_experiment.apply(lambda row: f"{row.experiment}_{row.sample_alias}", axis=1)
ena_experiment = ena_experiment.assign(alias=alias.values)
ena_experiment["study_alias"] = conf["ena_study"]["alias"]
with ChainedAssignment():
    for k, v in experiment_conf.items():
        ena_experiment.loc[:, k] = v
experiments = pd.concat([pd.DataFrame(experiment_head), ena_experiment], sort=False)[
    ena_experiment_cols
]

# ENA_run
ena_run_cols = ["alias", "experiment_alias", "file_name", "file_format"]
ena_run_comments = [
    "Unique identificator for each run.",
    "from experiment_metadata",
    "The name or relative pathname of a run data file.",
    "The run data file model.",
]

run_head = dict(zip(ena_run_cols, [[i] for i in ena_run_comments]))
run_conf = conf["ena_run"]
ena_run = sample_metadata[["run", "file", "experiment", "alias"]]
ena_run.rename(columns={"alias": "sample_alias"}, inplace=True)
pair = ena_run.apply(lambda row: re.search("(?<=R)[1,2]", row.file).group(0), axis=1)
ena_run = ena_run.assign(pair=pair.values)
fn = lambda row: os.path.basename(file_name_template.format(run=row.run, pair=row.pair))
file_name = ena_run.apply(fn, axis=1)
ena_run = ena_run.assign(file_name=file_name.values)
experiment_alias = ena_run.apply(
    lambda row: f"{row.experiment}_{row.sample_alias}", axis=1
)
ena_run = ena_run.assign(experiment_alias=experiment_alias.values)
with ChainedAssignment():
    for k, v in run_conf.items():
        ena_run.loc[:, k] = v
ena_run.rename(columns={"run": "alias"}, inplace=True)
runs = pd.concat([pd.DataFrame(run_head), ena_run], sort=False)[ena_run_cols]

with pd.ExcelWriter(snakemake.output.xlsx) as writer:
    study.to_excel(writer, sheet_name="ENA_study", index=False)
    samples.to_excel(writer, sheet_name="ENA_sample", index=False)
    experiments.to_excel(writer, sheet_name="ENA_experiment", index=False)
    runs.to_excel(writer, sheet_name="ENA_run", index=False)
