import argparse
import pathlib
import sys
import pandas as pd


FILE_FORMAT = "fastq"


def read_ena_xlsx_sheet(xlsx_path, sheet_name):
    file_path = pathlib.Path(xlsx_path)
    file_extension = file_path.suffix.lower()[1:]
    if file_extension == "xlsx":
        engine = "openpyxl"
    else:
        engine = None
    df = (
        pd.read_excel(
            xlsx_path,
            sheet_name=sheet_name,
            engine=engine,
            converters={"collection date": str},
        )
        .dropna(axis=0, how="all")
        .dropna(axis=1, how="all")
    )
    assert not df.empty, f"Sheet '{sheet_name}' is empty in {xlsx_path}"
    return df


def extract_data(xl_sheet, expected_columns, unique_key="alias"):
    if any(xl_sheet.columns.value_counts() > 1):
        sys.exit("Duplicated columns")
    for col in expected_columns:
        assert col in xl_sheet.columns, f"Expected metadata column {col} not found"
    assert not any(
        xl_sheet.duplicated(subset=[unique_key])
    ), f"{unique_key} identificators not unique"
    xl_sheet = xl_sheet.set_index(unique_key)
    return xl_sheet.to_dict("index")


def main(xlsx_path, out_path, action, viral_submission=False):

    # PARSE STUDIES
    #################
    xl_sheet = read_ena_xlsx_sheet(xlsx_path, sheet_name="ENA_study")
    if xl_sheet.shape[0] < 1:
        raise ValueError("No entries found in studies sheet")
    studies_col = ["alias", "title", "study_type", "study_abstract"]
    try:
        studies_dict = extract_data(xl_sheet, studies_col)
    except AssertionError as e:
        print("Sheet ENA_study: ", e)
        raise

    # PARSE SAMPLES
    #################
    xl_sheet = read_ena_xlsx_sheet(xlsx_path, sheet_name="ENA_sample")
    if xl_sheet.shape[0] < 1:
        raise ValueError("No entries found in samples")
    if viral_submission:
        samples_cols = [
            "alias",
            "title",
            "scientific_name",
            "sample_description",
            "geographic location (country and/or sea)",
            "host common name",
            "host health state",
            "host sex",
            "host scientific name",
            "collector name",
            "collection date",
            "collecting institution",
            "isolate",
        ]
    else:
        samples_cols = ["alias", "title", "scientific_name", "sample_description"]
    try:
        samples_dict = extract_data(xl_sheet, samples_cols)
    except AssertionError as e:
        print("Sheet ENA_sample: ", e)
        raise

    # PARSE EXPERIMENTS
    #################
    xl_sheet = read_ena_xlsx_sheet(xlsx_path, sheet_name="ENA_experiment")
    if xl_sheet.shape[0] < 1:
        raise ValueError("No experiments found in experiments sheet")
    exp_columns = [
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
    try:
        experiments_dict = extract_data(xl_sheet, exp_columns)
    except AssertionError as e:
        print("Sheet ENA_experiment: ", e)
        raise

    # PARSE RUNS SHEET
    #################
    xl_sheet = read_ena_xlsx_sheet(xlsx_path, sheet_name="ENA_run")
    if xl_sheet.shape[0] < 1:
        raise ValueError("No entries found in runs sheet")
    run_cols = ["alias", "experiment_alias", "file_name", "file_format"]
    try:
        runs_dict = extract_data(xl_sheet, run_cols, unique_key="file_name")
    except AssertionError as e:
        print("Sheet ENA_run: ", e)
        raise

    # DROP COMMENTS
    ###############
    studies_dict = {
        k: v
        for k, v in studies_dict.items()
        if k in set([v["study_alias"] for k, v in experiments_dict.items()])
    }
    assert bool(studies_dict), "No entries found in studies"
    experiments_dict = {
        k: v
        for k, v in experiments_dict.items()
        if v["study_alias"] in studies_dict.keys()
    }
    assert bool(experiments_dict), "No entries found in experiments"
    samples_dict = {
        k: v
        for k, v in samples_dict.items()
        if k in set([v["sample_alias"] for k, v in experiments_dict.items()])
    }
    assert bool(samples_dict), "No entries found in samples"
    runs_dict = {
        k: v
        for k, v in runs_dict.items()
        if v["experiment_alias"] in experiments_dict.keys()
    }
    assert bool(runs_dict), "No entries found in runs"

    # WRITE HEADERS TO TABLES
    studies_table = open(pathlib.Path(out_path) / "studies.tsv", "w")
    studies_table.write(
        "\t".join(
            [
                "alias",
                "status",
                "accession",
                "title",
                "study_type",
                "study_abstract",
                "pubmed_id",
                "submission_date",
            ]
        )
        + "\n"
    )
    samples_table = open(pathlib.Path(out_path) / "samples.tsv", "w")
    if viral_submission:
        optional_keys = list(samples_dict.values())[0]
        samples_table.write(
            "\t".join(
                [
                    i
                    for i in [
                        "alias",
                        "status",
                        "accession",
                        "title",
                        "scientific_name",
                        "taxon_id",
                        "sample_description",
                        "collection_date",
                        "geographic_location",
                        "geographic_location_region"
                        if optional_keys.get(
                            "geographic location (region and locality)"
                        )
                        else None,
                        "geographic_location_longitude"
                        if optional_keys.get("geographic location (longitude)")
                        else None,
                        "geographic_location_latitude"
                        if optional_keys.get("geographic location (latitude)")
                        else None,
                        "host_common_name",
                        "host_subject_id",
                        "host_health_state",
                        "host_sex",
                        "host_age" if optional_keys.get("host age") else None,
                        "host_scientific_name",
                        "collector_name",
                        "collecting_institution",
                        "isolate",
                        "submission_date",
                    ]
                    if i
                ]
            )
            + "\n"
        )
    else:
        samples_table.write(
            "\t".join(
                [
                    "alias",
                    "status",
                    "accession",
                    "title",
                    "scientific_name",
                    "taxon_id",
                    "sample_description",
                    "submission_date",
                ]
            )
            + "\n"
        )

    experiments_table = open(pathlib.Path(out_path) / "experiments.tsv", "w")
    experiments_table.write(
        "\t".join(
            [
                "alias",
                "status",
                "accession",
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
                "submission_date",
            ]
        )
        + "\n"
    )

    runs_table = open(pathlib.Path(out_path) / "runs.tsv", "w")
    runs_table.write(
        "\t".join(
            [
                "alias",
                "status",
                "accession",
                "experiment_alias",
                "file_name",
                "file_format",
                "file_checksum",
                "submission_date",
            ]
        )
        + "\n"
    )
    action = action

    # WRITE  DICTIONARIES TO TABLE FILES

    # ADD A TIMESTAMP TO THE ALIAS? SEEMS LIKE ENA REQUIRES ALL ENTRIES FOR A WEBIN TO HAVE UNIQUE IDS?
    # dt_oobj = datetime.now(tz=None)
    # timestamp = dt_oobj.strftime("%Y%m%d_%H:%M:%S")
    for study_alias, study in studies_dict.items():
        # study_alias = study_alias + '_' + timestamp
        studies_table.write(
            "\t".join(
                [
                    study_alias,
                    action,
                    "ENA_accession",
                    study["title"],
                    study["study_type"],
                    study["study_abstract"],
                    "",
                    "ENA_submission_data",
                ]
            )
            + "\n"
        )  # assuming no pubmed_id
    for sample_alias, sample in samples_dict.items():
        # sample_alias = sample_alias + '_' + timestamp
        if viral_submission:
            if sample["collector name"] == "":
                sample["collector name"] = "unknown"
            samples_table.write(
                "\t".join(
                    [
                        i
                        for i in [
                            sample_alias,
                            action,
                            "ena_accession",
                            sample["title"],
                            sample["scientific_name"],
                            "tax_id_updated_by_ENA",
                            sample["sample_description"],
                            sample["collection date"],
                            sample["geographic location (country and/or sea)"],
                            sample.get("geographic location (region and locality)"),
                            sample.get("geographic location (longitude)"),
                            sample.get("geographic location (latitude)"),
                            sample["host common name"],
                            "host subject id",
                            sample["host health state"],
                            sample["host sex"],
                            sample.get("host age"),
                            sample["host scientific name"],
                            sample["collector name"],
                            sample["collecting institution"],
                            sample["isolate"],
                            "ENA_submission_date",
                        ]
                        if i
                    ]
                )
                + "\n"
            )
        else:
            samples_table.write(
                "\t".join(
                    [
                        sample_alias,
                        action,
                        "ena_accession",
                        sample["title"],
                        sample["scientific_name"],
                        "tax_id_updated_by_ENA",
                        sample["sample_description"],
                    ]
                )
                + "\n"
            )
        for exp_alias, exp in experiments_dict.items():
            # should I check here if any experiment has a study or sample alias that is incorrect?
            # (not listed in the samples or study dict)
            # process the experiments for this sample
            if exp["sample_alias"] == sample_alias:
                if pd.isnull(exp["library_name"]):
                    if exp["sample_alias"] in exp_alias:
                        lib_alias = exp_alias
                    else:
                        lib_alias = exp_alias + "_" + exp["sample_alias"]
                else:
                    lib_alias = exp["library_name"]
                experiments_table.write(
                    "\t".join(
                        [
                            exp_alias,
                            action,
                            "ena_accession",
                            exp["title"],
                            exp["study_alias"],
                            sample_alias,
                            exp["design_description"],
                            lib_alias,
                            exp["library_strategy"],
                            exp["library_source"],
                            exp["library_selection"],
                            exp["library_layout"].lower(),
                            str(int(exp["insert_size"])),
                            exp["library_construction_protocol"],
                            exp["platform"],
                            exp["instrument_model"],
                            "submission_date_ENA",
                        ]
                    )
                    + "\n"
                )
                for file_name, run in runs_dict.items():
                    if run["experiment_alias"] == exp_alias:
                        runs_table.write(
                            "\t".join(
                                [
                                    run["alias"],
                                    action,
                                    "ena_run_accession",
                                    exp_alias,
                                    file_name,
                                    FILE_FORMAT,
                                    "file_checksum",
                                    "submission_date_ENA",
                                ]
                            )
                            + "\n"
                        )
    studies_table.close()
    samples_table.close()
    experiments_table.close()
    runs_table.close()


if __name__ == "__main__":

    # parser = argparse.ArgumentParser()
    # parser.add_argument("--form", dest="xlsx_path", required=True)
    # parser.add_argument("--out_dir", dest="out_path", required=True)
    # parser.add_argument("--action", dest="action", required=True)
    # parser.add_argument(
    #     "--vir", dest="viral_submission", required=False, action="store_true"
    # )
    # args = parser.parse_args()

    main(
        xlsx_path=snakemake.input[0],
        out_path=snakemake.params["out_dir"],
        action=snakemake.params["action"],
        viral_submission=snakemake.params["viral_submission"],
    )
