#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(lubridate)
library(glue)
library(writexl)

samples <- read_csv(args[1], col_types = list(collection_date = col_character()))
batches <- read_lines(args[2]) %>% 
  str_split("/") %>% 
  map(~.x[6:7]) %>% 
  do.call(rbind, .) %>% 
  as_tibble(.name_repair = "minimal") %>% 
  set_names(c("library_name", "file_name")) %>% 
  mutate(
    alias = str_extract(file_name, "^[A-Za-z0-9-]+"),
    alias = str_replace(alias, "-[^V]+", "")
  )

batches$path <- read_lines(args[2])

sample_metadata <- samples %>% 
  left_join(batches)

sample_metadata$path %>% 
  paste(collapse = " ") %>% 
  write_lines("results/data.txt")

standard_sheet <- function(vars, comments) {
  names(comments) <- vars
  return(as_tibble(rbind(comments)))
}

# ENA_study
ena_study_cols <- c(
  "alias",	
  "title",	
  "study_type",	
  "study_abstract"
  )
ena_study_comments <- c(
  "Unique identificator for a study. This is used to link experiments to the study.",
  "Title of the study as would be used in a publication.",
  "The STUDY_TYPE presents a controlled vocabulary for expressing the overall purpose of the study.",	
  "Briefly describes the goals, purpose, and scope of the Study.  This need not be listed if it can be inherited from a referenced publication."
)
study_head <- standard_sheet(ena_study_cols, ena_study_comments)

study_alias <- "KoroGeno-EST"
study_title <- "Whole genome sequencing of SARS-CoV-2 from Covid-19 patients from Estonia"
study_type <- "Whole Genome Sequencing"
study_abstract <- "Whole genome sequences of SARS-CoV-2 from naso-/oropharyngeal swabs obtained from Estonian Covid-19 patients."
ena_study <- tibble(
  alias = study_alias,
  title = study_title,
  study_type,
  study_abstract
)

ena_study <- bind_rows(
  study_head,
  ena_study
)

# ENA_sample
ena_sample_cols <- c(
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
  "isolate")
ena_sample_comments <- c(
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
  "individual isolate from which the sample was obtained")
sample_head <- standard_sheet(ena_sample_cols, ena_sample_comments)
ena_sample <- sample_metadata %>% 
  select(alias, `collection date` = collection_date)
ena_sample$title <- "PCR tiled amplicon WGS of SARS-CoV-2"
ena_sample$scientific_name	<- "Severe acute respiratory syndrome coronavirus 2"
ena_sample$sample_description	<- "Nasopharyngeal/oropharyngeal swab from Estonian Covid-19 patient."
ena_sample$`geographic location (country and/or sea)` <- "Estonia"
ena_sample$`host common name` <- "human"
ena_sample$`host health state`	<- "unknown"
ena_sample$`host sex` <- "unknown"
ena_sample$`host scientific name` <- "Homo sapiens"
ena_sample$`collector name` <- "unknown"
ena_sample$`collecting institution` <- "unknown"
ena_sample <- ena_sample %>% 
  mutate(isolate = glue("SARS-CoV-2/human/{`geographic location (country and/or sea)`}/{alias}/{year(`collection date`)}"))
ena_sample <- bind_rows(
  sample_head,
  ena_sample
)

# ENA_experiment
ena_experiment_cols <- c(
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
  "instrument_model")
ena_experiment_comments <- c(
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
  "Model of the sequencing instrument.")
experiment_head <- standard_sheet(ena_experiment_cols, ena_experiment_comments)
ena_experiment <- sample_metadata %>% 
  select(sample_alias = alias, library_name) %>% 
  mutate(alias = glue("{library_name}_{sample_alias}"))
ena_experiment$study_alias <- study_alias
ena_experiment$title <- "Illumina MiSeq sequencing of amplicons"
ena_experiment$design_description <- "Tiling amplicon PCR with primers from doi: 10.1093/ve/veaa027"
ena_experiment$library_strategy <- "AMPLICON"
ena_experiment$library_source <- "VIRAL RNA"	
ena_experiment$library_selection <- "RT-PCR"
ena_experiment$library_layout <- "PAIRED"
ena_experiment$insert_size	<- "250"
ena_experiment$library_construction_protocol <- "Illumina Nextera XT"
ena_experiment$platform <- "ILLUMINA"
ena_experiment$instrument_model <- "Illumina MiSeq"
ena_experiment <- ena_experiment %>% 
  select(alias,	title,	study_alias,	sample_alias, design_description, everything())
ena_experiment <- bind_rows(
  experiment_head,
  ena_experiment
)

# ENA_run
ena_run_cols <- c(
  "alias", 
  "experiment_alias",	
  "file_name",	
  "file_format"
)
ena_run_comments <- c(
  "Unique identificator for each run.",
  "from experiment_metadata",
  "The name or relative pathname of a run data file.",
  "The run data file model."
)
run_head <- standard_sheet(ena_run_cols, ena_run_comments)
ena_run <- sample_metadata %>% 
  select(file_name, library_name, sample_alias = alias) %>% 
  mutate(alias = str_remove(file_name, "\\.fastq.*"), 
         experiment_alias = glue("{library_name}_{sample_alias}"),
         file_format = "FASTQ") %>% 
  select(alias, experiment_alias, file_name, file_format)
ena_run <- bind_rows(
  run_head,
  ena_run
)

metadata <- list(
  ENA_study = ena_study,
  ENA_sample = ena_sample,
  ENA_experiment = ena_experiment,
  ENA_run = ena_run
) %>% 
  map(distinct)

write_xlsx(metadata, glue("results/samples_{Sys.Date()}_metadata.xlsx"))
