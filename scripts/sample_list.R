library(readxl)
library(dplyr)
library(readr)
library(glue)
library(lubridate)
samples <- read_excel("data/esimesed_12_üles_Taavile.xlsx", 
           sheet = 1, 
           range = "D3:E14", 
           col_names = c("alias", "collection_date"))
samples %>% 
  mutate_at("collection_date", ymd) %>% 
  write_csv(glue("results/samples_{Sys.Date()}.csv"))
