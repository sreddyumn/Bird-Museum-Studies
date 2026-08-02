# parsing data from Web of Science records with museum specimens AND birds OR avian OR aves
# written by Sushma Reddy, 2026-7-04 with the help of claude.ai

# required packages
install.packages("readr")
install.packages("dplyr")
install.packages("stringr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("purrr")
install.packages("countrycode")
# get world map
install.packages("sf")
install.packages("rnaturalearth")
# get plot pattern
install.packages("ggpattern")

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(purrr)
library(countrycode)
library(sf)
library(rnaturalearth)
library(ggpattern)


#need to remove before uplaod
setwd("~/phylogenyprograms/Rdir/MuseumStudies")

# Import tab-delimited data
# from Web of Science download [limit is 1000 records so had to download as 2 files]
df1 <- read_tsv("20260612_WoS1-1000.txt")
df2 <- read_tsv("20260612_WoS1001-1043.txt")
df<- rbind(df1, df2)

# Inspect date field PY
head(df$PY)
min(df$PY)
max(df$PY)

# create a variable with just publication year
PubsByPY <- data.frame(year = df$PY)

# plot number of publications by year
# Figure 1
ggplot(PubsByPY, aes(x = year)) +
  geom_histogram(binwidth = 2, fill = "#4e79a7", color = "white", alpha = 0.85) +
  theme_bw(base_size = 14) +
  labs(
    #title = "Distribution of Publications by Year",
    x = "Year",
    y = "Number of Publications"
  ) +
  scale_x_continuous(breaks = seq(from = 1930, to = 2026, by = 10)) +
  coord_cartesian(xlim = c(1930, NA)) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text()
  )

ggsave("publications_by_year1.png", width = 3.5, height = 7, units = "in", dpi = 300)

#------
# Extract Address information from each record
# Inspect the relevant column ("C1" =Addresses)
head(df1$C1)

# add another column for record number
df <- df %>% mutate(record_id = row_number())

# Updated Step 1: Split multi-author/institution records into separate rows 
# WoS often packs multiple addresses into one cell, separated by ;. Use separate_rows() to break these out:
# format is [Author1] institution, address, country; [Author2; Author3] institution, address, country

# need to redo for three formats
# multi-author in C1 [Author1] institution, address, country; [Author2; Author3] institution, address, country
# single author in C1 - no brackets
# C1 missing but some author info in RP [need to clean up the ~30 records that have RP but not C1]

#function to find author blocks
extract_blocks <- function(addr, rp = NA) {
  
  # If C1 is empty, fall back to RP field
  if (is.na(addr) || addr == "") {
    if (is.na(rp) || rp == "") return(character(0))
    addr <- rp
  }
  
  if (str_detect(addr, "\\[")) {
    # Bracketed format: [Author(s)] address, country
    str_extract_all(addr, "\\[[^\\]]+\\][^\\[]*")[[1]]
  } else {
    # No brackets or single address: split on ";" directly
    str_split(addr, ";\\s*")[[1]]
  }
}

df_blocks <- df %>%
  mutate(block = map2(C1, RP, extract_blocks)) %>%
  select(record_id, block) %>%
  unnest(block) %>%
  mutate(block = str_trim(block)) %>%
  filter(block != "")

#check if all records were parsed 
#extract all unique record_ids from new blocks
uni<-unique(df_blocks$record_id)
# [1043-1025=] 18 records have empty C1, RP fields; skip for now

#Updated Step 2: for records with multiple-authors per institution, extract authors/address/country, accounting for missing brackets
df_blocks <- df_blocks %>%
  mutate(
    has_authors = str_detect(block, "^\\["),
    authors_str = if_else(
      has_authors,
      str_extract(block, "(?<=\\[)[^\\]]+(?=\\])"),
      NA_character_
    ),
    address = if_else(
      has_authors,
      str_remove(block, "^\\[[^\\]]+\\]\\s*"),
      block
    ),
    address = str_remove(address, ";\\s*$"),
    country = str_trim(str_extract(address, "[^,]+$"))
  )

#Updated Step 3: expand to one row per author (NAs pass through safely)
df_authors <- df_blocks %>%
  separate_rows(authors_str, sep = ";\\s*") %>%
  rename(author = authors_str) %>%
  mutate(author = str_trim(author))

# Fix USA
df_authors <- df_authors %>%
  mutate(
    country = case_when(
      str_detect(country, "\\bUSA\\b") ~ "USA",
      TRUE ~ country
    )
  )

#clean up common quirks
# note: these were discovered by examining data that did not match using the functions above
df_authors <- df_authors %>%
  mutate(
    country = str_remove(country, "\\.$"),
    # USA addresses often end in "State ZIP" rather than "Country" ??? 
    # WoS usually still appends "USA" at the very end, but double-check
    # Standardize common variants
    country = case_when(
      str_detect(country, "^(ENGLAND|SCOTLAND|England|Scotland|Wales|N Ireland)$") ~ "United Kingdom",
      str_detect(country, "^(Peoples R China|Peoples R China.)$") ~ "China",
      str_detect(country, "^(Turkiye)$") ~ "Turkey",
      str_detect(country, "^USA$") ~ "USA",
      str_detect(country, "^WA 98195|CA 92093|CA 94720|HI 96718|AR 72701|CA 92112|CM 96950|FL 32611|GU 96910|HI 96817|HI 96822|MN 55108|MN 55455|VA 22630|VT 05091|USA$") ~ "USA",
      str_detect(country, "^Papua N Guinea$") ~ "Papua New Guinea",
      str_detect(country, "^North Ireland$") ~ "United Kingdom",
      str_detect(country, "^USSR$") ~ "Russia",
      str_detect(country, "^Ussr$") ~ "Russia",
      TRUE ~ country
    )
  )

#check for making sure there is country in each record
which(is.na(df_authors$country))
#if 0 then fine

# convert countries with all caps to title case
df_authors <- df_authors %>%
  mutate(country = tools::toTitleCase(tolower(country)))

#1. Get country counts from your data
# one count per paper per country (not per author)
rcountry_perpaper_counts <- df_authors %>%
  distinct(record_id, country) %>%
  count(country, sort = TRUE, name = "n_records")

# one count per author [better measure of usage]
rcountry_allauthors_counts <- df_authors %>%
  count(country, sort = TRUE, name = "n_records")

#check that only one name per country matched
country_counts <- rcountry_allauthors_counts %>%
  arrange(country)

# match countries to country code
rcountry_allauthors_counts <- rcountry_allauthors_counts %>%
  mutate(
    iso3 = countrycode(country, origin = "country.name", destination = "iso3c")
  )

#total number of countries = 236
#number of countries with less than 10 users [88-42=46]


# Check for unmatched countries
fixcountry<-rcountry_allauthors_counts %>% filter(is.na(iso3))
country_counts <- rcountry_allauthors_counts %>%
  arrange(iso3)

#create plot
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica")
coord_sf(crs = "+proj=robin")

# join your data to the map

map_data <- world %>%
  left_join(rcountry_allauthors_counts, by = c("iso_a3_eh" = "iso3"))

# check that all countries match map_data country codes
rcountry_allauthors_counts %>%
  filter(!iso3 %in% map_data$iso_a3_eh) %>%
  select(country, iso3, n_records)

#plot heatmap

library(ggplot2)

ggplot(map_data) +
  geom_sf(aes(fill = n_records), color = "gray40", linewidth = 0.1) +
  scale_fill_viridis_c(
    option = "cividis",
    #name = "Number of\nrecords",
    na.value = "grey90",
    trans = "log10"  # consider log scale if counts are highly skewed
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )  #+
 # labs(title = "Number of records by country")

ggsave("globalheatmap_perauthor1.png")


# plot without USA
usa_outline <- world %>% filter(iso_a3 == "USA")

ggplot(map_data %>% filter(iso_a3 != "USA")) +
  geom_sf(aes(fill = n_records), color = "gray40", linewidth = 0.1) +
  geom_sf_pattern(data = usa_outline,
                  pattern = "stripe",
                  fill = "white",
                  color = "gray40",
                  linewidth = 0.1,
                  pattern_colour = "gray60",
                  pattern_angle = 45,
                  pattern_density = 0.3,
                  pattern_spacing = 0.02) +
  scale_fill_viridis_c(
    option = "cividis",
    #name = "Number of\nrecords",
    na.value = "grey90",
    trans = "log10",  # consider log scale if counts are highly skewed
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )  #+
# labs(title = "Number of records by country without USA")

ggsave("globalheatmap_perauthor_woUSA1.png")

