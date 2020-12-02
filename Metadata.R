library(ggplot2)
library(tidyverse)
library(reshape2)
library(lubridate)
library(countrycode)


Metadata <- read.csv(file = 'GISAID20201201/metadata.tsv',sep = '\t')
nrow(Metadata)


# Preliminary filtering for Human-associated only
Human_associated <- Metadata %>% 
  filter(host == 'Human', ignore.case = TRUE) %>% 
  filter(length >= 29000) %>% 
  filter(!is.na(parse_date_time(date,orders="ymd")))

Human_associated <- Human_associated %>%
  mutate(sample_date = day(date), sample_month = month(date), sample_year = year(date))

Human_associated <- Human_associated %>%
  mutate(iso_week_no = isoweek(date))

Human_associated$continent <- countrycode(sourcevar = Human_associated[, "country"],
                                  origin = "country.name",
                                  destination = "continent")




# Metadata grouping by a country and getting a count
GISAID_country_top5 <- Human_associated %>% 
  group_by(country, continent) %>% 
  tally(sort = T) %>% 
  head(n=5)

ggplot(GISAID_country_top5) +
  geom_point(aes(x = n, y = country, color = continent)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 08, colour = 'black'),
        axis.text.y = element_text(hjust = 0.5, size = 08, colour = 'black'), 
        axis.line = element_line(colour = "blue"),
        panel.background = element_blank())

GISAID_country <- Human_associated %>% 
  group_by(country, continent) %>% 
  tally(sort = T) %>% 
  tail(n=-5)

ggplot(GISAID_country) +
  geom_point(aes(x = n, y = country, color = continent)) +
  theme_bw()











p <- ggplot(Metadata, aes(x=country)) + 
  geom_bar(aes(color = division ), fill = "white",
           position = "identity")
p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 08, colour = 'black'),
          axis.text.y = element_text(hjust = 0.5, size = 12, colour = 'black'), 
          axis.line = element_line(colour = "blue"),
          panel.background = element_blank())

