library(ggplot2)
library(tidyverse)
library(reshape2)
library(lubridate)
library(viridis)
library(hrbrthemes)
library(readxl)


Metadata <-
  read.csv(file = 'GISAID20210104/metadata.tsv', sep = '\t')

nrow(Metadata)

# Preliminary filtering for Human-associated only
Human_associated <- Metadata %>%
  filter(
  host == 'Human', ignore.case = TRUE,
  (length >= 29000),
  (!is.na(parse_date_time(date, orders = "ymd")))
  ) 

Human_associated <- Human_associated %>%
  mutate(
    sample_date = day(date),
    sample_month = month(date),
    sample_year = year(date)
  ) %>% 
  filter(
    !sample_year < 2019
  )

Human_associated <- Human_associated %>%
  mutate(iso_week_no = isoweek(date))

# Metadata grouping by a country and getting a count
# Top 05 contributing countries
GISAID_country_top05 <- Human_associated %>%
  group_by(country, region) %>%
  tally(sort = T) %>%
  head(n = 5)

ggplot(GISAID_country_top05, aes(x = n, y = country, color = region)) +
  geom_point(alpha=1.0, shape=21, stroke = 1.5) +
  geom_text(aes(label=n),hjust=0, vjust=2, show.legend = FALSE, size=03, color = "Black" )+
  xlab("GISAID submissions")+
  theme_ipsum()


# Metadata grouping by a country and getting a count
# Top 06 - 30 contributing countries
GISAID_country_top_05_30 <- Human_associated %>%
  group_by(country, region) %>%
  tally(sort = T) %>%
  head(n = 30) %>% 
  tail(n = -5)

ggplot(GISAID_country_top_05_30, aes(x = n, y = country, color = region)) +
  geom_point(alpha=1.0, shape=21, stroke = 1.5) +
  geom_text(aes(label=n),hjust=0, vjust=2, show.legend = FALSE, size=03, color = "Black" )+
  xlab("GISAID submissions")+
  theme_ipsum()



years <- c(2019, 2020, 2021)


# Top 30 contributing countries; per week distribution in 2019, 2020, 2021


GISAID_country_top30 <- Human_associated %>%
  group_by(country, region) %>%
  tally(sort = T) %>%
  head(n = 30)

Top30countries_weekly <- as.data.frame(Human_associated %>%
                                          filter(country %in% GISAID_country_top30$country) %>% 
                                          group_by(country, region, iso_week_no, sample_year) %>% 
                                          tally(sort = T)) 

Top30countries_weekly <- Top30countries_weekly[order(Top30countries_weekly$iso_week_no),]


for(year in years){
  if (year %in% Top30countries_weekly$sample_year ){
    print(ggplot(Top30countries_weekly %>% filter(sample_year == year), aes(x=iso_week_no, y=factor(country, level = rev(GISAID_country_top30$country)), color= region, size = n, label=n)) +
            geom_point(alpha=1.0, shape=21, stroke=2) +
            geom_text(aes(label=n),hjust=0, vjust=2, show.legend = FALSE, size=03, color = "Black" )+
            scale_size( range = c(1,10), name="Count per week") +
            guides(size=guide_legend())+
            scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A")+
            scale_x_continuous(breaks=seq(1, 54, 2))+
            theme_ipsum()+
            ylab("Countries") +
            xlab(paste("ISO Week Number", year))+
            ggtitle("SARS-CoV-2 Sampling distribution by ISO Week Number in GISAID (Top 30 Contributors)"))
  }
  }


Top30countries_weekly <- Top30countries_weekly %>% 
  rename(
    count = n
  )
Top30countries_weekly <- Top30countries_weekly[order(Top30countries_weekly$country, Top30countries_weekly$iso_week_no),]
write.table(Top30countries_weekly, file='~/GitHub/GISAID_CIHR/GISAID20210104/Top30_weekly_counts.tsv', quote=FALSE, sep='\t', row.names = FALSE)



# Divisions per week distribution in 2019, 2020, 2021
division_interest <- c('USA', 'United Kingdom', 'Canada')

for(div in division_interest){
    division_weekly <- as.data.frame(Human_associated %>%
                                             filter(country == div, ignore.case = TRUE) %>%
                                             filter(!division == div, ignore.case = TRUE) %>%
                                             group_by(country, division, iso_week_no, sample_year) %>% 
                                             tally(sort = T)) 
  
    division_weekly <- division_weekly[order(division_weekly$iso_week_no),]
  for(year in years){
    if (year %in% division_weekly$sample_year ){
      print(ggplot(division_weekly %>% filter(sample_year == year), aes(x=iso_week_no, y=division, size = n, color = country, label=n)) +
              geom_point(alpha=1.0, shape=21, stroke =2, color = 'blue') +
              geom_text(aes(label=n),hjust=0, vjust=2, show.legend = FALSE, size=03, color = "Black" )+
              scale_size( range = c(1,10), name="Count per week") +
              guides(size=guide_legend())+
              scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A")+
              scale_x_continuous(breaks=seq(1, 54, 2))+
              theme_ipsum()+
              ylab("Divisions/States/Provinces") +
              xlab(paste("ISO Week Number", year)) +
              ggtitle(paste("SARS-CoV-2 Sampling distribution by ISO Week Number",year,"within",div)))
      }
    }
    }




# Group of Interest per week distribution in 2019, 2020, 2021

Interest_group <- as.data.frame(Human_associated %>%
  filter(country %in% c('United Kingdom', 'USA', 'Hong Kong', 'Canada', 'Australia', 'China'), ignore.case = TRUE) %>%
  #filter(region == 'North America' ) %>%
  group_by(country, region, iso_week_no, sample_year) %>% 
  tally(sort = T)) 

Interest_group <- Interest_group[order(Interest_group$iso_week_no),]
for(year in years){
  if (year %in% Interest_group$sample_year ){
    print(ggplot(Interest_group %>% filter(sample_year == year), aes(x=iso_week_no, y=country, size = n, color = region, label=n)) +
            geom_point(alpha=0.7, shape=21, stroke =2) +
            geom_text(aes(label=n),hjust=0, vjust=2, show.legend = FALSE, size=03, color = 'black' )+
            scale_size( range = c(1,10), name="Count per week") +
            guides(size=guide_legend())+
            scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A")+
            scale_x_continuous(breaks=seq(1, 54, 2))+
            theme_ipsum()+
            ylab("Countries") +
            xlab(paste("ISO Week Number", year))+
            ggtitle(paste("SARS-CoV-2 Sampling distribution by ISO Week Number",year,"within",div)))
  }
  }


for(reg in unique(Human_associated$region)){
  region_weekly <- as.data.frame(Human_associated %>%
                                     filter(region == reg, ignore.case = TRUE) %>%
                                     #filter(!division == div, ignore.case = TRUE) %>%
                                     group_by(region, country, iso_week_no, sample_year) %>% 
                                     tally(sort = T)) 
  
  region_weekly <- region_weekly[order(region_weekly$iso_week_no),]
  
  for(year in years){
    if (year %in% region_weekly$sample_year ){
      print(ggplot(region_weekly %>% filter(sample_year == year), aes(x=iso_week_no, y=country, size = n, label=n)) +
              geom_point(alpha=0.7, shape=21, stroke =2, color = "blue") +
              geom_text(aes(label=n),hjust=0, vjust=2, show.legend = FALSE, size=03, color = 'black' )+
              scale_size( range = c(1,10), name="Count per week") +
              guides(size=guide_legend())+
              scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A")+
              scale_x_continuous(breaks=seq(1, 54, 2))+
              theme_ipsum()+
              ylab("Countries") +
              xlab(paste("ISO Week Number", year))+
              ggtitle(paste("SARS-CoV-2 Sampling distribution by ISO Week Number",year,"within",reg)))
    }
  }
  }
