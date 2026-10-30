## Lower Tuolumne River HS3 Power Analysis
## E. Suenaga 9/9/2026 - script adapted from Kiera McNeely 5/26/2026 [Cramer Fish Sciences]


# Libraries ----
library(effsize) 
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggpubr)
library(sf)
library(lmerTest)
library(tidyverse)
library(stringr) 


########################################## 
## Read in redd data and spatial Project riffle polygon data. 
##########################################
redd_data <- st_read("LTR_Redds_11-25new")
                     

# Subset redd data to only represent the values needed for this analysis: sample year is 2019-2025, only include species_si "Chinook" (omit O. mykiss and spring-run Chinook)
# subset redd data to year>2019, Species_si=="Chinook", ReddVerif=="New" for this analysis
redd_sf <- redd_data %>%
  filter(SampleYear > 2018) %>%
  filter(ReddVerif=="New") %>%
  mutate(GPS_Date = as.factor(GPS_Date)) %>%
  filter(Species_si == "Chinook") 

head(redd_sf)

# Export redd data as shapefile to be loaded into QGIS

st_write(redd_sf, "LTR_redds19-25.shp", delete_layer = TRUE)


## View riffle polygon data shapefile in R or upload direcly to QGIS:
spatial_data <- st_read("LTR_USOLGB_HS3_polygons2")
                        
View(spatial_data)

spatial_data <- as.data.frame(spatial_data) ## turning shapefile into dataframe


########################
### Use QGIS to spatially join (join attributes by location) redd data to the polygon data:


## 1. Load vector data for redds (LTR_Redds_11-25.shp) and polygon data (LTR_USOLGB_HS3_polygons2.shp)
## 2. In Processing Toolbox, use the "Clip" tool to clip redd data to only redds that are within the LTR_USOLGB_HS3_polygons2
## 3. Join the features in "Clipped" layer that intersect/are within by comparing to the polygon layer
## 4. Export the new joined feature as a .csv file (in this case, named: "LTR_USOLGB_HS3redddat.csv")


##############################
### Upload newly joined data and format for analysis

hs_dat <- read.csv("LTR_USOLGB_HS3redddat.csv") ## Load
head(hs_dat)
str(hs_dat)
View(hs_dat)

## Keep rows where SampleYear (esc year redds sampled) and esc.year from polygon data match (the spatial join in QGIS joined every redd across years to each esc.year of the polygons)
hs_dat$SampleYear <- as.integer(hs_dat$SampleYear) ## make sure in the format

hs_dat <- hs_dat %>%
  filter(SampleYear == esc_year)

## format data and calculate redd density and density of superimposed redds (imposed_redd_density) in each sub-site across years
redd_data <- hs_dat %>%
  group_by(treatment, riffle_ID, sub_id, esc_year, Species_si, site, area_m2, method, dam_rkm, tributary) %>%
  summarise(Total_redds = sum(ReddCount), .groups = "drop",
            Total_imposed = sum(ImpVerif == "Yes", na.rm = TRUE)) %>% 
  mutate(redd_density = ifelse(is.na(Total_redds / area_m2), 0, Total_redds / area_m2),
         imposed_redd_density = ifelse(is.na(Total_imposed / area_m2), 0, Total_imposed / area_m2))

View(redd_data) ### Only sub_ids that have redds are joined in QGIS, so we need to include the sub_ids with 0 redds for each year
head(redd_data)

### Join redd_data to the spatial_data dataframe
head(spatial_data)

sp_redd_data <- spatial_data %>%
  full_join(redd_data, by = c("esc_year","riffle_ID", "sub_id")) %>%
  mutate(
    site = coalesce(site.x, site.y),
    treatment= coalesce(treatment.x,treatment.y),
    method= coalesce(method.x, method.y),
    area_m2= coalesce(area_m2.x, area_m2.y),
    dat_rkm= coalesce(dam_rkm.x, dam_rkm.y),
    across(c(Total_redds,Total_imposed, redd_density,imposed_redd_density), ~ ifelse(is.na(.), 0, .))
  ) %>%
  dplyr::select(-site.x, -site.y, -treatment.x,-treatment.y,-area_m2.x, -area_m2.y,-method.x,-method.y,-tributary.x, -tributary.y,-dam_rkm.x,-dam_rkm.y) ## removing duplicate columns


View(sp_redd_data)

# ── Enhanced Riffle BACI Analysis ────────────────────────────────────────────────────────────────────
## Using the LMER 

dat_enhcd <- subset(sp_redd_data, treatment %in% c("Pre-project_enhanced","Post-project_enhanced","Pre-control_enhanced","Post-control_enhanced"))
View(dat_enhcd)

## Create rows for missing sub_ids


## Look at data plot

## formatting plots
font_size <- 15
ARS_theme <- theme_classic() +
  theme(
    text            = element_text(size = font_size),
    axis.text.x     = element_text(size = font_size, color = "black", angle = 45, vjust = 0.8),
    axis.text.y     = element_text(size = font_size, color = "black"),
    panel.grid.major = element_line(color = "grey", linetype = "dotted"),
    legend.text     = element_text(size = font_size),
    legend.title    = element_blank(),
    strip.text      = element_text(size = font_size),
    strip.background = element_blank(),
    panel.border    = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
    legend.position = "none"
  )



## Re-ordering treatment factor for plotting
dat_enhcd$treatment <- factor(dat_enhcd$treatment, 
                              levels = c("Pre-control_enhanced", "Post-control_enhanced",
                                         "Pre-project_enhanced", "Post-project_enhanced"),
                              labels = c("Control Before", "Control After", "Enhanced Project Before", "Enhanced Project After"))


### Plotting means and standard errors by treatment type (project/control) and time (pre-/post-restoration) --- **switch y axis to density of superimposed redds by using (imposed_redd_density)


RD_plot <- dat_enhcd %>%
  group_by(treatment) %>%
  summarise(
    mean_density = mean(redd_density, na.rm = TRUE),
    se_density   = sd(redd_density,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = treatment, y = mean_density, fill = treatment)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_density - se_density,
                    ymax = mean_density + se_density),
                width = 0.25) +
  scale_fill_manual(values = c("grey75", "grey34",
                               "darkolivegreen2", "darkolivegreen4",
                               "dodgerblue2", "dodgerblue4")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Site x Treatment", y = "Redd Density (redds/m²)") +
  theme_bw() + ARS_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)))

RD_plot

### Save plot if desired
ggsave(RD_plot,file ="BACI_LTR.USOLGB_Enhanced_ReddDens.png", width = 8, height = 6, units = "in") # change file name as needed

ggsave(RD_plot,file ="BACI_LTR.USOLGB_Enhanced_ImposedReddDens.png", width = 8, height = 6, units = "in") 

##### RUnning LMER to look at effect of restoration on enhanced riffles
## Look at interaction term between time (pre/post restoration) and site (impact [or project]/control), so need to split the 'treatment' variable up into two columns

dat_enhcd <- dat_enhcd %>%
  mutate(
    time = ifelse(esc_year < 2024, "pre", "post"),
    site_trt = case_when(
      grepl("control", treatment, ignore.case = TRUE) ~ "control",
      grepl("project", treatment, ignore.case = TRUE) ~ "impact",
      TRUE ~ NA_character_
    )
  )
View(dat_enhcd)

## Run LMER and anova to determine if interaction between site treatment and time (pre/post) is significant (p<0.05)
mod <- lmerTest::lmer(redd_density ~ site_trt*time + (1|sub_id), data=dat_enhcd) # Random effect for site controls for site-level differences. type*time interaction tests whether the change over time differs between control and impact sites.
summary(mod) 
anova(mod) # significant site_trt:time interaction (p<0.05) indicates a restoration effect.

## for superimposed redds
mod_imposed <- lmerTest::lmer(imposed_redd_density ~ site_trt*time + (1|sub_id), data=dat_enhcd) # Random effect for site controls for site-level differences. type*time interaction tests whether the change over time differs between control and impact sites.
summary(mod_imposed) 
anova(mod_imposed) # significant site_trt:time interaction (p<0.05) indicates a restoration effect.




# ── Created Riffle BACI Analysis ────────────────────────────────────────────────────────────────────

## Subset redd_data to include the created treatments, and also the enhanced control treatment that will be used a benchmark for visualization purposes
dat_created <- subset(sp_redd_data, treatment %in% c("Pre-project_created","Post-project_created","Pre-control_created","Post-control_created", "Pre-control_enhanced","Post-control_enhanced")) ## only rows with Post-project_created exist, since there were no redds to join in the other categories
View(dat_created)

## Look at data plot

## Re-ordering treatment factor for plotting
dat_created$treatment <- factor(dat_created$treatment, 
                              levels = c("Pre-control_created", "Post-control_created",
                                         "Pre-project_created", "Post-project_created",
                                         "Pre-control_enhanced", "Post-control_enhanced"),
                              labels = c("Control Before", "Control After", "Created Project Before", "Created Project After","Benchmark Before","Benchmark After"))


### Plotting means and standard errors by treatment type (project/control) and time (pre-/post-restoration) --- **switch y axis to density of superimposed redds by using (imposed_redd_density)


RD_plot <- dat_created %>%
  group_by(treatment) %>%
  summarise(
    mean_density = mean(redd_density, na.rm = TRUE),
    se_density   = sd(redd_density,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = treatment, y = mean_density, fill = treatment)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_density - se_density,
                    ymax = mean_density + se_density),
                width = 0.25) +
  scale_fill_manual(values = c("grey75", "grey34",
                               "darkolivegreen2", "darkolivegreen4",
                               "dodgerblue2", "dodgerblue4")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Site x Treatment", y = "Redd Density (redds/m²)") +
  theme_bw() + ARS_theme +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)))

RD_plot

### Save plot if desired
ggsave(RD_plot,file ="BACI_LTR.USOLGB_Created_ReddDens.png", width = 8, height = 6, units = "in")


##### RUnning LMER to look at effect of restoration on enhanced riffles
## Look at interaction term between time (pre/post restoration) and site (impact [or project]/control), so need to split the 'treatment' variable up into two columns

dat_created <- dat_created %>%
  mutate(
    time = ifelse(esc_year < 2024, "pre", "post"),
    site_trt = case_when(
      grepl("control", treatment, ignore.case = TRUE) ~ "control",
      grepl("project", treatment, ignore.case = TRUE) ~ "impact",
      TRUE ~ NA_character_
    )
  )
View(dat_created)

## remove benchmark site for analysis
dat_created_test <- na.omit(dat_created)

## Run LMER and anova to determine if interaction between site treatment and time (pre/post) is significant (p<0.05)
mod <- lmerTest::lmer(redd_density ~ site_trt*time + (1|sub_id), data=dat_created_test) # Random effect for site controls for site-level differences. type*time interaction tests whether the change over time differs between control and impact sites.
summary(mod) 
anova(mod) # significant site_trt:time interaction (p<0.05) indicates a restoration effect.

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------
#########################################
#### POWER ANALYSIS #####################

## Enhanced Site

head(dat_enhcd)

pow.dat <- dat_enhcd[,c(1:3,6:11)]
pow.dat <- pow.dat %>%
  mutate(time = case_when(
    str_detect(treatment, "^Pre")  ~ "Pre",
    str_detect(treatment, "^Post") ~ "Post",
    TRUE ~ NA_character_
  )) %>%
  mutate(type= case_when(
    str_detect(treatment, "project") ~"Impact",
    str_detect(treatment, "control") ~ "Control",
               TRUE ~ NA_character_))
pow.dat


pow.dat <- pow.dat %>%
  mutate(
    binary_type = if_else(type == "Impact", 1, 0),
    binary_time = if_else(time == "Post", 1, 0)
  )

#create site_id for the impact (1-15) and control (16-30) sites sub ids
pow.dat$site_id <-c(seq(1, 15, by = 1), seq(1, 15, by = 1),seq(1, 15, by = 1),seq(1, 15, by = 1),seq(1, 15, by = 1),seq(1, 15, by = 1),seq(1, 15, by = 1), seq(16, 30, by = 1), seq(16, 30, by = 1),seq(16, 30, by = 1),seq(16, 30, by = 1),seq(16, 30, by = 1),seq(16, 30, by = 1),seq(16, 30, by = 1))


pow.dat <- as.data.frame(pow.dat)
View(pow.dat)

## Can only compare two years (one pre, one post)
pow.dat_sub <- pow.dat %>% filter(esc_year %in% c(2025, 2019))

View(pow.dat_sub)

anyDuplicated(names(pow.dat))

# Plots -----
ggplot(pow.dat_sub %>% mutate(time = factor(time, levels = c("Pre", "Post"))))+
  geom_boxplot(aes(x = type, y = redd_density, fill = type))+
  facet_wrap(~time) +
  theme_bw() +
  labs(x = "Site Type", y = "Redd Density", title = "Test Redd Data")
# from this plot we would expect a significant interaction due to the increase in impact from before and after and the consistency of the control from before and after.


# Bootstrap Power Analysis -----
# need to subset by impact sites and control sites, then pull all data from those sites

impact_sites <- unique(pow.dat_sub$site_id[pow.dat_sub$type =="Impact"])  # List of impact sites to pull from for bootstrap
control_sites <- unique(pow.dat_sub$site_id[pow.dat_sub$type =="Control"]) # List of control sites to pull from for bootstrap

n_sites <- c(5,8,10,12,15,20,50,80) # samples to pull from each treatment
boot_itr <- 1000 # number of iterations

boot_baci <- function(impact_sites, control_sites, boot_itr, n_sites, data){ # need list of impact sites, control sites, data to filter from,
  
  # create results table
  total_rows <- length(n_sites)*boot_itr # pre-allocate results table
  boot_res <- data.frame(n = numeric(total_rows),         
                         effect_size = numeric(total_rows),
                         pvalue = numeric(total_rows))
  counter <- 1 # counter for storing results
  
  
  for(n_site in n_sites){# repeat bootstrap for each n value
    for(i in 1:boot_itr){
      
      # resample data: Get the sites to use in the bootstrap, need to pull data from those sites for all years.
      imp <- sample(impact_sites, size = n_site, replace = T) 
      con <- sample(control_sites, size = n_site, replace = T)
      all_sites <- c(imp,con)
      
      mod_data <- NULL
      for(s in all_sites){
        keep <- data[data$site_id == s,] # get the data for this site for all the years to use in the BACI analysis
        mod_data <- rbind(mod_data,keep)
      }
      
      mod_data <- as.data.frame(mod_data)
      # BACI
      mod <- tryCatch(
        lmerTest::lmer(redd_density ~ binary_time*binary_type + (1|site_id), data = mod_data),
        error = function(e) NULL
      )
      
      if (is.null(mod)) {
        p_value <- NA
        eff_size <- NA
      } else {
        mod_sum <- summary(mod)
        coef_names <- rownames(mod_sum$coefficients)
        interaction_term <- "binary_time:binary_type"
        
        if (interaction_term %in% coef_names) {
          p_value <- mod_sum$coefficients[interaction_term, "Pr(>|t|)"]
          eff_size <- ifelse(p_value <= 0.05, mod_sum$coefficients[interaction_term, "Estimate"], NA)
        } else {
          p_value <- NA
          eff_size <- NA
        }
      }
      
      # Save results
      res <- c(n_site, eff_size, p_value)
      boot_res[counter,] <- res
      counter <- counter +1
      
    } #/bootstrap
  } #/sample size
  
  # summary of results for each n (sample size group)
  summary <- boot_res %>% group_by(n) %>% summarize(mean_eff_size = mean(effect_size, na.rm=T),
                                                    mean_power = mean(pvalue < 0.05, na.rm=T)) #Accept Alternative if p-value is less than 0.05.
  return(summary)
}

## Results ----
set.seed(2026) # set before calling the function for reproducibility
baci_power <- boot_baci(impact_sites, control_sites, boot_itr, n_sites, data = pow.dat_sub) ## singular boundary is from sampling with replacement, this is ok.

baci_power
# Effect size is the difference of the difference, which is the effect of the restoration. It's how much the impact group changed related to the control group.

