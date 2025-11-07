
# scrap

# Check split patient overlap between Tacrolimus and Prednisolone

```{r}
df <- tacro_all_matched %>% 
  dplyr::select(Record.ID, Group, patient_class_median) %>% 
  distinct() %>%
  mutate(Immunosuppression = "Tacrolimus") %>%
  full_join(pred_all_matched %>% 
              dplyr::select(Record.ID, Group, patient_class_median) %>% 
              distinct() %>%
              mutate(Immunosuppression = "Prednisolone")
  ) %>%
  na.omit()

d <- df %>%
  group_by(Record.ID) %>%
  dplyr::count(patient_class_median) 

consistent <- d %>% dplyr::count(Record.ID) %>% filter(n == 1) %>% pull(Record.ID) # 31 stay consistent
consistently_high <- d %>% filter(Record.ID %in% consistent ) %>% filter(patient_class_median == "Higher") %>% pull(Record.ID) #16
consistently_low <- d %>% filter(Record.ID %in% consistent ) %>% filter(patient_class_median == "Lower") %>% pull(Record.ID) #15

df <- df %>% 
  mutate(patient_class = ifelse(Record.ID %in% consistently_high, "Higher", 
                                ifelse(Record.ID %in% consistently_low, "Lower", "Variable"))) 

p <- 
  ggplot(df, aes(x = Immunosuppression, y = Record.ID)) +
  geom_point(aes(fill = patient_class_median), shape = 21) +
  geom_path(aes(group = Record.ID)) +
  labs(title = "Patient consistency within first 6 months", fill = "Median level status") +
  theme(axis.text.y = element_blank(), plot.title = element_text(face = "bold")) + #, strip.text.y = element_text(angle = 0)
  facet_grid(Group+patient_class~., scales = "free", space = "free")

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/patient-comparison-tacro-pred.jpeg"), p, width = 12, height = 14, units = 'cm', dpi = 600)

df %>% 
  dplyr::select(Record.ID, Group, patient_class) %>% 
  distinct() %>% 
  group_by(Group) %>% 
  dplyr::count(patient_class) %>%
  mutate(percentage = round(100*n/sum(n)) )
```

```{r}
patient_groups <- df %>% 
  dplyr::select(Record.ID, patient_class) %>% 
  distinct()
```



# interpolation
#%>%
#group_by(Record.ID) %>%
#mutate(tacrolimus_interp = round(na.approx(tacrolimus_level, maxgap = 6, rule = 2, na.rm = FALSE), 2 )) %>%
#mutate(tacrolimus_interp = ifelse(Days == min(Days) & is.na(tacrolimus_level), NA, tacrolimus_interp),
#       tacrolimus_interp = ifelse(Days == max(Days) & is.na(tacrolimus_level), NA, tacrolimus_interp)) # remove values for lower and higher end of interpolation 

# in range calculation: is a patient more often in range, higher or lower
in_range_count <- tacro %>%
  group_by(Record.ID, Immunosuppression_interval_unique) %>%
  dplyr::count(within_range) %>%
  mutate(percentage_in_range_interval = n/sum(n)) %>%
  arrange(Record.ID, Immunosuppression_interval_unique, desc(percentage_in_range_interval), factor(within_range, levels = c("In range", "Above range", "Below range")) # tie-break: In range first
  ) %>%
  slice_head(n = 1) %>%
  dplyr::rename(within_range_interval = within_range) %>%
  dplyr::select(!n) 

#in range numeric score: avg_score <0 = mostly below, 0 = balanced/in range, >0 = mostly above.
tacro_score <- tacro %>%
  mutate(score = case_when(
    within_range == "Below range" ~ -1,
    within_range == "In range" ~ 0,
    within_range == "Above range" ~ 1,
    is.na(within_range) ~ NA
  )) %>%
  group_by(Record.ID, Immunosuppression_interval_unique) %>%
  summarise(avg_score_interval = mean(score, na.rm = TRUE)) %>%
  mutate(avg_score_interval_cat = case_when(avg_score_interval < 0 ~ "Below",
                                            avg_score_interval == 0 ~ "In range",
                                            avg_score_interval > 0 ~ "Above",
                                            is.na(avg_score_interval) ~ NA))

# Immune cells in the blood post-transplant

```{r}
immune_cells_blood <- master_metadata$bronch_metadata %>% 
  dplyr::select(Record.ID, Days, Patient_days, WCC, Lymphocytes, Neutrophils) %>% 
  filter(Record.ID %in% rna_patients) %>% 
  left_join(all_rna_metadata[, c("Record.ID", "Group")] %>% distinct()) %>%
  mutate(Months = round(Days/30.417, 1))

immune_cells_blood_average <- immune_cells_blood %>%
  group_by(Group, Record.ID, Months) %>%
  summarise(lymphocytes_average_levels = mean(Lymphocytes)) %>%
  na.omit()


ggplot(immune_cells_blood, aes(x = Months, y = Lymphocytes)) +
  #geom_path(size = 0.3) +
  geom_point(size = 1) +
  geom_smooth(se = FALSE) +
  labs(title = "Lymphocytes post-transplant", y = "Percentage", 
       x = "Time post-transplant (months)") +
  theme +
  theme(legend.position = "bottom")  +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 18, 24)) +
  scale_y_continuous(limits = c(0, 4)) 
```

# within range for interval
d <- tacro_all_matched %>% 
  dplyr::select(Immunosuppression_interval_unique, Record.ID, avg_score_interval_cat, Group) %>% 
  filter(Immunosuppression_interval_unique %in% c("<6 months")) %>%
  distinct() %>%
  group_by(Record.ID) %>%
  dplyr::count(avg_score_interval_cat) 

consistent <- d %>% dplyr::count(Record.ID) %>% filter(n == 1) %>% pull(Record.ID) # 35 stay either higher or lower or in range
consistently_high <- d %>% filter(Record.ID %in% consistent ) %>% filter(avg_score_interval_cat == "Above") %>% pull(Record.ID) #3
consistently_range <- d %>% filter(Record.ID %in% consistent ) %>% filter(avg_score_interval_cat == "In range") %>% pull(Record.ID) 
consistently_low <- d %>% filter(Record.ID %in% consistent ) %>% filter(avg_score_interval_cat == "Below") %>% pull(Record.ID) #32
variable <- d %>% dplyr::count(Record.ID) %>% filter(n > 1) #19 patients switch

tacro_all_matched <- tacro_all_matched %>% 
  mutate(patient_class_range = ifelse(Record.ID %in% consistently_high, "Higher", 
                                      ifelse(Record.ID %in% consistently_range, "In range", 
                                             ifelse(Record.ID %in% consistently_low, "Lower", "Variable"))))

# with high iqr for interval
d <- tacro_all_matched %>% 
  dplyr::select(Immunosuppression_interval_unique, Record.ID, tacro_iqr_3m, Group) %>% 
  filter(Immunosuppression_interval_unique %in% c("<6 months")) %>%
  distinct() %>%
  group_by(Record.ID) %>%
  dplyr::count(avg_score_interval_cat) 

consistent <- d %>% dplyr::count(Record.ID) %>% filter(n == 1) %>% pull(Record.ID) # 35 stay either higher or lower or in range
consistently_high <- d %>% filter(Record.ID %in% consistent ) %>% filter(avg_score_interval_cat == "Above") %>% pull(Record.ID) #3
consistently_range <- d %>% filter(Record.ID %in% consistent ) %>% filter(avg_score_interval_cat == "In range") %>% pull(Record.ID) 
consistently_low <- d %>% filter(Record.ID %in% consistent ) %>% filter(avg_score_interval_cat == "Below") %>% pull(Record.ID) #32
variable <- d %>% dplyr::count(Record.ID) %>% filter(n > 1) #19 patients switch

tacro_all_matched <- tacro_all_matched %>% 
  mutate(patient_class_iqr = ifelse(Record.ID %in% consistently_high, "Higher", 
                                    ifelse(Record.ID %in% consistently_range, "In range", 
                                           ifelse(Record.ID %in% consistently_low, "Lower", "Variable"))))

## Linear mixed model changes along immunosuppression dosage
### BAL
#### Tacrolimus

```{r}
bal_tacro_linear_test <- list()
for (i in 1:length(names(cell_set))) {
  
  bal_tacro_linear_test[[as.character(names(cell_set)[i])]] <- summary(lmerTest::lmer(TotalScore ~ Time_interval*tacrolimus_interp + (1|Record.ID), 
                                                                                      data = bal_tacro_matched %>% filter(Module == names(cell_set)[i]) ))$coefficients %>% 
    as.data.frame() %>% 
    mutate(Module = names(cell_set)[i]) 
}

bal_tacro_linear_test_df <- bind_rows(bal_tacro_linear_test) %>%
  rownames_to_column("Test") %>% 
  filter(!grepl("Intercept", Test)) %>%
  filter(`Pr(>|t|)` < 0.05) %>%
  filter(grepl("months:tac", Test)) 

sig <- bal_tacro_linear_test_df %>% pull(Module)

p <- 
  ggplot(bal_tacro_matched %>% filter(Module %in% sig) %>% arrange(Record.ID, Days), aes(x = tacrolimus_interp, y = TotalScore)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme +
  labs(title = "Modules with significant interaction with Tacrolimus", y = "Score", x = "Tacrolimus (trough)") +
  facet_wrap(Time_interval~Module, scales = "free")

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/bal-singscore-tacro.jpeg"), p, width = 15, height = 6, units = 'cm', dpi = 600)
```

#### Prednisolone

```{r}
bal_pred_linear_test <- list()
for (i in 1:length(names(cell_set))) {
  
  bal_pred_linear_test[[as.character(names(cell_set)[i])]] <- summary(lmerTest::lmer(TotalScore ~ Time_interval*pred_interp + (1|Record.ID), 
                                                                                     data = bal_pred_matched %>% filter(Module == names(cell_set)[i]) ))$coefficients %>% 
    as.data.frame() %>% 
    mutate(Module = names(cell_set)[i]) 
  
}

bal_pred_linear_test_df <- bind_rows(bal_pred_linear_test) %>% 
  rownames_to_column("Test") %>% 
  filter(!grepl("Intercept", Test)) %>%
  filter(`Pr(>|t|)` < 0.05) %>%
  filter(grepl("months:pred", Test)) 

sig <- bal_pred_linear_test_df %>% pull(Module)

p <- 
  ggplot(bal_pred_matched %>% filter(Module %in% sig, Time_interval == "6-12 months") %>% arrange(Record.ID, Days), aes(x = pred_interp, y = TotalScore)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme +
  labs(title = "Modules with significant interaction with Prednisolone", y = "Score", x = "Prednisolone (mg/kg)") +
  facet_wrap(Time_interval~Module, scales = "free")

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/bal-singscore-pred.jpeg"), p, width = 15, height = 12, units = 'cm', dpi = 600)
```





```{r}
# in range
tacro_all_matched %>% 
  filter(!is.na(within_range)) %>%
  ggplot(aes(x = Months, y = Record.ID)) +
  geom_tile(aes(fill = within_range))  +
  theme_bw() +
  theme(axis.text.y = element_blank() ) +
  labs(y = "Patient", title = "Tacrolimus trough levels", fill = "Alfred protocol\nrange") +
  facet_grid(patient_class_range~Immunosuppression_interval, scales = "free", space = "free")
```

```{r, fig.height=3, fig.width=4}
# tacro levels
ggplot() +
  #geom_point(data = tacro_all_matched %>% filter(Record.ID %in% c(12)), aes(x = Months, y = tacrolimus_interp), colour = "red", shape = 21) +
  geom_point(data = tacro_all_matched %>% filter(Record.ID %in% c(12)), aes(x = Months, y = tacrolimus_level))  +
  labs(title = "Patient 12", y = "Tacrolimus trough") +
  theme

# plot standard deviation
ggplot(tacro_all_matched %>% filter(Record.ID %in% c(2)), aes(x = Months, y = tacro_iqr_3m)) +
  geom_point() +
  geom_path()

ggplot() +
  geom_rect(data = expected_tacro %>% filter(Months <= 27), aes(xmin = Months-0.1, xmax = Months+0.1, ymin = exp_min, ymax = exp_max ), alpha = 0.2) +
  geom_point(data = tacro_all_matched %>% filter(Record.ID %in% c(1)), aes(x = Months, y = tacrolimus_level, fill = within_range), shape = 21) 
```

## Match to immunosuppression by interpolating

```{r, fig.height=3, fig.width=4}
pred_all_matched <- all_rna_metadata[, c("Patient_days", "Record.ID", "Days", "Months", "Group")] %>%
  mutate(Record.ID = factor(Record.ID)) %>%
  full_join(pred_all) %>% # do not remove pulses or we might interpolate wrong
  arrange(Record.ID, Days) %>%
  group_by(Record.ID) %>%
  mutate(pred_interp = round(na.approx(prednisolone_pulse_dose_mgd_per_kg, maxgap = 6, rule = 2, na.rm = FALSE), 2 )) %>%
  mutate(pred_interp = ifelse(Days == min(Days) & is.na(prednisolone_pulse_dose_mgd_per_kg), NA, pred_interp)) %>% # remove lower end of interpolation, not needed for higher end because they are all the same after 12M
  ungroup() %>%
  filter(Patient_days %in% all_rna_metadata$Patient_days) %>%
  filter(is.na(pulse) | pulse == "No")

ggplot() +
  geom_point(data = pred_all_matched %>% filter(Record.ID %in% c(12)), aes(x = Months, y = pred_interp ), colour = "red", shape = 21) +
  geom_point(data = pred_all_matched %>% filter(Record.ID %in% c(12)), aes(x = Months, y = prednisolone_pulse_dose_mgd_per_kg)) +
  labs(title = "Patient 12", y = "Prednisolone mg/kg") +
  theme +
  scale_y_log10()

length(pred_all_matched$Patient_days) # 252
```

# SD calculation
# use matched table to calculate standard deviation for 3months prior to bronchoscopy day using tacrolimus level (ignore NA values)
# at the end there will be SD values even if there is not tacro level values at that bronch date, there NA are there for SD is because all tacro levels were NA (e.g. tacro levels for first 2 weeks were not collected)
t_p <- list()

for (i in unique(tacro$Record.ID)){
  
  t_p[[i]] <- tacro %>% filter(Record.ID == i, !is.na(Group)) %>% dplyr::select(Days)
  
  t_sd <- list()
  
  for (e in unique(t_p[[i]]$Days)){
    
    t_sd[[e]] <- tacro %>% filter(Record.ID == i, Days <= e, Days > e - 90) %>% pull(tacrolimus_level) %>% sd(., na.rm = TRUE)
    
  }
  
  t_p[[i]]$tacro_sd_3m <- round(unlist(t_sd), 2)
  t_p[[i]] <- t_p[[i]] %>% mutate(Record.ID = i)
  
}

tacro_sd <- bind_rows(t_p)

#mean
p <- 
  
  ggplot() +
  geom_rect(data = expected_tacro %>% filter(Months <= 27), aes(xmin = Months-0.1, xmax = Months+0.1, ymin = exp_min, ymax = exp_max ), alpha = 0.2) +
  geom_ribbon(data = tacro_all_average, aes(x = Months, ymin = low_mean, ymax = high_mean), alpha = 0.3, fill = "blue", colour = "black", size = 0.2) +
  
  geom_path(data = tacro_all_average, aes(x = Months, y = tacro_average_levels), size = 0.3) +
  geom_point(data = tacro_all_average, aes(x = Months, y = tacro_average_levels), shape = 21, fill = "blue") +
  
  labs(title = "Tacrolimus post-transplant", y = "Trough level (average per month)", 
       x = "Time post-transplant (months)") +
  theme +
  theme(legend.position = "bottom")  +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 18, 24, 27))

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/tacrolimus-range-comparison-months-avg.jpeg"), p, width = 10, height = 10, units = 'cm', dpi = 600)


# in range
tacro_all_range <- tacro_all_matched %>%
  mutate(Months = round(Months)) %>%
  group_by(patient_class_range, Months) %>%
  summarise(tacro_median_levels = median(tacrolimus_level, na.rm = TRUE),
            high_median = median(tacrolimus_level, na.rm = TRUE) + iqr(tacrolimus_level, na.rm = TRUE),
            low_median  = median(tacrolimus_level, na.rm = TRUE) - iqr(tacrolimus_level, na.rm = TRUE)
  )

p <- 
  
  ggplot() +
  geom_rect(data = expected_tacro %>% filter(Months <= 27), aes(xmin = Months-0.1, xmax = Months+0.1, ymin = exp_min, ymax = exp_max ), alpha = 0.2) +
  #geom_ribbon(data = tacro_all_range, aes(x = Months, ymin = low_mean, ymax = high_mean), alpha = 0.3, fill = "blue", colour = "black", size = 0.2) +
  
  geom_path(data = tacro_all_range, aes(x = Months, y = tacro_median_levels, group = patient_class_range), size = 0.3) +
  geom_point(data = tacro_all_range, aes(x = Months, y = tacro_median_levels, fill = patient_class_range), shape = 21) +
  
  #geom_smooth(data = tacro_all_range, aes(x = Months, y = tacro_median_levels, group = patient_class_range), se = FALSE) +
  
  labs(title = "Tacrolimus post-transplant", y = "Trough level (median per month)", 
       x = "Time post-transplant (months)") +
  theme +
  theme(legend.position = "bottom")  +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 18, 24, 27)) +
  scale_fill_manual(values = c("High" = "red", "Low" = "yellow", "Variable" = "orange"))

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/tacrolimus-range-comparison-months-range-groups.jpeg"), p, width = 10, height = 10, units = 'cm', dpi = 600)


```{r}
tacro_all_matched %>% 
  dplyr::select(Immunosuppression_interval, Record.ID, patient_class_range, Group) %>% 
  distinct() %>%
  filter(Immunosuppression_interval %in% c("<3 months", "3-6 months")) %>%
  group_by(patient_class_range) %>%
  dplyr::count(Group) %>%
  mutate(Percentage = n/sum(n)*100) 

glm(Group ~ patient_class_median, data = tacro_all_matched, family = binomial) %>% summary()

# kaplan meyer curve
library(survival)
library(survminer)

tacro_all_matched$clad_status <- ifelse(is.na(tacro_all_matched$clad_status) | tacro_all_matched$clad_status == "Before", 0, 1 )

surv_obj <- Surv(time = tacro_all_matched$Days, event = tacro_all_matched$clad_status)
fit <- survfit(surv_obj ~ patient_class_range, data = tacro_all_matched)

ggsurvplot(
  fit,
  data = tacro_all_matched,
  risk.table = TRUE,       # show number at risk
  pval = TRUE,             # log-rank test p-value
  conf.int = TRUE,         # confidence intervals
  #legend.labs = c("Above median", "Below median"),
  legend.title = "Tacrolimus (0–6m)",
  xlab = "Time since transplant (months)",
  ylab = "CLAD-free survival probability",
  palette = c("#1f78b4", "#e31a1c", "green")
)

table(tacro_all_matched$clad_status) # 6 after CLAD
```


```{r}
# correlate
correlation_list <- list()

df1 <- t_p_average %>% filter(rownames(.) %in% colnames(rna_combined_heat_intervals))
df2 <- rna_combined_heat_intervals #%>% filter(rownames(.) %in% unique_features$Molecule)

correlation_list[["Genes"]] <- correlate_datasets(df1 = df1, 
                                                  df2 = df2, 
                                                  df1_type = "Immunosuppression", df2_type = "Genes", 
                                                  plots = FALSE, 
                                                  coef_threshold = 0, correlation = "pearson", pval_threshold = 1)
```

```{r}
#molecules_combined_heat <- as.data.frame(molecules_list$`group+months`$GroupCLAD$heat_sig_hits) 

molecules_combined_heat_intervals <- stable_small_molecules[molecules_list$`group+months`$GroupCLAD$sig_hits$Molecule,] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Patient_days") %>%
  left_join(stable_small_molecules_metadata[, c("Patient_days", "Month", "Record.ID", "Group")]) %>%
  dplyr::select(!c("Patient_days")) %>%
  group_by(Group, Record.ID, Month) %>%
  summarise(across(where(is.numeric), mean)) %>%
  ungroup() %>%
  mutate(Record_month = paste(Record.ID, Month, sep = "_")) %>%
  dplyr::select(!c("Group", "Record.ID", "Month")) %>%
  column_to_rownames("Record_month") %>%
  t() %>%
  as.data.frame()

df1 <- t_p_average %>% filter(rownames(.) %in% colnames(molecules_combined_heat_intervals))
df2 <- molecules_combined_heat_intervals #%>% filter(rownames(.) %in% unique_features$Molecule)

correlation_list[["Molecules"]] <- correlate_datasets(df1 = df1, 
                                                      df2 = df2, 
                                                      df1_type = "Immunosuppression", df2_type = "Molecules", 
                                                      plots = FALSE, 
                                                      coef_threshold = 0, correlation = "pearson", pval_threshold = 1)
```

```{r}
a <- correlation_list$Genes$coef_df %>% 
  dplyr::rename(Feature = df2_name) 

b <- correlation_list$Molecules$coef_df %>% 
  dplyr::rename(Feature = df2_name) 

corr_table <- rbind(a, b) %>% 
  dplyr::rename(Immunosuppression = df1_name) %>%
  mutate(Correlation2 = ifelse(pval >= 0.05, 0, Correlation),
         Sign = case_when(Correlation2 > 0 ~ "+", 
                          Correlation2 == 0 ~ "", 
                          Correlation2 < 0 ~ "-"), 
         Significance = case_when(pval >= 0.05 ~ "",
                                  pval < 0.05 & pval >= 0.01 ~ "*",
                                  pval < 0.01 & pval >= 0.001 ~ "**",
                                  pval < 0.001 & pval >= 0.0001 ~ "***",
                                  pval < 0.0001 ~ "****") )

p <- 
  ggplot(corr_table %>% filter(pval < 0.05), aes(x = Immunosuppression, y = Feature)) +
  geom_point(aes(size = -log(pval)), colour = "grey", stroke = 1) +
  geom_point(aes(colour = Correlation2, size = -log(pval))) +
  geom_text(aes(size = -log(pval), label = Sign), nudge_y = 0.1) +
  geom_text(aes(label = Significance), nudge_x = 0.15, nudge_y = 0.1) +
  facet_grid(df2_type~., scales = "free", space = "free") +
  labs(title = "Immunosuppression correlations", x = "", colour = "Pearson (r)", y = "Metabolite", size = "-Log(p-val)") +
  theme +
  theme(strip.text.y = element_text(angle = 0, size = 6),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  guides(fill = "none") +
  scale_colour_gradient2(low = "grey" , mid = "white", high = "orange")

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/summary-correlation.jpeg"), p, width = 13, height = 12, units = 'cm', dpi = 600) 
```

```{r}
data <- rna_combined_heat_intervals %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Record_month") %>%
  pivot_longer(!Record_month, names_to = "Gene", values_to = "Value") %>%
  left_join(t_p_average %>% rownames_to_column("Record_month")) %>%
  left_join(biomarker_rna_metadata %>% 
              mutate(Record_month = paste(Record.ID, Month, sep = "_")) %>% 
              dplyr::select(Record.ID, Record_month, Group) %>% distinct() )

sig_features <- corr_table %>%
  filter(pval < 0.05, Immunosuppression == "Tacrolimus") %>%
  pull(Feature) %>%
  unique()

p <-
  
  ggplot(data %>% filter(Gene %in% sig_features), aes(x = Tacrolimus, y = Value)) +
  geom_smooth(aes(group = Gene), method = "lm", colour = "darkgrey", lty = "dashed") +
  geom_point(aes(fill = Group), shape = 21, size = 1) +
  theme +
  labs(title = "Tacrolimus correlation", x = "Trough level", y = "Gene intensity", fill = "Group", 
       caption = "Each dot point is the average per patient per month") +
  facet_wrap(~Gene, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = c("CLAD" = "orange", "CLAD-free" = "darkgrey")) 

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/tacro-genes-correlation.jpeg"), p, width = 22, height = 8, units = 'cm', dpi = 600) 


sig_features <- corr_table %>%
  filter(pval < 0.05, Immunosuppression == "Prednisolone") %>%
  pull(Feature) %>%
  unique()

p <- 
  ggplot(data %>% filter(Gene %in% sig_features), aes(x = Prednisolone, y = Value)) +
  geom_smooth(aes(group = Gene), method = "lm", colour = "darkgrey", lty = "dashed") +
  geom_point(aes(fill = Group), shape = 21, size = 1) +
  theme +
  labs(title = "Prednisolone correlation", x = "mg (daily)", y = "Gene intensity", fill = "Group", 
       caption = "Each dot point is the average per patient per month") +
  facet_wrap(~Gene, ncol = 5, scales = "free_y") +
  scale_fill_manual(values = c("CLAD" = "orange", "CLAD-free" = "darkgrey")) 

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/pred-genes-correlation.jpeg"), p, width = 22, height = 12, units = 'cm', dpi = 600) 
```

```{r}
data <- molecules_combined_heat_intervals %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Record_month") %>%
  pivot_longer(!Record_month, names_to = "Molecule", values_to = "Value") %>%
  left_join(t_p_average %>% rownames_to_column("Record_month")) %>%
  left_join(biomarker_small_molecules_metadata %>% 
              mutate(Record_month = paste(Record.ID, Month, sep = "_")) %>% 
              dplyr::select(Record.ID, Record_month, Group) %>% distinct() ) 
unique(data$Molecule) # 9

sig_features <- corr_table %>%
  filter(pval < 0.05, Immunosuppression == "Tacrolimus") %>%
  pull(Feature) %>%
  unique()

p <- 
  ggplot(data %>% filter(Molecule %in% sig_features), aes(x = Tacrolimus, y = Value)) +
  geom_smooth(aes(group = Molecule), method = "lm", colour = "darkgrey", lty = "dashed") +
  geom_point(aes(fill = Group), shape = 21, size = 1) +
  theme +
  labs(title = "Tacrolimus correlation", x = "Trough level", y = "Metabolite intensity", fill = "Group", 
       caption = "Each dot point is the average per patient per month") +
  facet_wrap(~Molecule, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("CLAD" = "orange", "CLAD-free" = "darkgrey")) 

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/tacro-molecules-correlation.jpeg"), p, width = 17, height = 5, units = 'cm', dpi = 600) 


sig_features <- corr_table %>%
  filter(pval < 0.05, Immunosuppression == "Prednisolone") %>%
  pull(Feature) %>%
  unique()

p <- 
  ggplot(data %>% filter(Molecule %in% sig_features), aes(x = Prednisolone, y = Value)) +
  geom_smooth(aes(group = Molecule), method = "lm", colour = "darkgrey", lty = "dashed") +
  geom_point(aes(fill = Group), shape = 21, size = 1) +
  theme +
  labs(title = "Prednisolone correlation", x = "mg (daily)", y = "Metabolite intensity", fill = "Group", 
       caption = "Each dot point is the average per patient per month") +
  facet_wrap(~Molecule, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("CLAD" = "orange", "CLAD-free" = "darkgrey")) 

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/pred-molecules-correlation.jpeg"), p, width = 17, height = 12 , units = 'cm', dpi = 600) 
```


### Gene expression modelling

```{r}
# average tacro per month
tacro_stable_average <- tacro_stable %>%
  group_by(Month, Record.ID) %>%
  arrange(Days) %>%
  summarise(tacro_average_levels = mean(tacrolimus_level, na.rm = TRUE)) %>%
  left_join(expected_tacro) %>%
  mutate(tacro_deviation = tacro_average_levels - exp_min) %>%
  group_by(Record.ID) %>%
  #Between: do patients who are consistently above/below the guideline have different average expression?
  #Within: when a patient is above/below the guideline more than usual, does their expression shift?
  mutate(deviation_between = mean(tacro_deviation, na.rm = TRUE),
         deviation_within  = tacro_deviation - deviation_between) %>%
  ungroup()

# merge with data
cluster_genes <- readRDS(paste0(path, "Biomarker\ 2/clad-bal-multiomics-bm2/analysis/transcriptomics-small-molecules/tmixclust/rna/stable/cluster_genes.rds"))
stable_rna_adaptive <- cluster_genes[[2]]
length(intersect(stable_rna_adaptive, biomarker_rna_combined)) # 44

rna_combined_heat <- as.data.frame(deseq_list$months$ns.Months..df...3.1$heat_sig) %>% rownames_to_column("Gene") %>%
  full_join(as.data.frame(deseq_list$months$ns.Months..df...3.2$heat_sig) %>% rownames_to_column("Gene")) %>%
  full_join(as.data.frame(deseq_list$months$ns.Months..df...3.3$heat_sig) %>% rownames_to_column("Gene")) %>%
  filter(Gene %in% stable_rna_adaptive) %>%
  column_to_rownames("Gene")

# using actual normalised counts, not scaled counts assay(deseq_list$months$vsd)
rna_combined_heat_intervals <- rna_combined_heat[stable_rna_adaptive, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Patient_days") %>%
  left_join(stable_rna_metadata[, c("Patient_days", "Month", "Record.ID")]) %>%
  dplyr::select(!c("Patient_days")) %>%
  group_by(Record.ID, Month) %>%
  summarise(across(where(is.numeric), mean)) %>%
  ungroup() %>%
  mutate(Record.ID = factor(Record.ID)) 
dim(rna_combined_heat_intervals)[1] #164 

stable_rna_merge <- tacro_stable_average %>%
  inner_join(rna_combined_heat_intervals) %>%
  pivot_longer(-all_of(colnames(tacro_stable_average)), names_to="Gene", values_to="Expression") %>%
  filter(deviation_within < 10) # no significant genes when removing outlier

tacro_stable_average %>% dplyr::count(Record.ID, Month) %>% filter(n>1)
rna_combined_heat_intervals %>% dplyr::count(Record.ID, Month) %>% filter(n>1)
```

```{r}
fit_one_gene <- function(g) {
  dat <- stable_rna_merge %>% filter(Gene == g)
  # Keep spline small (df/k) given N=10x5
  m <- lmer(Expression ~ ns(Month, df = 3) + tacro_deviation + (1 | Record.ID), data = dat, REML = TRUE)
  out <- broom.mixed::tidy(m, effects = "fixed") %>% mutate(gene = g)
  list(model = m, tidy = out)
}

fits  <- lapply(stable_rna_adaptive, fit_one_gene)
results <- dplyr::bind_rows(lapply(fits, `[[`, "tidy"))

# Rank genes by |deviation_within| (strongest within-patient response to deviation from guideline)
ranked <- results %>%
  filter(term == "deviation_within") %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(desc(abs(estimate))) %>%
  filter(p_adj < 0.1)

ggplot(stable_rna_merge %>% filter(Gene == "JAK3"), aes(x = Month, y = Expression)) +
  geom_path(aes(group = Record.ID)) +
  geom_point() +
  geom_smooth()

ggplot(stable_rna_merge %>% filter(Gene == "ADORA2A"), aes(x = tacro_deviation, y = Expression)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)
```



### MCP

```{r}
# Example: 1 change point model
model_1cp <- list(
  TotalScore ~ 1,             # segment 1: intercept (flat)
  ~ 0 + tacro_interp_scaled     # segment 2: slope after cp_1
)

# Example: 2 change points (flat -> slope -> flat)
model_2cp <- list(
  TotalScore ~ 1,             # segment 1: intercept
  ~ 0 + tacro_interp_scaled,    # segment 2: slope
  ~ 1                         # segment 3: new plateau
)

# Fit for one module
fit <- mcp(
  model_2cp,
  data = bal_tacro_matched %>% filter(Module == "TCD8_Cytotoxic"),
  family = gaussian()
)

plot(fit)       # visualize fitted segments + change points
summary(fit)    # posterior estimates for slopes and cps
```

```{r}
fits <- list()
summaries <- list()

for (m in names(t_cell_set)) {
  dat <- bal_tacro_matched %>% filter(Module == m)
  
  fit <- mcp(model_2cp, data = dat, family = gaussian())
  fits[[m]] <- fit
  summaries[[m]] <- summary(fit)$summary %>%
    as.data.frame() %>%
    mutate(Module = m)
}

df <- bind_rows(summaries)
```


# GSVA

```{r}
# transcriptomics
gsvapar <- gsvaParam(all_rna_common_scaled, t_cell_set, kcdf = "Gaussian")
#gsvapar <- ssgseaParam(all_rna, t_cell_set)
bal_gsva_scores <- gsva(gsvapar)

# microarray
gsvapar <- gsvaParam(as.matrix(microarray_data_common_scaled), t_cell_set, kcdf = "Gaussian")
#gsvapar <- ssgseaParam(microarray.data, t_cell_set)
tissue_gsva_scores <- gsva(gsvapar)
```

## GSVA scores over time
## Transcriptomics

```{r}
d <- bal_gsva_scores %>% 
  as.data.frame() %>%
  rownames_to_column("module") %>%
  pivot_longer(!module, names_to = "Patient_days", values_to = "gsva_value") %>%
  right_join(all_rna_metadata[, c("Patient_days", "Record.ID", "Days", "Months", "Group", "days_to_clad")]) %>%
  mutate(Record.ID = factor(Record.ID)) 
# before CLAD onset or within a month of a CLAD diagnosis

ggplot(d, aes(x = Months, y = gsva_value)) +
  geom_point() +
  geom_smooth() + # line for all is the same as CLAD-free line
  #geom_smooth(aes(group = Group, colour = Group)) +
  facet_wrap(~module, scales = "free") +
  theme +
  #scale_colour_manual(values = c("CLAD" = "orange", "CLAD-free" = "darkgrey")) +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 18, 24))
```

```{r}
# Plot along immunosuppression
d <- bal_gsva_scores %>% 
  as.data.frame() %>%
  rownames_to_column("module") %>%
  pivot_longer(!module, names_to = "Patient_days", values_to = "gsva_value") %>%
  right_join(tacro_all_matched) %>%
  mutate(Record.ID = factor(Record.ID)) %>%
  filter(tacrolimus_level > 3 & tacrolimus_level < 15 ) #%>%
#mutate(tacrolimus_interp_scaled = scale(tacrolimus_interp))

d %>% 
  filter(module == "TCD8_TRM_exhausted_scRNAseq") %>%
  
  ggplot(aes(x = tacrolimus_level, y = gsva_value)) +
  geom_point(aes(colour = within_range)) +
  geom_smooth() +
  theme #+
#scale_colour_manual(values = c("CLAD" = "orange", "CLAD-free" = "darkgrey")) +
#facet_wrap(~Group)

unique(d$module)
```

```{r}
d <- gsva_scores %>% 
  as.data.frame() %>%
  rownames_to_column("module") %>%
  pivot_longer(!module, names_to = "Patient_days", values_to = "gsva_value") %>%
  right_join(pred_all_matched %>% filter(pulse == "No" | is.na(pulse))) %>%
  mutate(Record.ID = factor(Record.ID)) %>%
  mutate(pred_interp_scaled = scale(pred_interp)) %>%
  filter(pred_interp < 1 )

d %>% 
  filter(module == "TCD8_TRM_exhausted_scRNAseq") %>%
  
  ggplot(aes(x = pred_interp, y = gsva_value)) +
  geom_point() +
  geom_smooth() +
  theme 
```
## Microarray

```{r}
d <- microarray_gsva_scores %>% 
  as.data.frame() %>%
  rownames_to_column("module") %>%
  pivot_longer(!module, names_to = "Patient_days", values_to = "gsva_value") %>%
  right_join(microarray.metadata[, c("Patient_days", "Record.ID", "Days", "Months", "Group", "Days_interval")]) %>%
  mutate(Record.ID = factor(Record.ID)) 

ggplot(d, aes(x = Days_interval, y = gsva_value)) +
  geom_boxplot() +
  geom_path(aes(group = Record.ID)) +
  geom_point(aes(colour = Group)) +
  facet_wrap(~module, scales = "free") +
  theme
```

### Tissue

```{r}
tissue_tacro_matched <- tissue_score_df %>%
  mutate(Record.ID = factor(Record.ID)) %>%
  inner_join(tacro_all_matched) %>%
  mutate(tacro_interp_scaled = rescale(tacrolimus_interp, to = c(-1, 1))) %>%
  group_by(Record.ID) %>%
  arrange(Record.ID, Days) %>%
  ungroup()

tissue_pred_matched <- tissue_score_df %>%
  mutate(Record.ID = factor(Record.ID)) %>%
  inner_join(pred_all_matched) %>%
  mutate(pred_interp_scaled = rescale(pred_interp, to = c(-1, 1))) %>%
  group_by(Record.ID) %>%
  arrange(Record.ID, Days) %>%
  ungroup()
```

### Tissue
#### Tacrolimus

```{r}
tissue_tacro_linear_test <- list()
for (i in 1:length(names(t_cell_set))) {
  
  tissue_tacro_linear_test[[as.character(names(t_cell_set)[i])]] <- summary(lmerTest::lmer(TotalScore ~ Time_interval*tacrolimus_interp + (1|Record.ID), 
                                                                                           data = tissue_tacro_matched %>% filter(Module == names(t_cell_set)[i]) ))$coefficients %>% 
    as.data.frame() %>% 
    mutate(Module = names(t_cell_set)[i]) 
}

tissue_tacro_linear_test_df <- bind_rows(tissue_tacro_linear_test) %>%
  rownames_to_column("Test") %>% 
  filter(!grepl("Intercept", Test)) %>%
  filter(`Pr(>|t|)` < 0.05) %>%
  filter(grepl("months:tac", Test)) 

sig <- tissue_tacro_linear_test_df %>% pull(Module)

p <- 
  ggplot(tissue_tacro_matched %>% filter(Module %in% sig) %>% arrange(Record.ID, Days), aes(x = tacrolimus_interp, y = TotalScore)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme +
  labs(title = "Modules with significant interaction with Tacrolimus", y = "Score", x = "Tacrolimus (trough)") +
  facet_wrap(Module~Time_interval, scales = "free")

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/tissue-singscore-tacro.jpeg"), p, width = 15, height = 13, units = 'cm', dpi = 600)
```


#### Prednisolone

```{r}
tissue_pred_linear_test <- list()
for (i in 1:length(names(t_cell_set))) {
  
  tissue_pred_linear_test[[as.character(names(t_cell_set)[i])]] <- summary(lmerTest::lmer(TotalScore ~ Time_interval*pred_interp + (1|Record.ID), 
                                                                                          data = tissue_pred_matched %>% filter(Module == names(t_cell_set)[i]) ))$coefficients %>% 
    as.data.frame() %>% 
    mutate(Module = names(t_cell_set)[i]) 
  
}

tissue_pred_linear_test_df <- bind_rows(tissue_pred_linear_test) %>% 
  rownames_to_column("Test") %>% 
  filter(!grepl("Intercept", Test)) %>%
  filter(`Pr(>|t|)` < 0.05) %>%
  filter(grepl("months:pred", Test)) 

sig <- tissue_pred_linear_test_df %>% pull(Module)

p <- 
  ggplot(tissue_pred_matched %>% filter(Module %in% sig[1]) %>% arrange(Record.ID, Days), aes(x = pred_interp, y = TotalScore)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme +
  labs(title = "Modules with significant interaction with Prednisolone", y = "Score", x = "Prednisolone (mg/kg)") +
  facet_wrap(Module~Time_interval, scales = "free", ncol = 3)

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/tissue-singscore-pred.jpeg"), p, width = 15, height = 12, units = 'cm', dpi = 600)
```


# Tissue over time DE analysis

```{r}
design <- model.matrix(~Days_interval2 + Transplant_indication, microarray.metadata) # same for tx indication
dupcor <- duplicateCorrelation(microarray_data_common, design = design, block = microarray.metadata$Record.ID)

microarray_dge <- DGEList(microarray_data_common, 
                          samples = microarray.metadata, 
                          remove.zeros = TRUE)
microarray_dge <- calcNormFactors(microarray_dge)
v <- voom(microarray_dge, design, plot=TRUE)
vfit <- lmFit(v, design, block = microarray.metadata$Record.ID, correlation = dupcor$consensus)
efit <- eBayes(vfit)
summary(decideTests(efit, lfc = 0.1, p.value = 0.1)) 

topTreat(efit, coef="Days_interval2>12 months", n=Inf) %>% filter(adj.P.Val < 0.01)

```



```{r}
p <- 
  ggplot(bal_score_df, aes(x = TotalScore, y = TotalDispersion)) +
  geom_point(aes(fill = Immunosuppression_interval), shape = 21) +
  theme_bw() +
  theme +
  labs(title = "Score vs Dispersion", fill = "Immunosuppression interval") +
  facet_wrap(.~Module, scales = "free")

ggsave(paste0(path, "Biomarker\ 2/immunosuppression-analysis/figures/all/bal-singscore-disp.jpeg"), p, width = 16, height = 10, units = 'cm', dpi = 600)
```


```{r}
scoredf1 <- bal_score_df %>% filter(Module == "TCD8_Cytotoxic") %>% column_to_rownames("Patient_days")
scoredf2 <- bal_score_df %>% filter(Module == "Neutrophils") %>% column_to_rownames("Patient_days")
identical(rownames(scoredf1), rownames(scoredf2)) #TRUE

p <- plotScoreLandscape(
  scoredf1 = scoredf1,
  scoredf2 = scoredf2,
  scorenames = c('TCD8 score', 'Neutrophil score')
)

projectScoreLandscape(
  p,
  scoredf1 = scoredf1,
  scoredf2 = scoredf2,
  annot = scoredf1$Immunosuppression_interval,
  annot_name = 'Immunosuppression_interval',
  isInteractive = FALSE)


```

