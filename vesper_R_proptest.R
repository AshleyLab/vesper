library(matrixStats)
library(tidyverse)
library (dplyr)
library(stringr)

ex <- 9


data_SNV_pairs <- read.csv(sprintf("/Users/kaiser19/Documents/RBM20/CX_debug_19Jan2026_MAF_iPSC_BC/Exon%d/sample_snv_pairs/exon%d_co_travelling_snv_pairs.csv", ex, ex))
data_snvs <- read.csv(sprintf("/Users/kaiser19/Documents/RBM20/CX_debug_19Jan2026_MAF_iPSC_BC/Exon%d/exon%d_reformatted_mpileup_MAF_iPSC_BC.csv", ex, ex))

colnames(data_SNV_pairs) <- c("Variant1","Variant2", "calculated_Variant1_mpileup_count_B", "calculated_Variant1_mpileup_count_C", "calculated_Variant1_mpileup_count_iPSC",
                            "calculated_Variant2_mpileup_count_B", "calculated_Variant2_mpileup_count_C", "calculated_Variant2_mpileup_count_iPSC")

all_RBM20_vars_consequence <- read.csv("/Users/kaiser19/Documents/RBM20/Data/RBM20_transcript_all_poss_vars_for_R_code.csv")


# Add REF and observed_nucleotide columns
data_snvs$REF <- sub(".*([A-Z])>.*", "\\1", data_snvs$Variant)
data_snvs$observed_nucleotide <- sub(".*>([A-Z])", "\\1", data_snvs$Variant)

# Add AA_consequence, Consequence, and Pathogenic columns
all_consq_clean <- all_RBM20_vars_consequence %>%
  filter(Consequence != "REFERENCE ALLELE") %>%              # keep only non-reference
  mutate(
    Variant = str_remove(cDNA_Consequence, "^c\\.")             # "c.1881A>C" -> "1881A>C"
  ) %>%
  # If a Variant appears multiple times (e.g., multiple transcripts),
  # keep one row per Variant. If you prefer to combine, see the note below.
  distinct(Variant, .keep_all = TRUE)

## 2) Map AA_consequence / Consequence / Pathogenic into data_snvs by Variant
data_snvs <- data_snvs %>%
  left_join(all_consq_clean %>% select(Variant, AA_consequence, Consequence, Pathogenic),
    by = "Variant"
  )

data_snvs <- data_snvs %>%
  mutate(
    Consequence_clean = Consequence %>%
      str_replace_all('"', '') %>%
      str_to_lower(),
    Pathogenicity = case_when(
      Pathogenic == "Yes" ~ "Pathogenic",
      Pathogenic == "No" & str_detect(Consequence_clean, "stop_gained") ~ "Stop_Gained",
      Pathogenic == "No" & str_detect(Consequence_clean, "synonymous")  ~ "Synonymous",
      TRUE ~ "VUS"
    )
  ) %>%
  select(-Consequence_clean)

colnames(data_snvs)[colnames(data_snvs) == "REF"] <- "reference_nucleotide"
colnames(data_snvs)[colnames(data_snvs) == "AA_consequence"] <- "Var_p"
colnames(data_snvs)[colnames(data_snvs) == "Consequence"] <- "Var_type"
colnames(data_snvs)[colnames(data_snvs) == "Pathogenic"] <- "P_LP"
colnames(data_snvs)[colnames(data_snvs) == "Pathogenicity"] <- "PLP_Stopgain_Syn_VUS"

proptest1_dir <- sprintf("/Users/kaiser19/Documents/RBM20/CX_debug_19Jan2026_MAF_iPSC_BC/Exon%d/Exon%d_R_Proptest_Results/Proptest1", ex, ex)

# Create the proptest1 directory if it doesn't already exist
if (!dir.exists(proptest1_dir)) {
  dir.create(proptest1_dir, recursive = TRUE)
}

# proptest iPSC vs B+C 
proptest_results_iPSC_vs_BC <- data.frame(matrix(ncol = 9, nrow = 0))

for (i in 1:nrow(data_snvs)) {
  enrichment <- matrix()
  data_snvs$number_obs_exp<- data_snvs$Total_Count_iPSC
  data_snvs$number_obs_con <- data_snvs$Total_Count_B+data_snvs$Total_Count_C
  data_snvs$depth_exp <-data_snvs$Total_Depth_iPSC
  data_snvs$depth_con <-data_snvs$Total_Depth_B+data_snvs$Total_Depth_C
  results <- prop.test(x=c(data_snvs$number_obs_exp[i],data_snvs$number_obs_con[i]),n=c(data_snvs$depth_exp[i],data_snvs$depth_con[i])) #apply prop.test to row
  single_result <- c(data_snvs$Variant[i],data_snvs$reference_nucleotide[i], data_snvs$observed_nucleotide[i], data_snvs$PLP_Stopgain_Syn_VUS[i], results$estimate, results$conf.int, results$p.value) #extract test results
  single_result_row <- t(single_result)
  proptest_results_iPSC_vs_BC <- rbind(proptest_results_iPSC_vs_BC,single_result_row) #add to proptest results df
}
x <- c("Variant", "reference", "observed_nucleotide", "PLP_Stopgain_Syn_VUS", "proportion_exp","proportion_con","CI95_lower","CI95_upper","p_value")
colnames(proptest_results_iPSC_vs_BC) <- x
class(proptest_results_iPSC_vs_BC$proportion_exp)="numeric"
class(proptest_results_iPSC_vs_BC$proportion_con)="numeric"
class(proptest_results_iPSC_vs_BC$CI95_upper)="numeric"
class(proptest_results_iPSC_vs_BC$CI95_lower)="numeric"
class(proptest_results_iPSC_vs_BC$p_value)="numeric"
proptest_results_iPSC_vs_BC$exp_enrichment <- (proptest_results_iPSC_vs_BC$proportion_exp/proptest_results_iPSC_vs_BC$proportion_con)
proptest_results_iPSC_vs_BC$neglogP <- -(log10(proptest_results_iPSC_vs_BC$p_value))
proptest_results_iPSC_vs_BC$Var_p <-data_snvs$Var_p
proptest_results_iPSC_vs_BC$Var_type <-data_snvs$Var_type
proptest_results_iPSC_vs_BC$P_LP <-data_snvs$P_LP
proptest_results_iPSC_vs_BC$PLP_Stopgain_Syn_VUS <- data_snvs$PLP_Stopgain_Syn_VUS
proptest_results_iPSC_vs_BC$log2_exp_enrichment <- log2(proptest_results_iPSC_vs_BC$exp_enrichment)
# make a copy of proptest_results_iPSC_vs_BC for further analysis
proptest_results_iPSC_vs_BC_copy <- proptest_results_iPSC_vs_BC
proptest_results_iPSC_vs_BC_copy$Total_Count_B_C <- data_snvs$Total_Count_B+data_snvs$Total_Count_C
write.csv(proptest_results_iPSC_vs_BC_copy, file = file.path(proptest1_dir, "proptest_results_iPSC_vs_BC.csv"), row.names = FALSE)



#######end proptests for different phenotypes, now individually subtract each pair 
# and re-do prop-test and visualize the spread of new values
# take the statistically significant enrichment value closest to 1 to correct for any effect of other variants
data_SNV_pairs$SNVcount_Var1_Count_iPSC <- ""
data_SNV_pairs$SNVcount_Var1_Count_B <- ""
data_SNV_pairs$SNVcount_Var1_Count_C <- ""

data_SNV_pairs$SNVdepth_Var1_Depth_iPSC <- ""
data_SNV_pairs$SNVdepth_Var1_Depth_B <- ""
data_SNV_pairs$SNVdepth_Var1_Depth_C <- ""


for (i in 1:nrow(data_SNV_pairs)) {
  df <- data_snvs[which(data_snvs$Variant==data_SNV_pairs$Variant1[i]),]
    if (nrow(df) == 0) {
      next
     }
  data_SNV_pairs$SNVcount_Var1_Count_iPSC[i] <- df$Total_Count_iPSC
  data_SNV_pairs$SNVdepth_Var1_Depth_iPSC[i] <- df$Total_Depth_iPSC
  data_SNV_pairs$SNVcount_Var1_Count_B[i] <- df$Total_Count_B
  data_SNV_pairs$SNVdepth_Var1_Depth_B[i] <- df$Total_Depth_B
  data_SNV_pairs$SNVcount_Var1_Count_C[i] <- df$Total_Count_C
  data_SNV_pairs$SNVdepth_Var1_Depth_C[i] <- df$Total_Depth_C
  
}

data_SNV_pairs$SNVcount_Var2_Count_iPSC <- ""
data_SNV_pairs$SNVcount_Var2_Count_B <- ""
data_SNV_pairs$SNVcount_Var2_Count_C <- ""  

data_SNV_pairs$SNVdepth_Var2_Depth_iPSC <- ""
data_SNV_pairs$SNVdepth_Var2_Depth_B <- ""
data_SNV_pairs$SNVdepth_Var2_Depth_C <- ""      

for (i in 1:nrow(data_SNV_pairs)) {
  df <- data_snvs[which(data_snvs$Variant==data_SNV_pairs$Variant2[i]),]
    if (nrow(df) == 0) {
      next
     }
  data_SNV_pairs$SNVcount_Var2_Count_iPSC[i] <- df$Total_Count_iPSC
  data_SNV_pairs$SNVdepth_Var2_Depth_iPSC[i] <- df$Total_Depth_iPSC
  data_SNV_pairs$SNVcount_Var2_Count_B[i] <- df$Total_Count_B
  data_SNV_pairs$SNVdepth_Var2_Depth_B[i] <- df$Total_Depth_B
  data_SNV_pairs$SNVcount_Var2_Count_C[i] <- df$Total_Count_C
  data_SNV_pairs$SNVdepth_Var2_Depth_C[i] <- df$Total_Depth_C
  
}


##Now recalculate the counts in each group per variant when counts of reads with variant pairs are subtracted
data_SNV_pairs$deltaSNVcount_Var1_iPSC <- (as.numeric(data_SNV_pairs$SNVcount_Var1_Count_iPSC) - data_SNV_pairs$calculated_Variant1_mpileup_count_iPSC)
data_SNV_pairs$deltaSNVdepth_Var1_iPSC <- (as.numeric(data_SNV_pairs$SNVdepth_Var1_Depth_iPSC) - data_SNV_pairs$calculated_Variant1_mpileup_count_iPSC)
data_SNV_pairs$deltaSNVcount_Var1_iPSC <- ifelse(data_SNV_pairs$deltaSNVcount_Var1_iPSC == 0, 1, data_SNV_pairs$deltaSNVcount_Var1_iPSC)

data_SNV_pairs$deltaSNVcount_Var1_B <- (as.numeric(data_SNV_pairs$SNVcount_Var1_Count_B) - data_SNV_pairs$calculated_Variant1_mpileup_count_B)
data_SNV_pairs$deltaSNVdepth_Var1_B <- (as.numeric(data_SNV_pairs$SNVdepth_Var1_Depth_B) - data_SNV_pairs$calculated_Variant1_mpileup_count_B)
data_SNV_pairs$deltaSNVcount_Var1_B <- ifelse(data_SNV_pairs$deltaSNVcount_Var1_B == 0, 1, data_SNV_pairs$deltaSNVcount_Var1_B)

data_SNV_pairs$deltaSNVcount_Var1_C <- (as.numeric(data_SNV_pairs$SNVcount_Var1_Count_C) - data_SNV_pairs$calculated_Variant1_mpileup_count_C)
data_SNV_pairs$deltaSNVdepth_Var1_C <- (as.numeric(data_SNV_pairs$SNVdepth_Var1_Depth_C) - data_SNV_pairs$calculated_Variant1_mpileup_count_C)
data_SNV_pairs$deltaSNVcount_Var1_C <- ifelse(data_SNV_pairs$deltaSNVcount_Var1_C == 0, 1, data_SNV_pairs$deltaSNVcount_Var1_C)

data_SNV_pairs$deltaSNVcount_Var2_iPSC <- (as.numeric(data_SNV_pairs$SNVcount_Var2_Count_iPSC) - data_SNV_pairs$calculated_Variant2_mpileup_count_iPSC)
data_SNV_pairs$deltaSNVdepth_Var2_iPSC <- (as.numeric(data_SNV_pairs$SNVdepth_Var2_Depth_iPSC) - data_SNV_pairs$calculated_Variant2_mpileup_count_iPSC)
data_SNV_pairs$deltaSNVcount_Var2_iPSC <- ifelse(data_SNV_pairs$deltaSNVcount_Var2_iPSC == 0, 1, data_SNV_pairs$deltaSNVcount_Var2_iPSC)

data_SNV_pairs$deltaSNVcount_Var2_B <- (as.numeric(data_SNV_pairs$SNVcount_Var2_Count_B) - data_SNV_pairs$calculated_Variant2_mpileup_count_B)
data_SNV_pairs$deltaSNVdepth_Var2_B <- (as.numeric(data_SNV_pairs$SNVdepth_Var2_Depth_B) - data_SNV_pairs$calculated_Variant2_mpileup_count_B)
data_SNV_pairs$deltaSNVcount_Var2_B <- ifelse(data_SNV_pairs$deltaSNVcount_Var2_B == 0, 1, data_SNV_pairs$deltaSNVcount_Var2_B)

data_SNV_pairs$deltaSNVcount_Var2_C <- (as.numeric(data_SNV_pairs$SNVcount_Var2_Count_C) - data_SNV_pairs$calculated_Variant2_mpileup_count_C)
data_SNV_pairs$deltaSNVdepth_Var2_C <- (as.numeric(data_SNV_pairs$SNVdepth_Var2_Depth_C) - data_SNV_pairs$calculated_Variant2_mpileup_count_C)
data_SNV_pairs$deltaSNVcount_Var2_C <- ifelse(data_SNV_pairs$deltaSNVcount_Var2_C == 0, 1, data_SNV_pairs$deltaSNVcount_Var2_C)



##now do proptests across variants with delta counts for each phenotype
proptest2_dir <- sprintf("/Users/kaiser19/Documents/RBM20/CX_debug_19Jan2026_MAF_iPSC_BC/Exon%d/Exon%d_R_Proptest_Results/Proptest2", ex, ex)

# Create the proptest2 directory if it doesn't already exist
if (!dir.exists(proptest2_dir)) {
  dir.create(proptest2_dir, recursive = TRUE)
}

# iPSC vs B+C
proptest_variant1_iPSC_vs_B_C <- data.frame(matrix(ncol = 15, nrow = 0))

for (i in 1:nrow(data_SNV_pairs)) {
  enrichment <- matrix()
  if (is.na(data_SNV_pairs$deltaSNVcount_Var1_iPSC[i])
      |is.na(data_SNV_pairs$deltaSNVcount_Var1_B[i])
      |is.na(data_SNV_pairs$deltaSNVcount_Var1_C[i])) {
    next
  }
  results <- prop.test(x=c(data_SNV_pairs$deltaSNVcount_Var1_iPSC[i], (data_SNV_pairs$deltaSNVcount_Var1_B[i]+data_SNV_pairs$deltaSNVcount_Var1_C[i])),
                       n=c(data_SNV_pairs$deltaSNVdepth_Var1_iPSC[i], (data_SNV_pairs$deltaSNVdepth_Var1_B[i]+data_SNV_pairs$deltaSNVdepth_Var1_C[i]))) #apply prop.test to row

  single_result <- c(data_SNV_pairs$Variant1[i],data_SNV_pairs$Variant2[i],results$estimate, results$conf.int, results$p.value) #extract test results
  single_result_row <- t(single_result)
  proptest_variant1_iPSC_vs_B_C <- rbind(proptest_variant1_iPSC_vs_B_C,single_result_row) #add to proptest results df
}

colnames(proptest_variant1_iPSC_vs_B_C) <- c("Variant1","Variant2","proportion_exp","proportion_con","CI95_lower","CI95_upper","p_value") 

proptest_variant1_iPSC_vs_B_C$proportion_exp <- as.numeric(proptest_variant1_iPSC_vs_B_C$proportion_exp)
proptest_variant1_iPSC_vs_B_C$proportion_con <- as.numeric(proptest_variant1_iPSC_vs_B_C$proportion_con)
proptest_variant1_iPSC_vs_B_C$exp_enrichment <- (proptest_variant1_iPSC_vs_B_C$proportion_exp/proptest_variant1_iPSC_vs_B_C$proportion_con)
proptest_variant1_iPSC_vs_B_C$log2_exp_enrichment<- log2(proptest_variant1_iPSC_vs_B_C$exp_enrichment)
proptest_variant1_iPSC_vs_B_C$CI95_upper <- as.numeric(proptest_variant1_iPSC_vs_B_C$CI95_upper)
proptest_variant1_iPSC_vs_B_C$CI95_lower <- as.numeric(proptest_variant1_iPSC_vs_B_C$CI95_lower)
proptest_variant1_iPSC_vs_B_C$p_value <- as.numeric(proptest_variant1_iPSC_vs_B_C$p_value)
proptest_variant1_iPSC_vs_B_C$neglogP <- -(log10(proptest_variant1_iPSC_vs_B_C$p_value))
write.csv(proptest_variant1_iPSC_vs_B_C, file.path(proptest2_dir, "proptest_variant1_iPSC_vs_B_Cs.csv"), row.names = FALSE)

# iPSC vs B+C
proptest_variant2_iPSC_vs_B_C <- data.frame(matrix(ncol = 15, nrow = 0)) 
for (i in 1:nrow(data_SNV_pairs)) {
  enrichment <- matrix()
  if (is.na(data_SNV_pairs$deltaSNVcount_Var2_iPSC[i])
      |is.na(data_SNV_pairs$deltaSNVcount_Var2_B[i])
      |is.na(data_SNV_pairs$deltaSNVcount_Var2_C[i])) {
    next
  }
  results <- prop.test(x=c(data_SNV_pairs$deltaSNVcount_Var2_iPSC[i],
                           (data_SNV_pairs$deltaSNVcount_Var2_B[i]+data_SNV_pairs$deltaSNVcount_Var2_C[i])),
                       n=c(data_SNV_pairs$deltaSNVdepth_Var2_iPSC[i],
                           (data_SNV_pairs$deltaSNVdepth_Var2_B[i]+data_SNV_pairs$deltaSNVdepth_Var2_C[i]))) #apply prop.test to row

  single_result <- c(data_SNV_pairs$Variant2[i],data_SNV_pairs$Variant1[i],results$estimate, results$conf.int, results$p.value) #extract test results
  single_result_row <- t(single_result)
  proptest_variant2_iPSC_vs_B_C <- rbind(proptest_variant2_iPSC_vs_B_C,single_result_row) #add to proptest results df
}
rm(results)
rm(enrichment)
rm(i)
rm(single_result)
rm(single_result_row)               

colnames(proptest_variant2_iPSC_vs_B_C) <- c("Variant1","Variant2","proportion_exp","proportion_con","CI95_lower","CI95_upper","p_value")
proptest_variant2_iPSC_vs_B_C$proportion_exp <- as.numeric(proptest_variant2_iPSC_vs_B_C$proportion_exp)
proptest_variant2_iPSC_vs_B_C$proportion_con <- as.numeric(proptest_variant2_iPSC_vs_B_C$proportion_con)
proptest_variant2_iPSC_vs_B_C$exp_enrichment <- (proptest_variant2_iPSC_vs_B_C$proportion_exp/proptest_variant2_iPSC_vs_B_C$proportion_con)
proptest_variant2_iPSC_vs_B_C$log2_exp_enrichment<- log2(proptest_variant2_iPSC_vs_B_C$exp_enrichment)
proptest_variant2_iPSC_vs_B_C$CI95_upper <- as.numeric(proptest_variant2_iPSC_vs_B_C$CI95_upper)
proptest_variant2_iPSC_vs_B_C$CI95_lower <- as.numeric(proptest_variant2_iPSC_vs_B_C$CI95_lower)
proptest_variant2_iPSC_vs_B_C$p_value <- as.numeric(proptest_variant2_iPSC_vs_B_C$p_value)
proptest_variant2_iPSC_vs_B_C$neglogP <- -(log10(proptest_variant2_iPSC_vs_B_C$p_value))
write.csv(proptest_variant2_iPSC_vs_B_C, file.path(proptest2_dir, "proptest_variant2_iPSC_vs_B_C.csv"), row.names = FALSE)

# merge for each phenotype, finding highest and lowest values and value closest to 1 + p-value below
merged_delta_proptest_dir <- sprintf("/Users/kaiser19/Documents/RBM20/CX_debug_19Jan2026_MAF_iPSC_BC/Exon%d/Exon%d_R_Proptest_Results/Merged_Delta_Proptest", ex, ex)

# Create the merged_delta_proptest_dir directory if it doesn't already exist
if (!dir.exists(merged_delta_proptest_dir)) {
  dir.create(merged_delta_proptest_dir, recursive = TRUE)
}

# iPSC_vs_B_C
Delta_proptests_combined_iPSC_vs_B_C <- rbind(proptest_variant1_iPSC_vs_B_C,proptest_variant2_iPSC_vs_B_C)
rm(proptest_variant1_iPSC_vs_B_C,proptest_variant2_iPSC_vs_B_C)

# Delta_proptests_combined_iPSC_vs_B_C <- Delta_proptests_combined_iPSC_vs_B_C[which(Delta_proptests_combined_iPSC_vs_B_C$Variant1 %in% data_snvs$Variant),]
Delta_proptests_combined_iPSC_vs_B_C <- Delta_proptests_combined_iPSC_vs_B_C[
  Delta_proptests_combined_iPSC_vs_B_C$Variant1 %in% data_snvs$Variant |
  Delta_proptests_combined_iPSC_vs_B_C$Variant2 %in% data_snvs$Variant, ]       

Merged_delta_proptest_by_SNV_iPSC_vs_B_C <- 
  Delta_proptests_combined_iPSC_vs_B_C %>%
  group_by(Variant1)  %>%
  summarize(delta_enrichment_min = min(exp_enrichment),
            delta_enrichment_min_p_value = p_value[which.min(exp_enrichment)],
            delta_enrichment_max = max(exp_enrichment),
            delta_enrichment_max_p_value = p_value[which.max(exp_enrichment)],
            delta_enrichment_mean = mean(exp_enrichment),
            delta_enrichment_median = median(exp_enrichment),
            delta_enrichment_SD = sd(exp_enrichment),
            delta_enrichment_closest_to_1 = exp_enrichment[which.min(abs(exp_enrichment-1))],
            delta_enrichment_closest_to_1_p_value = p_value[which.min(abs(exp_enrichment-1))]
  ) 

#annotate by adding back to data_snvs
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$reference_nucleotide <- ""
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$observed_nucleotide <- ""
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$Var_p <- ""
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$Var_type <- "" 
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$PLP_Stopgain_Syn_VUS <-""
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$P_LP <- ""
for (i in 1:nrow(Merged_delta_proptest_by_SNV_iPSC_vs_B_C)) {
  df <- proptest_results_iPSC_vs_BC[which(proptest_results_iPSC_vs_BC$Variant==Merged_delta_proptest_by_SNV_iPSC_vs_B_C$Variant1[i]),]
  if (nrow(df)==0) {
    next }
  Merged_delta_proptest_by_SNV_iPSC_vs_B_C$reference_nucleotide[i] <- df$reference
  Merged_delta_proptest_by_SNV_iPSC_vs_B_C$observed_nucleotide[i] <- df$observed_nucleotide
  Merged_delta_proptest_by_SNV_iPSC_vs_B_C$Var_p[i] <- df$Var_p
  Merged_delta_proptest_by_SNV_iPSC_vs_B_C$Var_type[i] <- df$Var_type
  Merged_delta_proptest_by_SNV_iPSC_vs_B_C$PLP_Stopgain_Syn_VUS[i] <-df$PLP_Stopgain_Syn_VUS
  Merged_delta_proptest_by_SNV_iPSC_vs_B_C$P_LP[i] <- df$P_LP
}

# #_#_#_#_#_#_#_#

# # #calculate FDR p_adj
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$p_adj_ct1 <- p.adjust(Merged_delta_proptest_by_SNV_iPSC_vs_B_C$delta_enrichment_closest_to_1_p_value,"fdr")
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$p_adj_min <- p.adjust(Merged_delta_proptest_by_SNV_iPSC_vs_B_C$delta_enrichment_min_p_value,"fdr")
Merged_delta_proptest_by_SNV_iPSC_vs_B_C$p_adj_max <- p.adjust(Merged_delta_proptest_by_SNV_iPSC_vs_B_C$delta_enrichment_max_p_value,"fdr")
write.csv(Merged_delta_proptest_by_SNV_iPSC_vs_B_C, file.path(merged_delta_proptest_dir, "Merged_delta_proptest_by_SNV_iPSC_vs_B_C.csv"), row.names = FALSE)
