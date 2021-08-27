library(data.table)
library(ggplot2)
library(stringr)

# Load configuration file
source("lib/loadConfig.R")
config = loadConfig("config.yaml")

# Set arguments if not passed in already as config
UNIQUE_VARIANTs_PATH = config$unique_variants_path
VARIANT_PREDICTORS = config$variant_predictors
OCCURANCE_CUTOFF = config$occurance_cutoff

# Load genes
geneList = fread("pvals.csv")
genes = unique(geneList$gene)

###
# Survey variant predictors in genes
###
# Collect variant predictor status
predStatus = lapply(genes, function(gene) {
  # Load filtered variants
  filePath = sprintf("output/%s/%s_variants_filtered.csv", toupper(gene), gene)
  if (!file.exists(filePath)) return(NULL)
  filteredVars = fread(filePath)
  filteredVars = unique(filteredVars$variant_name)
  
  # Load variant predictors
  filePath = paste("input", toupper(gene), UNIQUE_VARIANTs_PATH, sep = "/")
  if (!file.exists(filePath)) return(NULL)
  varPredictors = fread(filePath, select = c("variant_name", "type", "hgvs_pro",
                                             unlist(VARIANT_PREDICTORS, use.names = F)))
  
  # Restrict to selected variants
  varPredictors = varPredictors[variant_name %in% filteredVars]
  
  # Categorize variant predictors:
  # "Pass" - predictors made predictions for at least the number of variants set by the OCCURANCE_CUTOFF
  # "Fail" - predictors made predictions for less than the number of variants set by the OCCURANCE_CUTOFF 
  predStatus = apply(varPredictors[, -c("variant_name", "type", "hgvs_pro"), with = F], 2, function(col) {
    if (length(na.omit(col)) >= OCCURANCE_CUTOFF)
      return(sprintf("Pass (at least %d variants)", OCCURANCE_CUTOFF))
    return("Fail")
  })
  predStatus = as.data.table(predStatus, keep.rownames = T)
  predStatus$gene_symbol = gene
  predStatus = predStatus[, .(gene_name = gene, predictor = rn, status = predStatus, n_muts = nrow(varPredictors))]
  
  return(predStatus)
})
predStatus = rbindlist(predStatus)

# Include missing gene-predictor combinations
fullGrid = expand.grid(gene_name = genes, predictor = unique(predStatus$predictor))
predStatus = as.data.table(merge(fullGrid, predStatus, all.x = T))
predStatus[is.na(status), status := "Fail"]

# Attach variant predictor names
predNames = as.data.table(unlist(VARIANT_PREDICTORS), keep.rownames = T)
colnames(predNames) = c("predictor_name", "predictor")
predStatus = merge(predStatus, predNames, by = "predictor")

# Plot gene-predictor status
plot = ggplot(predStatus) +
  geom_tile(aes(x = toupper(gene_name), y = predictor_name, fill = status)) +
  labs(x = "Burden test-associated human gene", y = "Variant effect predictor") +
  scale_fill_manual(name = "Variant Prediction Coverage Cutoff",
                    values = c("#edf6f9", "#006d77")) +
  ggtitle(sprintf("Number of genes = %d", length(genes))) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="bottom")
ggsave("variant_predictor_survey.png", plot, height = 7.5, width = 22,
       units = "in", dpi = 300)
