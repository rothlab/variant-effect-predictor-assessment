library(data.table)

# Set arguments if not passed in already as config
if (exists("config")) {
  GENE = tolower(config$gene)
} else {
  GENE = "abca1"
}

# Load filtered phenotypes and variants
phenotypes = fread(sprintf("output/%s/%s_phenotypes_filtered.csv", toupper(GENE), GENE))
variants = fread(sprintf("output/%s/%s_variants_filtered.csv", toupper(GENE), GENE))

# Merge variants and phenotypes
variants = merge(variants, phenotypes, by = "eid")

# Select quantitative phenotypes
columns = colnames(variants)
columns = columns[!columns %in% bPhenotypes[field_type == "categorical", field_id]]
quantPhenotype = variants[, ..columns]
quantPhenotype[quantPhenotype == ""] = NA

# Save quantitative and all phenotypes separately
fwrite(quantPhenotype, sprintf("output/%s/%s_quant-variants_measures.csv", toupper(GENE), GENE))
fwrite(variants, sprintf("output/%s/%s_variants_phenotypes.csv", toupper(GENE), GENE))