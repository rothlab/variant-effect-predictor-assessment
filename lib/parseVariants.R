library(data.table)
library(stringr)

# Set arguments if not passed in already as config
GENE = tolower(config$gene)
ALL_VARIANTs_PATH = config$all_variants_path
UNIQUE_VARIANTs_PATH = config$unique_variants_path
WITHDRAW_EIDs_PATH = config$withdraw_eids_path

# Load variants related files
# 1) All variants with HGVS
# 2) Unique vairants with weights
GENE = tolower(GENE)
allVariants = fread(paste("input", GENE, ALL_VARIANTs_PATH, sep = "/"), drop = 1)
uniqueVariants = fread(paste("input", GENE, UNIQUE_VARIANTs_PATH, sep = "/"))

# Select unique variants:
# 1) gnomAD MAF < 0.1% (0.001)
# 2) Missense
# 3) Not unknown (?)
uniqueVariants = uniqueVariants[filter_0.001 == "pass" &
                                  type %in% c("missense") &
                                  !str_detect(hgvs_pro, fixed("?"))]

# Subset all variants
mergedVariants = merge(allVariants, uniqueVariants[, .(variant_name, type)],
                       by = "variant_name")

# Select unique entry ID (eid)
eids = sort(unique(mergedVariants$eid))

# Remove withdrawn participants
eids = eids[eids > 0] # Negative EID are assigned to withdrawn participants
withdraws = readLines(WITHDRAW_EIDs_PATH)
withdraws = as.numeric(withdraws)
eids = eids[!eids %in% withdraws]

# Remove empty columns
emptycols = sapply(mergedVariants, function (k) all(is.na(k)))
mergedVariants = mergedVariants[, !emptycols, with = F]

# Save unique entry ID and the merged variants
outputPrefix = sprintf("output/%s/%s_", toupper(GENE), GENE)
writeLines(as.character(eids), paste0(outputPrefix, "unique_eids_filtered.txt"))
fwrite(mergedVariants, paste0(outputPrefix, "variants_filtered.csv"))
