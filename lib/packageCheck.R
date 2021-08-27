LIBRARIES = c("yaml", "data.table", "dplyr", "tidyr", "stringr", "tools",
              "ggtext", "glue", "forcats", "gridExtra", "fs", "Rmisc", "odbc",
              "RMariaDB", "reshape2", "logr", "yogiroc")

# Install needed libraries
cat("Checking missing libraries\n")
missingLibs = LIBRARIES[!LIBRARIES %in% installed.packages()[, "Package"]]
if (length(missingLibs)) {
  cat("Installing missing libraries:", paste0(missingLibs, collapse = ", "), "\n")
  install.packages(missingLibs)
} else {
  cat("All necessary libraries found.\n")
}
