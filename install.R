

cat("=== NGS Bioinformatics Project - Package Installation ===\n")
cat("Installing required R packages...\n\n")

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      install.packages(pkg, dependencies = TRUE)
      if (require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("✓", pkg, "installed successfully\n")
      } else {
        cat("✗ Failed to install", pkg, "\n")
      }
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }
}

# Function to install Bioconductor packages
install_bioc_if_missing <- function(packages) {
  # Install BiocManager if not available
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing Bioconductor package", pkg, "...\n")
      BiocManager::install(pkg, dependencies = TRUE)
      if (require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("✓", pkg, "installed successfully\n")
      } else {
        cat("✗ Failed to install", pkg, "\n")
      }
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }
}

# CRAN packages
cran_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "data.table",
  "ggplot2",
  "pheatmap",
  "VennDiagram",
  "RColorBrewer",
  "corrplot",
  "pcaMethods",
  "ggrepel",
  "gridExtra",
  "cowplot",
  "knitr",
  "kableExtra",
  "stringr",
  "magrittr",
  "glue"
)

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",
  "Limma",
  "Rsubread",
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

# Install CRAN packages
cat("Installing CRAN packages...\n")
install_if_missing(cran_packages)

cat("\nInstalling Bioconductor packages...\n")
install_bioc_if_missing(bioc_packages)

# Check installation status
cat("\n=== Installation Status Check ===\n")

all_packages <- c(cran_packages, bioc_packages)
installed_packages <- c()
missing_packages <- c()

for (pkg in all_packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    installed_packages <- c(installed_packages, pkg)
  } else {
    missing_packages <- c(missing_packages, pkg)
  }
}

cat("Successfully installed:", length(installed_packages), "packages\n")
cat("Missing packages:", length(missing_packages), "packages\n")

if (length(missing_packages) > 0) {
  cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Please try installing them manually:\n")
  for (pkg in missing_packages) {
    if (pkg %in% bioc_packages) {
      cat("BiocManager::install('", pkg, "')\n", sep = "")
    } else {
      cat("install.packages('", pkg, "')\n", sep = "")
    }
  }
}

# Create package info file
cat("\nCreating package information file...\n")
package_info <- data.frame(
  Package = all_packages,
  Type = ifelse(all_packages %in% bioc_packages, "Bioconductor", "CRAN"),
  Installed = all_packages %in% installed_packages,
  Version = sapply(all_packages, function(pkg) {
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      as.character(packageVersion(pkg))
    } else {
      "Not installed"
    }
  })
)

# Save package info
write.csv(package_info, "package_installation_status.csv", row.names = FALSE)
cat("Package information saved to: package_installation_status.csv\n")

# Final status
cat("\n=== Final Installation Status ===\n")
if (length(missing_packages) == 0) {
  cat("✓ All packages installed successfully!\n")
  cat("You can now run the analysis scripts.\n")
} else {
  cat("⚠ Some packages could not be installed.\n")
  cat("Please check the missing packages and install them manually.\n")
}

cat("\nInstallation script completed.\n")
cat("Check 'package_installation_status.csv' for detailed information.\n")
