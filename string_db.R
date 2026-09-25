#!/usr/bin/env Rscript

# STRING network-expansion module
#
# STRING is a protein-protein interaction resource, not a direct
# phenotype-to-gene database.
#
# Workflow:
# phenotype
#   -> curated phenotype-associated seed genes
#      from ClinVar / HPO / OMIM
#   -> STRING protein resolution
#   -> STRING interaction partners
#
# Usage:
#   Rscript string_db.R migraine
#
# Outputs:
#   AllPackagesGenes/<phenotype>_string_db.csv
#   AllPackagesGenes/<phenotype>_string_db_genes.csv

required_packages <- c("httr", "jsonlite")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package '", pkg,
      "' is not installed. Please install it before running this script."
    )
  }
}

suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(jsonlite))


clean_filename <- function(x) {
  x <- gsub("[^a-zA-Z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)

  if (!nzchar(x)) {
    x <- "phenotype"
  }

  x
}


empty_result <- function() {
  data.frame(
    Phenotype = character(),
    Gene = character(),
    Seed_Gene = character(),
    Seed_Sources = character(),
    Seed_Source_Count = integer(),
    Interaction_Score = numeric(),
    STRING_ID = character(),
    Source = character(),
    Evidence_Class = character(),
    Association_Type = character(),
    Query_Method = character(),
    stringsAsFactors = FALSE
  )
}


read_gene_file <- function(path, source_name) {

  if (!file.exists(path)) {
    cat("   ", source_name, ": file not found\n", sep = "")
    return(data.frame())
  }

  dat <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE),
    error = function(e) {
      cat(
        "   ",
        source_name,
        ": failed to read file - ",
        conditionMessage(e),
        "\n",
        sep = ""
      )
      NULL
    }
  )

  if (is.null(dat) || nrow(dat) == 0) {
    cat("   ", source_name, ": 0 genes\n", sep = "")
    return(data.frame())
  }

  candidate_columns <- c(
    "Gene",
    "GeneSymbol",
    "Gene_Symbol",
    "gene_symbol",
    "symbol"
  )

  gene_col <- candidate_columns[
    candidate_columns %in% colnames(dat)
  ]

  if (length(gene_col) == 0) {
    cat(
      "   ",
      source_name,
      ": no recognised gene column\n",
      sep = ""
    )
    return(data.frame())
  }

  genes <- trimws(
    as.character(
      dat[[gene_col[1]]]
    )
  )

  genes <- genes[
    !is.na(genes) &
    nzchar(genes)
  ]

  genes <- unique(genes)

  cat(
    "   ",
    source_name,
    ": ",
    length(genes),
    " genes\n",
    sep = ""
  )

  if (length(genes) == 0) {
    return(data.frame())
  }

  data.frame(
    Gene = genes,
    Seed_Source = source_name,
    stringsAsFactors = FALSE
  )
}


get_seed_genes <- function(
  phenotype,
  output_dir = "AllPackagesGenes",
  max_seeds = 10
) {

  clean_phenotype <- clean_filename(phenotype)

  source_files <- list(
    ClinVar = file.path(
      output_dir,
      paste0(clean_phenotype, "_clinvar_genes.csv")
    ),
    HPO = file.path(
      output_dir,
      paste0(clean_phenotype, "_hpo_genes.csv")
    ),
    OMIM = file.path(
      output_dir,
      paste0(clean_phenotype, "_omim_genes.csv")
    )
  )

  cat("🌱 Collecting curated phenotype-associated seed genes...\n")

  all_rows <- list()

  for (source_name in names(source_files)) {

    one <- read_gene_file(
      source_files[[source_name]],
      source_name
    )

    if (nrow(one) > 0) {
      all_rows[[length(all_rows) + 1]] <- one
    }
  }

  if (length(all_rows) == 0) {
    cat("   No seed genes available.\n")
    return(data.frame())
  }

  combined <- do.call(
    rbind,
    all_rows
  )

  all_genes <- sort(
    unique(combined$Gene)
  )

  ranked_list <- lapply(
    all_genes,
    function(gene) {

      sources <- sort(
        unique(
          combined$Seed_Source[
            combined$Gene == gene
          ]
        )
      )

      data.frame(
        Gene = gene,
        Seed_Sources = paste(
          sources,
          collapse = ";"
        ),
        Seed_Source_Count = length(sources),
        stringsAsFactors = FALSE
      )
    }
  )

  ranked <- do.call(
    rbind,
    ranked_list
  )

  ranked <- ranked[
    order(
      -ranked$Seed_Source_Count,
      ranked$Gene
    ),
    ,
    drop = FALSE
  ]

  if (nrow(ranked) > max_seeds) {
    ranked <- ranked[
      seq_len(max_seeds),
      ,
      drop = FALSE
    ]
  }

  rownames(ranked) <- NULL

  cat(
    "   ✅ Selected ",
    nrow(ranked),
    " seed genes\n",
    sep = ""
  )

  cat(
    "   Seeds: ",
    paste(
      ranked$Gene,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )

  ranked
}


resolve_string_seed <- function(seed_gene) {

  response <- GET(
    "https://string-db.org/api/json/get_string_ids",
    query = list(
      identifiers = seed_gene,
      species = 9606,
      limit = 1,
      caller_identity = "PhenotypeToGeneDownloaderR"
    ),
    timeout(30)
  )

  if (status_code(response) != 200) {
    return(NULL)
  }

  resolved <- fromJSON(
    content(
      response,
      "text",
      encoding = "UTF-8"
    ),
    simplifyVector = TRUE
  )

  if (
    length(resolved) == 0 ||
    !is.data.frame(resolved) ||
    nrow(resolved) == 0
  ) {
    return(NULL)
  }

  resolved[1, , drop = FALSE]
}


get_string_partners <- function(
  phenotype,
  seed_gene,
  seed_sources,
  seed_source_count,
  partner_limit = 50,
  required_score = 400
) {

  cat(
    "   Resolving seed protein:",
    seed_gene,
    "\n"
  )

  resolved <- tryCatch(
    resolve_string_seed(seed_gene),
    error = function(e) {
      cat(
        "      ⚠️ Resolution error:",
        conditionMessage(e),
        "\n"
      )
      NULL
    }
  )

  if (is.null(resolved)) {
    cat(
      "      ⚠️ Could not resolve seed\n"
    )
    return(empty_result())
  }

  resolved_name <- as.character(
    resolved$preferredName[1]
  )

  cat(
    "      Resolved as:",
    resolved_name,
    "\n"
  )

  response <- tryCatch(
    GET(
      "https://string-db.org/api/json/interaction_partners",
      query = list(
        identifiers = resolved_name,
        species = 9606,
        limit = partner_limit,
        required_score = required_score,
        caller_identity = "PhenotypeToGeneDownloaderR"
      ),
      timeout(30)
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(response) ||
    status_code(response) != 200
  ) {
    cat(
      "      ⚠️ Interaction request failed\n"
    )
    return(empty_result())
  }

  interactions <- fromJSON(
    content(
      response,
      "text",
      encoding = "UTF-8"
    ),
    simplifyVector = TRUE
  )

  if (
    length(interactions) == 0 ||
    !is.data.frame(interactions) ||
    nrow(interactions) == 0
  ) {
    cat(
      "      No interaction partners found\n"
    )
    return(empty_result())
  }

  valid <- (
    !is.na(interactions$preferredName_B) &
    interactions$preferredName_B != ""
  )

  if (!any(valid)) {
    return(empty_result())
  }

  batch <- data.frame(
    Phenotype = phenotype,
    Gene = as.character(
      interactions$preferredName_B[valid]
    ),
    Seed_Gene = seed_gene,
    Seed_Sources = seed_sources,
    Seed_Source_Count = as.integer(
      seed_source_count
    ),
    Interaction_Score = as.numeric(
      interactions$score[valid]
    ),
    STRING_ID = as.character(
      interactions$stringId_B[valid]
    ),
    Source = "STRING",
    Evidence_Class = "network_expansion",
    Association_Type = "indirect",
    Query_Method = "interaction_partners_from_curated_seed_gene",
    stringsAsFactors = FALSE
  )

  cat(
    "      Found ",
    nrow(batch),
    " partners\n",
    sep = ""
  )

  batch
}


download_string_genes <- function(
  phenotype,
  seed_table
) {

  cat(
    "\n🔗 STRING network expansion for phenotype:",
    phenotype,
    "\n"
  )

  if (
    is.null(seed_table) ||
    nrow(seed_table) == 0
  ) {
    cat(
      "   No phenotype-associated seed genes available\n"
    )
    return(empty_result())
  }

  all_results <- list()

  for (i in seq_len(nrow(seed_table))) {

    one <- get_string_partners(
      phenotype = phenotype,
      seed_gene = seed_table$Gene[i],
      seed_sources = seed_table$Seed_Sources[i],
      seed_source_count = seed_table$Seed_Source_Count[i],
      partner_limit = 50,
      required_score = 400
    )

    if (nrow(one) > 0) {
      all_results[[length(all_results) + 1]] <- one
    }

    Sys.sleep(0.5)
  }

  if (length(all_results) == 0) {
    return(empty_result())
  }

  result <- do.call(
    rbind,
    all_results
  )

  result <- result[
    order(
      -result$Interaction_Score,
      -result$Seed_Source_Count,
      result$Gene
    ),
    ,
    drop = FALSE
  ]

  result <- result[
    !duplicated(result$Gene),
    ,
    drop = FALSE
  ]

  rownames(result) <- NULL

  cat(
    "\n   ✅ Total unique STRING network genes:",
    nrow(result),
    "\n"
  )

  result
}


main <- function() {

  args <- commandArgs(
    trailingOnly = TRUE
  )

  if (length(args) < 1) {

    cat("STRING Network Expansion Module\n")
    cat("\n")
    cat("Usage:\n")
    cat("  Rscript string_db.R <phenotype>\n")
    cat("\n")
    cat("Example:\n")
    cat("  Rscript string_db.R migraine\n")

    quit(status = 1)
  }

  phenotype <- args[1]

  output_dir <- "AllPackagesGenes"

  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }

  clean_phenotype <- clean_filename(
    phenotype
  )

  output_file <- file.path(
    output_dir,
    paste0(
      clean_phenotype,
      "_string_db.csv"
    )
  )

  genes_only_file <- file.path(
    output_dir,
    paste0(
      clean_phenotype,
      "_string_db_genes.csv"
    )
  )

  cat(
    "🎯 Phenotype:         ",
    phenotype,
    "\n"
  )

  cat(
    "📁 Output directory:  ",
    output_dir,
    "\n"
  )

  cat(
    "🧬 STRING role:       network expansion / indirect evidence\n"
  )

  cat(
    "🌱 Maximum seeds:     10\n"
  )

  cat(
    "⏰ Start time:        ",
    format(Sys.time()),
    "\n\n"
  )

  seed_table <- get_seed_genes(
    phenotype = phenotype,
    output_dir = output_dir,
    max_seeds = 10
  )

  if (
    is.null(seed_table) ||
    nrow(seed_table) == 0
  ) {

    cat(
      "\n⚠️ SUCCESSFUL_NO_RESULTS\n"
    )

    cat(
      "No curated phenotype-associated seed genes were available.\n"
    )

    cat(
      "Phenotype text was NOT submitted directly to STRING.\n"
    )

    cat(
      "\n⏰ End time: ",
      format(Sys.time()),
      "\n"
    )

    quit(status = 0)
  }

  results <- download_string_genes(
    phenotype = phenotype,
    seed_table = seed_table
  )

  if (nrow(results) > 0) {

    write.csv(
      results,
      output_file,
      row.names = FALSE
    )

    genes_only <- data.frame(
      Gene = sort(
        unique(results$Gene)
      ),
      stringsAsFactors = FALSE
    )

    write.csv(
      genes_only,
      genes_only_file,
      row.names = FALSE
    )

    cat(
      "\n✅ SUCCESSFUL_WITH_RESULTS\n"
    )

    cat(
      "📊 Network records:  ",
      nrow(results),
      "\n"
    )

    cat(
      "🧬 Unique genes:     ",
      nrow(genes_only),
      "\n"
    )

    cat(
      "🌱 Seed genes used:  ",
      nrow(seed_table),
      "\n"
    )

    cat(
      "💾 Full results:     ",
      output_file,
      "\n"
    )

    cat(
      "🧬 Genes-only file:  ",
      genes_only_file,
      "\n"
    )

  } else {

    cat(
      "\n⚠️ SUCCESSFUL_NO_RESULTS\n"
    )

    cat(
      "Seed genes were available, but STRING returned no interaction partners.\n"
    )
  }

  cat(
    "\n⏰ End time: ",
    format(Sys.time()),
    "\n"
  )
}


if (!interactive()) {
  main()
}
