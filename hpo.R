 
#!/usr/bin/env Rscript

# ============================================================================
# HPO direct phenotype -> gene downloader
#
# FINAL BENCHMARK VERSION
#
# Official HPO sources:
#   - hp.obo
#   - genes_to_phenotype.txt
#
# Resolution:
#   1. exact HPO ID
#   2. exact normalized HPO primary term
#   3. exact normalized official HPO synonym
#
# No:
#   - hardcoded phenotype synonyms
#   - hardcoded genes
#   - disease-name expansion
#   - phenotype.hpoa -> many HPO features -> gene expansion
#
# Outputs:
#   AllPackagesGenes/<phenotype>_hpo.csv
#   AllPackagesGenes/<phenotype>_hpo_genes.csv
#   AllPackagesGenes/<phenotype>_hpo_resolution.csv
#
# Usage:
#   Rscript hpo.R migraine
#   Rscript hpo.R "gastro-oesophageal reflux"
#   Rscript hpo.R HP:0002076
# ============================================================================


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

OUTPUT_DIR <- "AllPackagesGenes"

CACHE_MAX_AGE_DAYS <- 7

HPO_GENE_URL <- paste0(
  "https://purl.obolibrary.org/obo/",
  "hp/hpoa/genes_to_phenotype.txt"
)

HPO_OBO_URL <- paste0(
  "https://purl.obolibrary.org/obo/",
  "hp.obo"
)

HPO_GENE_FILE <- "hpo_genes_to_phenotype.txt"
HPO_OBO_FILE  <- "hpo_ontology.obo"


# ----------------------------------------------------------------------------
# Packages
# ----------------------------------------------------------------------------

if (!requireNamespace("httr", quietly = TRUE)) {
  stop(
    "Required R package 'httr' is not installed."
  )
}

suppressPackageStartupMessages(
  library(httr)
)


# ----------------------------------------------------------------------------
# Basic helpers
# ----------------------------------------------------------------------------

safe_trim <- function(x) {

  x <- as.character(x)

  x[is.na(x)] <- ""

  trimws(x)
}


clean_filename <- function(x) {

  x <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    x
  )

  x <- gsub(
    "_+",
    "_",
    x
  )

  x <- gsub(
    "^_|_$",
    "",
    x
  )

  if (!nzchar(x)) {
    x <- "phenotype"
  }

  x
}


normalize_text <- function(x) {

  x <- safe_trim(x)

  x <- tolower(x)

  # Generic spelling normalization only.
  # These are not phenotype-specific aliases.

  x <- gsub(
    "oesoph",
    "esoph",
    x
  )

  x <- gsub(
    "tumour",
    "tumor",
    x
  )

  x <- gsub(
    "[^a-z0-9]+",
    " ",
    x
  )

  x <- gsub(
    "\\s+",
    " ",
    x
  )

  trimws(x)
}


file_is_fresh <- function(
  path,
  max_age_days = CACHE_MAX_AGE_DAYS
) {

  if (!file.exists(path)) {
    return(FALSE)
  }

  age <- as.numeric(
    difftime(
      Sys.time(),
      file.mtime(path),
      units = "days"
    )
  )

  (
    !is.na(age) &&
    age <= max_age_days
  )
}


# ----------------------------------------------------------------------------
# Download source files
# ----------------------------------------------------------------------------

download_file_if_needed <- function(
  url,
  destination
) {

  if (
    file_is_fresh(destination)
  ) {

    cat(
      "Using cached file:",
      destination,
      "\n"
    )

    return(destination)
  }


  cat(
    "Downloading:",
    url,
    "\n"
  )


  response <- tryCatch(

    httr::GET(
      url,
      httr::add_headers(
        `User-Agent` =
          "PhenotypeToGeneDownloaderR-HPO/3.0"
      ),
      httr::timeout(300)
    ),

    error = function(e) {

      cat(
        "access_failure:",
        conditionMessage(e),
        "\n"
      )

      NULL
    }
  )


  if (is.null(response)) {
    return(NULL)
  }


  if (
    httr::status_code(response) != 200
  ) {

    cat(
      "access_failure: HTTP",
      httr::status_code(response),
      "for",
      url,
      "\n"
    )

    return(NULL)
  }


  raw_data <- httr::content(
    response,
    as = "raw"
  )


  writeBin(
    raw_data,
    destination
  )


  if (
    !file.exists(destination) ||
    file.info(destination)$size < 100
  ) {

    cat(
      "access_failure: invalid downloaded file:",
      destination,
      "\n"
    )

    return(NULL)
  }


  cat(
    "Downloaded:",
    destination,
    "(",
    file.info(destination)$size,
    "bytes )\n"
  )


  destination
}


prepare_hpo_files <- function() {

  cat(
    "Checking official HPO files...\n"
  )


  gene_file <- download_file_if_needed(
    HPO_GENE_URL,
    HPO_GENE_FILE
  )


  ontology_file <- download_file_if_needed(
    HPO_OBO_URL,
    HPO_OBO_FILE
  )


  if (
    is.null(gene_file) ||
    !file.exists(gene_file)
  ) {

    stop(
      "access_failure: HPO genes_to_phenotype.txt unavailable."
    )
  }


  if (
    is.null(ontology_file) ||
    !file.exists(ontology_file)
  ) {

    stop(
      "access_failure: HPO ontology unavailable."
    )
  }


  list(
    genes = gene_file,
    ontology = ontology_file
  )
}


# ----------------------------------------------------------------------------
# Read genes_to_phenotype
# ----------------------------------------------------------------------------

load_gene_annotations <- function(path) {

  cat(
    "Loading HPO gene annotations...\n"
  )


  df <- tryCatch(

    read.delim(
      path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE
    ),

    error = function(e) {

      cat(
        "parsing_failure:",
        conditionMessage(e),
        "\n"
      )

      NULL
    }
  )


  if (
    is.null(df) ||
    nrow(df) == 0
  ) {

    stop(
      "parsing_failure: genes_to_phenotype contains no readable rows."
    )
  }


  colnames(df) <- sub(
    "^#",
    "",
    colnames(df)
  )


  required <- c(
    "ncbi_gene_id",
    "gene_symbol",
    "hpo_id",
    "hpo_name"
  )


  missing <- setdiff(
    required,
    colnames(df)
  )


  if (
    length(missing) > 0
  ) {

    stop(
      paste0(
        "parsing_failure: missing HPO gene columns: ",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }


  out <- data.frame(

    Gene_ID =
      safe_trim(
        df$ncbi_gene_id
      ),

    Gene =
      safe_trim(
        df$gene_symbol
      ),

    HPO_ID =
      safe_trim(
        df$hpo_id
      ),

    HPO_Term =
      safe_trim(
        df$hpo_name
      ),

    stringsAsFactors = FALSE
  )


  keep <- (
    nzchar(out$Gene) &
    !grepl(
      "^[0-9]+$",
      out$Gene
    ) &
    grepl(
      "^HP:[0-9]{7}$",
      out$HPO_ID
    ) &
    nzchar(out$HPO_Term)
  )


  out <- out[
    keep,
    ,
    drop = FALSE
  ]


  out <- out[
    !duplicated(
      paste(
        out$Gene,
        out$HPO_ID,
        sep = "|"
      )
    ),
    ,
    drop = FALSE
  ]


  rownames(out) <- NULL


  cat(
    "Gene-HPO associations:",
    nrow(out),
    "\n"
  )

  cat(
    "Unique genes:",
    length(
      unique(out$Gene)
    ),
    "\n"
  )

  cat(
    "Unique HPO terms:",
    length(
      unique(out$HPO_ID)
    ),
    "\n"
  )


  out
}


# ----------------------------------------------------------------------------
# Parse official HPO ontology
# ----------------------------------------------------------------------------

extract_synonym_text <- function(line) {

  m <- regexec(
    '^synonym:[[:space:]]+"([^"]+)"',
    line
  )

  hit <- regmatches(
    line,
    m
  )[[1]]


  if (
    length(hit) >= 2
  ) {

    return(
      hit[2]
    )
  }


  ""
}


parse_hpo_ontology <- function(path) {

  cat(
    "Parsing HPO ontology labels and synonyms...\n"
  )


  lines <- tryCatch(

    readLines(
      path,
      warn = FALSE,
      encoding = "UTF-8"
    ),

    error = function(e) {

      cat(
        "parsing_failure:",
        conditionMessage(e),
        "\n"
      )

      character()
    }
  )


  if (
    length(lines) == 0
  ) {

    stop(
      "parsing_failure: HPO ontology is empty."
    )
  }


  term_starts <- which(
    trimws(lines) == "[Term]"
  )


  if (
    length(term_starts) == 0
  ) {

    stop(
      "parsing_failure: no [Term] entries found in HPO ontology."
    )
  }


  primary_rows <- list()
  synonym_rows <- list()
  alt_rows <- list()

  p_idx <- 1
  s_idx <- 1
  a_idx <- 1


  for (
    i in seq_along(
      term_starts
    )
  ) {

    start <- term_starts[i] + 1


    if (
      i < length(term_starts)
    ) {

      end <- term_starts[i + 1] - 1

    } else {

      end <- length(lines)
    }


    block <- lines[
      start:end
    ]


    next_section <- which(
      grepl(
        "^\\[",
        trimws(block)
      )
    )


    if (
      length(next_section) > 0
    ) {

      first_section <- next_section[1]

      if (
        first_section > 1
      ) {

        block <- block[
          seq_len(
            first_section - 1
          )
        ]

      } else {

        next
      }
    }


    obsolete <- any(
      grepl(
        "^is_obsolete:[[:space:]]*true",
        block
      )
    )


    if (obsolete) {
      next
    }


    id_line <- grep(
      "^id:[[:space:]]*HP:[0-9]{7}",
      block,
      value = TRUE
    )


    name_line <- grep(
      "^name:[[:space:]]*",
      block,
      value = TRUE
    )


    if (
      length(id_line) == 0 ||
      length(name_line) == 0
    ) {

      next
    }


    hpo_id <- sub(
      "^id:[[:space:]]*",
      "",
      id_line[1]
    )


    hpo_name <- sub(
      "^name:[[:space:]]*",
      "",
      name_line[1]
    )


    hpo_id <- safe_trim(
      hpo_id
    )

    hpo_name <- safe_trim(
      hpo_name
    )


    if (
      !grepl(
        "^HP:[0-9]{7}$",
        hpo_id
      ) ||
      !nzchar(
        hpo_name
      )
    ) {

      next
    }


    primary_rows[[p_idx]] <- data.frame(

      HPO_ID = hpo_id,

      HPO_Name = hpo_name,

      Normalized_Name =
        normalize_text(
          hpo_name
        ),

      stringsAsFactors = FALSE
    )

    p_idx <- p_idx + 1


    # ------------------------------------------------------------------------
    # Official HPO synonyms
    # ------------------------------------------------------------------------

    syn_lines <- grep(
      "^synonym:[[:space:]]*",
      block,
      value = TRUE
    )


    if (
      length(syn_lines) > 0
    ) {

      for (
        syn_line in syn_lines
      ) {

        synonym <- extract_synonym_text(
          syn_line
        )


        if (
          nzchar(
            synonym
          )
        ) {

          synonym_rows[[s_idx]] <- data.frame(

            HPO_ID = hpo_id,

            HPO_Name = hpo_name,

            Synonym = synonym,

            Normalized_Synonym =
              normalize_text(
                synonym
              ),

            stringsAsFactors = FALSE
          )

          s_idx <- s_idx + 1
        }
      }
    }


    # ------------------------------------------------------------------------
    # Alternate HPO identifiers
    # ------------------------------------------------------------------------

    alt_lines <- grep(
      "^alt_id:[[:space:]]*HP:[0-9]{7}",
      block,
      value = TRUE
    )


    if (
      length(alt_lines) > 0
    ) {

      for (
        alt_line in alt_lines
      ) {

        alt_id <- sub(
          "^alt_id:[[:space:]]*",
          "",
          alt_line
        )


        alt_id <- safe_trim(
          alt_id
        )


        if (
          grepl(
            "^HP:[0-9]{7}$",
            alt_id
          )
        ) {

          alt_rows[[a_idx]] <- data.frame(

            Alt_ID = alt_id,

            HPO_ID = hpo_id,

            HPO_Name = hpo_name,

            stringsAsFactors = FALSE
          )

          a_idx <- a_idx + 1
        }
      }
    }
  }


  if (
    length(primary_rows) == 0
  ) {

    stop(
      "parsing_failure: HPO ontology produced no usable terms."
    )
  }


  primary <- do.call(
    rbind,
    primary_rows
  )


  synonyms <- if (
    length(synonym_rows) > 0
  ) {

    do.call(
      rbind,
      synonym_rows
    )

  } else {

    data.frame(
      HPO_ID = character(),
      HPO_Name = character(),
      Synonym = character(),
      Normalized_Synonym = character(),
      stringsAsFactors = FALSE
    )
  }


  alt_ids <- if (
    length(alt_rows) > 0
  ) {

    do.call(
      rbind,
      alt_rows
    )

  } else {

    data.frame(
      Alt_ID = character(),
      HPO_ID = character(),
      HPO_Name = character(),
      stringsAsFactors = FALSE
    )
  }


  primary <- primary[
    !duplicated(
      primary$HPO_ID
    ),
    ,
    drop = FALSE
  ]


  synonyms <- synonyms[
    !duplicated(
      paste(
        synonyms$HPO_ID,
        synonyms$Normalized_Synonym,
        sep = "|"
      )
    ),
    ,
    drop = FALSE
  ]


  alt_ids <- alt_ids[
    !duplicated(
      alt_ids$Alt_ID
    ),
    ,
    drop = FALSE
  ]


  rownames(primary) <- NULL
  rownames(synonyms) <- NULL
  rownames(alt_ids) <- NULL


  cat(
    "HPO primary terms:",
    nrow(primary),
    "\n"
  )

  cat(
    "HPO official synonyms:",
    nrow(synonyms),
    "\n"
  )

  cat(
    "HPO alternate IDs:",
    nrow(alt_ids),
    "\n"
  )


  list(
    primary = primary,
    synonyms = synonyms,
    alt_ids = alt_ids
  )
}


# ----------------------------------------------------------------------------
# Safe HPO term resolution
# ----------------------------------------------------------------------------

resolution_result <- function(
  accepted,
  reason,
  hpo_id = "",
  hpo_name = "",
  matched_text = ""
) {

  list(
    accepted = accepted,
    reason = reason,
    HPO_ID = hpo_id,
    HPO_Name = hpo_name,
    Matched_Text = matched_text
  )
}


resolve_hpo_term <- function(
  phenotype,
  ontology
) {

  raw_query <- safe_trim(
    phenotype
  )


  normalized_query <- normalize_text(
    phenotype
  )


  # --------------------------------------------------------------------------
  # 1. Exact current HPO ID
  # --------------------------------------------------------------------------

  if (
    grepl(
      "^HP:[0-9]{7}$",
      toupper(
        raw_query
      )
    )
  ) {

    query_id <- toupper(
      raw_query
    )


    current <- ontology$primary[
      ontology$primary$HPO_ID ==
        query_id,
      ,
      drop = FALSE
    ]


    if (
      nrow(current) == 1
    ) {

      return(
        resolution_result(
          TRUE,
          "HPO_ID_Exact",
          current$HPO_ID[1],
          current$HPO_Name[1],
          query_id
        )
      )
    }


    alt <- ontology$alt_ids[
      ontology$alt_ids$Alt_ID ==
        query_id,
      ,
      drop = FALSE
    ]


    if (
      nrow(alt) == 1
    ) {

      return(
        resolution_result(
          TRUE,
          "HPO_Alt_ID_Exact",
          alt$HPO_ID[1],
          alt$HPO_Name[1],
          query_id
        )
      )
    }


    return(
      resolution_result(
        FALSE,
        "UNKNOWN_HPO_ID"
      )
    )
  }


  # --------------------------------------------------------------------------
  # 2. Exact normalized primary HPO term
  # --------------------------------------------------------------------------

  exact_primary <- ontology$primary[
    ontology$primary$Normalized_Name ==
      normalized_query,
    ,
    drop = FALSE
  ]


  exact_primary <- exact_primary[
    !duplicated(
      exact_primary$HPO_ID
    ),
    ,
    drop = FALSE
  ]


  if (
    nrow(exact_primary) == 1
  ) {

    return(
      resolution_result(
        TRUE,
        "HPO_Term_Exact",
        exact_primary$HPO_ID[1],
        exact_primary$HPO_Name[1],
        exact_primary$HPO_Name[1]
      )
    )
  }


  if (
    nrow(exact_primary) > 1
  ) {

    return(
      resolution_result(
        FALSE,
        "AMBIGUOUS_HPO_PRIMARY_TERM"
      )
    )
  }


  # --------------------------------------------------------------------------
  # 3. Exact normalized official HPO synonym
  # --------------------------------------------------------------------------

  exact_synonym <- ontology$synonyms[
    ontology$synonyms$Normalized_Synonym ==
      normalized_query,
    ,
    drop = FALSE
  ]


  exact_synonym <- exact_synonym[
    !duplicated(
      exact_synonym$HPO_ID
    ),
    ,
    drop = FALSE
  ]


  if (
    nrow(exact_synonym) == 1
  ) {

    return(
      resolution_result(
        TRUE,
        "HPO_Synonym_Exact",
        exact_synonym$HPO_ID[1],
        exact_synonym$HPO_Name[1],
        exact_synonym$Synonym[1]
      )
    )
  }


  if (
    nrow(exact_synonym) > 1
  ) {

    return(
      resolution_result(
        FALSE,
        "AMBIGUOUS_HPO_SYNONYM"
      )
    )
  }


  # --------------------------------------------------------------------------
  # IMPORTANT:
  #
  # Do NOT automatically use substring/pattern/disease expansion here.
  #
  # A failure to map exactly is recorded as unsupported_query.
  # --------------------------------------------------------------------------

  resolution_result(
    FALSE,
    "NO_EXACT_HPO_TERM_OR_SYNONYM"
  )
}


# ----------------------------------------------------------------------------
# Resolution audit
# ----------------------------------------------------------------------------

save_resolution <- function(
  phenotype,
  resolution,
  path
) {

  audit <- data.frame(

    Input_Term =
      phenotype,

    Accepted =
      isTRUE(
        resolution$accepted
      ),

    Resolution_Method =
      resolution$reason,

    Matched_HPO_ID =
      resolution$HPO_ID,

    Matched_HPO_Name =
      resolution$HPO_Name,

    Matched_Text =
      resolution$Matched_Text,

    Retrieval_Timestamp =
      format(
        Sys.time(),
        "%Y-%m-%dT%H:%M:%S%z"
      ),

    stringsAsFactors = FALSE
  )


  write.csv(
    audit,
    path,
    row.names = FALSE
  )
}


# ----------------------------------------------------------------------------
# Retrieve direct gene-HPO annotations
# ----------------------------------------------------------------------------

retrieve_genes_for_hpo <- function(
  phenotype,
  resolution,
  gene_annotations
) {

  if (
    !isTRUE(
      resolution$accepted
    )
  ) {

    return(
      data.frame()
    )
  }


  result <- gene_annotations[
    gene_annotations$HPO_ID ==
      resolution$HPO_ID,
    ,
    drop = FALSE
  ]


  if (
    nrow(result) == 0
  ) {

    return(
      data.frame()
    )
  }


  # --------------------------------------------------------------------------
  # Every returned row is a DIRECT gene-HPO annotation.
  # --------------------------------------------------------------------------

  result$Source <- (
    resolution$reason
  )


  result$Match_Pattern <- (
    resolution$Matched_Text
  )


  result$Input_Term <- (
    phenotype
  )


  result$Matched_Term <- (
    resolution$HPO_Name
  )


  result$Matched_Identifier <- (
    resolution$HPO_ID
  )


  result$Matched_Ontology <- (
    "HPO"
  )


  result$Evidence_Class <- (
    "curated_clinical"
  )


  result$Association_Type <- (
    "gene_hpo_annotation"
  )


  result$Query_Method <- (
    resolution$reason
  )


  result$Source_Rank <- seq_len(
    nrow(result)
  )


  result$Retrieval_Timestamp <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z"
  )


  result <- result[
    ,
    c(
      "Gene",
      "HPO_ID",
      "HPO_Term",
      "Gene_ID",
      "Source",
      "Match_Pattern",
      "Input_Term",
      "Matched_Term",
      "Matched_Identifier",
      "Matched_Ontology",
      "Evidence_Class",
      "Association_Type",
      "Query_Method",
      "Source_Rank",
      "Retrieval_Timestamp"
    ),
    drop = FALSE
  ]


  result <- result[
    !duplicated(
      paste(
        result$Gene,
        result$HPO_ID,
        sep = "|"
      )
    ),
    ,
    drop = FALSE
  ]


  result <- result[
    order(
      result$Gene
    ),
    ,
    drop = FALSE
  ]


  result$Source_Rank <- seq_len(
    nrow(result)
  )


  rownames(result) <- NULL


  result
}


# ----------------------------------------------------------------------------
# Stale output removal
# ----------------------------------------------------------------------------

remove_stale_outputs <- function(
  phenotype
) {

  clean <- clean_filename(
    phenotype
  )


  files <- c(

    file.path(
      OUTPUT_DIR,
      paste0(
        clean,
        "_hpo.csv"
      )
    ),

    file.path(
      OUTPUT_DIR,
      paste0(
        clean,
        "_hpo_genes.csv"
      )
    ),

    file.path(
      OUTPUT_DIR,
      paste0(
        clean,
        "_hpo_resolution.csv"
      )
    )
  )


  for (
    f in files
  ) {

    if (
      file.exists(f)
    ) {

      file.remove(f)

      cat(
        "Removed stale HPO output:",
        f,
        "\n"
      )
    }
  }
}


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

main <- function() {

  args <- commandArgs(
    trailingOnly = TRUE
  )


  if (
    length(args) < 1
  ) {

    cat(
      "HPO Gene Downloader\n"
    )

    cat(
      "Usage: Rscript hpo.R <phenotype>\n"
    )

    cat(
      "Example: Rscript hpo.R migraine\n"
    )

    quit(
      status = 1
    )
  }


  phenotype <- args[1]


  if (
    !dir.exists(
      OUTPUT_DIR
    )
  ) {

    dir.create(
      OUTPUT_DIR,
      recursive = TRUE
    )
  }


  clean <- clean_filename(
    phenotype
  )


  full_file <- file.path(
    OUTPUT_DIR,
    paste0(
      clean,
      "_hpo.csv"
    )
  )


  genes_file <- file.path(
    OUTPUT_DIR,
    paste0(
      clean,
      "_hpo_genes.csv"
    )
  )


  resolution_file <- file.path(
    OUTPUT_DIR,
    paste0(
      clean,
      "_hpo_resolution.csv"
    )
  )


  cat(
    "\n============================================================\n"
  )

  cat(
    "HPO DIRECT STRUCTURED RETRIEVAL\n"
  )

  cat(
    "============================================================\n"
  )

  cat(
    "Input term:",
    phenotype,
    "\n"
  )

  cat(
    "Start time:",
    format(
      Sys.time()
    ),
    "\n\n"
  )


  remove_stale_outputs(
    phenotype
  )


  files <- tryCatch(

    prepare_hpo_files(),

    error = function(e) {

      cat(
        conditionMessage(e),
        "\n"
      )

      NULL
    }
  )


  if (
    is.null(files)
  ) {

    quit(
      status = 2
    )
  }


  gene_annotations <- tryCatch(

    load_gene_annotations(
      files$genes
    ),

    error = function(e) {

      cat(
        conditionMessage(e),
        "\n"
      )

      NULL
    }
  )


  if (
    is.null(gene_annotations)
  ) {

    quit(
      status = 3
    )
  }


  ontology <- tryCatch(

    parse_hpo_ontology(
      files$ontology
    ),

    error = function(e) {

      cat(
        conditionMessage(e),
        "\n"
      )

      NULL
    }
  )


  if (
    is.null(ontology)
  ) {

    quit(
      status = 3
    )
  }


  # --------------------------------------------------------------------------
  # Resolve query
  # --------------------------------------------------------------------------

  resolution <- resolve_hpo_term(
    phenotype,
    ontology
  )


  save_resolution(
    phenotype,
    resolution,
    resolution_file
  )


  cat(
    "\nResolution:",
    resolution$reason,
    "\n"
  )


  if (
    isTRUE(
      resolution$accepted
    )
  ) {

    cat(
      "Matched HPO term:",
      resolution$HPO_Name,
      "\n"
    )

    cat(
      "Matched HPO ID:",
      resolution$HPO_ID,
      "\n"
    )

    cat(
      "Matched text:",
      resolution$Matched_Text,
      "\n"
    )
  }


  if (
    !isTRUE(
      resolution$accepted
    )
  ) {

    cat(
      "\nunsupported_query: input could not be resolved to an exact ",
      "HPO term, official HPO synonym or HPO identifier.\n",
      sep = ""
    )

    cat(
      "Resolution audit:",
      resolution_file,
      "\n"
    )

    cat(
      "End time:",
      format(
        Sys.time()
      ),
      "\n"
    )

    quit(
      status = 0
    )
  }


  # --------------------------------------------------------------------------
  # Direct gene-HPO associations
  # --------------------------------------------------------------------------

  results <- retrieve_genes_for_hpo(
    phenotype,
    resolution,
    gene_annotations
  )


  if (
    nrow(results) == 0
  ) {

    cat(
      "\nsuccessful_no_results: HPO concept resolved successfully ",
      "but has no gene annotations.\n",
      sep = ""
    )

    cat(
      "Resolution audit:",
      resolution_file,
      "\n"
    )

    cat(
      "End time:",
      format(
        Sys.time()
      ),
      "\n"
    )

    quit(
      status = 0
    )
  }


  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------

  write.csv(
    results,
    full_file,
    row.names = FALSE
  )


  genes_only <- data.frame(

    Gene = sort(
      unique(
        results$Gene
      )
    ),

    stringsAsFactors = FALSE
  )


  write.csv(
    genes_only,
    genes_file,
    row.names = FALSE
  )


  cat(
    "\n============================================================\n"
  )

  cat(
    "SUCCESS\n"
  )

  cat(
    "============================================================\n"
  )

  cat(
    "Resolution method:",
    resolution$reason,
    "\n"
  )

  cat(
    "HPO term:",
    resolution$HPO_Name,
    "\n"
  )

  cat(
    "HPO ID:",
    resolution$HPO_ID,
    "\n"
  )

  cat(
    "Direct associations:",
    nrow(results),
    "\n"
  )

  cat(
    "Unique genes:",
    nrow(genes_only),
    "\n"
  )

  cat(
    "Full output:",
    full_file,
    "\n"
  )

  cat(
    "Genes output:",
    genes_file,
    "\n"
  )

  cat(
    "Resolution audit:",
    resolution_file,
    "\n"
  )


  cat(
    "\nGenes:\n"
  )


  preview <- head(
    genes_only$Gene,
    100
  )


  cat(
    paste(
      preview,
      collapse = ", "
    ),
    "\n"
  )


  cat(
    "\nEnd time:",
    format(
      Sys.time()
    ),
    "\n"
  )
}


if (
  !interactive()
) {

  main()
}
 