
#!/usr/bin/env Rscript

# ============================================================================
# OMIM / NCBI MedGen structured gene downloader
#
# Preferred path:
#   1. OMIM API, if an OMIM_API_KEY is available
#
# Structured fallback:
#   2. NCBI mim2gene_medgen
#   3. NCBI MedGen_HPO_OMIM_Mapping.txt.gz
#   4. NCBI Homo_sapiens.gene_info.gz
#
# IMPORTANT:
#   This script DOES NOT scrape OMIM HTML.
#
# The NCBI fallback is labelled:
#   OMIM_DERIVED_MEDGEN_GENE_MAP
#
# because the records are accessed through NCBI MedGen/Gene rather than
# directly through the OMIM API.
#
# Usage:
#   Rscript omim.R migraine
#
# Optional:
#   export OMIM_API_KEY="..."
#   Rscript omim.R migraine
#
# Force refresh of NCBI files:
#   Rscript omim.R migraine --refresh
# ============================================================================


required_packages <- c(
  "httr",
  "jsonlite"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages(
  library(httr)
)

suppressPackageStartupMessages(
  library(jsonlite)
)


`%||%` <- function(x, y) {
  if (
    is.null(x) ||
    length(x) == 0 ||
    all(is.na(x))
  ) {
    y
  } else {
    x
  }
}


# ============================================================================
# Configuration
# ============================================================================

OMIM_API_BASE <- "https://api.omim.org/api"

MIM2GENE_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/",
  "gene/DATA/mim2gene_medgen"
)

OMIM_HPO_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/",
  "pub/medgen/MedGen_HPO_OMIM_Mapping.txt.gz"
)

GENE_INFO_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/",
  "gene/DATA/GENE_INFO/Mammalia/",
  "Homo_sapiens.gene_info.gz"
)

CACHE_DIR <- "omim_cache"

if (!dir.exists(CACHE_DIR)) {
  dir.create(
    CACHE_DIR,
    recursive = TRUE
  )
}

CACHE_MAX_AGE_DAYS <- 7


# ============================================================================
# CLI
# ============================================================================

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) < 1) {

  cat(
    paste0(
      "\nOMIM structured gene downloader\n\n",
      "Usage:\n",
      "  Rscript omim.R <phenotype>\n",
      "  Rscript omim.R <phenotype> --refresh\n\n",
      "Optional OMIM API:\n",
      "  export OMIM_API_KEY='your_key'\n",
      "  Rscript omim.R migraine\n\n"
    )
  )

  quit(
    status = 1
  )
}


flags <- args[
  grepl(
    "^--",
    args
  )
]

positional <- args[
  !grepl(
    "^--",
    args
  )
]

phenotype <- positional[1]

refresh_requested <- (
  "--refresh"
  %in%
  flags
)

api_key <- Sys.getenv(
  "OMIM_API_KEY",
  unset = ""
)

if (
  length(positional) >= 2 &&
  nzchar(positional[2])
) {
  api_key <- positional[2]
}


# ============================================================================
# Utilities
# ============================================================================

normalise_text <- function(x) {

  x <- as.character(x)

  x[is.na(x)] <- ""

  x <- tolower(x)

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


clean_phenotype_filename <- function(x) {

  x <- gsub(
    "[^[:alnum:]_ -]",
    "",
    x
  )

  x <- gsub(
    "\\s+",
    "_",
    trimws(x)
  )

  x
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


download_file_safe <- function(
  url,
  destination,
  refresh = FALSE
) {

  if (
    !refresh &&
    file_is_fresh(destination)
  ) {

    cat(
      "   Using cached:",
      destination,
      "\n"
    )

    return(destination)
  }

  cat(
    "   Downloading:",
    url,
    "\n"
  )

  response <- tryCatch(
    GET(
      url,
      add_headers(
        `User-Agent` =
          "PhenotypeToGeneDownloaderR/OMIM-MedGen"
      ),
      timeout(300)
    ),
    error = function(e) {
      cat(
        "   Download error:",
        conditionMessage(e),
        "\n"
      )
      NULL
    }
  )

  if (is.null(response)) {
    return(NULL)
  }

  if (status_code(response) != 200) {

    cat(
      "   HTTP status:",
      status_code(response),
      "\n"
    )

    return(NULL)
  }

  raw_data <- content(
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
      "   Invalid downloaded file:",
      destination,
      "\n"
    )

    return(NULL)
  }

  cat(
    "   Saved:",
    destination,
    "(",
    file.info(destination)$size,
    "bytes )\n"
  )

  destination
}


find_existing_or_cache <- function(
  root_name,
  cache_name
) {

  candidates <- c(
    root_name,
    file.path(
      CACHE_DIR,
      cache_name
    )
  )

  for (candidate in candidates) {

    if (
      file.exists(candidate) &&
      file.info(candidate)$size > 100
    ) {

      return(candidate)
    }
  }

  NULL
}


# ============================================================================
# Bulk data acquisition
# ============================================================================

prepare_bulk_files <- function(
  refresh = FALSE
) {

  cat(
    "\nPreparing structured OMIM/MedGen files...\n"
  )

  # --------------------------------------------------------------------------
  # mim2gene_medgen
  # --------------------------------------------------------------------------

  root_mim <- (
    if (file.exists("mim2gene_medgen"))
      "mim2gene_medgen"
    else
      NULL
  )

  mim_file <- root_mim

  if (
    is.null(mim_file) ||
    refresh
  ) {

    mim_file <- download_file_safe(
      MIM2GENE_URL,
      file.path(
        CACHE_DIR,
        "mim2gene_medgen"
      ),
      refresh = refresh
    )
  }

  if (
    is.null(mim_file) ||
    !file.exists(mim_file)
  ) {

    stop(
      "mim2gene_medgen is unavailable."
    )
  }


  # --------------------------------------------------------------------------
  # OMIM-HPO mapping
  # --------------------------------------------------------------------------

  root_map <- (
    if (
      file.exists(
        "MedGen_HPO_OMIM_Mapping.txt.gz"
      )
    )
      "MedGen_HPO_OMIM_Mapping.txt.gz"
    else
      NULL
  )

  mapping_file <- root_map

  if (
    is.null(mapping_file) ||
    refresh
  ) {

    mapping_file <- download_file_safe(
      OMIM_HPO_URL,
      file.path(
        CACHE_DIR,
        "MedGen_HPO_OMIM_Mapping.txt.gz"
      ),
      refresh = refresh
    )
  }

  if (
    is.null(mapping_file) ||
    !file.exists(mapping_file)
  ) {

    stop(
      paste(
        "MedGen_HPO_OMIM_Mapping.txt.gz",
        "is unavailable."
      )
    )
  }


  # --------------------------------------------------------------------------
  # Human gene_info
  # --------------------------------------------------------------------------

  root_gene_info <- (
    if (
      file.exists(
        "Homo_sapiens.gene_info.gz"
      )
    )
      "Homo_sapiens.gene_info.gz"
    else
      NULL
  )

  gene_info_file <- root_gene_info

  if (
    is.null(gene_info_file) ||
    refresh
  ) {

    gene_info_file <- download_file_safe(
      GENE_INFO_URL,
      file.path(
        CACHE_DIR,
        "Homo_sapiens.gene_info.gz"
      ),
      refresh = refresh
    )
  }

  if (
    is.null(gene_info_file) ||
    !file.exists(gene_info_file)
  ) {

    stop(
      "Homo_sapiens.gene_info.gz is unavailable."
    )
  }


  list(
    mim2gene = mim_file,
    mapping = mapping_file,
    gene_info = gene_info_file
  )
}


# ============================================================================
# Parsers
# ============================================================================

load_mim2gene <- function(path) {

  cat(
    "Loading mim2gene_medgen...\n"
  )

  df <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE
  )

  if (ncol(df) < 5) {

    stop(
      "Unexpected mim2gene_medgen format."
    )
  }

  colnames(df)[1] <- "MIM_Number"

  expected <- c(
    "MIM_Number",
    "GeneID",
    "type",
    "Source",
    "MedGenCUI"
  )

  missing <- setdiff(
    expected,
    colnames(df)
  )

  if (length(missing) > 0) {

    stop(
      paste(
        "mim2gene_medgen missing columns:",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }


  df$MIM_Number <- trimws(
    as.character(
      df$MIM_Number
    )
  )

  df$GeneID <- trimws(
    as.character(
      df$GeneID
    )
  )

  df$type <- trimws(
    as.character(
      df$type
    )
  )

  df$Source <- trimws(
    as.character(
      df$Source
    )
  )


  # We want phenotype -> gene relationships.
  #
  # GeneMap is retained deliberately. This avoids treating ordinary
  # OMIM gene records as disease associations.

  keep <- (
    tolower(df$type) ==
      "phenotype"
  ) &
    df$GeneID != "-" &
    nzchar(df$GeneID) &
    df$Source == "GeneMap"


  df <- df[
    keep,
    ,
    drop = FALSE
  ]


  df <- df[
    !duplicated(
      paste(
        df$MIM_Number,
        df$GeneID,
        sep = "|"
      )
    ),
    ,
    drop = FALSE
  ]


  cat(
    "   Structured OMIM phenotype-gene relationships:",
    nrow(df),
    "\n"
  )

  df
}


load_omim_mapping <- function(path) {

  cat(
    "Loading MedGen HPO/OMIM mapping...\n"
  )

  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  df <- read.delim(
    con,
    header = TRUE,
    sep = "|",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE
  )


  # First field starts with '#'.
  if (
    "#OMIM_CUI"
    %in%
    colnames(df)
  ) {

    colnames(df)[
      colnames(df) ==
        "#OMIM_CUI"
    ] <- "OMIM_CUI"
  }


  required <- c(
    "OMIM_CUI",
    "MIM_number",
    "OMIM_name",
    "relationship",
    "HPO_ID",
    "HPO_name"
  )


  missing <- setdiff(
    required,
    colnames(df)
  )

  if (length(missing) > 0) {

    stop(
      paste(
        "OMIM mapping missing columns:",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }


  df$MIM_number <- trimws(
    as.character(
      df$MIM_number
    )
  )

  df$OMIM_name <- trimws(
    as.character(
      df$OMIM_name
    )
  )

  df$HPO_name <- trimws(
    as.character(
      df$HPO_name
    )
  )

  df$OMIM_name_norm <- normalise_text(
    df$OMIM_name
  )

  df$HPO_name_norm <- normalise_text(
    df$HPO_name
  )


  cat(
    "   OMIM/HPO mapping rows:",
    nrow(df),
    "\n"
  )

  df
}


load_gene_info <- function(path) {

  cat(
    "Loading NCBI human gene_info...\n"
  )

  con <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(con),
    add = TRUE
  )

  df <- read.delim(
    con,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE
  )


  if (
    "#tax_id"
    %in%
    colnames(df)
  ) {

    colnames(df)[
      colnames(df) ==
        "#tax_id"
    ] <- "tax_id"
  }


  required <- c(
    "GeneID",
    "Symbol"
  )

  missing <- setdiff(
    required,
    colnames(df)
  )

  if (length(missing) > 0) {

    stop(
      paste(
        "gene_info missing columns:",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }


  if (
    "tax_id"
    %in%
    colnames(df)
  ) {

    df <- df[
      as.character(df$tax_id) ==
        "9606",
      ,
      drop = FALSE
    ]
  }


  df$GeneID <- as.character(
    df$GeneID
  )

  df$Symbol <- trimws(
    as.character(
      df$Symbol
    )
  )


  keep_columns <- intersect(
    c(
      "GeneID",
      "Symbol",
      "description",
      "Synonyms"
    ),
    colnames(df)
  )


  df <- df[
    ,
    keep_columns,
    drop = FALSE
  ]


  df <- df[
    !duplicated(
      df$GeneID
    ),
    ,
    drop = FALSE
  ]


  cat(
    "   Human genes loaded:",
    nrow(df),
    "\n"
  )

  df
}


# ============================================================================
# Structured phenotype matching
# ============================================================================

find_direct_omim_matches <- function(
  phenotype,
  mapping
) {

  query <- normalise_text(
    phenotype
  )


  if (!nzchar(query)) {

    return(
      data.frame()
    )
  }


  # --------------------------------------------------------------------------
  # Tier 1: exact OMIM disease-name match
  # --------------------------------------------------------------------------

  exact <- mapping[
    mapping$OMIM_name_norm ==
      query,
    ,
    drop = FALSE
  ]


  if (nrow(exact) > 0) {

    exact$Match_Type <- (
      "OMIM_DISEASE_NAME_EXACT"
    )

    exact$Match_Score <- 1.00
  }


  # --------------------------------------------------------------------------
  # Tier 2: OMIM disease name contains complete query phrase.
  #
  # Example:
  # migraine
  # -> FAMILIAL HEMIPLEGIC MIGRAINE
  #
  # This is still disease-name matching, not HPO-feature expansion.
  # --------------------------------------------------------------------------

  query_pattern <- paste0(
    "(^| )",
    gsub(
      "([][{}()+*^$|\\\\?.])",
      "\\\\\\1",
      query
    ),
    "( |$)"
  )


  contains_idx <- grepl(
    query_pattern,
    mapping$OMIM_name_norm,
    perl = TRUE
  )


  contains <- mapping[
    contains_idx,
    ,
    drop = FALSE
  ]


  if (nrow(contains) > 0) {

    contains$Match_Type <- (
      "OMIM_DISEASE_NAME_CONTAINS"
    )

    contains$Match_Score <- 0.90
  }


  # --------------------------------------------------------------------------
  # Combine exact + contains only.
  #
  # We intentionally DO NOT use an HPO feature match in the default OMIM
  # genes list because that would make the OMIM module partly dependent
  # on HPO annotations and would contaminate leave-source-out validation.
  # --------------------------------------------------------------------------

  pieces <- list()

  if (nrow(exact) > 0) {
    pieces[[length(pieces) + 1]] <- exact
  }

  if (nrow(contains) > 0) {
    pieces[[length(pieces) + 1]] <- contains
  }


  if (length(pieces) == 0) {

    return(
      data.frame()
    )
  }


  out <- do.call(
    rbind,
    pieces
  )


  out <- out[
    order(
      -out$Match_Score,
      out$OMIM_name,
      out$MIM_number
    ),
    ,
    drop = FALSE
  ]


  # If same MIM appears through exact and contains,
  # retain strongest match.

  out <- out[
    !duplicated(
      out$MIM_number
    ),
    ,
    drop = FALSE
  ]


  rownames(out) <- NULL

  out
}


# ============================================================================
# Join OMIM phenotype MIM numbers to GeneMap relationships
# ============================================================================

build_medgen_results <- function(
  phenotype,
  mapping,
  mim2gene,
  gene_info
) {

  matches <- find_direct_omim_matches(
    phenotype,
    mapping
  )


  if (nrow(matches) == 0) {

    cat(
      "   No direct OMIM disease-name matches.\n"
    )

    return(
      data.frame()
    )
  }


  cat(
    "   Matched OMIM disorders:",
    nrow(matches),
    "\n"
  )


  preview_n <- min(
    10,
    nrow(matches)
  )


  for (
    i in seq_len(preview_n)
  ) {

    cat(
      "      ",
      matches$MIM_number[i],
      " | ",
      matches$OMIM_name[i],
      " | ",
      matches$Match_Type[i],
      "\n",
      sep = ""
    )
  }


  relationships <- merge(
    matches,
    mim2gene,
    by.x = "MIM_number",
    by.y = "MIM_Number",
    all = FALSE
  )


  if (nrow(relationships) == 0) {

    cat(
      "   Matched OMIM disorders had no GeneMap relationships.\n"
    )

    return(
      data.frame()
    )
  }


  relationships$GeneID <- as.character(
    relationships$GeneID
  )


  relationships <- merge(
    relationships,
    gene_info,
    by = "GeneID",
    all.x = TRUE
  )


  relationships <- relationships[
    !is.na(
      relationships$Symbol
    ) &
      nzchar(
        relationships$Symbol
      ),
    ,
    drop = FALSE
  ]


  if (nrow(relationships) == 0) {

    cat(
      "   Gene IDs could not be mapped to human symbols.\n"
    )

    return(
      data.frame()
    )
  }


  retrieval_time <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z"
  )


  gene_name <- rep(
    "",
    nrow(relationships)
  )


  if (
    "description"
    %in%
    colnames(relationships)
  ) {

    gene_name <- as.character(
      relationships$description
    )
  }


  result <- data.frame(

    Gene_Symbol =
      relationships$Symbol,

    Gene_Name =
      gene_name,

    MIM_Number =
      relationships$MIM_number,

    Entry_Title =
      relationships$OMIM_name,

    Search_Term =
      phenotype,

    Phenotype_Map =
      relationships$OMIM_name,

    Inheritance =
      "",

    Entrez_GeneIDs =
      relationships$GeneID,

    Ensembl_IDs =
      "",

    Source =
      "OMIM_DERIVED_MEDGEN_GENE_MAP",

    input_term =
      phenotype,

    matched_term =
      relationships$OMIM_name,

    matched_identifier =
      paste0(
        "OMIM:",
        relationships$MIM_number
      ),

    matched_ontology =
      "OMIM",

    source_release =
      paste0(
        "NCBI_MedGen_",
        format(
          Sys.Date(),
          "%Y-%m-%d"
        )
      ),

    source_record_id =
      relationships$MIM_number,

    gene_symbol_original =
      relationships$Symbol,

    gene_symbol_current =
      relationships$Symbol,

    gene_identifier =
      paste0(
        "NCBIGene:",
        relationships$GeneID
      ),

    evidence_class =
      "curated_clinical",

    association_type =
      "OMIM_gene_disease_map",

    source_score =
      relationships$Match_Score,

    source_rank =
      NA_integer_,

    query_method =
      relationships$Match_Type,

    retrieval_timestamp =
      retrieval_time,

    MedGen_CUI =
      relationships$MedGenCUI,

    OMIM_CUI =
      relationships$OMIM_CUI,

    HPO_ID_context =
      relationships$HPO_ID,

    HPO_name_context =
      relationships$HPO_name,

    stringsAsFactors = FALSE
  )


  result <- result[
    order(
      -result$source_score,
      result$Entry_Title,
      result$Gene_Symbol
    ),
    ,
    drop = FALSE
  ]


  result <- result[
    !duplicated(
      paste(
        result$Gene_Symbol,
        result$MIM_Number,
        sep = "|"
      )
    ),
    ,
    drop = FALSE
  ]


  # Source rank after deterministic ordering.

  result$source_rank <- seq_len(
    nrow(result)
  )


  rownames(result) <- NULL


  cat(
    "   Structured gene-disease rows:",
    nrow(result),
    "\n"
  )

  cat(
    "   Unique genes:",
    length(
      unique(
        result$Gene_Symbol
      )
    ),
    "\n"
  )


  result
}


# ============================================================================
# Optional OMIM API path
# ============================================================================

safe_text <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0
  ) {
    return("")
  }

  as.character(x)
}


omim_get <- function(
  path,
  query = list(),
  api_key
) {

  url <- paste0(
    OMIM_API_BASE,
    path
  )


  response <- GET(
    url,
    query = c(
      query,
      list(
        format = "json"
      )
    ),
    add_headers(
      ApiKey = api_key,
      `Accept-Encoding` = "gzip",
      `User-Agent` =
        "PhenotypeToGeneDownloaderR-OMIM/2.0"
    ),
    timeout(60)
  )


  if (
    status_code(response) != 200
  ) {

    stop(
      "OMIM API HTTP ",
      status_code(response)
    )
  }


  content(
    response,
    as = "parsed",
    type = "application/json",
    simplifyVector = FALSE
  )
}


extract_api_rows <- function(
  wrapper,
  phenotype
) {

  entry <- wrapper$entry %||% NULL

  if (is.null(entry)) {
    return(list())
  }


  mim <- entry$mimNumber %||% NA

  title <- (
    entry$titles$preferredTitle
    %||%
    ""
  )


  maps <- (
    entry$geneMapList
    %||%
    list()
  )


  rows <- list()


  for (gm_wrapper in maps) {

    gm <- (
      gm_wrapper$geneMap
      %||%
      NULL
    )

    if (is.null(gm)) {
      next
    }


    symbols <- safe_text(
      gm$approvedGeneSymbols
    )


    if (!nzchar(symbols)) {

      symbols <- safe_text(
        gm$geneSymbols
      )
    }


    if (!nzchar(symbols)) {
      next
    }


    symbols <- unique(
      trimws(
        unlist(
          strsplit(
            symbols,
            ","
          )
        )
      )
    )


    symbols <- symbols[
      nzchar(symbols)
    ]


    for (gene in symbols) {

      rows[[length(rows) + 1]] <- data.frame(

        Gene_Symbol =
          gene,

        Gene_Name =
          safe_text(
            gm$geneName
          ),

        MIM_Number =
          as.character(mim),

        Entry_Title =
          as.character(title),

        Search_Term =
          phenotype,

        Phenotype_Map =
          as.character(title),

        Inheritance =
          "",

        Entrez_GeneIDs =
          safe_text(
            gm$geneIDs
          ),

        Ensembl_IDs =
          safe_text(
            gm$ensemblIDs
          ),

        Source =
          "OMIM_API",

        input_term =
          phenotype,

        matched_term =
          as.character(title),

        matched_identifier =
          paste0(
            "OMIM:",
            mim
          ),

        matched_ontology =
          "OMIM",

        source_release =
          paste0(
            "OMIM_API_",
            format(
              Sys.Date(),
              "%Y-%m-%d"
            )
          ),

        source_record_id =
          as.character(mim),

        gene_symbol_original =
          gene,

        gene_symbol_current =
          gene,

        gene_identifier =
          safe_text(
            gm$geneIDs
          ),

        evidence_class =
          "curated_clinical",

        association_type =
          "OMIM_gene_disease_map",

        source_score =
          1,

        source_rank =
          NA_integer_,

        query_method =
          "OMIM_API_entry_search",

        retrieval_timestamp =
          format(
            Sys.time(),
            "%Y-%m-%dT%H:%M:%S%z"
          ),

        MedGen_CUI =
          "",

        OMIM_CUI =
          "",

        HPO_ID_context =
          "",

        HPO_name_context =
          "",

        stringsAsFactors = FALSE
      )
    }
  }


  rows
}


run_api_mode <- function(
  phenotype,
  api_key
) {

  if (!nzchar(api_key)) {

    cat(
      "No OMIM API key available; using structured NCBI fallback.\n"
    )

    return(
      data.frame()
    )
  }


  cat(
    "Trying OMIM API first...\n"
  )


  result <- tryCatch({

    response <- omim_get(
      "/entry/search",
      query = list(
        search = phenotype,
        include = "geneMap",
        sort = "score desc",
        start = 0,
        limit = 100
      ),
      api_key = api_key
    )


    entries <- (
      response$
        omim$
        searchResponse$
        entryList
      %||%
      list()
    )


    rows <- list()


    for (entry in entries) {

      extracted <- extract_api_rows(
        entry,
        phenotype
      )

      if (length(extracted) > 0) {

        rows <- c(
          rows,
          extracted
        )
      }
    }


    if (length(rows) == 0) {

      data.frame()

    } else {

      out <- do.call(
        rbind,
        rows
      )


      out <- out[
        !duplicated(
          paste(
            out$Gene_Symbol,
            out$MIM_Number,
            sep = "|"
          )
        ),
        ,
        drop = FALSE
      ]


      out$source_rank <- seq_len(
        nrow(out)
      )


      out
    }

  }, error = function(e) {

    cat(
      "OMIM API failed:",
      conditionMessage(e),
      "\n"
    )

    data.frame()
  })


  result
}


# ============================================================================
# Save
# ============================================================================

save_results <- function(
  df,
  phenotype
) {

  output_dir <- "AllPackagesGenes"

  if (!dir.exists(output_dir)) {

    dir.create(
      output_dir,
      recursive = TRUE
    )
  }


  clean <- clean_phenotype_filename(
    phenotype
  )


  full_file <- file.path(
    output_dir,
    paste0(
      clean,
      "_omim.csv"
    )
  )


  genes_file <- file.path(
    output_dir,
    paste0(
      clean,
      "_omim_genes.csv"
    )
  )


  write.csv(
    df,
    full_file,
    row.names = FALSE
  )


  genes <- sort(
    unique(
      df$Gene_Symbol
    )
  )


  gene_df <- data.frame(
    Gene = genes,
    stringsAsFactors = FALSE
  )


  write.csv(
    gene_df,
    genes_file,
    row.names = FALSE
  )


  list(
    full_file = full_file,
    genes_file = genes_file,
    genes = gene_df
  )
}


remove_stale_outputs <- function(
  phenotype
) {

  clean <- clean_phenotype_filename(
    phenotype
  )


  targets <- c(

    file.path(
      "AllPackagesGenes",
      paste0(
        clean,
        "_omim.csv"
      )
    ),

    file.path(
      "AllPackagesGenes",
      paste0(
        clean,
        "_omim_genes.csv"
      )
    )
  )


  for (path in targets) {

    if (file.exists(path)) {

      file.remove(path)

      cat(
        "Removed stale OMIM output:",
        path,
        "\n"
      )
    }
  }
}


# ============================================================================
# Main
# ============================================================================

main <- function() {

  cat(
    "\n============================================================\n"
  )

  cat(
    "OMIM STRUCTURED GENE RETRIEVAL\n"
  )

  cat(
    "============================================================\n"
  )

  cat(
    "Input term:       ",
    phenotype,
    "\n"
  )

  cat(
    "OMIM API key:     ",
    ifelse(
      nzchar(api_key),
      "available",
      "not available"
    ),
    "\n"
  )

  cat(
    "Refresh bulk data:",
    refresh_requested,
    "\n"
  )

  cat(
    "Start time:       ",
    format(Sys.time()),
    "\n\n"
  )


  # Important: old HTML-scraper outputs must not survive.

  remove_stale_outputs(
    phenotype
  )


  # --------------------------------------------------------------------------
  # 1. API
  # --------------------------------------------------------------------------

  result <- run_api_mode(
    phenotype,
    api_key
  )


  # --------------------------------------------------------------------------
  # 2. Structured NCBI fallback
  # --------------------------------------------------------------------------

  if (nrow(result) == 0) {

    cat(
      "\nUsing NCBI structured OMIM/MedGen fallback.\n"
    )


    files <- prepare_bulk_files(
      refresh =
        refresh_requested
    )


    mim2gene <- load_mim2gene(
      files$mim2gene
    )


    mapping <- load_omim_mapping(
      files$mapping
    )


    gene_info <- load_gene_info(
      files$gene_info
    )


    result <- build_medgen_results(
      phenotype,
      mapping,
      mim2gene,
      gene_info
    )
  }


  # --------------------------------------------------------------------------
  # Outcome
  # --------------------------------------------------------------------------

  if (nrow(result) == 0) {

    cat(
      "\nNo structured OMIM gene-disease association ",
      "was found for: ",
      phenotype,
      "\n",
      sep = ""
    )

    cat(
      "This is a successful no-result query, ",
      "not an access or parsing failure.\n"
    )

    cat(
      "\nEnd time:",
      format(Sys.time()),
      "\n"
    )

    quit(
      status = 0
    )
  }


  saved <- save_results(
    result,
    phenotype
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
    "Rows:         ",
    nrow(result),
    "\n"
  )

  cat(
    "Unique genes: ",
    nrow(saved$genes),
    "\n"
  )

  cat(
    "Source:       ",
    paste(
      unique(
        result$Source
      ),
      collapse = ", "
    ),
    "\n"
  )

  cat(
    "Full output:  ",
    saved$full_file,
    "\n"
  )

  cat(
    "Genes output: ",
    saved$genes_file,
    "\n"
  )


  cat(
    "\nTop genes:\n"
  )


  preview <- head(
    saved$genes$Gene,
    30
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
    format(Sys.time()),
    "\n"
  )
}


main()
 