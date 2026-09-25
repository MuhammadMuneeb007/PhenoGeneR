 
#!/usr/bin/env Rscript

# ============================================================================
# Open Targets Platform phenotype/disease/trait -> gene downloader
#
# Improvements for final benchmark:
#   - does NOT blindly accept the first Open Targets search hit
#   - performs explicit term-resolution validation
#   - rejects ambiguous or semantically weak matches
#   - does NOT manually map medication terms to diseases
#   - retains matched concept ID/name and provenance
#   - writes a query-resolution audit table
#   - distinguishes unsupported_query from successful_no_results
#   - removes stale outputs before execution
#
# Usage:
#   Rscript opentargets.R migraine
#   Rscript opentargets.R "blood pressure medication"
#
# Outputs:
#   AllPackagesGenes/<phenotype>_opentargets.csv
#   AllPackagesGenes/<phenotype>_opentargets_genes.csv
#   AllPackagesGenes/<phenotype>_opentargets_resolution.csv
# ============================================================================


# ----------------------------------------------------------------------------
# Packages
# ----------------------------------------------------------------------------

required_packages <- c(
  "httr",
  "jsonlite"
)

for (pkg in required_packages) {

  if (
    !requireNamespace(
      pkg,
      quietly = TRUE
    )
  ) {

    stop(
      "Required package '",
      pkg,
      "' is not installed."
    )
  }
}

suppressPackageStartupMessages(
  library(httr)
)

suppressPackageStartupMessages(
  library(jsonlite)
)


# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

OT_URL <- (
  "https://api.platform.opentargets.org/api/v4/graphql"
)

OUTPUT_DIR <- "AllPackagesGenes"


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

safe_text <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0 ||
    all(is.na(x))
  ) {

    return("")
  }

  trimws(
    as.character(x)[1]
  )
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

  x <- safe_text(x)

  x <- tolower(x)

  # Generic British/American spelling normalization.
  x <- gsub(
    "oesoph",
    "esoph",
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


comparison_text <- function(x) {

  x <- normalize_text(x)

  # Generic ontology-label suffixes should not determine whether
  # two concepts are equivalent.
  generic_words <- c(
    "disease",
    "disorder",
    "syndrome",
    "condition",
    "phenotype",
    "trait"
  )

  words <- unlist(
    strsplit(
      x,
      "\\s+"
    )
  )

  words <- words[
    nzchar(words)
  ]

  words <- words[
    !words %in% generic_words
  ]

  paste(
    words,
    collapse = " "
  )
}


tokenize <- function(x) {

  x <- comparison_text(x)

  if (!nzchar(x)) {
    return(character())
  }

  unique(
    unlist(
      strsplit(
        x,
        "\\s+"
      )
    )
  )
}


character_similarity <- function(
  a,
  b
) {

  a <- comparison_text(a)
  b <- comparison_text(b)

  if (
    !nzchar(a) ||
    !nzchar(b)
  ) {

    return(0)
  }

  if (a == b) {
    return(1)
  }

  distance <- as.numeric(
    adist(
      a,
      b
    )[1]
  )

  denominator <- max(
    nchar(a),
    nchar(b)
  )

  if (denominator == 0) {
    return(0)
  }

  max(
    0,
    1 - distance / denominator
  )
}


token_jaccard <- function(
  a,
  b
) {

  aa <- tokenize(a)
  bb <- tokenize(b)

  if (
    length(aa) == 0 ||
    length(bb) == 0
  ) {

    return(0)
  }

  union_set <- union(
    aa,
    bb
  )

  intersection_set <- intersect(
    aa,
    bb
  )

  length(intersection_set) /
    length(union_set)
}


query_token_coverage <- function(
  query,
  candidate
) {

  q <- tokenize(query)
  c <- tokenize(candidate)

  if (length(q) == 0) {
    return(0)
  }

  length(
    intersect(
      q,
      c
    )
  ) /
    length(q)
}


# ----------------------------------------------------------------------------
# GraphQL
# ----------------------------------------------------------------------------

run_query <- function(
  query,
  variables = NULL
) {

  response <- tryCatch(

    httr::POST(
      OT_URL,
      body = list(
        query = query,
        variables = variables
      ),
      encode = "json",
      httr::content_type_json(),
      httr::add_headers(
        `User-Agent` =
          "PhenotypeToGeneDownloaderR-OpenTargets/2.0"
      ),
      httr::timeout(90)
    ),

    error = function(e) {

      cat(
        "access_failure: Open Targets request failed:",
        conditionMessage(e),
        "\n"
      )

      return(NULL)
    }
  )

  if (is.null(response)) {
    return(NULL)
  }

  status <- httr::status_code(
    response
  )

  if (status != 200) {

    cat(
      "access_failure: Open Targets returned HTTP",
      status,
      "\n"
    )

    return(NULL)
  }

  parsed <- tryCatch(

    httr::content(
      response,
      as = "parsed",
      type = "application/json",
      simplifyVector = FALSE
    ),

    error = function(e) {

      cat(
        "parsing_failure: Open Targets JSON parsing failed:",
        conditionMessage(e),
        "\n"
      )

      NULL
    }
  )

  if (is.null(parsed)) {
    return(NULL)
  }

  if (
    !is.null(parsed$errors) &&
    length(parsed$errors) > 0
  ) {

    messages <- vapply(
      parsed$errors,
      function(x) {
        safe_text(x$message)
      },
      character(1)
    )

    cat(
      "execution_failure: Open Targets GraphQL error:",
      paste(
        messages,
        collapse = " | "
      ),
      "\n"
    )

    return(NULL)
  }

  parsed
}


# ----------------------------------------------------------------------------
# Search Open Targets concepts
# ----------------------------------------------------------------------------

search_disease_candidates <- function(
  phenotype
) {

  cat(
    "Searching Open Targets for:",
    phenotype,
    "\n"
  )

  query <- '
    query searchDisease($term: String!) {
      search(
        queryString: $term,
        entityNames: ["disease"]
      ) {
        hits {
          id
          name
        }
      }
    }
  '

  result <- run_query(
    query,
    list(
      term = phenotype
    )
  )

  if (is.null(result)) {
    return(NULL)
  }

  hits <- tryCatch(
    result$data$search$hits,
    error = function(e) NULL
  )

  if (
    is.null(hits) ||
    length(hits) == 0
  ) {

    return(
      data.frame(
        Candidate_Rank = integer(),
        Candidate_ID = character(),
        Candidate_Name = character(),
        Exact_Match = logical(),
        Character_Similarity = numeric(),
        Token_Jaccard = numeric(),
        Query_Token_Coverage = numeric(),
        Resolution_Score = numeric(),
        stringsAsFactors = FALSE
      )
    )
  }

  rows <- list()

  for (
    i in seq_along(hits)
  ) {

    hit <- hits[[i]]

    id <- safe_text(
      hit$id
    )

    name <- safe_text(
      hit$name
    )

    if (
      !nzchar(id) ||
      !nzchar(name)
    ) {

      next
    }

    exact <- (
      comparison_text(phenotype) ==
      comparison_text(name)
    )

    char_sim <- character_similarity(
      phenotype,
      name
    )

    jaccard <- token_jaccard(
      phenotype,
      name
    )

    coverage <- query_token_coverage(
      phenotype,
      name
    )

    score <- max(
      char_sim,
      jaccard
    )

    if (exact) {
      score <- 1
    }

    rows[[
      length(rows) + 1
    ]] <- data.frame(

      Candidate_Rank = i,

      Candidate_ID = id,

      Candidate_Name = name,

      Exact_Match = exact,

      Character_Similarity =
        round(
          char_sim,
          4
        ),

      Token_Jaccard =
        round(
          jaccard,
          4
        ),

      Query_Token_Coverage =
        round(
          coverage,
          4
        ),

      Resolution_Score =
        round(
          score,
          4
        ),

      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {

    return(
      data.frame()
    )
  }

  candidates <- do.call(
    rbind,
    rows
  )

  candidates <- candidates[
    order(
      -as.integer(
        candidates$Exact_Match
      ),
      -candidates$Resolution_Score,
      candidates$Candidate_Rank
    ),
    ,
    drop = FALSE
  ]

  rownames(
    candidates
  ) <- NULL

  candidates
}


# ----------------------------------------------------------------------------
# Safe concept resolution
# ----------------------------------------------------------------------------

resolve_disease_concept <- function(
  phenotype,
  candidates
) {

  if (
    is.null(candidates) ||
    nrow(candidates) == 0
  ) {

    return(
      list(
        accepted = FALSE,
        reason = "NO_SEARCH_HITS",
        candidate = NULL
      )
    )
  }


  # --------------------------------------------------------------------------
  # 1. Exact normalized concept match.
  # --------------------------------------------------------------------------

  exact <- candidates[
    candidates$Exact_Match,
    ,
    drop = FALSE
  ]

  if (nrow(exact) == 1) {

    return(
      list(
        accepted = TRUE,
        reason = "EXACT_NORMALIZED_MATCH",
        candidate = exact[1, , drop = FALSE]
      )
    )
  }


  if (nrow(exact) > 1) {

    return(
      list(
        accepted = FALSE,
        reason = "AMBIGUOUS_EXACT_MATCH",
        candidate = exact[1, , drop = FALSE]
      )
    )
  }


  # --------------------------------------------------------------------------
  # 2. Strong non-exact lexical match.
  #
  # Conditions deliberately conservative:
  #
  #   - character similarity >= 0.88
  #       OR
  #   - token Jaccard >= 0.80 AND all query tokens are represented
  #
  # We also require at least two meaningful query tokens for a non-exact
  # token-subset match. This prevents a generic one-word query such as
  # "migraine" from silently resolving to a specific migraine subtype when
  # no exact migraine concept exists.
  # --------------------------------------------------------------------------

  q_tokens <- tokenize(
    phenotype
  )

  acceptable <- candidates[
    (
      candidates$Character_Similarity >= 0.88
    )
    |
    (
      length(q_tokens) >= 2
      &
      candidates$Token_Jaccard >= 0.80
      &
      candidates$Query_Token_Coverage >= 1.00
    ),
    ,
    drop = FALSE
  ]


  if (nrow(acceptable) == 0) {

    return(
      list(
        accepted = FALSE,
        reason = "NO_SAFE_CONCEPT_MATCH",
        candidate = candidates[1, , drop = FALSE]
      )
    )
  }


  acceptable <- acceptable[
    order(
      -acceptable$Resolution_Score,
      acceptable$Candidate_Rank
    ),
    ,
    drop = FALSE
  ]


  best <- acceptable[
    1,
    ,
    drop = FALSE
  ]


  # --------------------------------------------------------------------------
  # Ambiguity check.
  #
  # If two non-exact concepts score almost identically, do not choose one
  # automatically.
  # --------------------------------------------------------------------------

  if (nrow(acceptable) >= 2) {

    second_score <- acceptable$
      Resolution_Score[2]

    best_score <- best$
      Resolution_Score[1]

    if (
      abs(
        best_score -
        second_score
      ) <= 0.03
    ) {

      return(
        list(
          accepted = FALSE,
          reason = "AMBIGUOUS_NONEXACT_MATCH",
          candidate = best
        )
      )
    }
  }


  list(
    accepted = TRUE,
    reason = "STRONG_LEXICAL_MATCH",
    candidate = best
  )
}


# ----------------------------------------------------------------------------
# Retrieve disease-target associations
# ----------------------------------------------------------------------------

get_genes <- function(
  disease_id,
  page_size = 1000
) {

  cat(
    "Retrieving Open Targets associated targets for:",
    disease_id,
    "\n"
  )

  all_rows <- list()

  page_index <- 0


  repeat {

    query <- '
      query diseaseAssociations(
        $efoId: String!,
        $index: Int!,
        $size: Int!
      ) {
        disease(efoId: $efoId) {
          associatedTargets(
            page: {
              index: $index,
              size: $size
            }
          ) {
            rows {
              target {
                id
                approvedSymbol
                approvedName
              }
              score
            }
          }
        }
      }
    '


    result <- run_query(
      query,
      list(
        efoId = disease_id,
        index = page_index,
        size = page_size
      )
    )


    if (is.null(result)) {

      return(NULL)
    }


    disease_object <- tryCatch(
      result$data$disease,
      error = function(e) NULL
    )


    if (is.null(disease_object)) {

      cat(
        "unsupported_query: resolved Open Targets ID has no disease object.\n"
      )

      return(
        data.frame()
      )
    }


    rows <- tryCatch(
      disease_object$
        associatedTargets$
        rows,
      error = function(e) NULL
    )


    if (
      is.null(rows) ||
      length(rows) == 0
    ) {

      break
    }


    all_rows <- c(
      all_rows,
      rows
    )


    cat(
      "  Page",
      page_index,
      "-",
      length(rows),
      "associations\n"
    )


    if (
      length(rows) <
      page_size
    ) {

      break
    }


    page_index <-
      page_index + 1


    Sys.sleep(
      0.25
    )
  }


  if (
    length(all_rows) == 0
  ) {

    return(
      data.frame()
    )
  }


  gene_rows <- list()


  for (
    i in seq_along(
      all_rows
    )
  ) {

    row <- all_rows[[i]]

    target <- row$target


    if (is.null(target)) {
      next
    }


    gene <- safe_text(
      target$approvedSymbol
    )

    target_id <- safe_text(
      target$id
    )

    gene_name <- safe_text(
      target$approvedName
    )


    score <- suppressWarnings(
      as.numeric(
        row$score
      )
    )


    if (!nzchar(gene)) {
      next
    }


    gene_rows[[
      length(gene_rows) + 1
    ]] <- data.frame(

      Gene = gene,

      TargetID = target_id,

      GeneName = gene_name,

      Score = score,

      stringsAsFactors = FALSE
    )
  }


  if (
    length(gene_rows) == 0
  ) {

    return(
      data.frame()
    )
  }


  genes <- do.call(
    rbind,
    gene_rows
  )


  # Keep strongest association if a target occurs more than once.

  genes <- genes[
    order(
      -genes$Score,
      genes$Gene
    ),
    ,
    drop = FALSE
  ]


  genes <- genes[
    !duplicated(
      genes$Gene
    ),
    ,
    drop = FALSE
  ]


  rownames(
    genes
  ) <- NULL


  genes
}


# ----------------------------------------------------------------------------
# Output helpers
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
        "_opentargets.csv"
      )
    ),

    file.path(
      OUTPUT_DIR,
      paste0(
        clean,
        "_opentargets_genes.csv"
      )
    ),

    file.path(
      OUTPUT_DIR,
      paste0(
        clean,
        "_opentargets_resolution.csv"
      )
    )
  )


  for (
    f in files
  ) {

    if (
      file.exists(f)
    ) {

      file.remove(
        f
      )

      cat(
        "Removed stale output:",
        f,
        "\n"
      )
    }
  }
}


save_resolution_audit <- function(
  phenotype,
  candidates,
  resolution,
  path
) {

  if (
    is.null(candidates) ||
    nrow(candidates) == 0
  ) {

    audit <- data.frame(

      Input_Term = phenotype,

      Candidate_Rank = NA_integer_,

      Candidate_ID = "",

      Candidate_Name = "",

      Exact_Match = FALSE,

      Character_Similarity = NA_real_,

      Token_Jaccard = NA_real_,

      Query_Token_Coverage = NA_real_,

      Resolution_Score = NA_real_,

      Accepted = FALSE,

      Resolution_Reason =
        resolution$reason,

      stringsAsFactors = FALSE
    )

  } else {

    audit <- candidates


    audit$Input_Term <- (
      phenotype
    )


    audit$Accepted <- FALSE


    if (
      isTRUE(
        resolution$accepted
      ) &&
      !is.null(
        resolution$candidate
      )
    ) {

      accepted_id <- (
        resolution$
          candidate$
          Candidate_ID[1]
      )


      accepted_name <- (
        resolution$
          candidate$
          Candidate_Name[1]
      )


      audit$Accepted <- (
        audit$Candidate_ID ==
          accepted_id
        &
        audit$Candidate_Name ==
          accepted_name
      )
    }


    audit$Resolution_Reason <- (
      resolution$reason
    )


    audit <- audit[
      ,
      c(
        "Input_Term",
        "Candidate_Rank",
        "Candidate_ID",
        "Candidate_Name",
        "Exact_Match",
        "Character_Similarity",
        "Token_Jaccard",
        "Query_Token_Coverage",
        "Resolution_Score",
        "Accepted",
        "Resolution_Reason"
      ),
      drop = FALSE
    ]
  }


  write.csv(
    audit,
    path,
    row.names = FALSE
  )
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
      "Open Targets Gene Downloader\n"
    )

    cat(
      "Usage: Rscript opentargets.R <phenotype>\n"
    )

    cat(
      "Example: Rscript opentargets.R migraine\n"
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
      "_opentargets.csv"
    )
  )


  genes_file <- file.path(
    OUTPUT_DIR,
    paste0(
      clean,
      "_opentargets_genes.csv"
    )
  )


  resolution_file <- file.path(
    OUTPUT_DIR,
    paste0(
      clean,
      "_opentargets_resolution.csv"
    )
  )


  cat(
    "\n============================================================\n"
  )

  cat(
    "OPEN TARGETS STRUCTURED RETRIEVAL\n"
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


  # --------------------------------------------------------------------------
  # Concept discovery
  # --------------------------------------------------------------------------

  candidates <- search_disease_candidates(
    phenotype
  )


  if (is.null(candidates)) {

    cat(
      "access_failure: Open Targets search could not be completed.\n"
    )

    quit(
      status = 2
    )
  }


  if (
    nrow(candidates) > 0
  ) {

    cat(
      "\nTop Open Targets candidates:\n"
    )


    preview_n <- min(
      10,
      nrow(candidates)
    )


    print(
      candidates[
        seq_len(preview_n),
        c(
          "Candidate_Rank",
          "Candidate_ID",
          "Candidate_Name",
          "Exact_Match",
          "Character_Similarity",
          "Token_Jaccard",
          "Query_Token_Coverage",
          "Resolution_Score"
        ),
        drop = FALSE
      ],
      row.names = FALSE
    )

  } else {

    cat(
      "No Open Targets search candidates returned.\n"
    )
  }


  # --------------------------------------------------------------------------
  # Safe resolution
  # --------------------------------------------------------------------------

  resolution <- resolve_disease_concept(
    phenotype,
    candidates
  )


  save_resolution_audit(
    phenotype,
    candidates,
    resolution,
    resolution_file
  )


  cat(
    "\nResolution:",
    resolution$reason,
    "\n"
  )


  if (
    !isTRUE(
      resolution$accepted
    )
  ) {

    cat(
      "unsupported_query: no sufficiently safe Open Targets ",
      "disease/trait concept could be assigned automatically.\n",
      sep = ""
    )


    if (
      !is.null(
        resolution$candidate
      )
    ) {

      cat(
        "Best candidate was:",
        resolution$
          candidate$
          Candidate_Name[1],
        "(",
        resolution$
          candidate$
          Candidate_ID[1],
        ")\n"
      )
    }


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


  matched <- (
    resolution$
      candidate
  )


  matched_id <- (
    matched$
      Candidate_ID[1]
  )


  matched_name <- (
    matched$
      Candidate_Name[1]
  )


  match_score <- (
    matched$
      Resolution_Score[1]
  )


  cat(
    "\nAccepted concept:",
    matched_name,
    "(",
    matched_id,
    ")\n"
  )


  cat(
    "Resolution method:",
    resolution$reason,
    "\n"
  )


  cat(
    "Resolution score:",
    match_score,
    "\n"
  )


  # --------------------------------------------------------------------------
  # Retrieve target associations
  # --------------------------------------------------------------------------

  genes <- get_genes(
    matched_id
  )


  if (is.null(genes)) {

    cat(
      "access_failure: target retrieval did not complete.\n"
    )

    quit(
      status = 2
    )
  }


  if (
    nrow(genes) == 0
  ) {

    cat(
      "\nsuccessful_no_results: Open Targets concept resolved ",
      "successfully but returned no associated genes.\n",
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
  # Add provenance
  # --------------------------------------------------------------------------

  genes$Source <- (
    "Open Targets"
  )


  genes$Input_Term <- (
    phenotype
  )


  genes$Matched_Term <- (
    matched_name
  )


  genes$Matched_Identifier <- (
    matched_id
  )


  id_prefix <- sub(
    "[:_].*$",
    "",
    matched_id
  )


  genes$Matched_Ontology <- (
    id_prefix
  )


  genes$Evidence_Class <- (
    "integrated_disease_association"
  )


  genes$Association_Type <- (
    "disease_target_association"
  )


  genes$Query_Method <- (
    resolution$reason
  )


  genes$Term_Resolution_Score <- (
    match_score
  )


  genes$Source_Rank <- seq_len(
    nrow(genes)
  )


  genes$Retrieval_Timestamp <- format(
    Sys.time(),
    "%Y-%m-%dT%H:%M:%S%z"
  )


  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------

  write.csv(
    genes,
    full_file,
    row.names = FALSE
  )


  genes_only <- data.frame(

    Gene = genes$Gene,

    stringsAsFactors = FALSE
  )


  genes_only <- genes_only[
    !duplicated(
      genes_only$Gene
    ),
    ,
    drop = FALSE
  ]


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
    "Matched concept:",
    matched_name,
    "\n"
  )

  cat(
    "Matched ID:",
    matched_id,
    "\n"
  )

  cat(
    "Total genes:",
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
    "\nTop 20 associated genes:\n"
  )


  preview <- head(
    genes[
      ,
      c(
        "Gene",
        "Score",
        "Source_Rank"
      ),
      drop = FALSE
    ],
    20
  )


  print(
    preview,
    row.names = FALSE
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
 