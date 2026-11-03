# Split supplied gene assignments into readable fields, retaining every feature.
# Field layout follows the manufacturer documentation and reviewed examples.
# Source: https://documents.thermofisher.com/TFS-Assets/LSG/manuals/tac_user_manual.pdf
amd_format_annotation <- function(annotation) {
  stopifnot(
    is.data.frame(annotation), nrow(annotation) > 0L,
    all(c("ID", "gene_assignment") %in% colnames(annotation))
  )

  truncation_marker <- "[WARNING: THIS FIELD TRUNCATED]"

  # Separate assignments; retain a row for features without an assignment.
  original_text <- as.character(annotation$gene_assignment)
  split_text <- original_text
  split_text[is.na(split_text) | trimws(split_text) == ""] <- "---"
  assignments <- strsplit(split_text, "///", fixed = TRUE)
  assignments <- lapply(assignments, function(x) if (length(x)) x else "---")
  source_row <- rep(seq_len(nrow(annotation)), lengths(assignments))
  records <- trimws(unlist(assignments, use.names = FALSE))
  source_truncated <- !is.na(original_text) &
    grepl(truncation_marker, original_text, fixed = TRUE)
  assignment_truncated <- grepl(truncation_marker, records, fixed = TRUE)

  # Split the five fields, converting placeholders to missing values.
  fields <- lapply(seq_along(records), function(i) {
    value <- records[i]
    if (is.na(value) || value %in% c("", "---")) {
      return(rep(NA_character_, 5L))
    }
    parts <- trimws(strsplit(value, "//", fixed = TRUE)[[1]])
    if (assignment_truncated[i]) {
      # Retain only complete fields before the warning; missing text is not inferred.
      warning_position <- which(grepl(truncation_marker, parts, fixed = TRUE))[1L]
      prefix <- parts[seq_len(warning_position - 1L)]
      if (length(prefix) > 5L) {
        stop(sprintf(
          "Feature %s has an unexpected truncated format. Review its original text.",
          annotation$ID[source_row[i]]
        ), call. = FALSE)
      }
      parts <- rep(NA_character_, 5L)
      parts[seq_along(prefix)] <- prefix
    } else if (length(parts) != 5L) {
      stop(sprintf(
        "Feature %s has an assignment with %d fields; expected 5. Review its original text.",
        annotation$ID[source_row[i]], length(parts)
      ), call. = FALSE)
    }
    parts[parts %in% c("", "---")] <- NA_character_
    parts
  })
  fields <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
  colnames(fields) <- c(
    "reference_id", "gene_symbol", "description", "cytoband", "entrez_id"
  )

  # Keep source links and original text alongside the readable columns.
  cbind(
    data.frame(
      ID = as.character(annotation$ID[source_row]),
      platform_row = source_row,
      assignment_index = sequence(lengths(assignments)),
      stringsAsFactors = FALSE
    ),
    fields,
    data.frame(
      annotation_truncated = source_truncated[source_row],
      assignment_truncated = assignment_truncated,
      gene_assignment_raw = original_text[source_row],
      assignment_raw = records,
      stringsAsFactors = FALSE
    )
  )
}
