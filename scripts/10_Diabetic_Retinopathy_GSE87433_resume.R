# Portable repository copy; start R in the repository root.
# This is a processed MICROARRAY precheck: no additional log2 transformation,
# CPM calculation or TMM normalization is applied to the supplied values.
# DR row labels are Entrez IDs; supplied platform annotation links them to gene symbols.
# Original working files are preserved; GitHub sync did not execute this code.
# Capture the project folder before the original workflow changes working directory.
dr_project_dir <- normalizePath(getwd(), mustWork = TRUE)
# Reload the DR checkpoint after its save block has returned TRUE.
# Use this script when returning instead of rerunning the full precheck.
dr_resume_file <- file.path(dr_project_dir, "results", "diabetic_retinopathy", "DR_GSE87433_precheck_checkpoint_2026-10-31.RData")
load(dr_resume_file)
setwd(dr_output_dir)

# Rebuild the XML objects from the extracted metadata file.
dr_metadata_xml <- xml2::read_xml(dr_metadata_file)
dr_sample_nodes <- xml2::xml_find_all(
  dr_metadata_xml, "/geo:MINiML/geo:Sample", ns = dr_ns
)
dr_platform_nodes <- xml2::xml_find_all(
  dr_metadata_xml, "/geo:MINiML/geo:Platform", ns = dr_ns
)
dr_characteristic_nodes <- xml2::xml_find_all(
  dr_sample_nodes, ".//geo:Characteristics", ns = dr_ns
)
dr_characteristic_tags <- xml2::xml_attr(dr_characteristic_nodes, "tag")

# Restore the familiar reader for individual characteristic fields.
dr_get_characteristic <- function(tag) {
  dr_path <- paste0(".//geo:Characteristics[@tag='", tag, "']")
  xml2::xml_text(
    xml2::xml_find_first(dr_sample_nodes, dr_path, ns = dr_ns),
    trim = TRUE
  )
}

# Restore the plotting helpers and display the aligned table dimensions.
source(file.path(dr_project_dir, "scripts", "11_DR_projection_style.R"))
dim(dr_expr_all)
dim(dr_metadata_matched)
dim(dr_gene_annotation)
