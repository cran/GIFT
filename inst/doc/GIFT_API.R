## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.width = 8, fig.height = 8)
# Packages --------------------------------------------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    library("GIFT")
    library("knitr")
    library("kableExtra")
  })
})

options(tinytex.verbose = TRUE)

## ----gottingen_logo, fig.show = 'hold', out.width = "20%", echo = FALSE-------
knitr::include_graphics("../man/figures/biodiv_gottingen_logo.png")
knitr::include_graphics("../man/figures/GIFT.png")

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  # GIFT_traits_raw => reference_traits traits_raw
#  # GIFT_traits_meta => traits_meta
#  # GIFT_traits => traits
#  # GIFT_taxonomy => taxonomy
#  # GIFT_taxgroup => taxonomy
#  # GIFT_species_lookup => names_matched names_matched_unique
#  # GIFT_species_distribution => species_distr overlap
#  # GIFT_species => species
#  # GIFT_richness => taxonomy, traits_cov, species_num
#  # GIFT_regions => regions
#  # GIFT_references => references
#  # GIFT_overlap => overlap_misc
#  # GIFT_no_overlap => overlap
#  # GIFT_lists => lists
#  # GIFT_env_meta_raster => env_raster, references_citavi
#  # GIFT_env_meta_misc => env_misc, references_citavi
#  # GIFT_env => geoentities_env_misc, geoentities_env_raster
#  # GIFT_checklists_raw => taxonomy, checklists
#  # GIFT_checklists_conditional => lists, taxonomy
#  # GIFT_checklists => lists, taxonomy, overlap
#  # GIFT_phylogeny => phylogeny

## ----echo = FALSE-------------------------------------------------------------
query_table <-
  data.frame(
    Query = c("checklists", "env_misc", "env_raster", "geoentities_env_misc",
              "geoentities_env_raster", "lists", "names_matched",
              "names_matched_unique", "overlap", "overlap_misc", "phylogeny",
              "references_citavi", "references", "reference_traits", "regions",
              "species", "species_distr", "species_num", "taxonomy", "traits",
              "traits_cov", "traits_meta", "traits_raw", "versions"),
    Arguments = c("listid (integer)\n, taxonid, namesmatched, filter", "", "",
                  "envvar",
                  "layername, sumstat", "", "genus, epithet", "genus, epithet",
                  "taxon, startat",
                  "layer", "taxon, startat, limit", "", "", "", "",
                  "startat, limit", "nameid",
                  "taxonid",
                  "", "traitid, biasref, biasderiv, startat, limit",
                  "traitid, taxonid", "", "traitid, deriv, biasderiv, refid",
                  ""
    ),
    R_function = c(
      "GIFT_checklists_raw()",  "GIFT_env_meta_misc()",
      "GIFT_env_meta_raster()", "GIFT_env()", "GIFT_env()",
      "GIFT_lists(), GIFT_checklists_conditional(), GIFT_checklists()",
      "GIFT_species_lookup()", "GIFT_species_lookup()",
      "GIFT_species_distribution(), GIFT_no_overlap(), GIFT_checklists()",
      "GIFT_overlap()",
      "GIFT_phylogeny()",
      "GIFT_env_meta_raster(), GIFT_env_meta_misc()",
      "GIFT_references()", "GIFT_traits_raw()", "GIFT_regions()",
      "GIFT_species()", "GIFT_species_distribution()", "GIFT_richness()",
      "GIFT_taxonomy(), GIFT_taxgroup(), GIFT_richness(),
      GIFT_checklists_raw(), GIFT_checklists(), GIFT_checklists_conditional()",
      "GIFT_traits()", "GIFT_richness()", "GIFT_traits_meta()",
      "GIFT_traits_raw()", "All GIFT functions"
    )
  )

kable(query_table, "html") %>%
  kable_styling(full_width = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  # Query for trait and chunks
#  # Default value for end is 10000
#  # https://gift.uni-goettingen.de/api/extended/index.php?query=traits&
#  # traitid=1.1.1&biasref=1&biasderiv=1&startat=100000&limit=10
#  
#  # To retrieve the geojson
#  paste0("https://gift.uni-goettingen.de/geojson/geojson_smaller",
#         ifelse(GIFT_version == "beta", "", GIFT_version), "/",
#         entity_ID[i], ".geojson")

## ----fig.cap = "", out.width = "100%",echo = FALSE----------------------------
knitr::include_graphics("../man/figures/GIFT_network_functions.png")

