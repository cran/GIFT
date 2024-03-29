#' Species list in GIFT
#'
#' Retrieve all name matching information for one taxonomic name. All results 
#' are returned, where the name is either found in the unstandardized or 
#' taxonomically standardized names.
#'
#' @param genus character string defining the genus name to be looked for.
#' 
#' @param epithet character string defining the specific epithet to be looked
#' for.
#' 
#' @param namesmatched Logical `FALSE` by default, set to `TRUE` if you want
#' to look for the species not only in the standardized species names but also 
#' in the original species names as they came in the original resources.
#' 
#' @template GIFT_version_api
#' 
#' @return
#' A data frame with 19 columns (or 24 if namesmatched = TRUE).
#'
#' @details Here is what each column refers to:
#' \emph{orig_ID} - Identification number of the species before taxonomic
#' harmonization\cr
#' \emph{orig_genus} - Genus before taxonomic harmonization\cr
#' \emph{name_ID} - Identification number of the genus before taxonomic
#' harmonization\cr
#' \emph{cf_genus}- Whether the genus name is uncertain\cr
#' \emph{genus}- Genus before taxonomic harmonization\cr
#' \emph{cf_species}- Whether the species' epithet is uncertain\cr
#' \emph{aff_species}- Species' epithet uncertain\cr
#' \emph{species_epithet}- Epithet of the species before taxonomic
#'  harmonization\cr
#' \emph{subtaxon}- Subtaxon of the species before taxonomic harmonization\cr
#' \emph{author}- Author who described the species (before taxonomic
#'  harmonization)\cr
#' \emph{matched}- Is the species matched in the taxonomic backbone\cr
#' \emph{epithetscore}- Matching score for the epithet\cr
#' \emph{overallscore}- Overall matching score for the species\cr
#' \emph{resolved}- Is the species name resolved in the taxonomic backbone\cr
#' \emph{synonym}- Is the species name a synonym in the taxonomic backbone\cr
#' \emph{matched_subtaxon}- Is the subtaxon matched in the taxonomic backbone\cr
#' \emph{accepted}- Is the species name accepted in the taxonomic backbone\cr
#' \emph{service}- Service use for the taxonomic harmonization\cr
#' \emph{work_ID}- Identification number of the species after taxonomic
#' harmonization\cr
#' \emph{taxon_ID}- Identification number of the taxonomic group\cr
#' \emph{work_genus}- Identification number of the genus after taxonomic
#' harmonization\cr
#' \emph{work_species_epithet}- Identification number of the species epithet
#' after taxonomic harmonization\cr
#' \emph{work_species} - Species name (after taxonomic harmonization)\cr
#' \emph{work_author}-  Author who described the species (after taxonomic
#' harmonization)
#'
#' @references
#'      Denelle, P., Weigelt, P., & Kreft, H. (2023). GIFT—An R package to
#'      access the Global Inventory of Floras and Traits. Methods in Ecology
#'      and Evolution, 14, 2738-2748.
#'      https://doi.org/10.1111/2041-210X.14213
#' 
#'      Weigelt, P, König, C, Kreft, H. GIFT – A Global Inventory of Floras and
#'      Traits for macroecology and biogeography. J Biogeogr. 2020; 47: 16– 43.
#'      https://doi.org/10.1111/jbi.13623
#'
#' @seealso [GIFT::GIFT_checklists()]
#'
#' @examples
#' \donttest{
#' ex <- GIFT_species_lookup(genus = "Fagus", epithet = "sylvatica")
#' }
#' 
#' @importFrom jsonlite read_json
#' @importFrom dplyr mutate_at
#' 
#' @export

GIFT_species_lookup <-
  function(genus = "", epithet = "", namesmatched = FALSE, 
           api = "https://gift.uni-goettingen.de/api/extended/",
           GIFT_version = "latest"){
    # 1. Controls ----
    api_check <- check_api(api)
    if(is.null(api_check)){
      return(NULL)
    } else{
      GIFT_version <- check_gift_version(GIFT_version)
      
      if(!is.character(genus)){
        stop("genus must be a character string indicating genus to look for.")
      }
      
      if(!is.character(epithet)){
        stop(
          "epithet must be a character string indicating the specific epithet
      to look for.")
      }
      
      if(length(namesmatched) != 1 || !is.logical(namesmatched) ||
         is.na(namesmatched)){
        stop("'namesmatched' must be a logical stating whether you only want to 
    look for the species not only in the standardized species names or also 
    in the original species names as they came in the original resources")
      }
      
      # 2. Function ----
      # Return the name matching information
      if(namesmatched){
        tmp <- jsonlite::read_json(paste0(
          api, "index", ifelse(GIFT_version == "beta", "", GIFT_version),
          ".php?query=names_matched&genus=", genus, "&epithet=", epithet
        ), simplifyVector = TRUE)
        
        if(length(tmp)>0){
          tmp <- dplyr::mutate_at(
            tmp, c("orig_ID", "name_ID", "cf_genus", "cf_species",
                   "aff_species", "matched", "epithetscore", "overallscore",
                   "resolved", "synonym", "matched_subtaxon", "accepted",
                   "work_ID", "taxon_ID"),
            as.numeric)
        } else {
          tmp <- data.frame(orig_ID = numeric(), orig_genus = character(), 
                            name_ID = numeric(), cf_genus = numeric(), 
                            genus = character(), cf_species = numeric(), 
                            aff_species = numeric(), 
                            species_epithet = character(),
                            subtaxon = character(), author = character(),
                            matched = numeric(), epithetscore = numeric(),
                            overallscore = numeric(), resolved = numeric(),
                            synonym = numeric(), matched_subtaxon = numeric(), 
                            accepted = numeric(),
                            service = character(), work_ID = numeric(),
                            taxon_ID = numeric(), work_genus = character(),
                            work_species_epithet = character(), 
                            work_species = character(),
                            work_author = character())
        }
      } else {
        tmp <- jsonlite::read_json(paste0(
          api, "index", ifelse(GIFT_version == "beta", "", GIFT_version),
          ".php?query=names_matched_unique&genus=", genus, "&epithet=", epithet
        ), simplifyVector = TRUE)
        
        if(length(tmp)>0){
          tmp <- dplyr::mutate_at(
            tmp, c("name_ID", "matched", "epithetscore", "overallscore", 
                   "resolved", "synonym", "matched_subtaxon", "accepted", 
                   "work_ID", "taxon_ID"),
            as.numeric)
        } else {
          tmp <- data.frame(name_ID = numeric(), genus = character(),
                            species_epithet = character(),
                            subtaxon = character(), author = character(),
                            matched = numeric(), epithetscore = numeric(),
                            overallscore = numeric(), resolved = numeric(),
                            synonym = numeric(), matched_subtaxon = numeric(), 
                            accepted = numeric(),
                            service = character(), work_ID = numeric(),
                            taxon_ID = numeric(), work_genus = character(),
                            work_species_epithet = character(), 
                            work_species = character(),
                            work_author = character())
        }
      }
      if(nrow(tmp)== 0){
        message("No species names found.")
      }
      return(tmp)
    }
  }
