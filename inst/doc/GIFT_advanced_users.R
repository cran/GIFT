## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.width = 8, fig.height = 8, fig.align = "center")
# Packages --------------------------------------------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    library("GIFT")
    library("knitr")
    library("kableExtra")
    library("ggplot2")
    library("sf")
  })
})

options(tinytex.verbose = TRUE)

## ----gottingen_logo, fig.show = 'hold', out.width = "20%", echo = FALSE-------
knitr::include_graphics("../man/figures/biodiv_gottingen_logo.png")
knitr::include_graphics("../man/figures/GIFT.png")

## -----------------------------------------------------------------------------
versions <- GIFT_versions()
kable(versions, "html") %>%
  kable_styling(full_width = FALSE)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  list_latest <- GIFT_lists(GIFT_version = "latest") # default value
#  list_1 <- GIFT_lists(GIFT_version = "1.0")

## -----------------------------------------------------------------------------
ref <- GIFT_references()
ref <- ref[which(ref$ref_ID %in% c(22, 10333, 10649)),
           c("ref_ID", "ref_long", "geo_entity_ref")]

# 3 first rows of that table
kable(ref, "html") %>%
  kable_styling(full_width = FALSE)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  listID_1 <- GIFT_checklists_raw(list_ID = c(11926))
#  listID_1_tax <- GIFT_checklists_raw(list_ID = c(11926), namesmatched = TRUE)
#  
#  ncol(listID_1) # 16 columns
#  ncol(listID_1_tax) # 33 columns
#  length(unique(listID_1$work_ID)); length(unique(listID_1_tax$orig_ID))

## -----------------------------------------------------------------------------
data("western_mediterranean")

## ---- fig.cap = "Figure 1. GIFT spatial", out.width = "50%", echo = FALSE-----
knitr::include_graphics("../man/figures/GIFT_spatial.svg")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  med_centroid_inside  <- GIFT_spatial(shp = western_mediterranean,
#                                       overlap = "centroid_inside")
#  med_extent_intersect <- GIFT_spatial(shp = western_mediterranean,
#                                       overlap = "extent_intersect")
#  med_shape_intersect <- GIFT_spatial(shp = western_mediterranean,
#                                      overlap = "shape_intersect")
#  med_shape_inside <- GIFT_spatial(shp = western_mediterranean,
#                                   overlap = "shape_inside")

## ---- echo = FALSE, eval = FALSE----------------------------------------------
#  med_shape_inside <- GIFT_spatial(shp = western_mediterranean,
#                                   overlap = "shape_inside")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  length(unique(med_extent_intersect$entity_ID))
#  length(unique(med_shape_intersect$entity_ID))
#  length(unique(med_centroid_inside$entity_ID))
#  length(unique(med_shape_inside$entity_ID))

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  geodata_extent_intersect <- GIFT_shapes(med_extent_intersect$entity_ID)
#  
#  geodata_shape_inside <-
#    geodata_extent_intersect[which(geodata_extent_intersect$entity_ID %in%
#                                     med_shape_inside$entity_ID), ]
#  geodata_centroid_inside <-
#    geodata_extent_intersect[which(geodata_extent_intersect$entity_ID %in%
#                                     med_centroid_inside$entity_ID), ]
#  geodata_shape_intersect <-
#    geodata_extent_intersect[which(geodata_extent_intersect$entity_ID %in%
#                                     med_shape_intersect$entity_ID), ]

## ---- echo = TRUE, eval=FALSE, fig.width = 8, fig.height = 4------------------
#  par_overlap <- par(mfrow = c(2, 2), mai = c(0, 0, 0.5, 0))
#  plot(sf::st_geometry(geodata_shape_inside),
#       col = geodata_shape_inside$entity_ID,
#       main = paste("shape inside\n",
#                    length(unique(med_shape_inside$entity_ID)),
#                    "polygons"))
#  plot(sf::st_geometry(western_mediterranean), lwd = 2, add = TRUE)
#  
#  plot(sf::st_geometry(geodata_centroid_inside),
#       col = geodata_centroid_inside$entity_ID,
#       main = paste("centroid inside\n",
#                    length(unique(med_centroid_inside$entity_ID)),
#                    "polygons"))
#  points(geodata_centroid_inside$point_x, geodata_centroid_inside$point_y)
#  plot(sf::st_geometry(western_mediterranean), lwd = 2, add = TRUE)
#  
#  plot(sf::st_geometry(geodata_shape_intersect),
#       col = geodata_shape_intersect$entity_ID,
#       main = paste("shape intersect\n",
#                    length(unique(med_shape_intersect$entity_ID)),
#                    "polygons"))
#  plot(sf::st_geometry(western_mediterranean), lwd = 2, add = TRUE)
#  
#  plot(sf::st_geometry(geodata_extent_intersect),
#       col = geodata_extent_intersect$entity_ID,
#       main = paste("extent intersect\n",
#                    length(unique(med_extent_intersect$entity_ID)),
#                    "polygons"))
#  plot(sf::st_geometry(western_mediterranean), lwd = 2, add = TRUE)
#  par(par_overlap)

## ---- fig.cap = "", out.width = "100%",echo = FALSE---------------------------
knitr::include_graphics("../man/figures/advanced_overlap.png")

## ---- echo = FALSE, eval = TRUE-----------------------------------------------
med_shape_inside <- data.frame(
  entity_ID = c(145, 146, 147, 148, 149, 150, 151, 414, 415, 416, 417, 547, 548,
                549, 550, 551, 552, 586, 591, 592, 736, 738, 739, 1033, 1036,
                10001, 10034, 10071, 10072, 10104, 10184, 10303, 10422, 10430,
                10751, 10860, 10978, 11028, 11029, 11030, 11031, 11033, 11035,
                11036, 11037, 11038, 11039, 11040, 11041, 11042, 11043, 11044,
                11045, 11046, 11434, 11455, 11461, 11474, 11477, 11503, 12065,
                12071, 12078, 12230, 12231, 12232, 12233, 12551, 12632, 12633,
                12634, 12635))

## ---- message=FALSE, fig.width = 10, fig.height = 6---------------------------
length(med_shape_inside$entity_ID)
length(GIFT_no_overlap(med_shape_inside$entity_ID, area_threshold_island = 0,
                       area_threshold_mainland = 100, overlap_threshold = 0.1))

# The following polygons are overlapping:
GIFT_no_overlap(med_shape_inside$entity_ID, area_threshold_island = 0,
                area_threshold_mainland = 100, overlap_threshold = 0.1)

## ---- eval=FALSE, echo = TRUE-------------------------------------------------
#  # Example of two overlapping polygons: Spain mainland and Andalusia
#  overlap_shape <- GIFT_shapes(entity_ID = c(10071, 12078))

## ---- include=FALSE, eval = TRUE----------------------------------------------
overlap_shape <- GIFT_shapes(entity_ID = c(10071, 12078))

## ---- message=FALSE, fig.width = 10, fig.height = 6---------------------------
par_overlap_shp <- par(mfrow = c(1, 1))
plot(sf::st_geometry(overlap_shape),
     col = c(rgb(red = 1, green = 0, blue = 0, alpha = 0.5),
             rgb(red = 0, green = 0, blue = 1, alpha = 0.3)),
     lwd = c(2, 1),
     main = "Overlapping polygons")
par(par_overlap_shp)

GIFT_no_overlap(c(10071, 12078), area_threshold_island = 0,
                area_threshold_mainland = 100, overlap_threshold = 0.1)
GIFT_no_overlap(c(10071, 12078), area_threshold_island = 0,
                area_threshold_mainland = 100000, overlap_threshold = 0.1)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  ex <- GIFT_checklists(taxon_name = "Tracheophyta", by_ref_ID = FALSE,
#                        list_set_only = TRUE)
#  ex2 <- GIFT_checklists(taxon_name = "Tracheophyta",
#                         remove_overlap = TRUE, by_ref_ID = TRUE,
#                         list_set_only = TRUE)
#  ex3 <- GIFT_checklists(taxon_name = "Tracheophyta",
#                         remove_overlap = TRUE, by_ref_ID = FALSE,
#                         list_set_only = TRUE)
#  
#  length(unique(ex$lists$ref_ID))
#  length(unique(ex2$lists$ref_ID))
#  length(unique(ex3$lists$ref_ID))

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  unique(ex2$lists$ref_ID)[!(unique(ex2$lists$ref_ID) %in%
#                               unique(ex3$lists$ref_ID))] # 25 references

## ---- include = FALSE, eval = TRUE--------------------------------------------
pilbara <- GIFT_shapes(entity_ID = c(10043, 12172, 11398, 11391, 10918))

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  # Pilbara region Australy and overlapping shapes
#  pilbara <- GIFT_shapes(entity_ID = c(10043, 12172, 11398, 11391, 10918))

## -----------------------------------------------------------------------------
ggplot(pilbara) +
  geom_sf(aes(fill = as.factor(entity_ID)), alpha = 0.5) +
  scale_fill_brewer("entity_ID", palette = "Set1")

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  species <- GIFT_species()

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  # Add Family
#  species$Family <- GIFT_taxgroup(
#    as.numeric(species$work_ID), taxon_lvl = "family", return_ID = FALSE,
#    species = species)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  GIFT_taxgroup(as.numeric(species$work_ID[1:5]), taxon_lvl = "order",
#                return_ID = FALSE)
#  GIFT_taxgroup(as.numeric(species$work_ID[1:5]),
#                taxon_lvl = "higher_lvl", return_ID = FALSE,
#                species = species)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  Fagus <- GIFT_species_lookup(genus = "Fagus", epithet = "sylvatica",
#                               namesmatched = TRUE)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  taxo <- GIFT_taxonomy()

## -----------------------------------------------------------------------------
glonaf <- GIFT_overlap(resource = "glonaf")

kable(glonaf[1:5, ], "html") %>%
  kable_styling(full_width = FALSE)

gmba <- GIFT_overlap(resource = "gmba")

kable(gmba[1:5, ], "html") %>%
  kable_styling(full_width = FALSE)

