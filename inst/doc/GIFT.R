## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.width = 4, fig.height = 4, fig.align = "center")
# Packages --------------------------------------------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    library("GIFT")
    library("dplyr")
    library("knitr")
    library("kableExtra")
    library("ggplot2")
    library("sf")
    library("rnaturalearth")
    library("rnaturalearthdata")
    library("tidyr")
    library("ape")
    library("patchwork")
  })
})

options(tinytex.verbose = TRUE)

## ----gottingen_logo, fig.show = 'hold', out.width = "20%", echo = FALSE-------
knitr::include_graphics("../man/figures/biodiv_gottingen_logo.png")
knitr::include_graphics("../man/figures/GIFT.png")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  library("GIFT")
#  library("dplyr")
#  library("knitr")
#  library("kableExtra")
#  library("ggplot2")
#  library("sf")
#  library("rnaturalearth")
#  library("rnaturalearthdata")
#  library("tidyr")
#  library("ape")
#  library("patchwork")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  options(timeout = max(1000, getOption("timeout")))

## ----fig.width = 3, fig.height = 3--------------------------------------------
data("western_mediterranean")

world <- ne_coastline(scale = "medium", returnclass = "sf")
world_countries <- ne_countries(scale = "medium", returnclass = "sf")
# Fixing polygons crossing dateline
world <- st_wrap_dateline(world)
world_countries <- st_wrap_dateline(world_countries)

# Eckert IV projection
eckertIV <-
  "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

ggplot(world) +
  geom_sf(color = "gray50") +
  geom_sf(data = western_mediterranean, fill = "darkblue", color = "black",
          alpha = 0.5, size = 1) +
  labs(title = "Western Mediterranean basin") +
  lims(x = c(-20, 20), y = c(24, 48)) +
  theme_void()

## ----fig.cap = "Figure 1. Principle of the complete_taxon argument", out.width = "50%",echo = FALSE----
knitr::include_graphics("../man/figures/GIFT_complete_taxon.svg")

## ----fig.cap = "Figure 2. Principle of the complete_floristic argument", out.width = "50%", echo = FALSE----
knitr::include_graphics("../man/figures/GIFT_complete_floristic.svg")

## ----fig.cap = "Figure 3. Principle of GIFT_spatial()", out.width = "50%", echo = FALSE----
knitr::include_graphics("../man/figures/GIFT_spatial.svg")

## ----fig.cap = "Figure 4. Removing overlapping polygons with remove_overlap argument", out.width = "50%", echo = FALSE----
knitr::include_graphics("../man/figures/GIFT_no_overlap.svg")

## ----fig.cap = "Figure 5. Principle of the overlap_th argument", out.width = "50%", echo = FALSE----
knitr::include_graphics("../man/figures/GIFT_overlap_th.svg")

## ----eval=FALSE---------------------------------------------------------------
#  ex_meta <- GIFT_checklists(taxon_name = "Angiospermae",
#                             shp = western_mediterranean,
#                             overlap = "centroid_inside",
#                             list_set_only = TRUE)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  medit <- GIFT_checklists(taxon_name = "Angiospermae",
#                           complete_taxon = TRUE,
#                           floristic_group = "native",
#                           complete_floristic = TRUE,
#                           geo_type = "All",
#                           shp = western_mediterranean,
#                           overlap = "centroid_inside",
#                           remove_overlap = FALSE,
#                           taxonomic_group = TRUE) # this argument adds two
#  # columns to the checklist: plant family and taxonomic group of each species

## ----eval=FALSE---------------------------------------------------------------
#  # Number of references covered
#  length(unique(medit[[2]]$ref_ID))
#  #   22 references
#  
#  # Number of checklists covered (one reference can have several lists inside)
#  length(unique(medit[[2]]$list_ID))
#  #   115 checklists
#  
#  # Number of species
#  length(unique(medit[[2]]$work_species))
#  #   12840 plant species

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  medit_no_overlap <- GIFT_checklists(shp = western_mediterranean,
#                                      overlap = "centroid_inside",
#                                      taxon_name = "Angiospermae",
#                                      remove_overlap = TRUE)
#  
#  # Number of references covered
#  length(unique(medit[[2]]$ref_ID)) # 23 references
#  length(unique(medit_no_overlap[[2]]$ref_ID)) # 22 references

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  custom_point <- cbind(9.9, 51) # coordinates of Göttingen
#  
#  got <- GIFT_checklists(coordinates = custom_point,
#                         overlap = "extent_intersect",
#                         taxon_name = "Angiospermae",
#                         remove_overlap = TRUE,
#                         list_set_only = TRUE)

## ----eval=FALSE, echo = TRUE--------------------------------------------------
#  gift_shapes <- GIFT_shapes() # retrieves all shapefiles by default

## ----echo = TRUE, eval = FALSE, message=FALSE, fig.width = 8, fig.height = 4----
#  angio_rich <- GIFT_richness(taxon_name = "Angiospermae")
#  
#  rich_map <- dplyr::left_join(gift_shapes, angio_rich, by = "entity_ID") %>%
#    dplyr::filter(stats::complete.cases(total))
#  
#  ggplot(world) +
#    geom_sf(color = "gray50") +
#    geom_sf(data = rich_map, aes(fill = total + 1)) +
#    scale_fill_viridis_c("Species number\n(log-transformed)", trans = "log10",
#                         labels = scales::number_format(accuracy = 1)) +
#    labs(title = "Angiosperms", subtitle = "Projection EckertIV") +
#    coord_sf(crs = eckertIV) +
#    theme_void()

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/angio_rich_map.png")

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/angio_rich_map2.png")

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  # Background box
#  xmin <- st_bbox(world)[["xmin"]]; xmax <- st_bbox(world)[["xmax"]]
#  ymin <- st_bbox(world)[["ymin"]]; ymax <- st_bbox(world)[["ymax"]]
#  bb <- sf::st_union(sf::st_make_grid(st_bbox(c(xmin = xmin,
#                                                xmax = xmax,
#                                                ymax = ymax,
#                                                ymin = ymin),
#                                              crs = st_crs(4326)),
#                                      n = 100))
#  
#  # Equator line
#  equator <- st_linestring(matrix(c(-180, 0, 180, 0), ncol = 2, byrow = TRUE))
#  equator <- st_sfc(equator, crs = st_crs(world))
#  
#  # Color code from Barthlott 2007
#  hexcode_barthlott2007 <- c("#fbf9ed", "#f3efcc", "#f6e39e", "#cbe784",
#                             "#65c66a", "#0e8d4a", "#4a6fbf",
#                             "#b877c2", "#f24dae", "#ed1c24")
#  
#  ggplot(world) +
#    geom_sf(data = bb, fill = "aliceblue") +
#    geom_sf(data = equator, color = "gray50", linetype = "dashed",
#            linewidth = 0.1) +
#    geom_sf(data = world_countries, fill = "antiquewhite1", color = NA) +
#    geom_sf(color = "gray50", linewidth = 0.1) +
#    geom_sf(data = bb, fill = NA) +
#    geom_sf(data = rich_map,
#            aes(fill = ifelse(rich_map$entity_class %in%
#                                c("Island/Mainland", "Mainland",
#                                  "Island Group", "Island Part"),
#                              total + 1, NA)),
#            size = 0.1) +
#    geom_point(data = rich_map,
#               aes(color = ifelse(rich_map$entity_class %in%
#                                    c("Island"),
#                                  total + 1, NA),
#                   geometry = geometry),
#               stat = "sf_coordinates", size = 1, stroke = 0.5) +
#    scale_color_gradientn(
#      "Species number", trans = "log10", limits = c(1, 40000),
#      colours = hexcode_barthlott2007,
#      breaks = c(1, 10, 100, 1000, 10000, 40000),
#      labels = c(1, 10, 100, 1000, 10000, 40000),
#      na.value = "transparent") +
#    scale_fill_gradientn(
#      "Species number", trans = "log10", limits = c(1, 40000),
#      colours = hexcode_barthlott2007,
#      breaks = c(1, 10, 100, 1000, 10000, 40000),
#      labels = c(1, 10, 100, 1000, 10000, 40000),
#      na.value = "transparent") +
#    labs(title = "Angiosperms", subtitle = "Projection EckertIV") +
#    coord_sf(crs = eckertIV) +
#    theme_void()

## ----echo = TRUE, eval = FALSE, message=FALSE, fig.width = 8, fig.height = 4----
#  med_shape <- gift_shapes[which(gift_shapes$entity_ID %in%
#                                   unique(medit[[2]]$entity_ID)), ]
#  
#  med_rich <- angio_rich[which(angio_rich$entity_ID %in%
#                                 unique(medit[[2]]$entity_ID)), ]
#  
#  med_map <- dplyr::left_join(med_shape, med_rich, by = "entity_ID") %>%
#    dplyr::filter(stats::complete.cases(total))
#  
#  ggplot(world) +
#    geom_sf(color = "gray50") +
#    geom_sf(data = western_mediterranean,
#            fill = "darkblue", color = "black", alpha = 0.1, size = 1) +
#    geom_sf(data = med_map, aes(fill = total)) +
#    scale_fill_viridis_c("Species number") +
#    labs(title = "Angiosperms in the Western Mediterranean basin") +
#    lims(x = c(-20, 20), y = c(24, 48)) +
#    theme_void()

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/angio_medit.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  all_sp <- GIFT_species()

## ----anemone_lookup-----------------------------------------------------------
anemone_lookup <- GIFT_species_lookup(genus = "Anemone", epithet = "nemorosa")

kable(anemone_lookup, "html") %>%
  kable_styling(full_width = FALSE)

## -----------------------------------------------------------------------------
anemone_distr <- GIFT_species_distribution(
  genus = "Anemone", epithet = "nemorosa", aggregation = TRUE,
  GIFT_version = "2.2")

anemone_statuses <- anemone_distr %>%
  mutate(native = ifelse(native == 1, "native", "non-native"),
         naturalized = ifelse(naturalized == 1, "naturalized",
                              "non-naturalized"),
         endemic_list = ifelse(endemic_list == 1, "endemic_list",
                               "non-endemic_list")) %>%
  dplyr::select(entity_ID, native, naturalized, endemic_list)

table(anemone_statuses$endemic_list)

## -----------------------------------------------------------------------------
table(paste(anemone_statuses$native, anemone_statuses$naturalized,
            sep = "_"))

## ----echo = TRUE, eval = FALSE, fig.width = 8, fig.height = 4-----------------
#  # We rename the statuses based on the distinct combinations
#  anemone_statuses <- anemone_statuses %>%
#    mutate(Status = case_when(
#      native == "native" & naturalized == "non-naturalized" ~ "native",
#      native == "native" & is.na(naturalized) ~ "native",
#      native == "non-native" & is.na(naturalized) ~ "non-native",
#      native == "non-native" & naturalized == "naturalized" ~ "naturalized",
#      native == "non-native" & naturalized == "non-naturalized" ~ "non-native",
#      is.na(native) & is.na(naturalized) ~ "unknown"
#    ))
#  
#  # Merge with the shapes
#  anemone_shape <- gift_shapes[which(gift_shapes$entity_ID %in%
#                                       unique(anemone_distr$entity_ID)), ]
#  anemone_map <- dplyr::left_join(anemone_shape, anemone_statuses,
#                                  by = "entity_ID")
#  
#  # Area of distribution with floristic status
#  ggplot(world) +
#    geom_sf(color = "gray70") +
#    geom_sf(data = anemone_map, color = "black", aes(fill = as.factor(Status))) +
#    scale_fill_brewer("Status", palette = "Set2") +
#    labs(title = expression(paste("Distribution map of ",
#                                  italic("Anemone nemorosa"))),
#         subtitle = "Unprojected (GCS: WGS84)") +
#    lims(x = c(-65, 170), y = c(-45, 70)) +
#    theme_void()

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/anemone_map.png")

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/anemone_map2.png")

## ----eval=FALSE---------------------------------------------------------------
#  anemone_map_plot_bg_parts <-
#    ggplot(world) +
#    geom_sf(data = bb, fill = "aliceblue", color = NA) +
#    geom_sf(data = equator, color = "gray50", linetype = "dashed",
#            linewidth = 0.1) +
#    geom_sf(data = world_countries, fill = "antiquewhite1", color = NA) +
#    geom_sf(color = "gray50", linewidth = 0.1) +
#    geom_sf(data = bb, fill = NA) +
#    geom_sf(data = anemone_map, color = "black", aes(fill = as.factor(Status))) +
#    scale_fill_manual("Status",
#                      values = c("native" = "#2c7bb6",
#                                 "naturalized" = "#d7191c",
#                                 "non-native" = "#fdae61",
#                                 "unknown" = "#abd9e9")) +
#    labs(title = expression(paste("b) Distribution map of ",
#                                  italic("Anemone nemorosa")))) +
#    theme_void() +
#    theme(axis.title = element_blank(),
#          axis.text = element_blank(),
#          axis.ticks = element_blank())
#  
#  (anemone_map_plot_bg_parts +
#      lims(x = c(-69, 61), y = c(37, 70)) + # Europe & Newfoundland
#      theme(panel.border = element_rect(fill = NA, linewidth = 1)) +
#      theme(legend.position = "bottom")
#    |
#      anemone_map_plot_bg_parts +
#      lims(x = c(165, 178), y = c(-47, -35)) + # new zealand
#      labs(title = "") +
#      guides(fill = "none") +
#      theme(panel.border = element_rect(fill = NA, linewidth = 1)))

## -----------------------------------------------------------------------------
trait_meta <- GIFT_traits_meta()
trait_meta[which(trait_meta$Trait2 == "Plant_height_max"), ]

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  height <- GIFT_traits(trait_IDs = c("1.6.2"), agreement = 0.66,
#                        bias_ref = FALSE, bias_deriv = FALSE)
#  
#  height_raw <- GIFT_traits_raw(trait_IDs = c("1.6.2"))
#  
#  # Raw values
#  as.numeric(height_raw[which(height_raw$work_species == "Fagus sylvatica"),
#                        "trait_value"])
#  
#  # Aggregated value
#  as.numeric(height[which(height$work_species == "Fagus sylvatica"),
#                    "trait_value_1.6.2"])

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  references <- GIFT_references(GIFT_version = "beta")
#  
#  unique(unlist(strsplit(height$references_1.6.2, ",")))
#  
#  references <- references[
#    which(references$ref_ID %in%
#            unique(unlist(strsplit(height$references_1.6.2, ",")))), ]
#  references[1:2, ]

## ----echo=TRUE, eval=FALSE----------------------------------------------------
#  trait_tax <- GIFT_traits_tax(trait_IDs = c("1.1.1", "1.2.1", "1.4.1"),
#                               bias_ref = FALSE, bias_deriv = FALSE)
#  
#  trait_tax[1:3, ]

## ----echo = TRUE, eval = FALSE, fig.width = 8, fig.height = 4-----------------
#  angio_height <- GIFT_coverage(what = "trait_coverage",
#                                taxon_name = "Angiospermae",
#                                trait_ID = "1.6.2")
#  
#  angio_height_shape <- gift_shapes[which(gift_shapes$entity_ID %in%
#                                            unique(angio_height$entity_ID)), ]
#  
#  angio_height_map <- dplyr::left_join(
#    angio_height_shape, angio_height, by = "entity_ID")
#  
#  angio_height_map <- angio_height_map[complete.cases(angio_height_map$native), ]
#  
#  ggplot(world) +
#    geom_sf(color = "gray50") +
#    geom_sf(data = angio_height_map[complete.cases(angio_height_map$native), ],
#            aes(fill = native)) +
#    scale_fill_viridis_c("Coverage (%)") +
#    labs(title = "Coverage for maximal vegetative height of Angiosperms",
#         subtitle = "Projection EckertIV") +
#    coord_sf(crs = eckertIV) +
#    theme_void()

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/angio_height_plot.png")

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/angio_height_plot2.png")

## ----eval=FALSE---------------------------------------------------------------
#  ggplot(world) +
#    geom_sf(data = bb, fill = "aliceblue") +
#    geom_sf(data = equator, color = "gray50", linetype = "dashed",
#            linewidth = 0.1) +
#    geom_sf(data = world_countries, fill = "antiquewhite1", color = NA) +
#    geom_sf(color = "gray50", linewidth = 0.1) +
#    geom_sf(data = bb, fill = NA) +
#    geom_sf(data = angio_height_map,
#            aes(fill = ifelse(angio_height_map$entity_class %in%
#                                c("Island/Mainland", "Mainland",
#                                  "Island Group", "Island Part"),
#                              100*native, NA)), size = 0.1) +
#    geom_point(data = angio_height_map,
#               aes(color = ifelse(angio_height_map$entity_class %in%
#                                    c("Island"),
#                                  100*native, NA),
#                   geometry = geometry),
#               stat = "sf_coordinates", size = 1, stroke = 0.5) +
#    scale_color_gradientn(
#      "Coverage (%)",
#      colours = rev(RColorBrewer::brewer.pal(9, name = "PuBuGn")),
#      limits = c(0, 100),
#      na.value = "transparent") +
#    scale_fill_gradientn(
#      "Coverage (%)",
#      colours = rev(RColorBrewer::brewer.pal(9, name = "PuBuGn")),
#      limits = c(0, 100),
#      na.value = "transparent") +
#    labs(title = "Coverage for maximal vegetative height of Angiosperms",
#         subtitle = "Projection EckertIV") +
#    coord_sf(crs = eckertIV) +
#    theme_void()

## ----env_meta-----------------------------------------------------------------
misc_env <- GIFT_env_meta_misc()
raster_env <- GIFT_env_meta_raster()

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  med_env <- GIFT_env(entity_ID = unique(medit[[2]]$entity_ID),
#                      miscellaneous = c("perimeter", "biome"),
#                      rasterlayer = c("mn30_grd", "wc2.0_bio_30s_01"),
#                      sumstat = list(c("mean", "med"), "max"))
#  
#  med_env[1, ]

## ----echo = TRUE, eval = FALSE, fig.width = 8, fig.height = 4-----------------
#  world_temp <- GIFT_env(entity_ID = unique(angio_rich$entity_ID),
#                         rasterlayer = c("wc2.0_bio_30s_01"),
#                         sumstat = c("mean"))
#  
#  temp_shape <- gift_shapes[which(gift_shapes$entity_ID %in%
#                                    unique(angio_rich$entity_ID)), ]
#  
#  temp_map <- dplyr::left_join(temp_shape, world_temp, by = "entity_ID")
#  
#  ggplot(world) +
#    geom_sf(color = "gray50") +
#    geom_sf(data = temp_map, aes(fill = mean_wc2.0_bio_30s_01)) +
#    scale_fill_viridis_c("Celsius degrees") +
#    labs(title = "Average temperature",
#         subtitle = "Projection EckertIV") +
#    coord_sf(crs = eckertIV) +
#    theme_void()

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/temp_plot.png")

## ----fig.cap = "", out.width = "90%",echo = FALSE-----------------------------
knitr::include_graphics("../man/figures/temp_plot2.png")

## ----eval=FALSE---------------------------------------------------------------
#  ggplot(world) +
#    geom_sf(data = bb, fill = "aliceblue") +
#    geom_sf(data = equator, color = "gray50", linetype = "dashed",
#            linewidth = 0.1) +
#    geom_sf(data = world_countries, fill = "antiquewhite1", color = NA) +
#    geom_sf(color = "gray50", linewidth = 0.1) +
#    geom_sf(data = bb, fill = NA) +
#    geom_sf(data = temp_map,
#            aes(fill = ifelse(temp_map$entity_class %in%
#                                c("Island/Mainland", "Mainland",
#                                  "Island Group", "Island Part"),
#                              mean_wc2.0_bio_30s_01, NA)), size = 0.1) +
#    geom_point(data = temp_map,
#               aes(color = ifelse(temp_map$entity_class %in%
#                                    c("Island"),
#                                  mean_wc2.0_bio_30s_01, NA),
#                   geometry = geometry),
#               stat = "sf_coordinates", size = 1, stroke = 0.5) +
#    scale_color_gradientn(
#      "°C",
#      colours = RColorBrewer::brewer.pal(9, name = "Reds"),
#      limits = c(-20, 30),
#      na.value = "transparent") +
#    scale_fill_gradientn(
#      "°C",
#      colours = RColorBrewer::brewer.pal(9, name = "Reds"),
#      limits = c(-20, 30),
#      na.value = "transparent") +
#    labs(title = "Average temperature",
#         subtitle = "Projection EckertIV") +
#    coord_sf(crs = eckertIV) +
#    theme_void()

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  # Retrieving phylogeny, taxonomy and species from GIFT
#  phy <- GIFT_phylogeny(clade = "Tracheophyta", GIFT_version = "beta")
#  tax <- GIFT_taxonomy(GIFT_version = "beta")
#  gift_sp <- GIFT_species(GIFT_version = "beta")
#  
#  gf <- GIFT_traits(trait_IDs = "1.2.1", agreement = 0.66, bias_ref = FALSE,
#                    bias_deriv = FALSE, GIFT_version = "beta")

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  # Replacing space with _ for the species names
#  gf$work_species <- gsub(" ", "_", gf$work_species, fixed = TRUE)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  # Retrieving family of each species
#  sp_fam <- GIFT_taxgroup(work_ID = unique(gift_sp$work_ID),
#                          taxon_lvl = "family", GIFT_version = "beta")
#  sp_genus_fam <- data.frame(
#    work_ID = unique(gift_sp$work_ID),
#    work_species = unique(gift_sp$work_species),
#    family = sp_fam)
#  sp_genus_fam <- left_join(sp_genus_fam,
#                            gift_sp[, c("work_ID", "work_genus")],
#                            by = "work_ID")
#  colnames(sp_genus_fam)[colnames(sp_genus_fam) == "work_genus"] <- "genus"
#  
#  # Problem with hybrid species on the tip labels of the phylo tree
#  phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"] <-
#    substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"],
#              3,
#              nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))
#  
#  phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"] <-
#    substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"],
#              3,
#              nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  sp_genus_fam <- left_join(sp_genus_fam,
#                            gf[, c("work_ID", "trait_value_1.2.1")],
#                            by = "work_ID")
#  
#  genus_gf <- sp_genus_fam %>%
#    group_by(genus) %>%
#    mutate(prop_gf = round(100*sum(is.na(trait_value_1.2.1))/n(), 2)) %>%
#    ungroup() %>%
#    dplyr::select(-work_ID, -work_species, -family, -trait_value_1.2.1) %>%
#    distinct(.keep_all = TRUE)
#  
#  fam_gf <- sp_genus_fam %>%
#    group_by(family) %>%
#    mutate(prop_gf = round(100*sum(is.na(trait_value_1.2.1))/n(), 2)) %>%
#    ungroup() %>%
#    dplyr::select(-work_ID, -work_species, -genus, -trait_value_1.2.1) %>%
#    distinct(.keep_all = TRUE)
#  
#  sp_genus_fam$species <- gsub("([[:punct:]])|\\s+", "_",
#                               sp_genus_fam$work_species)
#  
#  # Keeping one species per genus only
#  one_sp_per_gen <- data.frame()
#  for(i in 1:n_distinct(sp_genus_fam$genus)){ # loop over genera
#    # Focal genus
#    focal_gen <- unique(sp_genus_fam$genus)[i]
#    # All species in that genus
#    gen_sp_i <- sp_genus_fam[which(sp_genus_fam$genus == focal_gen),
#                             "species"]
#    # Species from the genus available in the phylogeny
#    gen_sp_i <- gen_sp_i[gen_sp_i %in% phy$tip.label]
#    # Taking the first one (if at least one is available)
#    gen_sp_i <- gen_sp_i[1]
#  
#    one_sp_per_gen <- rbind(one_sp_per_gen,
#                            data.frame(species = gen_sp_i,
#                                       genus = focal_gen))
#  }
#  
#  # Adding the trait coverage per genus
#  one_sp_per_gen <- left_join(one_sp_per_gen, genus_gf, by = "genus")
#  
#  # Adding the trait coverage per family
#  one_sp_per_gen <- left_join(one_sp_per_gen,
#                              sp_genus_fam[!duplicated(sp_genus_fam$genus),
#                                           c("genus", "family")],
#                              by = "genus")
#  colnames(one_sp_per_gen)[colnames(one_sp_per_gen) == "prop_gf"] <-
#    "prop_gf_gen"
#  one_sp_per_gen <- left_join(one_sp_per_gen, fam_gf, by = "family")
#  colnames(one_sp_per_gen)[colnames(one_sp_per_gen) == "prop_gf"] <-
#    "prop_gf_fam"

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  phy_gen <- ape::keep.tip(
#    phy = phy,
#    tip = one_sp_per_gen[complete.cases(one_sp_per_gen$species), "species"])

## ----echo=TRUE, eval = FALSE--------------------------------------------------
#  library("BiocManager")
#  install("ggtree")
#  library("ggtree")
#  library("tidytree")
#  install("ggtreeExtra")
#  library("ggtreeExtra")

## ----echo = TRUE, eval = FALSE, fig.width = 6, fig.height = 6-----------------
#  ggtree(phy_gen, color = "grey70", layout = "circular") %<+% one_sp_per_gen +
#    geom_fruit(geom = geom_tile,
#               mapping = aes(fill = prop_gf_gen),
#               width = 50,
#               offset = 0.1) +
#    geom_fruit(geom = geom_tile,
#               mapping = aes(color = prop_gf_fam, fill = prop_gf_fam),
#               width = 50,
#               offset = 0.1,
#               show.legend = FALSE) +
#    scale_color_viridis_c() +
#    scale_fill_viridis_c("Growth form availability per genus (%)") +
#    theme(legend.position = "bottom")

## ----fig.cap = "", out.width = "100%",echo = FALSE----------------------------
knitr::include_graphics("../man/figures/phy_tree_plot.png")

