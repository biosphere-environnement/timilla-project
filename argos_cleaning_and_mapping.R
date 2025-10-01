# ============================================================
#
# SCRIPT: ARGOS_CLEANING_AND_MAPPING.R
#
#
# AUTHOR: Raphaël Musseau — BioSphère Environnement
# YEAR: 2025
# LICENSE: MIT
# SPDX-License-Identifier: MIT
#
# DESCRIPTION:
#   Processes and filters ARGOS tracking data and produces maps
#   of individual trajectories (LC 1–3 kept; LC 0/A/B/C kept if smaj ≤ 6 km).
#
# CITATION:
#   If you use this script in scientific work, please cite:
#   Musseau, R. (2025). ARGOS data filtering and mapping script (v1.0.0).
#   BioSphère Environnement. GitHub. https://github.com/biosphere-environnement/timilla-project.git
#
# CONTACT: r.musseau@biosphere-environnement.com  ORCID: 0000-0003-3825-6418
#
# ============================================================
##
## ============================================================
##
## STRUCTURE DU SCRIPT & LOGIQUE DE TRAITEMENT
##
## 0) PACKAGES
##    - Charge les bibliothèques nécessaires pour lire, filtrer, cartographier.
##
## 1) Source des données (URL ou local)
##    - Lit un CSV (séparateur auto-détecté “;” ou “,”).
##    - Entrée attendue : colonnes ARGOS standard (N° ID, Loc idx, Date de loc., Longitude, Latitude, Qualité loc., Demi-grand axe, Date Msg [optionnel]).
##    - Sortie : data.frame brut `raw`.
##
## 2) PARAMETRES (A ADAPTER)
##    - Seuil ellipse pour LC faibles (smaj ≤ 6 km).
##    - Règles d’affichage : ruptures visuelles (Δt > 12 h, saut > 150 km).
##    - Ordre de qualité LC (3 > 2 > 1 > 0 > A > B > C).
##    - Libellés des balises (283694/95/96/97).
##
## 3) UTILITAIRES
##    - `to_num()` : convertit virgule → point et chaîne → numérique.
##    - `to_time()` : parse robustement les dates (plusieurs formats).
##
## 4) PARSE & NORMALISATION
##    - Convertit/typage des champs clés (id, loc_idx, date_loc, lon, lat, lc, smaj).
##    - Détecte l’unité du demi-grand axe (m vs km) et met en km.
##    - Ne garde que les balises d’intérêt (283694–97) et LC autorisées.
##    - Ajoute `lc_rank` (rang de qualité) et `date_loc_s` (arrondi à la seconde).
##    - Sortie : tibble propre `d`.
##
## 5) DEDOUBLONNAGE PAR (ID, DATE_LOC_S, LOC_IDX)
##    - Objectif : ne garder qu’UNE ligne “représentative” quand plusieurs lignes
##      décrivent la même localisation opérationnelle.
##    - Règle de tri (priorité décroissante) :
##        (1) meilleure LC (3>2>1>0>A>B>C),
##        (2) ellipse plus petite (smaj_km ; NA en dernier),
##        (3) dernier message (`Date Msg` si présent, sinon `Date de loc.`).
##    - Sortie : `d1`.
##
## 6) FILTRE QUALITE (ELLIPSE + LC)
##    - Garde LC 1/2/3 systématiquement.
##    - Pour LC 0/A/B/C : garde si `smaj_km` est défini et ≤ seuil (6 km par défaut).
##    - Sortie : `d_final` (données retenues pour analyse/carte).
##
## 7) COMPTES LC (POUR LA LEGENDE)
##    - Calcule n par LC et génère des étiquettes “LC (n = …)”.
##    - Détermine l’ordre d’affichage en 2 lignes : “1,2,3” puis “0,A,B”.
##
## 8) RUPTURES VISUELLES (CARTE UNIQUEMENT)
##    - Marque des segments pointillés quand Δt > 12 h ou saut > 150 km.
##    - Ne modifie PAS `d_final` (effet graphique uniquement).
##    - Sorties : `traj_cont` (segments continus), `traj_gap` (ponts pointillés).
##
## 9) FOND CARTOGRAPHIQUE
##    - Récupère les polygones pays (incluant union Maroc / Sahara occidental).
##    - Coastlines pour un trait de côte net.
##
## 10) COULEURS & LABELS
##    - Palette par balise ; remappe les ID en “283694 (NAJAT)”, etc.
##
## 11) CARTE
##    - Tracé des trajectoires (lignes), ruptures (pointillés) et points (formes par LC).
##    - Deux légendes empilées : balises (sans titre), puis “Classes de localisation”.
##    - Échelle style “ticks” (trait fin) et flèche nord.
##
## 12) EXPORTS (OPTIONNELS)
##    - `ggsave()` pour la carte PNG.
##    - `write.csv()` pour le tableau filtré (si souhaité).
##
## *** RESUME DE LA SELECTION ***
##
##  - Détection de doublons opérationnels via (id, date_loc au sec, loc_idx).
##  - Choix de la meilleure ligne par priorité LC, puis ellipse (plus petite), puis dernier message.
##  - Conservation des LC élevées (1/2/3) sans restriction d’ellipse.
##  - Pour LC faibles (0/A/B/C), conservation si incertitude compatible (smaj ≤ 6 km),
##    ce qui élimine la majorité des localisations grossières ou aberrantes.
##  - Aucune contrainte de vitesse/angle n’est appliquée (ne biaise pas la migration),
##    les “sauts” n’affectent que la représentation (pointillés), pas la sélection.
##
## *** REPRODUCTIBILITE ***
##
##  - Paramètres regroupés en tête de script.
##  - Fonctions utilitaires isolées.
##  - Ordre des opérations identique à chaque exécution.
##
## ============================================================

## 0) PACKAGES ----
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(ggspatial)
  library(sf)
  library(geosphere)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(scales)
})

## 1) PARAMÈTRES ----
csv_url       <- "https://raw.githubusercontent.com/biosphere-environnement/timilla-project/main/argos_data_turtledove_morocco_2025.csv"

ellipse_max_km <- 6      # seuil demi-grand axe (km) pour garder 0/A/B/C
gap_dt_h       <- 12     # rupture visuelle si Δt > 12 h
gap_dist_km    <- 150    # rupture visuelle si saut > 150 km

# Classes ARGOS considérées et ordre de qualité
lc_keep  <- c("1","2","3","0","A","B","C")
lc_order <- c("3","2","1","0","A","B","C")
lc_rank  <- setNames(seq_along(lc_order), lc_order)

# Noms “humains” pour les 4 balises
balise_labels <- c(
  "283694" = "283694 (NAJAT)",
  "283695" = "283695 (INES)",
  "283696" = "283696 (NAYLE)",
  "283697" = "283697 (BENOIT)"
)

## 2) FONCTIONS UTILITAIRES ----
to_num  <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
to_time <- function(x) suppressWarnings(as.POSIXct(
  as.character(x), tz = "UTC",
  tryFormats = c("%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M",
                 "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M"))
)

## 3) LECTURE DU CSV (séparateur ;) ----
raw <- read.csv(
  csv_url,
  sep = ";",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)

# Vérif colonnes minimales (exactement les libellés ARGOS FR)
required <- c("N° ID","Loc idx","Date de loc.","Longitude","Latitude",
              "Qualité loc.","Demi-grand axe")
if (!all(required %in% names(raw))) {
  stop("Colonnes requises manquantes. Trouvées :\n",
       paste(names(raw), collapse = " | "))
}

has_date_msg <- "Date Msg" %in% names(raw)

## 4) PRÉPARATION DES VARIABLES (simple & lisible) ----
# 4.1 Détection unité du demi-grand axe et conversion en km
smaj_raw <- to_num(raw[["Demi-grand axe"]])
smaj_km  <- if (isTRUE(stats::median(smaj_raw, na.rm = TRUE) > 50)) smaj_raw/1000 else smaj_raw

# 4.2 Construction du data.frame de travail (typé)
d <- tibble(
  id       = as.character(raw[["N° ID"]]),
  loc_idx  = suppressWarnings(as.integer(to_num(raw[["Loc idx"]]))),   # peut être NA
  date_loc = to_time(raw[["Date de loc."]]),
  lon      = to_num(raw[["Longitude"]]),
  lat      = to_num(raw[["Latitude"]]),
  lc       = toupper(trimws(as.character(raw[["Qualité loc."]]))),
  smaj_km  = smaj_km
)

# 4.3 Ajout date_msg si dispo (facultatif mais utile pour départager)
if (has_date_msg) {
  d$date_msg <- to_time(raw[["Date Msg"]])
} else {
  d$date_msg <- as.POSIXct(NA, tz="UTC")
}

# 4.4 Filtres de base & colonnes auxiliaires
d <- d %>%
  filter(!is.na(id), !is.na(date_loc), is.finite(lon), is.finite(lat)) %>%
  filter(id %in% names(balise_labels)) %>%     # ne garder que les 4 balises du projet
  filter(lc %in% lc_keep) %>%
  mutate(
    lc_rank    = lc_rank[lc],
    date_loc_s = floor_date(date_loc, "1 second")  # arrondi à la seconde
  ) %>%
  arrange(id, date_loc_s, loc_idx)

## 5) DÉDOUBLONNAGE (id, date_loc_s, loc_idx) ----
# Règle de tri décroissante :
#  (1) meilleure LC (3>2>1>0>A>B>C)
#  (2) ellipse plus petite (NA en dernier)
#  (3) dernier message (Date Msg si présent, sinon Date de loc.)
d1 <- d %>%
  mutate(loc_idx_chr = ifelse(is.na(loc_idx), "__NA__", as.character(loc_idx))) %>%
  group_by(id, date_loc_s, loc_idx_chr) %>%
  arrange(
    lc_rank,
    is.na(smaj_km), smaj_km,
    dplyr::desc(coalesce(date_msg, date_loc)),
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(-lc_rank, -loc_idx_chr) %>%
  arrange(id, date_loc)

## 6) FILTRE QUALITÉ (ellipse + LC) ----
d_final <- d1 %>%
  filter(
    lc %in% c("1","2","3") |
      (lc %in% c("0","A","B","C") & is.finite(smaj_km) & smaj_km <= ellipse_max_km)
  )

## 7) COMPTES LC (pour la légende) ----
lc_show <- c("1","2","3","0","A","B","C")      # ordre global
d_final$lc <- factor(d_final$lc, levels = lc_show)

tab_lc <- as.data.frame(table(factor(d_final$lc, levels = lc_show)))
tab_lc <- tab_lc[tab_lc$Freq > 0, ]
lab_lc <- setNames(paste0(tab_lc$Var1, " (n=", tab_lc$Freq, ")"), tab_lc$Var1)

# On n’affiche que 1,2,3 puis 0,A,B (C retirée de la légende)
lc_breaks <- c("1","2","3","0","A","B")
lc_breaks <- lc_breaks[lc_breaks %in% names(lab_lc)]

## 8) RUPTURES VISUELLES (pour la carte uniquement) ----
d_plot <- d_final %>%
  group_by(id) %>%
  arrange(date_loc, .by_group = TRUE) %>%
  mutate(
    dt_h = as.numeric(difftime(date_loc, lag(date_loc), units = "hours")),
    ok   = is.finite(lon) & is.finite(lat) & is.finite(lag(lon)) & is.finite(lag(lat)),
    dist_km = ifelse(ok,
                     geosphere::distHaversine(cbind(lon, lat), cbind(lag(lon), lag(lat)))/1000,
                     NA_real_),
    is_break = is.na(dt_h) | dt_h > gap_dt_h | dist_km > gap_dist_km,
    seg_id   = cumsum(replace(is_break, is.na(is_break), FALSE))
  ) %>% ungroup() %>% select(-ok)

traj_gap <- d_plot %>% group_by(id) %>% arrange(date_loc, .by_group = TRUE) %>%
  mutate(lon_lag = lag(lon), lat_lag = lag(lat)) %>% ungroup() %>%
  filter(is_break & is.finite(lon_lag) & is.finite(lat_lag))

traj_cont <- d_plot %>% filter(!is_break)

## 9) FOND CARTOGRAPHIQUE ----
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
maroc_claimed <- world %>%
  dplyr::filter(admin %in% c("Morocco","Western Sahara")) %>%
  dplyr::mutate(admin = "Maroc") %>%
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
world_custom <- world %>%
  dplyr::filter(!(admin %in% c("Morocco","Western Sahara"))) %>%
  dplyr::bind_rows(maroc_claimed)
coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")

# Cadre automatique (avec petite marge)
xpad <- max(diff(range(d_plot$lon, na.rm = TRUE))*0.06, 0.3)
ypad <- max(diff(range(d_plot$lat, na.rm = TRUE))*0.06, 0.3)
xlim_vals <- range(d_plot$lon, na.rm = TRUE) + c(-xpad, xpad)
ylim_vals <- range(d_plot$lat, na.rm = TRUE) + c(-ypad, ypad)

## 10) COULEURS & LABELS BALISES ----
pal_id <- scales::hue_pal()(length(unique(d_plot$id)))
names(pal_id) <- sort(unique(d_plot$id))
balise_labels_present <- balise_labels[names(pal_id)]

## 11) CARTE ----
p <- ggplot() +
  geom_sf(data = world_custom, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = coast, color = "grey60", linewidth = 0.2) +
  geom_path(data = traj_cont,
            aes(lon, lat, group = interaction(id, seg_id), colour = id),
            linewidth = 1, alpha = 0.9) +
  geom_segment(data = traj_gap,
               aes(x = lon_lag, y = lat_lag, xend = lon, yend = lat, colour = id),
               linewidth = 0.8, linetype = "22", alpha = 0.9) +
  geom_point(data = d_plot,
             aes(lon, lat, colour = id, shape = lc), size = 2.2, stroke = 0.2) +
  coord_sf(xlim = xlim_vals, ylim = ylim_vals, expand = FALSE) +

  # Légende 1 (balises) — sans titre
  scale_colour_manual(
    name   = NULL,
    values = pal_id,
    labels = balise_labels_present,
    guide  = guide_legend(nrow = 4, byrow = TRUE, order = 1)  # 4 lignes = 4 balises
  ) +

  # Légende 2 (LC) — titre sur sa propre ligne, 2 lignes de symboles
  scale_shape_manual(
    name   = "Classes de localisation\n",
    values = c("1"=16, "2"=17, "3"=15, "0"=3, "A"=8, "B"=13, "C"=7),
    breaks = lc_breaks,
    labels = lab_lc[lc_breaks],
    guide  = guide_legend(order = 2, nrow = 3, byrow = TRUE, title.position = "top")
  ) +

  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 9) +
  theme(
    legend.position  = "right",     # <<< LÉGENDE À DROITE
    legend.box       = "vertical",  # empile les deux blocs
    legend.justification = "top",
    legend.title     = element_text(size = 10, hjust = 0),
    legend.text      = element_text(size = 10),
    legend.key.size  = unit(0.9, "lines"),
    legend.spacing.y = unit(2, "mm"),
    legend.margin    = margin(4, 4, 4, 4, "pt")
  ) +

  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1.2, "cm"), width = unit(1.2, "cm")) +
  annotation_scale(location = "br", width_hint = 0.3,
                 text_cex = 0.7, style = "ticks", line_width = 0.2)
print(p)
## 12) (OPTION) EXPORTS ----
# ggsave("argos_map.png", p, width = 8, height = 9, dpi = 300)
# write.csv(d_final, "argos_clean.csv", row.names = FALSE, fileEncoding = "UTF-8")
