# Timilla Project — Analysis of ARGOS Data
R scripts, test datasets, and tutorials for bird movements monitored with ARGOS devices, applied to Turtle Dove movement ecology (Timilla Project).

Espace collaboratif pour le partage de **scripts R**, **datasets tests** et **tutoriels** destinés à l’analyse des données **ARGOS** appliquées à l’écologie du mouvement de la **Tourterelle des bois (Streptopelia turtur)**.

Ce dépôt s’inscrit dans le cadre du **projet Timilla**, développé en partenariat avec :  
- **BioSphère Environnement** (France)  
- **Office Français de la Biodiversité (OFB)**  
- **Institut Scientifique, Université Mohammed V de Rabat (IS-UM5)** (Maroc)  
- **GREPOM / BirdLife Maroc** (Maroc)  

---

## 🎯 Objectif du projet
Le projet Timilla vise à documenter les **stratégies de migration et d’hivernage** des populations de Tourterelles des bois se reproduisant au Maroc (sous-espèce *arenicola*) grâce à la pose de balises **ARGOS**, afin de :  
- identifier les sites clés d’escale migratoire et d’hivernage de l’espèce sur le continent africain,  
- renforcer les compétences locales en matière de développement de projets similaires et d’analyse de données.  

---

## 📂 Contenu du dépôt
- `scripts/` : scripts R pour le nettoyage et l’analyse de données ARGOS (cartographie et analyse des espaces exploités par l’espèce)  
- `data_example/` : jeux de données anonymisés pour tests  
- `docs/` : tutoriels Quarto ou RMarkdown  

---

## 🔧 Prérequis techniques
- R ≥ 4.3  
- Packages recommandés :  
  ```r
  install.packages(c("tidyverse","sf","ggplot2","lubridate","janitor",
                     "rnaturalearth","rnaturalearthdata","geosphere"))
