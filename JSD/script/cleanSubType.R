clean.subtypes <- function(name.vec){
  cleaner <- c(
    "NotSunExposedSuprapubic" = "Suprapubic (Not Sun Exposed)",
    "CerebellarHemisphere" = "Cerebellar Hemisphere",
    "Spinalcordcervicalc_1" = "Spinal Cord (Cervical, C1)",
    "Nucleusaccumbensbasalganglia" = "Basal Ganglia (Nucleus Accumbens)",
    "SunExposedLowerleg" = "Lower Leg (Sun Exposed)",
    "EBV_transformedlymphocytes" = "Lymphocytes (EBV-transformed)",
    "GastroesophagealJunction" = "Gastroesoph. Junct.",
    "Substantianigra" = "Substantia Nigra",
    "Putamenbasalganglia" = "Basal Ganglia (Putamen)",
    "FrontalCortexBA9" = "Frontal Cortex (Brodmann A. 9)",
    "MammaryTissue" = "Mammary Tissue",
    "VisceralOmentum" = "Visceral Omentum",
    "TerminalIleum" = "Terminal Ileum",
    "AnteriorcingulatecortexBA24" = "Ant. Cing. Cortex (Brodmann A. 24)",
    "Caudatebasalganglia" = "Basal Ganglia (Caudate)",
    "AtrialAppendage" = "Atrial Appendage",
    "MinorSalivaryGland" = "Minor Salivary Gland",
    "WholeBlood" = "Whole Blood",
    "AdrenalGland" = "Adrenal Gland",
    "LeftVentricle" = "Left Ventricle",
    "SmallIntestine" = "Small Intestine",
    
    "Not_Sun_Exposed_Suprapubic" = "Suprapubic (Not Sun Exposed)",
    "Sun_Exposed_Lower_leg" = "Lower Leg (Sun Exposed)",
    "Gastroesophageal_Junction" = "Gastroesoph. Junct.",
    "Mammary_Tissue" = "Mammary Tissue",
    "Visceral_Omentum" = "Visceral Omentum",
    "Atrial_Appendage" = "Atrial Appendage",
    "MinorSalivaryGland" = "Minor Salivary Gland",
    "WholeBlood" = "Whole Blood",
    "AdrenalGland" = "Adrenal Gland",
    "Left_Ventricle" = "Left Ventricle",
    "SmallIntestine" = "Small Intestine",
    "Cultured_fibroblasts" = "Cultured fibroblasts"
  )
  
  name.vec %>% 
    sapply(
      .,
      function(x){
        ifelse(x %in% names(cleaner),
               cleaner[[x]],
               x
               )
      }
    )
}

