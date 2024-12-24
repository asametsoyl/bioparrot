# 3D Molek??ler G??rselle??tirme Mod??l?? UI
StructureViewerUI <- function(id) {
  ns <- NS(id)  # Namespace tan??mlamas??
  tagList(
    fileInput(ns("pdbFile"), "Upload PDB File", accept = c(".pdb")),  # PDB dosyas??n?? y??klemek i??in
    selectInput(ns("viewStyle"), "Select View Style", 
                choices = c("Cartoon", "Stick", "Sphere", "Surface", "Line"),
                selected = "Cartoon"),
    selectInput(ns("colorBy"), "Color By", 
                choices = c("Atom Type", "Residue"),
                selected = "Atom Type"),
    r3dmolOutput(ns("structure3D"), height = "600px"),  # 3D g??rselle??tirme
    downloadButton(ns("downloadModel"), "Download Model")
  )
}
