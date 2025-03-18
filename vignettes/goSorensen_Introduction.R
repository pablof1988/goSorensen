## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----css, echo=FALSE, results='asis'------------------------------------------
cat("
<style>
/* Ocultar inicialmente el contenido colapsable */
.collapsible-content {
    display: none;
}

/* Cambiar el puntero del ratón en los encabezados colapsables */
.collapsible-header {
    cursor: pointer;
}

.collapsible-header:hover {
    text-decoration: underline;
}

/* Título de Nivel 1 */
h1 {
    font-size: 1.2em; 
    color: #0e5775; 
    font-weight: bold;
}

/* Título de Nivel 2 */
h2.collapsible-header {
    font-size: 1.1em; 
    color: #12769f; 
    margin-left: 20px; 
    font-weight: bold;
}

/* Título de Nivel 3 */
h3.collapsible-header {
    font-size: 1.0em; 
    color: #16a3dc; 
    margin-left: 40px; 
    font-weight: bold;
}

/* Mantén el contenido alineado con los encabezados */
.collapsible-content {
    margin-left: inherit;
}

#appendix-appendix {
    display: none;
}

/* Estilo para el contenedor de referencias */
#refs {
    list-style-type: none; 
    margin-left: 20px;
    padding-left: 20px; 
    text-align: justify;
}


/* Estilo para la viñeta */
#refs > div::before {
    content: '•'; /* Define la viñeta como un círculo */
    font-size: 1.2em; /* Tamaño de la viñeta */
    color: #000; /* Color de la viñeta */
    margin-right: 10px; /* Espacio entre la viñeta y el texto */
    flex-shrink: 0; /* Asegura que la viñeta no se reduce en líneas largas */
}

.sidenote {
    text-align: justify; 
    float: right; /* Posiciona el pie de página en el lado derecho */
    width: 28%; 
    max-width: 28%; /* Ajusta el ancho relativo al contenedor principal */
    padding-left: 20px; /* Añade margen interno dentro del pie de página */
    box-sizing: border-box; /* Asegura que padding no desborde el ancho del contenedor */
}

</style>
")

## ----js, echo=FALSE, results='asis'-------------------------------------------
cat("
<script>
document.addEventListener('DOMContentLoaded', function () {
    const excludedSections = ['Abstract', 'Package', 'Author Information', 'Date'];

    const headers = document.querySelectorAll('h2, h3');
    headers.forEach(header => {
        const headerText = header.textContent.trim();

        // Excluir encabezados específicos
        if (excludedSections.some(section => headerText.includes(section))) {
            return;
        }

        // Encuentra todo el contenido hasta el próximo encabezado
        const content = [];
        let next = header.nextElementSibling;
        while (next && !/^H[1-6]$/.test(next.tagName)) {
            content.push(next);
            next = next.nextElementSibling;
        }

        // Crear contenedor colapsable
        const contentWrapper = document.createElement('div');
        contentWrapper.className = 'collapsible-content';
        content.forEach(node => contentWrapper.appendChild(node));

        // Agregar funcionalidad de colapsar/expandir
        header.classList.add('collapsible-header');
        header.addEventListener('click', () => {
            const isVisible = contentWrapper.style.display === 'block';
            contentWrapper.style.display = isVisible ? 'none' : 'block';
        });

        // Insertar el contenido después del encabezado
        header.parentNode.insertBefore(contentWrapper, next);
    });
});
</script>
")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  comment = ""
)

## ----env, message = FALSE, warning = FALSE, echo = TRUE-----------------------
library(goSorensen)

## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("goSorensen", quietly = TRUE)) {
#     BiocManager::install("goSorensen")
# }
# library(goSorensen)

## -----------------------------------------------------------------------------
data("allOncoGeneLists")

## -----------------------------------------------------------------------------
sapply(allOncoGeneLists, length)

# First 15 gene identifiers of gene lists atlas and sanger:
allOncoGeneLists[["atlas"]][1:15]
allOncoGeneLists[["sanger"]][1:15]

## ----message = FALSE, warning = FALSE, eval = FALSE---------------------------
# if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#     BiocManager::install("org.Hs.eg.db")
# }

## ----message = FALSE, warning = FALSE, eval = TRUE----------------------------
library(org.Hs.eg.db)

## ----message = FALSE, warning = FALSE, echo = TRUE----------------------------
humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")

## ----echo=FALSE---------------------------------------------------------------
options(max.print = 50)

## ----eval=TRUE----------------------------------------------------------------
enrichedAtlas <- enrichedIn(allOncoGeneLists[["atlas"]],
           geneUniverse = humanEntrezIDs, 
           orgPackg = "org.Hs.eg.db",
           onto = "BP", GOLevel = 4)
enrichedAtlas

## -----------------------------------------------------------------------------
fullEnrichedAtlas <- enrichedIn(allOncoGeneLists[["atlas"]],
           geneUniverse = humanEntrezIDs, 
           orgPackg = "org.Hs.eg.db",
           onto = "BP", GOLevel = 4, 
           onlyEnriched = FALSE)
fullEnrichedAtlas

## -----------------------------------------------------------------------------
# number of GO terms in enrichedAtlas
length(enrichedAtlas)

# number of GO terms in fullEnrichedAtlas
length(fullEnrichedAtlas)

## ----echo=FALSE---------------------------------------------------------------
options(max.print = 100)

## ----eval=FALSE---------------------------------------------------------------
# enrichedInBP4 <- enrichedIn(allOncoGeneLists,
#                                geneUniverse = humanEntrezIDs,
#                                orgPackg = "org.Hs.eg.db",
#                                onto = "BP", GOLevel = 4)
# enrichedInBP4

## ----echo=FALSE---------------------------------------------------------------
data("enrichedInBP4")
enrichedInBP4

## ----eval=FALSE---------------------------------------------------------------
# fullEnrichedInBP4 <- enrichedIn(allOncoGeneLists,
#                                geneUniverse = humanEntrezIDs,
#                                orgPackg = "org.Hs.eg.db",
#                                onto = "BP", GOLevel = 4,
#                                onlyEnriched = FALSE)
# fullEnrichedInBP4

## ----echo=FALSE---------------------------------------------------------------
data("fullEnrichedInBP4")
fullEnrichedInBP4

## -----------------------------------------------------------------------------
# number of GO terms (rows) in enrichedInBP4
nrow(enrichedInBP4)

# number of GO terms (rows) in fullEnrichedInBP4
nrow(fullEnrichedInBP4)

## ----eval=FALSE---------------------------------------------------------------
# cont_atlas.sanger_BP4 <- buildEnrichTable(allOncoGeneLists$atlas,
#                                           allOncoGeneLists$sanger,
#                                           listNames = c("atlas", "sanger"),
#                                           geneUniverse = humanEntrezIDs,
#                                           orgPackg = "org.Hs.eg.db",
#                                           onto = "BP",
#                                           GOLevel = 4)
# cont_atlas.sanger_BP4

## ----echo=FALSE---------------------------------------------------------------
data("cont_atlas.sanger_BP4")
cont_atlas.sanger_BP4

## ----eval=FALSE---------------------------------------------------------------
# cont_all_BP4 <- buildEnrichTable(allOncoGeneLists,
#                                  geneUniverse = humanEntrezIDs,
#                                  orgPackg = "org.Hs.eg.db",
#                                  onto = "BP",
#                                  GOLevel = 4)

## ----eval=FALSE---------------------------------------------------------------
# allContTabs <- allBuildEnrichTable(allOncoGeneLists,
#                                    geneUniverse = humanEntrezIDs,
#                                    orgPackg = "org.Hs.eg.db",
#                                    ontos = c("BP", "CC", "MF"),
#                                    GOLevels = 3:10)

## ----eval=FALSE---------------------------------------------------------------
# eqTest_atlas.sanger_BP4 <- equivTestSorensen(allOncoGeneLists$atlas,
#                                              allOncoGeneLists$sanger,
#                                              listNames = c("atlas", "sanger"),
#                                              geneUniverse = humanEntrezIDs,
#                                              orgPackg = "org.Hs.eg.db",
#                                              onto = "BP", GOLevel = 4,
#                                              d0 = 0.4444,
#                                              conf.level = 0.95)
# eqTest_atlas.sanger_BP4

## ----echo=FALSE---------------------------------------------------------------
data("eqTest_atlas.sanger_BP4")
eqTest_atlas.sanger_BP4

## -----------------------------------------------------------------------------
equivTestSorensen(cont_atlas.sanger_BP4, 
                  d0 = 0.4444, 
                  conf.level = 0.95)

## -----------------------------------------------------------------------------
upgrade(eqTest_atlas.sanger_BP4, d0 = 0.2857, 
        conf.level = 0.99, boot = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# eqTest_all_BP4 <- equivTestSorensen(allOncoGeneLists,
#                                     geneUniverse = humanEntrezIDs,
#                                     orgPackg = "org.Hs.eg.db",
#                                     onto = "BP",
#                                     GOLevel = 4,
#                                     d0 = 0.4444,
#                                     conf.level = 0.95)

## ----eval=FALSE---------------------------------------------------------------
# eqTest_all_BP4 <- equivTestSorensen(cont_all_BP4,
#                                     d0 = 0.4444,
#                                     conf.level = 0.95)

## ----echo=FALSE---------------------------------------------------------------
data(eqTest_all_BP4)

## ----echo=FALSE---------------------------------------------------------------
options(digits = 4)

## -----------------------------------------------------------------------------
getDissimilarity(eqTest_all_BP4, simplify = FALSE)

## -----------------------------------------------------------------------------
getPvalue(eqTest_all_BP4, simplify = FALSE)

## ----eval=FALSE---------------------------------------------------------------
# allEqTests <- allEquivTestSorensen(allOncoGeneLists,
#                                    geneUniverse = humanEntrezIDs,
#                                    orgPackg = "org.Hs.eg.db",
#                                    ontos = c("BP", "CC", "MF"),
#                                    GOLevels = 3:10,
#                                    d0 = 0.4444,
#                                    conf.level = 0.95)

## ----eval=FALSE---------------------------------------------------------------
# allEqTests <- allEquivTestSorensen(allContTabs,
#                                    d0 = 0.4444,
#                                    conf.level = 0.95)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

