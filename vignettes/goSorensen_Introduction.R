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

## ----lib, eval=FALSE----------------------------------------------------------
# if (!requireNamespace("goSorensen", quietly = TRUE)) {
#   BiocManager::install("goSorensen")
# }
# library(goSorensen)

## ----data---------------------------------------------------------------------
data("allOncoGeneLists")

## ----allonco------------------------------------------------------------------
sapply(allOncoGeneLists, length)

# First 15 gene identifiers of gene lists atlas and sanger:
allOncoGeneLists[["atlas"]][1:15]
allOncoGeneLists[["sanger"]][1:15]

## ----bioc, message = FALSE, warning = FALSE, eval = FALSE---------------------
# if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#   BiocManager::install("org.Hs.eg.db")
# }

## ----org, message = FALSE, warning = FALSE, eval = TRUE-----------------------
library(org.Hs.eg.db)

## ----human, message = FALSE, warning = FALSE, echo = TRUE---------------------
humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")

## ----op, echo=FALSE-----------------------------------------------------------
options(max.print = 50)

## ----en-atlas-----------------------------------------------------------------
enrichedAtlas <- enrichedIn(allOncoGeneLists[["atlas"]],
  geneUniverse = humanEntrezIDs,
  orgPackg = "org.Hs.eg.db",
  onto = "BP", GOLevel = 4
)
enrichedAtlas

## ----full-atlas---------------------------------------------------------------
fullEnrichedAtlas <- enrichedIn(allOncoGeneLists[["atlas"]],
  geneUniverse = humanEntrezIDs,
  orgPackg = "org.Hs.eg.db",
  onto = "BP", GOLevel = 4,
  onlyEnriched = FALSE
)
fullEnrichedAtlas

## ----en-atlas-BP--------------------------------------------------------------
enrichedAtlasBP <- enrichedIn(allOncoGeneLists[["atlas"]],
  geneUniverse = humanEntrezIDs,
  orgPackg = "org.Hs.eg.db",
  onto = "BP", GOLevel = NULL
)
enrichedAtlasBP

## ----atr-atlas----------------------------------------------------------------
attr(enrichedAtlasBP, "nTerms")

## ----n-terms------------------------------------------------------------------
# number of GO terms in enrichedAtlas
length(enrichedAtlas)

# number of GO terms in fullEnrichedAtlas
length(fullEnrichedAtlas)

## ----symbol,message=FALSE-----------------------------------------------------
library(AnnotationDbi)

# Convert all gene lists from ENTREZID to SYMBOL
allOncoGeneListsSYMBOL <- lapply(allOncoGeneLists, function(geneList) {
  symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
    keys = geneList,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  unique(na.omit(symbols))
})

# Obtain the human gene universe using SYMBOL identifiers
humanSymbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

## ----en-atlas-symbol----------------------------------------------------------
enrichedAtlasSymbolBP4 <- enrichedIn(allOncoGeneListsSYMBOL$atlas,
  geneUniverse = humanSymbols,
  orgPackg = "org.Hs.eg.db",
  keyType = "SYMBOL",
  onto = "BP", GOLevel = 4
)
enrichedAtlasSymbolBP4

## ----op-100, echo=FALSE-------------------------------------------------------
options(max.print = 100)

## ----en-BP4, eval=FALSE-------------------------------------------------------
# enrichedInBP4 <- enrichedIn(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP", GOLevel = 4
# )
# enrichedInBP4

## ----en-b-BP4, echo=FALSE-----------------------------------------------------
data("enrichedInBP4")
enrichedInBP4

## ----f-en-BP4, eval=FALSE-----------------------------------------------------
# fullEnrichedInBP4 <- enrichedIn(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP", GOLevel = 4,
#   onlyEnriched = FALSE
# )
# fullEnrichedInBP4

## ----f-en-b-BP4, echo=FALSE---------------------------------------------------
data("fullEnrichedInBP4")
fullEnrichedInBP4

## ----en-BP, eval=FALSE--------------------------------------------------------
# enrichedInBP <- enrichedIn(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP", GOLevel = NULL
# )
# enrichedInBP

## ----en-b-BP,echo=FALSE-------------------------------------------------------
data("enrichedInBP")
enrichedInBP

## ----nrow-BP4-----------------------------------------------------------------
# number of GO terms (rows) in enrichedInBP4
nrow(enrichedInBP4)

# number of GO terms (rows) in fullEnrichedInBP4
nrow(fullEnrichedInBP4)

## ----cont-atlas-sanger-BP4, eval=FALSE----------------------------------------
# cont_atlas.sanger_BP4 <- buildEnrichTable(allOncoGeneLists$atlas,
#   allOncoGeneLists$sanger,
#   listNames = c("atlas", "sanger"),
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP",
#   GOLevel = 4
# )
# cont_atlas.sanger_BP4

## ----cont-b-atlas-sanger-BP4, echo=FALSE--------------------------------------
data("cont_atlas.sanger_BP4")
cont_atlas.sanger_BP4

## ----cont-atlas-sanger-BP,eval=FALSE------------------------------------------
# cont_atlas.sanger_BP <- buildEnrichTable(allOncoGeneLists$atlas,
#   allOncoGeneLists$sanger,
#   listNames = c("atlas", "sanger"),
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP",
#   GOLevel = NULL
# )
# cont_atlas.sanger_BP

## ----cont-b-atlas-sanger-BP,echo=FALSE----------------------------------------
data("cont_atlas.sanger_BP")
cont_atlas.sanger_BP

## ----cont-all-BP4, eval=FALSE-------------------------------------------------
# cont_all_BP4 <- buildEnrichTable(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP",
#   GOLevel = 4
# )

## ----cont-b-all-BP4,echo=FALSE, eval=FALSE------------------------------------
# data("cont_all_BP4")
# cont_all_BP4

## ----cont-all-BP, eval=FALSE--------------------------------------------------
# cont_all_BP <- buildEnrichTable(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP",
#   GOLevel = NULL
# )

## ----cont-b-all-BP,echo=FALSE,eval=FALSE--------------------------------------
# data("cont_all_BP")
# cont_all_BP

## ----allcontTabs, eval=FALSE--------------------------------------------------
# allContTabs <- allBuildEnrichTable(allOncoGeneLists,
#   geneUnierse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   ontos = c("BP", "CC", "MF"),
#   GOLevels = 3:10
# )

## ----allcontTabsNoLevel,eval=FALSE--------------------------------------------
# allContTabsNoLevel <- allBuildEnrichTable(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   ontos = c("BP", "CC", "MF"),
#   GOLevels = NULL
# )

## ----eqTest-atlas-sanger-BP4, eval=FALSE--------------------------------------
# eqTest_atlas.sanger_BP4 <- equivTestSorensen(allOncoGeneLists$atlas,
#   allOncoGeneLists$sanger,
#   listNames = c("atlas", "sanger"),
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP", GOLevel = 4,
#   d0 = 0.4444,
#   conf.level = 0.95
# )
# eqTest_atlas.sanger_BP4

## ----eqTest-b-atlas-sanger-BP4, echo=FALSE------------------------------------
data("eqTest_atlas.sanger_BP4")
eqTest_atlas.sanger_BP4

## ----eqTest-atlas-sanger-BP,eval=FALSE----------------------------------------
# eqTest_atlas.sanger_BP <- equivTestSorensen(allOncoGeneLists$atlas,
#   allOncoGeneLists$sanger,
#   listNames = c("atlas", "sanger"),
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP", GOLevel = NULL,
#   d0 = 0.4444,
#   conf.level = 0.95
# )
# eqTest_atlas.sanger_BP

## ----eqTest-b-atlas-sanger-BP,echo=FALSE--------------------------------------
data("eqTest_atlas.sanger_BP")
eqTest_atlas.sanger_BP

## ----equivtestSorensen--------------------------------------------------------
equivTestSorensen(cont_atlas.sanger_BP4,
  d0 = 0.4444,
  conf.level = 0.95
)

## ----upgrade-t----------------------------------------------------------------
upgrade(eqTest_atlas.sanger_BP4,
  d0 = 0.2857,
  conf.level = 0.99, boot = TRUE
)

## ----eqTest-all-BP4, eval=FALSE-----------------------------------------------
# eqTest_all_BP4 <- equivTestSorensen(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP",
#   GOLevel = 4,
#   d0 = 0.4444,
#   conf.level = 0.95
# )

## ----alleqTestsBP,eval=FALSE--------------------------------------------------
# allEqTestsBP <- equivTestSorensen(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   onto = "BP",
#   GOLevel = NULL,
#   d0 = 0.4444,
#   conf.level = 0.95
# )

## ----eqTest-b-all-BP4, eval=FALSE---------------------------------------------
# eqTest_all_BP4 <- equivTestSorensen(cont_all_BP4,
#   d0 = 0.4444,
#   conf.level = 0.95
# )

## ----eqTest-db-all-BP4, echo=FALSE--------------------------------------------
data(eqTest_all_BP4)

## ----op-d, echo=FALSE---------------------------------------------------------
options(digits = 4)

## ----getdissimilarity---------------------------------------------------------
getDissimilarity(eqTest_all_BP4, simplify = FALSE)

## ----getpvalue----------------------------------------------------------------
getPvalue(eqTest_all_BP4, simplify = FALSE)

## ----alleqTests, eval=FALSE---------------------------------------------------
# allEqTests <- allEquivTestSorensen(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   ontos = c("BP", "CC", "MF"),
#   GOLevels = 3:10,
#   d0 = 0.4444,
#   conf.level = 0.95
# )

## ----alleqTestsNoLevel,eval=FALSE---------------------------------------------
# allEqTestsNoLevel <- allEquivTestSorensen(allOncoGeneLists,
#   geneUniverse = humanEntrezIDs,
#   orgPackg = "org.Hs.eg.db",
#   ontos = c("BP", "CC", "MF"),
#   GOLevels = NULL,
#   d0 = 0.4444,
#   conf.level = 0.95
# )

## ----all-b-eqTests, eval=FALSE------------------------------------------------
# allEqTests <- allEquivTestSorensen(allContTabs,
#   d0 = 0.4444,
#   conf.level = 0.95
# )

## ----session-info-------------------------------------------------------------
sessionInfo()

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

