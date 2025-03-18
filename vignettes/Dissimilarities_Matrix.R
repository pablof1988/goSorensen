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

## -----------------------------------------------------------------------------
data("allOncoGeneLists")

## -----------------------------------------------------------------------------
# name and length of the gene lists:
sapply(allOncoGeneLists, length)

## ----warning=FALSE, message=FALSE, eval=FALSE---------------------------------
# # Previously load the genomic annotation package for the studied specie:
# library(org.Hs.eg.db)
# humanEntrezIDs <- keys(org.Hs.eg.db, keytype = "ENTREZID")
# 
# # Computing the irrelevance-threshold matrix of dissimilarities
# dissMatrx_BP4 <- sorenThreshold(allOncoGeneLists,
#                             geneUniverse = humanEntrezIDs,
#                             orgPackg = "org.Hs.eg.db",
#                             onto = "BP",
#                             GOLevel = 4,
#                             trace = FALSE)
# dissMatrx_BP4

## ----echo=FALSE---------------------------------------------------------------
options(digits = 4)
data("dissMatrx_BP4")
data("cont_all_BP4")
dissMatrx_BP4

## ----warning=FALSE, message=FALSE, comment=NA, eval=FALSE---------------------
# # Previously compute the enrichment contingency tables:
# cont_all_BP4 <- buildEnrichTable(allOncoGeneLists,
#                                  geneUniverse = humanEntrezIDs,
#                                  orgPackg = "org.Hs.eg.db",
#                                  onto = "BP",
#                                  GOLevel = 4)

## -----------------------------------------------------------------------------
# Computing the irrelevance-threshold matrix of dissimilarities from the 
# enrichment contingency table "cont_all_BP4":
dissMatrx_BP4 <- sorenThreshold(cont_all_BP4, 
                                trace = FALSE)
dissMatrx_BP4

## ----eval=FALSE---------------------------------------------------------------
# boot_dissMatrx_BP4 <- sorenThreshold(allOncoGeneLists,
#                             geneUniverse = humanEntrezIDs,
#                             orgPackg = "org.Hs.eg.db",
#                             onto = "BP",
#                             GOLevel = 4,
#                             boot = TRUE, # use bootstrap distribution
#                             trace = FALSE)
# boot_dissMatrx_BP4

## ----echo=FALSE---------------------------------------------------------------
sorenThreshold(cont_all_BP4, boot = TRUE, trace = FALSE)

## -----------------------------------------------------------------------------
boot_dissMatrx_BP4 <- sorenThreshold(cont_all_BP4, 
                                     boot = TRUE, 
                                     trace = FALSE)
boot_dissMatrx_BP4

## ----warning=FALSE, message=FALSE, comment=NA, eval=FALSE---------------------
# # For example, for GO levels 3 to 10 and for the ontologies BP, CC and MF:
# allDissMatrx <- allSorenThreshold(allOncoGeneLists,
#                                   geneUniverse = humanEntrezIDs,
#                                   orgPackg = "org.Hs.eg.db",
#                                   ontos = c("BP", "CC", "MF"),
#                                   GOLevels = 3:10,
#                                   trace = FALSE)

## ----warning=FALSE, message=FALSE, comment=NA, eval=FALSE---------------------
# # Previously compute the enrichment contingency tables:
# allContTabs <- allBuildEnrichTable(allOncoGeneLists,
#                                    geneUniverse = humanEntrezIDs,
#                                    orgPackg = "org.Hs.eg.db",
#                                    ontos = c("BP", "CC", "MF"),
#                                    GOLevels = 3:10)
# 
# # Computing the irrelevance-threshold matrix of dissimilarities from the
# # enrichment contingency tables "allContTabs":
# allDissMatrx <- allSorenThreshold(allContTabs,
#                                   trace = FALSE)

## ----warning=FALSE, message=FALSE---------------------------------------------
clust.threshold <- hclustThreshold(dissMatrx_BP4)
plot(clust.threshold)

## ----warning=FALSE, message=FALSE---------------------------------------------
# multidimensional scaling analysis:
mds <- as.data.frame(cmdscale(dissMatrx_BP4, k = 2))

## ----warning=FALSE, message=FALSE, fig.align='left', fig.height=3, fig.width=5----
library(ggplot2)
library(ggrepel)
graph <- ggplot() +
 geom_point(aes(mds[,1], mds[,2]), color = "blue", size = 3) +
 geom_text_repel(aes(mds[,1], mds[,2], label = attr(dissMatrx_BP4, "Labels")),
                 color = "black", size = 3) +
 xlab("Dim 1") +
 ylab("Dim 2") +
 theme_light()
graph

## -----------------------------------------------------------------------------
# Split the axis 20% to the left, 60% to the middle and 20% to the right:
prop <- c(0.2, 0.6, 0.2) 

# Sort according  dimension 1:
sorted <- mds[order(mds[, 1]), ] 

# Determine the range for dimension 1.
range <- sorted[, 1][c(1, nrow(mds))] 

# Find the cut-points to split the axis:
cutpoints <- (cumsum(prop)[1:2] * diff(range)) + range[1]

# Identify lists to the left:  
lleft <- rownames(sorted[sorted[, 1] < cutpoints[1], ])
lleft

# Identify lists to the right
lright <- rownames(sorted[sorted[, 1] > cutpoints[2], ])
lright

## ----warning=FALSE, message=FALSE, fig.align='left', fig.height=3, fig.width=5----
graph +
  geom_vline(xintercept = cutpoints, color = "red", 
             linetype = "dashed", linewidth = 0.75)

## ----echo=FALSE---------------------------------------------------------------
options(max.print = 16)

## -----------------------------------------------------------------------------
# enrichment contingency tables:
contTablesBP4 <- attr(dissMatrx_BP4, "all2x2Tables")

# matrix of GO term enrichment for the lists located at the extreme left
tableleft <- attr(contTablesBP4, "enriched")[, lleft, drop = FALSE]
tableleft

## ----echo=FALSE---------------------------------------------------------------
options(max.print = 50)

## -----------------------------------------------------------------------------
# matrix of GO term enrichment for the lists located at the extreme right
tableright <- attr(contTablesBP4, "enriched")[, lright, drop = FALSE]
tableright

## ----warning=FALSE, message=FALSE---------------------------------------------
# function to compute mean and standard deviation:
mean_sd <- function(x){
  c("mean" = mean(x), "sd" = ifelse(length(x) == 1, 0, sd(x)))
}

# mean and sd of the lists to the extreme left:
lmnsd <- apply(tableleft, 1, mean_sd)

# mean and sd of the lists to the extreme right:
rmnsd <- apply(tableright, 1, mean_sd)

## ----warning=FALSE, message=FALSE---------------------------------------------
nl <- ncol(tableleft) 
nr <- ncol(tableright) 
t_stat <- abs(lmnsd[1, ] - rmnsd[1, ]) / 
  sqrt((((lmnsd[2, ] / nl) + (rmnsd[2, ] / nr))) + 0.00000001)

## ----warning=FALSE, message=FALSE---------------------------------------------
result <- t_stat[t_stat == max(t_stat)]
result

## ----warning=FALSE, message=FALSE---------------------------------------------
library(GO.db)
library(DT)
# Previous function to get the description of the identified GO terms:
get_go_description <- function(go_id) {
  go_term <- Term(GOTERM[[go_id]])
  return(go_term)
}

# GO terms description:
datatable(data.frame(GO_term = names(result),
                     Description = sapply(names(result), get_go_description, 
                                           USE.NAMES = TRUE)), 
          filter = "top", rownames = FALSE)

