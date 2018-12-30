
# Constants ----

.geneSetNameInput <- "geneSetName"
.shinyLabelsPlotChoices <- c(
    "Barplot (#)"="barplotPredictionCount",
    "Barplot (%)"="barplotPredictionProportion")

# Plotting functions ----

#' Generate the panels in the app body
#'
#' Constructs the panels in the main body of the app to examine and rename gene signatures.
#'
#' @rdname INTERNAL_panelGeneration
#'
#' @param genesets A set of signatures of class inheriting from "\code{\link{tbl_geneset}}".
#'
#' @return
#' A HTML tag object containing the UI elements for the main body of the app.
#' This includes the output plots/tables as well as UI elements to control them.
#' @importFrom shiny plotOutput
.panelGeneration <- function(genesets) {
    NSETS <- nlevels(genesets$set)
    panelList <- list()
    for (id in seq_len(NSETS)) {
        id0 <- id
        geneSetName0 <- levels(genesets$set)[id0]
        geneIds0 <- genesets[genesets$set == geneSetName0, "gene", drop=TRUE]
        geneIdText <- paste(geneIds0, collapse=", ")
        plotName0 <- paste0("plot", id0)
        signaturePlot <- plotOutput(plotName0, height="400px")
        panelList[[id0]] <- box(
            textInput(
                inputId=paste0(.geneSetNameInput, id0),
                label=paste("Gene set name", id0),
                value=geneSetName0, width="100%",
                placeholder=paste("Gene set name", id0)
            ),
            column(signaturePlot, width=12),
            HTML(sprintf("<strong>Features:</strong><p>%s<p>", geneIdText)),
            width=12, title=paste("Gene set", id0)
        )
    }
    outValue <- do.call(fluidRow, panelList)

    return(outValue)
}

# App ----

#nocov start
#' Interactively Inspect and Name Gene Signatures
#'
#' This function launches a Shiny app to inspect the gene signatures defined in a \code{tbl_geneset}
#' for the purpose of (re-)naming those signatures interactively.
#' The app returns the updated \code{tbl_geneset} when closed using the \code{"Done"} button.
#'
#' @param gs A set of gene signatures inheriting from \code{\link{tbl_geneset}}.
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#'
#' @return The updated set of gene signatures as a \code{\link{tbl_geneset}}.
#' @export
#' @importFrom methods is
#' @importFrom shiny shinyApp reactiveValues observeEvent stopApp isolate
#' fluidRow column icon textInput actionButton renderUI renderPlot
#' HTML uiOutput
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#' dashboardBody box menuItem
#'
#' @examples
#' # Example data ----
#' library(GeneSet)
#' tgs <- tbl_geneset(
#'     "Cell type 1"=c("Gene001", "Gene002"),
#'     "Cell type 2"=c("Gene002", "Gene003", "Gene004")
#' )
#'
#' library(SummarizedExperiment)
#' nsamples <- 100
#' u <- matrix(rpois(20000, 2), ncol=nsamples)
#' rownames(u) <- paste0("Gene", sprintf("%03d", seq_len(nrow(u))))
#' colnames(u) <- paste0("Cell", sprintf("%03d", seq_len(ncol(u))))
#' se <- SummarizedExperiment(assays=list(counts=u))
#' colData(se)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se), replace=TRUE))
#'
#' se1 <- predict(tgs, se, method="ProportionPositive", cluster.col="cluster")
#'
#' # Example usage ----
#' if (interactive()){
#'     x <- runApp(shinyLabels(tgs, se1))
#' }
shinyLabels <- function(gs, se) {

    stopifnot(is(gs, "tbl_geneset"))
    stopifnot(identical(
        levels(gs$set),
        levels(se$Hancock$prediction)
    ))

    REACTIVE <- reactiveValues(
        GS=gs,
        SE=se
    )

    app_ui <- dashboardPage(
        dashboardHeader(),
        dashboardSidebar(
            actionButton(inputId="Done", label="Done", icon=icon("sign-out"), width="50%"),
            selectizeInput(inputId="Plot", label="Plot", choices=.shinyLabelsPlotChoices, selected="barplotPredictions")
        ),
        dashboardBody(
            uiOutput("mainPanels")
        ),
        title="Hancock: Label signatures"
    )

    app_server <- function(input, output, session) {

        # App-specific constants ----

        NSETS <- nlevels(gs$set)

        # Main panels ----

        panelList <- list()
        for (id in seq_len(NSETS)) {
            local({
                id0 <- id
                plotName0 <- paste0("plot", id0)
                output[[plotName0]] <- renderPlot({
                    plotFUN <- get(input[["Plot"]])
                    geneSetName0 <- levels(REACTIVE$GS$set)[id0]
                    plotFUN(REACTIVE$SE, geneSetName0)
                })
            })
        }

        output$mainPanels <- renderUI({
            .panelGeneration(REACTIVE$GS)
        })

        # Observer for the gene set names ----

        for (id in seq_len(NSETS)) {
            local({
                id0 <- id
                inputId0 <- paste0(.geneSetNameInput, id0)
                observeEvent(input[[inputId0]], {
                    newValue <- input[[inputId0]]
                    levels(REACTIVE$GS$set)[id0] <- as.character(newValue)
                    levels(REACTIVE$SE$Hancock$prediction)[id0] <- as.character(newValue)
                })
            })
        }

        # Observer for closing the app and returning the updated object ----

        observeEvent(input$Done, {
            stopApp(invisible(REACTIVE$GS))
        })

    }

    shinyApp(ui=app_ui, server=app_server)
}
#nocov end
