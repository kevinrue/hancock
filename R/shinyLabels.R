
# Constants ----

.doneInput <- "Done"
.resetInput <- "Reset"
.tourInput <- "Tour"

.plotOutput <- "Plot"

.geneSetNameInput <- "GeneSetName"
.plotFunction <- "PlotType"
.redDimTypeInput <- "RedDimType"
.xAxisInput <- "XAxis"
.yAxisInput <- "YAxis"
.showLabelsInput <- "ShowLabels"

.shinyLabelsPlotChoices <- c(
    "Barplot (#)"="barplotPredictionCount",
    "Barplot (%)"="barplotPredictionProportion")

.actionbutton_biocstyle <- "color: #ffffff; background-color: #0092AC; border-color: #2e6da4"


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
        plotName0 <- paste0(.plotOutput, id0)
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
            width=12, title=paste("Gene set", id0)#, collapsible=TRUE, collapsed=(id0 > 1L)
            # TODO: remember whether each box should be open/closed at each refresh
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
#' @importFrom methods is as
#' @importFrom shiny shinyApp reactiveValues observeEvent stopApp isolate
#' fluidRow column icon textInput actionButton renderUI renderPlot
#' HTML uiOutput checkboxInput selectizeInput conditionalPanel
#' updateSelectizeInput
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#' dashboardBody box menuItem dropdownMenu notificationItem
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDimNames
#' @importFrom rintrojs introjsUI introjs
#' @importFrom utils read.delim
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
#' # Add reduced dimension results to enable app features
#' library(SingleCellExperiment)
#' library(Rtsne)
#' sce1 <- as(se1, "SingleCellExperiment")
#' reducedDim(sce1, "PCA") <- prcomp(t(assay(sce1)))$x
#' reducedDim(sce1, "TSNE") <- Rtsne(X=reducedDim(sce1, "PCA"))$Y
#'
#' # Example usage ----
#' if (interactive()){
#'     library(shiny)
#'     x <- runApp(shinyLabels(tgs, sce1))
#' }
shinyLabels <- function(gs, se) {

    # Initialization of reactive values ----

    stopifnot(is(gs, "tbl_geneset"))
    stopifnot(identical(
        levels(gs$set),
        levels(se$Hancock$prediction)
    ))

    se <- as(se, "SingleCellExperiment")

    # Storage for all the reactive objects
    REACTIVE <- reactiveValues(
        GS=gs,
        SE=se
    )

    # Initialization of non-reactive values ----

    NSETS <- nlevels(gs$set)

    # Storage for persistent non-reactive objects.
    pObjects <- new.env()

    # TODO: "expand all"/"collapse all" button
    # pObjects$collapsedBox <- c(FALSE, rep(TRUE, NSETS-1))

    # Initial state for reduced dimension inputs ----
    xAxisMax <- 0L
    yAxisMax <- 0L
    if (!is.null(reducedDim(se))) {
        .shinyLabelsPlotChoices <- c(
            .shinyLabelsPlotChoices,
            "reducedDim"="reducedDimPrediction")
        xAxisMax <- yAxisMax <- ncol(reducedDim(se, 1L))
    }

    app_ui <- dashboardPage(
        dashboardHeader(
            dropdownMenu(
                type="tasks",
                icon=icon("question-circle fa-1g"),
                badgeStatus=NULL,
                headerText="Documentation",
                notificationItem(
                    text=actionButton(
                        .tourInput, "Interactive tour",
                        icon("hand-o-right"),
                        style=.actionbutton_biocstyle
                    ),
                    icon=icon(""), # tricking it to not have additional icon
                    status="primary"
                )
            )
        ),
        dashboardSidebar(
            actionButton(inputId=.doneInput, label="Done", icon=icon("sign-out"), width="50%"),
            actionButton(inputId=.resetInput, label="Reset", icon=icon("undo"), width="50%"),
            selectizeInput(inputId=.plotFunction, label="Plot type:", choices=.shinyLabelsPlotChoices, selected="barplotPredictions"),
            checkboxInput(.showLabelsInput, "Show labels", TRUE),
            conditionalPanel(
                condition = sprintf("input.%s == 'reducedDimPrediction'", .plotFunction),
                selectizeInput(
                    .redDimTypeInput, "Type",
                    reducedDimNames(se)
                ),
                selectizeInput(
                    .xAxisInput, "x-axis",
                    seq_len(xAxisMax),
                    min(1L, xAxisMax, na.rm=TRUE)
                ),
                selectizeInput(
                    .yAxisInput, "y-axis",
                    seq_len(yAxisMax),
                    min(2L, yAxisMax, na.rm=TRUE)
                )
            )
        ),
        dashboardBody(
            introjsUI(), # must be included in UI
            uiOutput("mainPanels")
        ),
        title="Hancock: Label signatures"
    )

    app_server <- function(input, output, session) {

        # Main panels ----

        panelList <- list()
        for (id in seq_len(NSETS)) {
            local({
                id0 <- id
                plotName0 <- paste0(.plotOutput, id0)
                output[[plotName0]] <- renderPlot({
                    geneSetName0 <- levels(REACTIVE$GS$set)[id0]
                    .plotWrapper(
                        REACTIVE$SE, geneSetName0,
                        plotType=input[[.plotFunction]],
                        labels=input[[.showLabelsInput]],
                        redDimType=input[[.redDimTypeInput]],
                        x=as.numeric(input[[.xAxisInput]]),
                        y=as.numeric(input[[.yAxisInput]])
                    )
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

        # Observer to update the available x/y-axis dropdown menus ----

        observeEvent(ignoreInit = TRUE, input[[.redDimTypeInput]], {
            redDimType <- input[[.redDimTypeInput]]
            updateSelectizeInput(
                session, .xAxisInput,
                choices=seq_len(ncol(reducedDim(se, redDimType))),
                selected=min(1L, ncol(reducedDim(se, redDimType)), na.rm=TRUE))
            updateSelectizeInput(
                session, .yAxisInput,
                choices=seq_len(ncol(reducedDim(se, redDimType))),
                selected=min(2L, ncol(reducedDim(se, redDimType)), na.rm=TRUE))
        })

        # Observer for the tour ----

        observeEvent(input[[.tourInput]], {
            tour <- read.delim(
                system.file("extdata", "intro_shinyLabels.txt", package="Hancock"),
                sep=";", stringsAsFactors=FALSE, row.names=NULL, quote="")
            introjs(session, options=list(steps=tour))
        })

        # Observer for resetting the gene signature object ----

        observeEvent(input[[.resetInput]], {
            REACTIVE$GS <- gs
        })

        # Observer for closing the app and returning the updated object ----

        observeEvent(input[[.doneInput]], {
            stopApp(invisible(REACTIVE$GS))
        })

    }

    shinyApp(ui=app_ui, server=app_server)
}
#nocov end
