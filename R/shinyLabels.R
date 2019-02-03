#nocov start
#' Interactively Inspect and Name Gene Signatures
#'
#' This function launches a Shiny app to inspect the gene signatures defined in a `BaseSets`
#' for the purpose of (re-)naming those signatures interactively.
#' The app returns the updated `BaseSets` when closed using the "Done" button.
#'
#' @param gs A set of gene signatures inheriting from [`BaseSets`].
#' @param se An object of class inheriting from [`SummarizedExperiment`][RangedSummarizedExperiment-class].
#'
#' @return The updated set of gene signatures as a [`BaseSets`].
#' @export
#' @importFrom methods is as
#' @importFrom shiny shinyApp reactiveValues observeEvent stopApp isolate
#' fluidRow column icon textInput actionButton renderUI renderPlot
#' HTML uiOutput plotOutput checkboxInput selectizeInput conditionalPanel
#' updateSelectizeInput updateTextInput
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#' dashboardBody box menuItem dropdownMenu notificationItem
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDimNames
#' @importFrom rintrojs introjsUI introjs
#' @importFrom utils read.delim
#'
#' @examples
#' # Example data ----
#'
#' bs <- BaseSets(
#'     relations=DataFrame(
#'         element = c("Gene001", "Gene002", "Gene003", "Gene004"),
#'         set     = c(rep("Cell type 1", 2), rep("Cell type 2", 2))
#'     )
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
#' se1 <- predict(bs, se, method="ProportionPositive", cluster.col="cluster")
#'
#' # Add reduced dimension results to enable app features
#' library(SingleCellExperiment)
#' library(Rtsne)
#' sce1 <- as(se1, "SingleCellExperiment")
#' reducedDim(sce1, "PCA") <- prcomp(t(assay(sce1)))$x
#' reducedDim(sce1, "TSNE") <- Rtsne(X=reducedDim(sce1, "PCA"))$Y
#'
#' # Example usage ----
#'
#' if (interactive()){
#'     library(shiny)
#'     x <- runApp(shinyLabels(bs, sce1))
#' }
shinyLabels <- function(gs, se) {

    stopifnot(is(gs, "BaseSets"))
    stopifnot(all(levels(colData(se)[[getPackageName()]][["prediction"]]) %in% setIds(gs)))

    se <- as(se, "SingleCellExperiment")

    # Initial state for reduced dimension inputs ----
    xAxisMax <- 0L
    yAxisMax <- 0L
    if (!is.null(reducedDim(se))) {
        .shinyLabelsPlotChoices <- c(
            .shinyLabelsPlotChoices,
            "Reduced dimension"="reducedDimPrediction")
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
            actionButton(inputId=.collapseAllInput, label="Collapse all", icon=icon("window-minimize"), width="50%"),
            actionButton(inputId=.expandAllInput, label="Expand all", icon=icon("window-maximize"), width="50%"),
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
        title="hancock: Label signatures"
    )

    app_server <- function(input, output, session) {

        # Initialization of non-reactive values ----

        # Storage for persistent non-reactive objects.
        pObjects <- new.env()

        # TODO: "expand all"/"collapse all" button
        pObjects[[.boxOpen]] <- c(TRUE, rep(FALSE, nSets(gs)-1))

        # Initialization of reactive values ----

        # Storage for all the reactive objects
        REACTIVE <- reactiveValues(
            GS=gs,
            SE=se,
            refresh=0L
        )

        # Main panels ----

        panelList <- list()
        for (id in seq_len(nSets(gs))) {
            local({
                id0 <- id
                plotName0 <- paste0(.plotOutput, id0)
                output[[plotName0]] <- renderPlot({
                    geneSetName0 <- setIds(REACTIVE$GS)[id0]
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
            force(REACTIVE$refresh)

            panelList <- list()
            setList <- as(gs, "list")
            for (id in seq_len(nSets(gs))) {
                id0 <- id
                geneSetName0 <- setIds(gs)[id0]
                geneIds0 <- ids(setList[[geneSetName0]])
                geneIdText <- head(paste(geneIds0, collapse=", "), 50) # TODO: change 50
                plotName0 <- paste0(.plotOutput, id0)
                signaturePlot <- plotOutput(plotName0, height="400px")
                panelList[[id0]] <- iSEE:::collapseBox(
                    id=paste0(.boxOpen, id0),
                    title=paste("Gene set", id0),
                    open=pObjects[[.boxOpen]][id0],
                    textInput(
                        inputId=paste0(.geneSetNameInput, id0),
                        label=paste("Gene set name", id0),
                        value=geneSetName0, width="100%",
                        placeholder=paste("Gene set name", id0)
                    ),
                    column(signaturePlot, width=12),
                    HTML(sprintf("<strong>Features:</strong><p>%s<p>", geneIdText))
                )
            }
            outValue <- do.call(fluidRow, panelList)

            outValue
        })

        # Observers for the gene set names ----

        for (id in seq_len(nSets(gs))) {
            local({
                id0 <- id
                inputId0 <- paste0(.geneSetNameInput, id0)
                observeEvent(input[[inputId0]], {
                    newValue <- input[[inputId0]]
                    setIds(REACTIVE$GS)[id0] <- as.character(newValue)
                    levels(colData(REACTIVE$SE)[[getPackageName()]][["prediction"]])[id0] <- as.character(newValue)
                })
            })
        }

        # Observers for the collapsible panels ----

        for (id in seq_len(nSets(gs))) {
            local({
                id0 <- id
                open_field <- paste0(.boxOpen, id0)
                observeEvent(input[[open_field]], {
                    pObjects[[.boxOpen]][id0] <- input[[open_field]]
                })
            })
        }

        observeEvent(input[[.collapseAllInput]], {
            pObjects[[.boxOpen]] <- rep(FALSE, nSets(gs))
            REACTIVE$refresh <- REACTIVE$refresh + 1
        })

        observeEvent(input[[.expandAllInput]], {
            pObjects[[.boxOpen]] <- rep(TRUE, nSets(gs))
            REACTIVE$refresh <- REACTIVE$refresh + 1
        })

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
                system.file("extdata", "intro_shinyLabels.txt", package=getPackageName()),
                sep=";", stringsAsFactors=FALSE, row.names=NULL, quote="")
            introjs(session, options=list(steps=tour))
        })

        # Observer for resetting the gene signature object ----

        observeEvent(input[[.resetInput]], {
            REACTIVE$GS <- gs
            REACTIVE$SE <- se
            for (id in seq_len(nSets(gs))) {
                local({
                    id0 <- id
                    inputId0 <- paste0(.geneSetNameInput, id0)
                    updateTextInput(session, inputId0, value=setIds(gs)[id0])
                })
            }
        })

        # Observer for closing the app and returning the updated object ----

        observeEvent(input[[.doneInput]], {
            stopApp(invisible(REACTIVE$GS))
        })

    }

    shinyApp(ui=app_ui, server=app_server)
}
#nocov end
