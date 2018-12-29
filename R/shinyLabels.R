
# Constants ----

.geneSetNameInput <- "geneSetName"

# App ----

#nocov start
#' Interactively Inspect and Name Gene Signatures
#'
#' This function launches a Shiny app to inspect the gene signatures defined in a \code{tbl_geneset}
#' for the purpose of (re-)naming those signatures interactively.
#' The app returns the updated \code{tbl_geneset} when closed using the \code{"Done"} button.
#'
#' @param gs A set of gene signatures inheriting from \code{\link{tbl_geneset}}.
#'
#' @return The update set of gene signatures as a \code{\link{tbl_geneset}}.
#' @export
#' @importFrom methods is
#' @importFrom shiny shinyApp reactiveValues observeEvent stopApp
#' fluidRow column icon textInput actionButton renderUI HTML uiOutput
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
#' # Example usage ----
#' if (interactive()){
#' x <- runApp(shinyLabels(tgs))
#' }
shinyLabels <- function(gs) {

    stopifnot(is(gs, "tbl_geneset"))

    REACTIVE <- reactiveValues(
        geneset=gs
    )

    app_ui <- dashboardPage(
        dashboardHeader(),
        dashboardSidebar(
            actionButton(inputId="Done", label="Done", icon=icon("sign-out"), width="50%")
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

        output$mainPanels <- renderUI({

            panelList <- list()
            for (id in seq_len(NSETS)) {
                id0 <- id
                geneSetName0 <- levels(gs$set)[id0]
                geneIds0 <- gs[gs$set == geneSetName0, "gene", drop=TRUE]
                geneIdText <- paste(geneIds0, collapse=", ")
                panelList[[id0]] <- box(
                    textInput(
                        inputId=paste0(.geneSetNameInput, id0),
                        label=paste("Gene set name", id0),
                        value=geneSetName0, width="50%",
                        placeholder=paste("Gene set name", id0)
                    ),
                    HTML(sprintf("<p>%s<p>", geneIdText)),
                    width=12, title=paste("Gene set name", id0)
                )
            }
            outValue <- do.call(fluidRow, panelList)

            return(outValue)
        })

        # Observer for the gene set names ----

        for (id in seq_len(NSETS)) {
            local({
                id0 <- id
                inputId0 <- paste0(.geneSetNameInput, id0)
                observeEvent(input[[inputId0]], {
                    newValue <- input[[inputId0]]
                    levels(REACTIVE$geneset$set)[id0] <- as.character(newValue)
                })
            })
        }

        # Observer for closing the app and returning the updated object ----

        observeEvent(input$Done, {
            stopApp(invisible(REACTIVE$geneset))
        })

    }

    shinyApp(ui=app_ui, server=app_server)
}
#nocov end
