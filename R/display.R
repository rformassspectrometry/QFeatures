##' @title Interactive MultiAssayExperiment Explorer
##'
##' @description
##'
##' A shiny app to browser and explore the assays in an
##' `MultiAssayExperiment` object. Each assay can be selected from the
##' dropdown meny in the side panel, and the quantitative data and row
##' metadata are displayed in the respective *Assay* and *Row data*
##' tabs. The *Heatmap* tab displays a heatmap of the assay. The
##' selection of rows in the *Row data* table is used to subset the
##' features displayed in the *Assay* table and the heatmap to those
##' currectly selected. See [QFeatures] for an example.
##'
##' @param object An instance inheriting from `MultiAssayExperiment`.
##'
##' @param n A `numeric(1)` indicating the maximum number of features
##'     (rows) to consider before disabling row clustering and
##'     displaying feature names for speed purposes. Default is 100.
##'
##' @param ...  Additional parameters (other than `Rowv` and `labRow`,
##'     which are set internally based on the value of `n`) passed to
##'     heatmap.
##'
##' @return Used for its side effect.
##'
##' @md
##'
##' @author Laurent Gatto
##'
##' @importFrom stats heatmap
##'
##' @export
display <- function(object, n = 100, ...) {
    stopifnot(inherits(object, "MultiAssayExperiment"))

    requireNamespace("DT")
    requireNamespace("shiny")
    requireNamespace("shinydashboard")

    ui <- shinydashboard::dashboardPage(
        shinydashboard::dashboardHeader(title = "MultiAssayExperiment Viewer"),
        shinydashboard::dashboardSidebar(
            shiny::actionButton("clear_rows", "Clear row selection"),
            shiny::selectInput("name", "Assay:", names(object))
        ),
        shinydashboard::dashboardBody(
            shinydashboard::tabBox(
                width = 12,
                shiny::tabPanel(title = "Heatmap",
                                shiny::plotOutput("heatmap")),
                shiny::tabPanel(title = "Assay",
                                DT::dataTableOutput("assay")),
                shiny::tabPanel(title = "Row data",
                                DT::dataTableOutput("rowdata"))
            )
        )
    )

    server <- function(input, output) {

        shiny::observeEvent(input$clear_rows, {
            DT::selectRows(proxy_rowdata, NULL)
        })

        output$heatmap <- shiny::renderPlot({
            .Rowv <- NULL
            .labRow <- NULL
            .assay <- assay(object[[input$name]])
            sel <- input$rowdata_rows_selected
            if (!length(sel)) sel <- TRUE
            if ((is.logical(sel) & sel) || (length(sel) > n)) {
                .Rowv <- NA
                .labRow <- NA
            }
            heatmap(.assay[sel, , drop = FALSE],
                    Rowv = .Rowv, labRow = .labRow, ...)
        })

        output$assay <- DT::renderDataTable({
            .assay <- assay(object[[input$name]])
            sel <- input$rowdata_rows_selected
            if (!length(sel)) sel <- TRUE
            data.frame(.assay[sel, , drop = FALSE])
        })

        output$rowdata <- DT::renderDataTable({
            .rowdata <- rowData(object[[input$name]])
            data.frame(.rowdata)
        }, selection = list(target = 'row'))

        proxy_rowdata <- DT::dataTableProxy('rowdata')
    }
    shiny::shinyApp(ui, server)
}
