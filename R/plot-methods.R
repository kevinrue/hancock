
#' @describeIn predictByProportionPositive
#' Returns a \code{Heatmap} displaying the proportion (on a scale from 0 to 100) of samples that are positive for each individual signature in each cluster.
#'
#' @aliases plotProportionPositive
#'
#' @param ... Additional arguments to be passed to methods.
#'
#' @export
#' @importFrom ComplexHeatmap Heatmap
plotProportionPositive <- function(
    se, ...
){
    ppbc <- metadata(se)[["Hancock"]][["ProportionPositiveByCluster"]]
    if (is.null(ppbc)) {
        stop("Method 'ProportionPositive' was not run yet.")
    }
    Heatmap(matrix=t(ppbc*100), name="Proportion (%)", ...)
}

#' @describeIn predictHancock Returns a \code{ggplot} bar plot displaying
#' the count of samples predicted for each gene signature.
#'
#' @param highlight Character vector indicating names of signatures to highlight.
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom BiocGenerics ncol
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes aes_string geom_bar guides labs rel
#' scale_fill_manual scale_x_discrete
#' theme element_blank
#' @importFrom cowplot theme_cowplot
barplotPredictionCount <- function(
    se, highlight=character(0), labels=TRUE
) {
    ggFrame <- as.data.frame(colData(se)[, "Hancock"], row.names=seq_len(ncol(se)))
    ggFrame$highlight <- FALSE
    if (length(highlight) > 0) {
        ggFrame[which(ggFrame$prediction %in% highlight), "highlight"] <- TRUE
    }
    gg <- ggplot(ggFrame, aes_string("prediction")) +
        geom_bar(aes_string(fill="highlight")) +
        scale_fill_manual(values=c("TRUE"="black", "FALSE"="grey")) +
        scale_x_discrete(drop=FALSE) +
        guides(fill="none") +
        labs(y="Count", x="Prediction") +
        theme_cowplot()
    if (!labels) {
        gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    }
    gg
}

#' @describeIn predictHancock Returns a \code{ggplot} bar plot displaying
#' the proportion of samples predicted for each gene signature.
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom BiocGenerics ncol
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes_string geom_col guides position_fill labs
#' scale_fill_manual scale_x_discrete scale_y_continuous
#' theme element_blank
#' @importFrom cowplot theme_cowplot
#' @importFrom scales percent
barplotPredictionProportion <- function(
    se, highlight=character(0), labels=TRUE
) {
    # TODO: refactor with barplotPredictionCount above
    ggFrame <- as.data.frame(table(colData(se)$Hancock$prediction))
    ggFrame$Proportion <- ggFrame$Freq / sum(ggFrame$Freq)
    ggFrame$highlight <- FALSE
    if (length(highlight) > 0) {
        ggFrame[which(ggFrame$Var1 %in% highlight), "highlight"] <- TRUE
    }
    gg <- ggplot(ggFrame, aes_string("Var1", "Proportion")) +
        geom_col(aes_string(fill="highlight")) +
        scale_fill_manual(values=c("TRUE"="black", "FALSE"="grey")) +
        scale_x_discrete(drop=FALSE) +
        scale_y_continuous(labels=scales::percent) +
        guides(fill="none") +
        labs(y="Proportion", x="Prediction") +
        theme_cowplot()
    if (!labels) {
        gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    }
    gg
}

#' @describeIn predictHancock Returns a \code{ggplot} bar plot displaying
#' the first reduced dimension result in \code{reducedDims(se)}.
#'
#' @param redDimType Name of the reduced dimension result type to display.
#' @param x Name of the covariate to display on the x-axis.
#' @param y Name of the covariate to display on the y-axis.
#' @param labels Logical value indicating whether to display labels.
#'
#' @return A \code{ggplot} object.
#' @export
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom ggplot2 ggplot aes_string geom_point geom_label geom_text
#' guides position_fill labs scale_color_manual
#' @importFrom cowplot theme_cowplot
reducedDimPrediction <- function(
    se, highlight=character(0), redDimType="PCA", x=1, y=2, labels=TRUE
) {
    # TODO: refactor with barplotPredictionCount above
    ggFrame <- as.data.frame(reducedDim(se, redDimType))
    ggFrame <- ggFrame[, c(x, y)]
    colnames(ggFrame) <- c("X", "Y")
    ggFrame$prediction <- se$Hancock$prediction
    ggFrame$highlight <- FALSE
    if (length(highlight) > 0) {
        ggFrame[which(ggFrame$prediction %in% highlight), "highlight"] <- TRUE
    }
    gg <- ggplot() +
        geom_point(aes_string("X", "Y"), subset(ggFrame, !highlight), color="grey") +
        geom_point(aes_string("X", "Y"), subset(ggFrame, highlight), color="red")
    if (labels) {
        ggLabels <- data.frame(
            prediction=levels(ggFrame$prediction),
            X=tapply(ggFrame[, "X"], ggFrame$prediction, FUN="mean"),
            Y=tapply(ggFrame[, "Y"], ggFrame$prediction, FUN="mean"),
            highlight=tapply(ggFrame$highlight, ggFrame$prediction, FUN="unique")
        )
        gg <- gg +
            geom_text(aes_string("X", "Y", label="prediction"), subset(ggLabels, !highlight), color="black", alpha=0.9, size=rel(4)) +
            geom_label(aes_string("X", "Y", label="prediction"), subset(ggLabels, highlight), color="black", alpha=0.8, size=rel(5))
    }
    gg <- gg +
        labs(y="Dimension 1", x="Dimension 2") +
        theme_cowplot()
    gg
}

#' Plotting wrapper
#'
#' Plotting wrapper that dispatches relevant arguments to individual functions.
#'
#' @rdname INTERNAL_plotWrapper
#'
#' @param se An object of class inheriting from "\code{\link{SummarizedExperiment}}".
#' @param highlight Character vector indicating names of signatures to highlight.
#' @param plotType Name of a plot type. See Details.
#' @param labels Logical value indicating whether to display prediction labels.
#' @param ... Additional argument passed to individual plotting functions.
#'
#' @details
#' The \code{plotType} argument should be the name of a plotting function:
#' one of \code{"barplotPredictionCount"}, \code{"barplotPredictionProportion"}, \code{"reducedDimPrediction"}.
#'
#' @return A \code{ggplot} object.
#' @author Kevin Rue-Albrecht
.plotWrapper <- function(se, highlight, plotType, labels=TRUE, ...) {
    dots <- list(...)
    extras <- ""
    if (identical(plotType, "reducedDimPrediction")) {
        redDimType <- dots[["redDimType"]]
        xAxis <- dots[["x"]]
        yAxis <- dots[["y"]]
        extras <- paste0(extras, sprintf(
            ", redDimType='%s', x=%i, y=%i",
            redDimType, xAxis, yAxis))
    }
    # Write explicit prediction(s) to highlight
    highlight <- deparse(highlight)
    # Assemble function call
    functionCall <- sprintf("%s(se, highlight=%s, labels=%s%s)", plotType, highlight, labels, extras)
    ggPlot <- eval(parse(text=functionCall))
    ggPlot
}
