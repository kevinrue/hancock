
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
#' @importFrom ggplot2 ggplot aes aes_string geom_bar guides labs
#' scale_fill_manual scale_x_discrete
#' @importFrom cowplot theme_cowplot
barplotPredictionCount <- function(se, highlight=character(0)) {
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
#' @importFrom cowplot theme_cowplot
#' @importFrom scales percent
barplotPredictionProportion <- function(se, highlight=character(0)) {
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
    gg
}
