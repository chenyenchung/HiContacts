#' @title Contact map insulation
#' @name getDiamondInsulation
#' @rdname getDiamondInsulation
#' @description
#' 
#' Computes diamond insulation score along the entire genome 
#'
#' @param x A `HiCExperiment` object over a full genome
#' @param window_size Which window size to use to compute diamond insulation 
#' score (default: 10 * resolution)
#' @param weak_threshold Less stringent cutoff to call borders in the 
#' diamond insulation score
#' @param strong_threshold More stringent cutoff to call borders in the 
#' diamond insulation score
#' @param BPPARAM BiocParallel parallelization settings 
#' @return a `HiCExperiment` object with additional `insulation` metadata, 
#' containing the diamond insulation score computed 
#'
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples 
#' library(HiContacts)
#' hic <- contacts_yeast() |> 
#'   refocus('II:1-300000') |> 
#'   zoom(1000)
#' diams <- getDiamondInsulation(hic)
#' getDiamondInsulation(diams)

getDiamondInsulation <- function(
    x, 
    window_size = NULL, 
    BPPARAM = BiocParallel::bpparam()
) {

    message( "Going through preflight checklist..." )
    # - Check resolutions 
    si <- seqinfo(x)
    resolution <- resolution(x)
    if(is.null(window_size)) window_size <- 8 * resolution
    ress <- resolutions(x)
    if (!resolution %in% ress) stop("Chosen resolution not available")
    if (!window_size %in% ress) stop("Chosen window_size not available")
    fact <- window_size / resolution
    cm <- gi2cm(interactions(x), 'balanced')
    cm_mat <- as.matrix(cm)

    # - Scan each bin individually
    message( "Scan each window and compute diamond insulation score..." )
    idx_ids <- seq(fact+1, nrow(cm)-fact-1)
    scores <- BiocParallel::bplapply(BPPARAM = BPPARAM, 
        idx_ids, 
        function(idx) {
            # Extract matrix from chr-level, centered on idx
            mat <- cm_mat[seq(idx-fact, idx+fact), seq(idx-fact, idx+fact)]
            # Extract top-right quadrant
            mat_trq <- mat[seq(1, fact) + 1 - 1, seq(fact, fact + fact - 1) + 1 + 1]
            sum(mat_trq,na.rm = TRUE)
        }
    ) |> unlist()
    fbins <- regions(cm)[idx_ids]
    fbins$score <- scores
    fbins$insulation <- log2( scores / mean(scores) )
    fbins <- fbins[!is.na(fbins$insulation)]
    fbins <- fbins[!is.nan(fbins$insulation)]
    fbins <- fbins[is.finite(fbins$insulation)]

    # - Annotate boundaries
    message( "Annotating diamond score prominence for each window..." )
    mins <- which(diff(sign(diff(fbins$insulation)))==+2)+1
    maxs <- which(diff(sign(diff(fbins$insulation)))==-2)+1

    # It is not safe to assume that mins and maxs vectors can only differ by 1
    # in a noisy environment.
    
    # if (length(maxs) < length(mins)) maxs[[length(maxs)+1]] <- length(fbins)
    # if (length(maxs) > length(mins)) maxs <- utils::tail(maxs, -1)
    fbins$min <- seq_along(fbins) %in% mins
    fbins$prominence <- NA
    
    fbins$prominence[fbins$min] <- sapply(
        mins, function(min_pos) {
            # Find closest maximum before this minimum
            left_maxs <- maxs[maxs < min_pos]
            # Find closest maximum after this minimum
            right_maxs <- maxs[maxs > min_pos]

            # Only calculate prominence if we have maxima on both sides
            if (length(left_maxs) > 0 && length(right_maxs) > 0) {
                left_max <- max(left_maxs)
                right_max <- min(right_maxs)
                
                # Calculate prominence as depth relative to lower of the two surrounding peaks
                left_height <- fbins$insulation[left_max]
                right_height <- fbins$insulation[right_max]
                reference_height <- min(left_height, right_height)
                
                # Prominence is how far the minimum dips below the reference height
                out <- -(fbins$insulation[min_pos] - reference_height)
            } else {
                out <- NA
            }
        }
    )

    # Update min indicator to only mark valid minima
    fbins$min[mins[is.na(fbins$prominence[mins])]] <- FALSE
    
    ## -- Return initial HiCExperiment object with insulation in metadata
    metadata(x)$insulation <- fbins
    return(x)

}

#' @rdname getDiamondInsulation
#' @export

getBorders <- function(x, weak_threshold = 0.2, strong_threshold = 0.5) {

    fbins <- metadata(x)$insulation
    fbins$border <- fbins$prominence >= weak_threshold & fbins$min
    fbins$type <- ifelse(fbins$prominence < weak_threshold, NA, ifelse(fbins$prominence < strong_threshold, 'weak', 'strong'))
    gr <- fbins[which(fbins$border)]
    gr$score <- gr$prominence
    names(gr) <- gr$type
    S4Vectors::mcols(gr) <- data.frame(score = gr$score)
    topologicalFeatures(x, 'borders') <- gr
    return(x)

}

.getDiamondScore <- function(cm, bin, fact) {
    ids <- seq(bin$bin_id - fact, bin$bin_id + fact)
    idx_ids <- seq(fact+1, nrow(cm_mat)-fact-1)
    balanced <- l[ids+1, ids+1]
    s <- sum(
        balanced[seq(1, fact) + 1 - 1, seq(fact, fact + fact - 1) + 1 + 1],
        na.rm = TRUE
    )
    return(s)
}
