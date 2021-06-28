#' @title HeatmapCNPairs()
#' @description Plot absolute copy number heatmap of a pair of samples from copy number call scores data.
#' @param mat A matrix with the call scores per bin across the genome (columns) and the pair of samples (rows). 
#' @param breakdir A character vector with the directory to the shared breakpoints exported with exportSharedBreaks() function.
#' @param brk A vector of break values to plot. By default, call score < 0 will be losses (blue) and > 0 will be gains (red).
#' @param col A character vector of colors for each sample to plot. By default, purple and orange.
#' @examples HeatmapCNPairs(mat, breakdir = ".")
#' @export
HeatmapCNPairs <- function(mat, breakdir, brk = c(-1,0,1), col = c("#8a4e97", "#ff762f")){
    if ( nrow(mat) != 2 ){
        stop("mat row size must be 2")
    }
    #read the breakpoints
    breaks <- dir(breakdir, "*.tsv", full.names = TRUE)
    pairs <- gsub(".tsv", "", dir(breakdir, "*.tsv", full.names = FALSE))
    
    breaks <- lapply(breaks, data.table::fread)
    names(breaks) <- pairs
    
    shared <-  breaks[names(breaks) == paste0(rownames(mat)[1], "_", rownames(mat)[2])][[1]]
    shared$seqnames <- gsub("chr", "", shared$seqnames)
    
    #build shared breakpoints and copy number annotation bar for the heatmap
    annot <- rep(NA, ncol(mat))
    m <- as.data.frame(t(mat))
    m.chr <- gsub("\\:.*","",rownames(m))
    m.start <- as.numeric(gsub('.*\\:',"",gsub("\\-.*","",rownames(m))))
    m.end <- as.numeric(gsub('.*\\-',"",rownames(m)))
    
    start.tab <- na.omit(shared[1:which(shared[,1] == "seqnames")-1,])
    start.tab.new <- start.tab
    if ( nrow(start.tab) != 0 ){
        for (i in 1:nrow(start.tab)){
            message("Line ", i, " out of ", nrow(start.tab))
            if ( start.tab$start[i] == 1 ){
                start.tab.new <- start.tab.new[-i,]
            }
        }
    }
    start.tab <- start.tab.new
    
    end.tab <- na.omit(shared[which(shared[,1] == "seqnames")+1:nrow(shared),])
    end.tab.new <- end.tab
    if ( nrow(end.tab) != 0 ){
        for (i in 1:nrow(end.tab)){
            mx <- max(m.end[which(m.chr == end.tab$seqnames[i])]) 
            if ( end.tab$end[i] == mx ){
                end.tab.new <- end.tab.new[-i,]
            }
        }
    }
    end.tab <- end.tab.new
    
    samples.ids <- colnames(m)
    
    mix <- rbind(start.tab, end.tab)
    
    shared.regions <- duplicated(mix)
    dup <- mix[shared.regions,]
    colnames(dup) <- c("Chr", "Start", "End", "width", "strand")
    
    for(i in 1:nrow(dup)){
        index <- which(m.chr == dup$Chr[i] & m.start >= as.numeric(dup$Start[i]) & m.end <= as.numeric(dup$End[i]))
        annot[index] <- "Shared region"
    }
    
    annot_brk <- rep(NA, ncol(mat))
    
    if ( nrow(start.tab) != 0 ){
        for (i in 1:nrow(start.tab)){
            a=which(colnames(mat) == paste0(start.tab$seqnames[i], ":", start.tab$start[i], "-", as.integer(as.integer(start.tab$start[i])+99999) ) )
            annot_brk[rep(a:(70+a))] <- "Shared breakpoint"
        }
    }
    
    if ( nrow(end.tab) != 0 ){
        for (i in 1:nrow(end.tab)){
            a=which(colnames(mat) == paste0(end.tab$seqnames[i], ":", as.integer(as.integer(end.tab$end[i])-99999), "-", end.tab$end[i] ) )
            annot_brk[rep((a-70):a)] <- "Shared breakpoint"
        }
    }
    
    require(ComplexHeatmap)
    ha = HeatmapAnnotation(Shared = annot, Shared_break = annot_brk,
                           col = list(Shared = c("Shared region" = "antiquewhite4"), Shared_break = c("Shared breakpoint" = "black")), 
                           na_col = "white", 
                           height = unit(0.4, "cm"), 
                           annotation_label = c(" "),
                           simple_anno_size_adjust = TRUE, 
                           show_legend = FALSE,
                           show_annotation_name = FALSE
    )
    
    #create chromosome info
    chromosome_info <- data.frame(chrom = gsub("\\:.*", "", colnames(mat)))
    
    #plot heatmap
    mat <- mat[c(1,1,2,2),]
    
    ht <- Heatmap(
        mat, 
        border = T,
        show_parent_dend_line = F,
        cluster_columns=F, 
        column_order = 1:ncol(mat),
        cluster_row_slices = F,
        row_order = 1:nrow(mat),
        row_split = c("A", "A", "B", "B"),
        column_title_side = "top",
        show_row_names = F,
        show_column_names = F,
        column_title_gp = gpar(fontsize = 7),
        column_title_rot = 90,
        column_split = factor(as.vector(unlist(chromosome_info)), levels = c(1:22, "X")),
        col = circlize::colorRamp2(breaks = brk, c("#91A5E1", "white", "#EF8C80")),
        bottom_annotation = ha,
        right_annotation = rowAnnotation(Case = c("s1", "s1", "s2", "s2"), 
                                         col = list(Case = c("s1" = col[1], "s2" = col[2])), 
                                         width = unit(bar.size-0.1, "cm"), 
                                         show_annotation_name = F,
                                         show_legend = F,
                                         simple_anno_size_adjust = T),
        use_raster = F,
        show_heatmap_legend = F,
        show_row_dend = F,
        row_title = " "
    )
    
    return(ht)
}