#' @title plotSegs()
#' @description Plot copy number profile 
#' @param copyNumberSegmented A QDNAseq object with segmented and copy number data
#' @param id A character vector with sample ID
#' @param color A character vector with bin color 
#' @param limits A numeric vector with y axis limits
#' @details
#' @examples plotSegs(copyNumberSegmented)
#' @export
plotSegs <- function(copyNumberSegmented, id, color = "#8a4e97", limits = c(-1.5, 2)){
    #build template
    segs <- log2(Biobase::assayDataElement(copyNumberSegmented, "segmented")) 
    copynumber <- log2(Biobase::assayDataElement(copyNumberSegmented, "copynumber"))
    
    template <- data.frame(region = rownames(segs), copynumber = as.numeric(copynumber[,id]), segs = as.numeric(segs[,id]))
    template$bin <- rep(1:nrow(template))
    template$chr <- factor(gsub("\\:.*", "", template$region))
    template <- template[template$chr != "Y",]
    template$start <- as.integer(gsub("\\-.*","",gsub(".*:","",template$region))) 
    template$end <- as.integer(gsub(".*-","",template$region))
    
    #define centromere limits
    rlechr <- rle(as.vector(template$chr))
    binchrend <- c()
    currentbin <- 0
    binchrmdl <- c()
    for (i in seq_along(rlechr$values)) {
        currentmiddle <- currentbin+rlechr$lengths[i]/2
        currentbin <- currentbin+rlechr$lengths[i]
        binchrend <- append(binchrend, currentbin) 
        binchrmdl <- append(binchrmdl, currentmiddle) 
    }
    
    cytobands <- data.table::fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz") #download cytoband details from UCSC
    cytobands[[1]] <- gsub("chr","",cytobands[[1]])
    cytobands <- cytobands[cytobands[[1]] %in% c(1:22, "X"),]
    cytobands[[4]] <- gsub("[[:punct:]]", "", gsub('[0-9]+', "",as.character(cytobands[[4]])))
    positions <- c()
    for ( i in rlechr$values ){
        tmp <- cytobands[cytobands[[1]] == i,]
        tmp.p <- tmp[tmp[[4]] == "p",]
        p = max(tmp.p[,3])
        positions <- c(positions, p)
    }
    
    loo <- template 
    loo$chr <- as.character(loo$chr)
    
    bincytoend <- c()
    for (i in 1:length(positions)){
        bincytoend <- c(bincytoend, which(loo$chr == i & loo$start <= positions[i] & loo$end >= positions[i])) 
    }
    
    #define y axis limits
    bottom <- limits[1]
    cap <- limits[2]
    
    #plot
    library(ggpubr)
    library(ggrastr)
    ggplot2::ggplot() +
        scale_y_continuous(name = "Log2 Ratio", limits = c(bottom,cap), breaks = seq(bottom, cap, 0.5), expand=c(0.01,0.01)) +
        scale_x_continuous(name = "Chromosome", limits = c(0,tail(binchrend,1)), breaks = binchrmdl, labels = rlechr$values, expand = c(0,0)) +
        
        geom_hline(yintercept = 0, color = '#666666', size = 0.3) +
        
        geom_hline(yintercept = 0.5, color = '#666666', linetype = "dotted",size=0.2) +
        geom_hline(yintercept = 1, color = '#666666', linetype = "dotted",size=0.2) +
        geom_hline(yintercept = 1.5, color = '#666666', linetype = "dotted",size=0.2) +
        geom_hline(yintercept = 2, color = '#666666', linetype = "dotted",size=0.2) +
        geom_hline(yintercept = -0.5, color = '#666666', linetype = "dotted",size=0.2) +
        geom_hline(yintercept = -1, color = '#666666', linetype = "dotted",size=0.2 ) +
        geom_hline(yintercept = -1.5, color = '#666666', linetype = "dotted",size=0.2) +
        
        geom_vline(xintercept = binchrend, color = "#666666", linetype = "solid",size=0.2) +
        
        geom_vline(xintercept = bincytoend, color = "#666666", linetype = "dashed" ,size=0.2 ) +
        
        geom_point_rast(aes(x=bin, y=copynumber), data = na.omit(template), color = color, shape = ".", raster.dpi=300) +
        geom_point_rast(aes(x=bin, y=segs), data = na.omit(template), color = "black", size=0.1, raster.dpi=300) +
        
        theme_classic() + 
        theme(
            axis.line = element_line(color='black'), 
            axis.ticks = element_line(color='black'), 
            axis.text.y = element_text(color='black', size = 7),
            axis.text.x = element_text(color="black", size = 7, angle = 90),                
            axis.title.y = element_text(color = "black", size = 7),               
            axis.title.x = element_text(color = "black", size = 7),
            panel.border = element_rect(color = "black", fill=NA, size = 0.5)
        ) 
}
