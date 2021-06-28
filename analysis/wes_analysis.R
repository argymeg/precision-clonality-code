# load libraries ----------------------------------------------------------------------------------------------------
library(pacman)
p_load(tidyverse, reshape2)
p_load(glue, janitor, params)
p_load(limma, broom)
p_load(forcats, effsize)
p_load(MultiAssayExperiment, SummarizedExperiment)
p_load(ComplexHeatmap, cowplot, ggsci, ggpubr, ggstatsplot)
p_load(LBSPR, magrittr,maftools, stargazer)
library(gridExtra);library(grid);library(rhdf5)
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(wranglr)

# read files ----------------------------------------------------------------------------------------------------
df3 = read.csv("mutations.csv")

# VAF plots  ----------------------------------------------------------------------------------------------------
p_load(ggrepel)

df_vaf1 = df3 %>% dcast(patient + gene + variant_classification + key + recurrence_type + aaannotation + shared ~ timepoint, value.var = "tumor_f")
df_vaf1[is.na(df_vaf1)] <- 0

df_vaf1.1 = df_vaf %>% dcast(final_clonality + patient + gene + variant_classification + key + recurrence_type + shared_hotspot_updated~ timepoint, value.var = "tumor_f")
df_vaf1.1[is.na(df_vaf1.1)] <- 0

df_vaf1.2 = df_vaf %>% dcast(final_clonality + patient + gene + variant_classification + key + recurrence_type + shared~ timepoint, value.var = "tumor_f")
df_vaf1.2[is.na(df_vaf1.2)] <- 0

#df_vaf2 = df_vaf1 %>% dplyr::filter(patient == "NKI1335") #NKI9127
col = c("Missense_Mutation" = "purple", "Nonsense_Mutation" = "gold", "Silent" = "tomato", "Nonstop_Mutation" = "darkgreen")
col2 = c("private" = "red", "shared_no_hotspot" = "chartreuse1", "shared_hotspot" = "blue")
col_fig = c("private-primary" = "steelblue", "private-recurrence" =  "powderblue", "shared" = "black")
col_fig = c("private-primary" = "lightsalmon4", "private-recurrence" =  "lightsalmon", "shared" = "black")
#col_fig = c("private-primary" = "darkslategray4", "private-recurrence" =  "darkseagreen2", "shared" = "black")

# VAF BOXPLOT - MISHA -------
df3_nonsync = df3 %>% dplyr::filter(!patient %in% duke_samples, !patient %in% c("9100205", "NKI2010", "DCIS-270")) %>% dplyr::filter(final_clonality == "clonal") %>% mutate(mut_type = paste0(shared,"_",timepoint ))
df_vaf_nonsync = df3_nonsync %>% dcast(final_clonality + patient + gene + variant_classification + key + recurrence_type + shared~ timepoint, value.var = "tumor_f")
df_vaf_nonsync[is.na(df_vaf_nonsync)] <- 0

my_comparisons = list(c("private_primary", "shared_primary"), c("private_recurrence", "shared_recurrence")) #, "recurrence"
p0 = ggboxplot(df3_nonsync, x = "mut_type", y = "tumor_f", add.params = list(size = 0.1, jitter = 0.2), 
               color = "mut_type", palette = c("steelblue","black", "powderblue", "black"),#fill = "patient"
               add = "jitter", facet.by = "recurrence_type") + scale_y_continuous(trans='log2') +
  xlab("Mutation Type") + ylab("log2(VAF)") + ggtitle("VAF for private vs shared mutations") +
  theme(title = element_text(size = 8, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt"))) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test");p0
p1 = ggboxplot(df3_sync, x = "mut_type", y = "tumor_f", add.params = list(size = 0.1, jitter = 0.2), order = c("private_primary", "shared_primary","private_recurrence", "shared_recurrence"),
               color = "mut_type", palette = c("lightsalmon4","black", "lightsalmon", "black"),#fill = "patient"
               add = "jitter", facet.by = "recurrence_type") + scale_y_continuous(trans='log2') +
  xlab("Mutation Type") + ylab("log2(VAF)") + ggtitle("VAF for private vs shared mutations") +
  theme(title = element_text(size = 8, face = "bold"),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt"))) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif");p1
p = cowplot::plot_grid(p0, p1, nrow = 1)
pdf(paste0("vaf_comparions_nonsync.pdf"), width=6, height=8)
print(p0)
dev.off()
pdf(paste0("vaf_comparions_sync.pdf"), width=4, height=8)
print(p1)
dev.off()

for(i in 1:length(unique(df_vaf$patient))){
  
  pt = unique(df_vaf$patient)[i]
  df_vaf2 = df_vaf1 %>% dplyr::filter(patient == unique(df_vaf$patient)[i]) #NKI9127
  df_vaf2.1 = df_vaf1.1 %>% dplyr::filter(patient == unique(df_vaf$patient)[i])
  df_vaf2.1 = df_vaf2.1 %>% dplyr::mutate(mut_status = case_when(primary > 0 & recurrence == 0 ~ "private-primary", 
                                                                 primary == 0 & recurrence > 0 ~ "private-recurrence",
                                                                 primary > 0 & recurrence > 0 ~ "shared"))
  pt_clonality = df3 %>% dplyr::filter(patient == unique(df_vaf$patient)[i]) %>% dplyr::select(final_clonality) %>% unique() %>% pull(final_clonality)
  message("working on- ", pt)
  
  gs1 =   ggplot(df_vaf2.1) +
    geom_point(aes(x = primary, y = recurrence, color = mut_status), shape = 1, size = 1) + #, position = "dodge"
    #scale_fill_manual(values = c("purple", "green", "hotpink")) + guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)))  +
    #scale_fill_calc() + 
    theme_cowplot() + 
    #geom_text_repel(min.segment.length = 0, seed = 42, box.padding = 0.5) + 
    
    geom_text_repel(data = df_vaf2.1 %>% dplyr::filter(gene %in% c("APC", "AKT1", "CHD2","NOS2","ERBB2","NF1","ATM","RUNX1","GATA3","ARID1A","CHEK2","ERBB2","CDH1",
                                                                   "FOXA3","KMT2C", "SLC26A7","TPR","USP15","SMARCA4", "PIK3CA", "FOXA1", "FLT3", "EGFR", "PTEN", "TP53")), aes(x = primary, y = recurrence,label = gene),size = 4,
                    # #         #label = gene, 
                    #         # nudge_x = 0.015, nudge_y = 0.05, 
                    #          check_overlap = T
                    min.segment.length = 0, seed = 42, box.padding = 0.5
    ) +
    scale_color_manual(values = col_fig) + 
    xlab("Primary") + ylab("Recurrence") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs1
  
  gs2 =   ggplot(df_vaf2.1) +
    geom_point(aes(x = primary, y = recurrence, color = mut_status), shape = 1, size = 2) + #, position = "dodge"
    #scale_fill_manual(values = c("purple", "green", "hotpink")) + guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)))  +
    #scale_fill_calc() + 
    theme_cowplot() + 
    # geom_text(data = df_vaf2.1 %>% dplyr::filter(gene %in% c("FOXA1", "FLT3", "EGFR", "PTEN", "TP53")), aes(x = primary, y = recurrence,label = gene),size = 4,
    # #         #label = gene, 
    #         # nudge_x = 0.015, nudge_y = 0.05, 
    #          check_overlap = T
    #  ) +
    scale_color_manual(values = col_fig) + 
    xlab("Primary") + ylab("Recurrence") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs2
  
  gs3 =   ggplot(df_vaf2.1) +
    geom_point(aes(x = primary, y = recurrence, color = mut_status), shape = 16, size = 2, stroke = 0.2) + #, position = "dodge"
    #scale_fill_manual(values = c("purple", "green", "hotpink")) + guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)))  +
    #scale_fill_calc() + 
    theme_cowplot() +
    #theme_void() +
    geom_text_repel(data = df_vaf2.1 %>% dplyr::filter(gene %in% c("APC", "AKT1", "CHD2","NOS2","ERBB2",
                                                                   "FOXA3","KMT2C", "SLC26A7","TPR","USP15","SMARCA4", "PIK3CA", "FOXA1", "FLT3", "EGFR", "PTEN", "TP53")), aes(x = primary, y = recurrence,label = gene),size = 4,
                    # #         #label = gene, 
                    #         # nudge_x = 0.015, nudge_y = 0.05, 
                    #          check_overlap = T
                    min.segment.length = 0, seed = 42, box.padding = 0.5
    ) +
    scale_color_manual(values = col_fig) + 
    xlab("Primary") + ylab("Recurrence") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs3
  
  p = cowplot::plot_grid(gs1, gs2, gs3, nrow = 1)
  pdf(paste0("vaf_",pt,"_", pt_clonality,"orange.pdf"), width=18, height=3.5)
  print(p)
  dev.off()
  
  p1 = gs3
  pdf(paste0("vaf_",pt,"_paper_", pt_clonality,"orange.pdf"), width=4.5, height=2.5)
  print(p1)
  dev.off()
  
}


# VAF Scatter all -------
p_load(ggrepel)

df_vaf1test = df_vaf1
df_vaf1test1 = df_vaf1test %>% dplyr::select(patient,gene) %>% unique() %>% dplyr::count(gene) %>% dplyr::arrange(-n) %>% top_n(50)
col_rec = c("dcis" = "seagreen2", "invasive" = "maroon", "dcis_invasive" = "blue")
genes_selected = intersect(df_vaf1test1$gene, mut_sig_genes_all)

gs0 =   ggplot(ms, aes(x=-log10(p.x), y=-log10(p.y), label = gene))  + #+ annotation_custom(grob) 
  geom_point(shape = 21, size = 2, 
             alpha = 0.7, stroke = 0.5) + #, position = "dodge"
  geom_label_repel(data = filter(ms, -log10(p.y)>3)) + 
 theme_cowplot() + 
  ggtitle("Top mutated genes in recurrence based on mutsigCV") + 
  scale_fill_manual(values = col) + 
  xlab("Primary log10(p)") + ylab("Recurrence log10(p)") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs0

df_vaf11 = df_vaf1 %>% dplyr::filter(gene %in% c("FOXA1", "SMAD4", "NOTCH2", "DNMT3A")) %>% #AKT1", "SMAD4", "DNMT3", "GCM2", "ADH6" 
  dplyr::select(patient, gene, aaannotation, primary, recurrence) %>% arrange(gene)
gs01 =  ggplot(df_vaf11, aes(x=primary, y=recurrence, fill = gene, label = paste0(gene, "-", aaannotation)))  + #+ annotation_custom(grob) 
  geom_point(shape = 21, size = 2, 
             alpha = 0.7, stroke = 0.5, aes(fill = gene)) + #, position = "dodge"
  geom_label_repel(data = filter(df_vaf1, gene %in% c("FOXA1", "SMAD4", "NOTCH2", "DNMT3A"))) + 
  theme_cowplot() + 
  ggtitle("Top mutated genes in recurrence based on mutsigCV") + 
  scale_fill_brewer(palette = "Set3") +
  #scale_fill_manual(values = col_rec) + 
  xlab("Primary") + ylab("Recurrence") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs01

gs1 =   ggplot(df_vaf1)  + #+ annotation_custom(grob) 
  geom_point(aes(x=primary, y=recurrence, fill = variant_classification), shape = 21, size = 2, 
             alpha = 0.7, stroke = 0.5) + #, position = "dodge"
  theme_cowplot() + 
  ggtitle("All mutations by type") + 
  scale_fill_manual(values = col) + 
  xlab("All Primary") + ylab("All Recurrence") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs1

gs2 =   ggplot(df_vaf1.1) +
  geom_point(aes(x = primary, y = recurrence, fill = shared_hotspot_updated), shape = 21, size = 2, 
             alpha = 0.8, stroke = 0.5) + #, position = "dodge"
  theme_cowplot() + 
  ggtitle("All mutations by hotspots") + 
  scale_fill_manual(values = col2) + 
  xlab("All Primary") + ylab("All Recurrence") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs2

gs3 =   ggplot(df_vaf1.1) +
  geom_point(aes(x = primary, y = recurrence, fill = recurrence_type), shape = 21, size = 2, 
             alpha = 0.8, stroke = 0.5) + #, position = "dodge"
  ggtitle("All mutations by recurrence") + 
  theme_cowplot() + 
  scale_fill_manual(values = col_rec) + 
  xlab("All Primary") + ylab("All Recurrence") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs3

gs4 =   df_vaf1 %>% dplyr::filter(recurrence_type == "dcis") %>% 
  ggplot() +
  geom_point(aes(x = primary, y = recurrence, fill = variant_classification), shape = 21, size = 2, 
             alpha = 0.8, stroke = 0.5) + #, position = "dodge"
  theme_cowplot() +
  ggtitle("DCIS only mutations by type") + 
  scale_fill_manual(values = col) + 
  xlab("DCIS Primary") + ylab("DCIS Recurrence") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs4

gs5 =   df_vaf1 %>% dplyr::filter(recurrence_type == "invasive") %>% 
  ggplot() +
  geom_point(aes(x = primary, y = recurrence, fill = variant_classification), shape = 21, size = 2, 
             alpha = 0.8, stroke = 0.5) + #, position = "dodge"
  theme_cowplot() + 
  ggtitle("Invasive Recurrence mutations by type") + 
  scale_fill_manual(values = col) + 
  xlab("DCIS Primary") + ylab("Inv Recurrence") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs5


p = cowplot::plot_grid(gs0, gs01, gs1, gs2, gs3,gs4, gs5, ncol = 7)
pdf(paste0("sync_vaf_all_patients.pdf"), width=42, height=4)
print(p)
dev.off()

# SHARED VAF heatmaps  ----

df_htvaf = df3
df_htvaf1 = df_htvaf %>% dcast(patient + gene + key + shared_hotspot_updated + recurrence_type ~ timepoint, value.var = "tumor_f")
df_htvaf1[is.na(df_htvaf1)] <- 0

df_htvaf1.1 = df_htvaf %>% dcast(patient + gene + variant_classification + key + recurrence_type + shared_hotspot_updated~ timepoint, value.var = "tumor_f")
df_htvaf1.1[is.na(df_htvaf1.1)] <- 0

#df_htvaf2 = df_htvaf1 %>% dplyr::filter(patient == "NKI1335") #NKI9127
col = c("Missense_Mutation" = "purple", "Nonsense_Mutation" = "gold", "Silent" = "tomato", "Nonstop_Mutation" = "darkgreen")
col2 = c("private" = "red", "shared_no_hotspot" = "chartreuse1", "shared_hotspot" = "blue")

for(i in 1:length(unique(df3$patient))){
  
  pt = unique(df3$patient)[i]
  df_htvaf2 = df_htvaf1 %>% dplyr::filter(patient == unique(df3$patient)[i]) #NKI9127
  df_htvaf2.1 = df_htvaf1.1 %>% dplyr::filter(patient == unique(df3$patient)[i])
  pt_clonality = df3 %>% dplyr::filter(patient == unique(df3$patient)[i]) %>% dplyr::select(final_clonality) %>% unique() %>% pull(final_clonality)
  message("working on- ", pt)
  
  gs1 =   ggplot(df_htvaf1)  + #+ annotation_custom(grob) 
    geom_col(aes(x=primary, y=recurrence, fill = variant_classification), shape = 21, size = 2, 
             alpha = 0.7, stroke = 0.5) + #, position = "dodge"
    #scale_fill_manual(values = c("purple", "green", "hotpink")) + guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)))  +
    #scale_fill_calc() + 
    theme_cowplot() + 
    scale_fill_manual(values = col) + 
    xlab("Primary") + ylab("Recurrence") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs1
  
  gs2 =   ggplot(df_htvaf2.1) +
    geom_point(aes(x = primary, y = recurrence, fill = shared_hotspot_updated), shape = 21, size = 2, 
               alpha = 0.8, stroke = 0.5) + #, position = "dodge"
    #scale_fill_manual(values = c("purple", "green", "hotpink")) + guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)))  +
    #scale_fill_calc() + 
    theme_cowplot() + 
    scale_fill_manual(values = col2) + 
    xlab("Primary") + ylab("Recurrence") +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = "black"));gs2
  
  p = cowplot::plot_grid(gs1, gs2)
  pdf(paste0("vaf_",pt,"_", pt_clonality,".pdf"), width=12, height=3.5)
  print(p)
  dev.off()
  
}

# oncoplots -----------------
# get the main genes from mutsigCN and panel seq

# get the main genes -----
genes_list = c("TP53", "PIK3CA", "DDX6", "NOTCH2","PLEC", "ERBB2", "KMT2C", "AKT1","TMEM132C", "ATM", "EGFR", "FLG2", "FOXA1",
               "TENM2", "ABCA13", "CDH23", "BTBD10", "XPO1", "DYNC1H1", "IL8", "DMGDH", "FBXW11", "SLC24A1", "MEF2A", "VHL", "SMAD4", "DDX")

# column annotations
sample_sequence = unique(df4_meta$patient)
coldata = data.frame(Recurrence = df4_meta$recurrence_type,
                     final_clonality = df4_meta$final_clonality,#minor_group_factor
                     #sample_id = df4_meta$sample_id,
                     #timepoint = df4_meta$timepoint,
                     patient_id = df4_meta$patient,
                     stringsAsFactors = F) #%>% dplyr::arrange(sample_sequence) #Recurrence, final_clonality
coldata1 = coldata[match(sample_sequence, coldata$patient_id),]

col_rec = c("dcis" = "seagreen2", "invasive" = "maroon", "dcis_invasive" = "blue")
col_clon = c("independent" = "#800000", "ambiguous" = "#FFB347", "clonal" = "#FDFD96")
col_time = c("primary" = "lightblue", "recurrence" = "lightpink")

# top annotation -----
ha_bottom = dplyr::select(coldata1, final_clonality, Recurrence) %>% 
  HeatmapAnnotation(df = ., show_annotation_name = FALSE,
                    annotation_legend_param = list(final_clonality = list(title = "Final Clonality")),
                    annotation_height = unit(0.55, "mm"), border = TRUE, 
                    #width = unit(0.1, "cm"),
                    simple_anno_size = unit(0.25, "cm"),
                    col = list(Recurrence = col_rec, final_clonality = col_clon)) #, timepoint = col_time
draw(ha_bottom)

# top annotation -----
df4_bot = df3 %>% dplyr::filter(patient %in% patient_seq)
df4_bot1  = df4_bot %>% group_by(patient, timepoint, shared) %>% dplyr::count()
df4_bot2 = df4_bot1 %>% dplyr::filter(patient %in% sample_sequence) %>% dplyr::mutate(num_mut = n, mut_count = if_else(shared == "shared", "shared", paste0("private-", timepoint))) %>% ungroup() %>% 
  dplyr::select(patient, num_mut, mut_count) %>% distinct()
df4_bot2 = data.frame(df4_bot2)
df4_bot3 = reshape2::dcast(df4_bot2, formula = mut_count ~ patient, value.var = "num_mut", fun.aggregate = sum)
mat4_bot3 = data.matrix(df4_bot3)
mat4_bot3 = mat4_bot3[,-1]
rownames(mat4_bot3) = df4_bot3$mut_count
#mat4_bot4 = t(mat4_bot3) 
mat4_bot4 = mat4_bot3[,sample_sequence]
mat4_bot5 = prop.table(t(mat4_bot4), margin = 1)
col_meno = c("private-primary" = "steelblue", "private-recurrence" =  "powderblue", "shared" = "black")
ha_top = HeatmapAnnotation(Mut_freq = anno_barplot(mat4_bot5, border = FALSE,show_annotation_name = FALSE, 
                                                   gp = gpar(fill = col_meno, border = NA, lty = "blank"),
                                                   axis = TRUE, size = unit(0.95, "mm"), bar_width = 0.85))

samps_sorted = sample_sequence

# rearrnge the matrix
#mat = df5_wide[, samps_sorted]
mat = tt_select
mat = as.matrix(tt_select[,-1])
rownames(mat) = tt_select$gene

# right row annotation 1 ----
df4_uniq = df4 %>% distinct(gene, patient, recurrence_type)
mat_right = table(df4_uniq$gene, df4_uniq$recurrence_type)
num_dcis = df4 %>% dplyr::filter(recurrence_type == "dcis") %>% pull(patient) %>% unique() 
num_dcis = length(num_dcis)
num_inv = df4 %>% dplyr::filter(recurrence_type == "invasive") %>% pull(patient) %>% unique() 
num_inv = length(num_inv)

mat_right1 = as.data.frame(mat_right)
mat_right1$Var1 = as.character(mat_right1$Var1)
mat_right11 = mat_right1 %>% dplyr::filter(Var1 %in% gene_seq)
mat_right2 = mat_right11 %>% dplyr::mutate(Proportion = if_else(Var2 == "dcis", round(-Freq*100/num_dcis), round(Freq*100/num_inv)))
mat_right2$Var1 = factor(mat_right2$Var1, levels = rev(gene_seq), ordered = TRUE)
brks <- seq(-100, 100, 50)
lbls = paste0(as.character(c(seq(100, 0, -50), seq(50, 100, 50))), "%")
p_load(ggthemes)
g1 = ggplot(mat_right2, aes(y = Var1, x = Proportion, fill = Var2)) +   # Fill column
  geom_bar(stat = "identity", width = .6, colour = "black", size = 0.06) +   # draw the bars
  scale_x_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  xlab("Frequency") + ylab("") + ggtitle("Freq of top mutated genes in \n DCIS vs Invasive samples \n Non-syncronous samples") +
  #scale_y_reverse() + 
  #theme_cowplot()+
  theme_tufte() +  # Tufte theme from ggfortify
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = col_rec)
#Bar chart showing the frequency of affected driver genes in metastatic breast cancer (dark blue) versus primary breast cancer 
pdf("mirror_plot_dcis_invasive_nonsync_withoutaxis.pdf", width = 4, height = 4)
g1
dev.off()

# right row annotation 2 ----
df42_uniq = df4 %>% distinct(gene, patient, final_clonality)
mat_right3 = table(df42_uniq$gene, df42_uniq$final_clonality)
num_clonal = df4 %>% dplyr::filter(final_clonality == "clonal") %>% pull(patient) %>% unique() 
num_clonal = length(num_clonal)
num_nc = df4 %>% dplyr::filter(final_clonality == "independent") %>% pull(patient) %>% unique() 
num_nc = length(num_nc)

mat_right3 = as.data.frame(mat_right3)
mat_right3$Var1 = as.character(mat_right3$Var1) 
mat_right3 = mat_right3 %>% dplyr::filter(Var1 %in% gene_seq)
mat_right4 = mat_right3 %>% dplyr::mutate(Proportion_clonality = if_else(Var2 == "clonal", round(-Freq*100/num_clonal), round(Freq*100/num_nc)))
mat_right4$Var1 = factor(mat_right4$Var1, levels = rev(gene_seq), ordered = TRUE)
brks <- seq(-100, 100, 50)
lbls = paste0(as.character(c(seq(100, 0, -50), seq(50, 100, 50))), "%")
p_load(ggthemes)
g2 = ggplot(mat_right4, aes(y = Var1, x = Proportion_clonality, fill = Var2)) +   # Fill column
  geom_bar(stat = "identity", width = .6,colour = "black", size = 0.06) +   # draw the bars
  scale_x_continuous(breaks = brks,   # Breaks
                     labels = lbls) + # Labels
  xlab("Frequency") + ylab("") + ggtitle("Freq of top mutated genes in \n Clonal vs Indep & Ambiguous samples \n Non-syncronous samples") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = .5), 
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = col_clon)
pdf("mirror_plot_cloanl_nonclonal_nonsync_withaxis.pdf", width = 5, height = 5)
g2
dev.off()

col_mut = c("Missense_Mutation" = "deeppink1", "Nonsense_Mutation" = "deepskyblue1", "Silent" = "tomato", "Multi_Hit" = "darkgreen")

library(ComplexHeatmap)
alter_fun = list(
  background = alter_graphic("rect", fill = "#e6e6e6"), ##CCCCCC
  #"0" = alter_graphic("rect", fill = "#CCCCCC"),
  Multi_Hit = alter_graphic("rect", width = 0.66, height = 0.66, fill = col_mut["Multi_Hit"]),
  Missense_Mutation = alter_graphic("rect", width = 0.66, height = 0.66, fill = col_mut["Missense_Mutation"]),
  Nonsense_Mutation = alter_graphic("rect", width = 0.66, height = 0.66, fill = col_mut["Nonsense_Mutation"]), 
  Silent = alter_graphic("rect", width = 0.66, fill = col_mut["Silent"]))
mat[is.na(mat)] <- ""
mat[mat == 0] <- ""

col = c("Missense_Mutation" = "purple", "Nonsense_Mutation" = "gold", "Silent" = "tomato", "Nonstop_Mutation" = "darkgreen")
col2 = c("private" = "red", "shared_no_hotspot" = "chartreuse1", "shared_hotspot" = "blue")


p1 = oncoPrint(mat, col = col_mut,
               column_order = sample_sequence,
               right_annotation = NULL,
               show_column_names = TRUE,
               heatmap_legend_param = list(
                 title = "Alterations",
                 legend_direction = "horizontal", 
                 legend_width = unit(2, "cm")
               ),
               top_annotation = c(ha_top,ha_bottom),
               #bottom_annotation = ha_top2,
               alter_fun = alter_fun, remove_empty_rows = FALSE,
               row_names_gp = gpar(fontsize = 8),pct_gp = gpar(fontsize = 5),
               height = unit(7, "cm"))
pdf("htmaps_non_syncronous_patient.pdf", width = 8, height = 12)
ht = draw(p1, heatmap_legend_side = "bottom")
dev.off()
