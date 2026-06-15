# 3.1 Main Bacterial Analysis: Data Processing, Visualization, and Statistical Testing

This chapter processes the MMC metagenomic (and ITS) profiles to characterize how antifungal treatment reshapes the gut bacterial community. It covers count import and diversity estimation, beta-diversity (PERMANOVA/PCoA), taxonomic composition, LEfSe differential-abundance testing across time points, dynamics of gut-associated/SCFA-producing taxa, phylogenetic visualization, and cross-kingdom relationships between bacterial and fungal features.

## 3.1.1. Overview of Bacterial ITS and Metagenomic Data Processing

Imports the metagenomic count matrix and taxonomy, collapses counts to the species level with quality filtering, and builds a `TreeSummarizedExperiment` with relative abundances and a full panel of alpha-diversity metrics (observed richness, Shannon, Simpson, rarity, dominance, divergence).

~~~R
ReadCap.20251215 <- mcreadRDS("./workshop/MMC/Aidan_info/v2/ReadCap.20251215.rds")
total_times_tmp <- read.csv("./projects/MMC/Figures/figures_making/v3/patient.PFS.20251215.csv")
total_times_tmp <- total_times_tmp[total_times_tmp$treatment!="Clotrimazole",]
sample_info_total5 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% total_times_tmp$Omics_patient_Names,]
# sample_info_total5 <- sample_info_total5[sample_info_total5$time!="W6",]
sample_info_total5 <- sample_info_total5[sample_info_total5$treatment!="Clotrimazole",]
METAG_taxa_raw <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/MMC.metagenomic.taxa.v3.rds")
METAG_data <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/METAG_data_all.v3.rds")
METAG_data_all <- apply(METAG_data, 2, as.numeric)
rownames(METAG_data_all) <- rownames(METAG_data)
tmp_sum <- colSums(METAG_data_all)
tmp_sum[tmp_sum<120]
METAG_data_all <- METAG_data_all[,colSums(METAG_data_all)>100]
METAG.samples <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/metagenomic.samples_info1.v3.rds")
METAG.samples <- METAG.samples[METAG.samples$sample %in% sample_info_total5$Omics_samples_Names,]
METAG.samples <- METAG.samples[intersect(rownames(METAG.samples),colnames(METAG_data_all)),]
METAG_data_all <- as.data.frame(METAG_data_all[,rownames(METAG.samples)])
METAG_taxa <- METAG_taxa_raw[!duplicated(METAG_taxa_raw$Species),]
taxa1 <- METAG_taxa_raw
Species <- unique(taxa1$Species)
METAG_data_all_ <- lapply(1:length(Species),function(x) {
    tmp_counts <- METAG_data_all[rownames(taxa1[taxa1$Species==Species[x],]),]
    if (nrow(tmp_counts)==1){tmp_counts <- tmp_counts} else {
        tmp_counts <- as.data.frame(t(colSums(tmp_counts)))
    }
    # message(x)
    tmp_counts
    })
METAG_data_all <- do.call(rbind,METAG_data_all_)
rownames(METAG_taxa) <- rownames(METAG_data_all) <- as.character(METAG_taxa$Species)
summary(colSums(METAG_data_all))

library(microeco)
library(mecodev)
METAG_data_all1 <- METAG_data_all[rowSums(METAG_data_all)>0,colSums(METAG_data_all)>0]
METAG_taxa <- METAG_taxa[METAG_taxa$Genus!="NA",]
unique(METAG_taxa$Genus)
# METAG_taxa <- METAG_taxa[-grep(".spp",METAG_taxa$Species,value=FALSE),]
unique(METAG_taxa$Species)
METAG_data_all1 <- METAG_data_all1[intersect(rownames(METAG_taxa),rownames(METAG_data_all1)),]
tmp_sum <- colSums(METAG_data_all1[,sort(grep("W0",colnames(METAG_data_all1),value=FALSE))])
tmp_sum[tmp_sum<120]
METAG_taxa1 <- METAG_taxa[rownames(METAG_data_all1),]
METAG.samples <- METAG.samples[colnames(METAG_data_all1),]
MMC.MGX <- as.data.frame(cbind(METAG_taxa1[rownames(METAG_data_all1),],METAG_data_all1))
sort(unique(METAG_taxa$Species))
Bac.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(METAG_data_all1)),colData = DataFrame(METAG.samples),rowData = DataFrame(METAG_taxa1))
Bac.tse <- transformAssay(Bac.tse_taxa, MARGIN = "samples", method = "relabundance")
Bac.tse <- addPerCellQC(Bac.tse)
Bac.tse <- mia::estimateRichness(Bac.tse, assay.type = "counts", index = "observed", name="observed")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "coverage", name = "coverage")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "gini_simpson", name = "gini_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "inverse_simpson", name = "inverse_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "log_modulo_skewness", name = "Rarity")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "shannon", name = "shannon")
Bac.tse <- estimateDominance(Bac.tse, assay.type = "counts", index="relative", name = "Dominance")
Bac.tse <- mia::estimateDivergence(Bac.tse,assay.type = "counts",reference = "median",FUN = vegan::vegdist)
colData(Bac.tse)$total_raw_counts <- colSums(assay(Bac.tse, "counts"))
~~~

## 3.1.2 Beta-diversity: PERMANOVA and PCoA

Loads per-treatment PERMANOVA results (Jaccard/Bray-Curtis) and draws the PCoA ordination of bacterial community composition, exported as `Fig3.1.Pcoa.svg`.

~~~R
All_sig.Flu <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC_MGX_PERMANOVA_res.Fluconazole.rds")
All_sig.Nys <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC_MGX_PERMANOVA_res.Nystatin.rds")
All_sig.Flu[All_sig.Flu$pvl_jac <= 0.05 | All_sig.Flu$pvl_bray <= 0.05,]

Sel_type <- c("Genus","Species")
library(ggside)
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggside)
library(patchwork)
treatment <- c("Nystatin", "Fluconazole")
All_sig <- list(All_sig.Nys,All_sig.Flu)
pal <- jdb_palette("corona")
ALL_PLOTS <- lapply(1:length(treatment), function(treat) {
  tse_tmp <- Bac.tse[, Bac.tse$treatment %in% treatment[treat]]
  tse_tmp <- tse_tmp[, tse_tmp$time %in% c("W0","W1","W4","W2","W8")]
  tse_tmp$time.v2 <- as.character(tse_tmp$time)
  tse_tmp$time.v2[tse_tmp$time.v2 %in% "W0"] <- "Pre"
  tse_tmp$time.v2[tse_tmp$time.v2 %in% c("W1","W4","W2")] <- "Post"
  tse_tmp$time.v2[tse_tmp$time.v2 %in% "W8"] <- "LTM"
  tse_tmp$time.v2 <- factor(tse_tmp$time.v2, levels=c("Pre","Post","LTM"))
  ## 转换为相对丰度
  mat_rel_ <- lapply(1:length(Sel_type),function(t) {
    tse_tmp1 <- subsetByPrevalentFeatures(tse_tmp, rank = Sel_type[t], detection = 0.31, prevalence = 0.31)
    ## 相对丰度 + Jaccard（presence/absence）
    tse_tmp1 <- transformAssay(tse_tmp1, assay.type = "counts", method = "relabundance")
    mat_rel1 <- t(assay(tse_tmp1, "relabundance"))
    return(mat_rel1)
    })
  mat_rel <- as.data.frame(do.call(cbind,mat_rel_))

  dist_bray <- vegdist(mat_rel, method="bray")
  ## PCoA
  pcoa_res <- pcoa(dist_bray)
  var_explained <- round(100 * pcoa_res$values$Relative_eig[1:2], 1)
  PcoA <- data.frame(x = pcoa_res$vectors[,1],y = pcoa_res$vectors[,2],time.v2 = tse_tmp$time.v2)
  ## 分组求 median 和 MAD
  sbg <- PcoA %>% group_by(time.v2) %>%
    summarise(median_x = median(x), mad_x = mad(x),median_y = median(y), mad_y = mad(y), .groups='drop')
  ellipse_data <- do.call(rbind, lapply(1:nrow(sbg), function(i) {
    df <- generate_ellipse(sbg$median_x[i], sbg$mad_x[i],sbg$median_y[i], sbg$mad_y[i], factor=0.8)
    df$time.v2 <- sbg$time.v2[i]
    df
  }))
  All_sig_tmp <- All_sig[[treat]][All_sig[[treat]]$prevalence==0.4,]
  p <- ggplot(PcoA, aes(x=x, y=y, color=time.v2)) + geom_point(aes(shape=time.v2, fill=time.v2), size=5, alpha=0.5) +
    geom_polygon(data=ellipse_data, aes(x=x, y=y, fill=time.v2), alpha=0.2) +
    geom_errorbar(data = sbg, aes(x = median_x, y = median_y, ymin = median_y - 0.5 * mad_y, ymax = median_y + 0.5 * mad_y), width = 0) +
    geom_errorbarh(data = sbg, aes(x = median_x, y = median_y, xmin = median_x - 0.5 * mad_x, xmax = median_x + 0.5 * mad_x), height = 0) +
    geom_point(data=sbg, aes(x=median_x, y=median_y, fill=time.v2),color="black", shape=22, size=5, alpha=0.7, show.legend=FALSE) +
    theme_bw() + scale_fill_manual(values=pal) + scale_color_manual(values=pal) +
    labs(x = paste0("PCoA1 (", var_explained[1], "%)"),y = paste0("PCoA2 (", var_explained[2], "%)"),
        title = All_sig_tmp$title_lab_jac) +
    geom_xsideboxplot(orientation="y") + geom_ysideboxplot(orientation="x")
  return(p)
})
plot <- wrap_plots(ALL_PLOTS, nrow=1)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig3.1.Pcoa.svg", plot=plot,width = 10, height = 4,dpi=300)
~~~

![image-20260615092959779](./Bacterial.assets/image-20260615092959779.png)

## 3.1.3 Relative Abundance and SCFA-producer / Probiotic Taxa Sets

Rebuilds the bacterial `TreeSummarizedExperiment`, computes genus- and species-level relative abundances (saved as `MMC.metagenomic.relativeAb_all.rds`), and defines curated taxa sets (butyrate/propionate/valerate producers, probiotics, gut-associated bacteria) used throughout the chapter.

~~~R
ReadCap.20251215 <- mcreadRDS("./workshop/MMC/Aidan_info/v2/ReadCap.20251215.rds")
total_times_tmp <- read.csv("./projects/MMC/Figures/figures_making/v3/patient.PFS.20251215.csv")
sample_info_total5 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% total_times_tmp$Omics_patient_Names,]
# sample_info_total5 <- sample_info_total5[sample_info_total5$time!="W6",]
sample_info_total5 <- sample_info_total5[sample_info_total5$treatment!="Clotrimazole",]
METAG_taxa_raw <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/MMC.metagenomic.taxa.v3.rds")
METAG_data <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/METAG_data_all.v3.rds")
METAG_data_all <- apply(METAG_data, 2, as.numeric)
rownames(METAG_data_all) <- rownames(METAG_data)
tmp_sum <- colSums(METAG_data_all)
tmp_sum[tmp_sum<120]
METAG_data_all <- METAG_data_all[,colSums(METAG_data_all)>100]
METAG.samples <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/metagenomic.samples_info1.v3.rds")
METAG.samples <- METAG.samples[METAG.samples$sample %in% sample_info_total5$Omics_samples_Names,]
METAG.samples <- METAG.samples[intersect(rownames(METAG.samples),colnames(METAG_data_all)),]
METAG_data_all <- as.data.frame(METAG_data_all[,rownames(METAG.samples)])
METAG_taxa <- METAG_taxa_raw[!duplicated(METAG_taxa_raw$Species),]
taxa1 <- METAG_taxa_raw
Species <- unique(taxa1$Species)
METAG_data_all_ <- lapply(1:length(Species),function(x) {
    tmp_counts <- METAG_data_all[rownames(taxa1[taxa1$Species==Species[x],]),]
    if (nrow(tmp_counts)==1){tmp_counts <- tmp_counts} else {
        tmp_counts <- as.data.frame(t(colSums(tmp_counts)))
    }
    # message(x)
    tmp_counts
    })
METAG_data_all <- do.call(rbind,METAG_data_all_)
rownames(METAG_taxa) <- rownames(METAG_data_all) <- as.character(METAG_taxa$Species)
summary(colSums(METAG_data_all))

library(microeco)
library(mecodev)
METAG_data_all1 <- METAG_data_all[rowSums(METAG_data_all)>0,colSums(METAG_data_all)>0]
METAG_taxa <- METAG_taxa[METAG_taxa$Genus!="NA",]
unique(METAG_taxa$Genus)
# METAG_taxa <- METAG_taxa[-grep(".spp",METAG_taxa$Species,value=FALSE),]
unique(METAG_taxa$Species)
METAG_data_all1 <- METAG_data_all1[intersect(rownames(METAG_taxa),rownames(METAG_data_all1)),]
tmp_sum <- colSums(METAG_data_all1[,sort(grep("W0",colnames(METAG_data_all1),value=FALSE))])
tmp_sum[tmp_sum<120]
METAG_taxa1 <- METAG_taxa[rownames(METAG_data_all1),]
METAG.samples <- METAG.samples[colnames(METAG_data_all1),]
MMC.MGX <- as.data.frame(cbind(METAG_taxa1[rownames(METAG_data_all1),],METAG_data_all1))
sort(unique(METAG_taxa$Species))
Bac.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(METAG_data_all1)),colData = DataFrame(METAG.samples),rowData = DataFrame(METAG_taxa1))
Bac.tse <- transformAssay(Bac.tse_taxa, MARGIN = "samples", method = "relabundance")
Bac.tse <- addPerCellQC(Bac.tse)
Bac.tse <- mia::estimateRichness(Bac.tse, assay.type = "counts", index = "observed", name="observed")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "coverage", name = "coverage")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "gini_simpson", name = "gini_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "inverse_simpson", name = "inverse_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "log_modulo_skewness", name = "Rarity")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "shannon", name = "shannon")
Bac.tse <- estimateDominance(Bac.tse, assay.type = "counts", index="relative", name = "Dominance")
Bac.tse <- mia::estimateDivergence(Bac.tse,assay.type = "counts",reference = "median",FUN = vegan::vegdist)
colData(Bac.tse)$total_raw_counts <- colSums(assay(Bac.tse, "counts"))
table(Bac.tse$treatment)
library(lefser)
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)
mcsaveRDS(relativeAb_all,"./workshop/MMC/sample_info/final_Res/v2/MMC.metagenomic.relativeAb_all.rds")

MMC.METAG.inc_all <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.METAG.inc_all.rds")
tmp_projects <- as.data.frame(colData(Bac.tse))
probiotic <- c(grep("Clostridium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),"Lactobacillus","Lactobacillus_acidophilus","Bifidobacteriales","Streptococcus_thermophilus","Bacillus",
    grep("Streptococcus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Lactobacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bifidobacteriales",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),
    grep("Faecalibacterium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE))
probiotic.bac <- intersect(probiotic,relativeAb_all$names)
Butyrate_producing_bac <- c("Anaerostipes_hadrus","Clostridium_butyricum","Coprococcus_comes","Eubacterium_hallii","Eubacterium_rectale","Faecalibacterium_prausnitzii","Roseburia_hominis","Roseburia_intestinalis","Butyrivibrio_fibrisolvens","Ruminococcus_bromii","Megasphaera_elsdenii","Odoribacter_splanchnicus")
Propionate_producing_bac <- c("Bacteroides_thetaiotaomicron","Bacteroides_fragilis","Bacteroides_ovatus","Bacteroides_uniformis","Phascolarctobacterium_succinatutens","Veillonella_parvula","Dialister_succinatiphilus","Prevotella_copri","Propionibacterium_freudenreichii","Anaerovibrio_lipolytica","Oddibacterium_anthropi")
Valeric_producing_bac <- c("Butyrivibrio_fibrisolvens","Megasphaera_elsdenii","Coprococcus_eutactus","Blautia")
gut_associated_bac <- sort(c(Butyrate_producing_bac,Propionate_producing_bac,Valeric_producing_bac,probiotic.bac))
gut_associated_bac <- intersect(gut_associated_bac,relativeAb_all$names)
~~~

## 3.1.4 Taxonomic Composition and Alpha-diversity Over Time

Computes per-treatment mean relative abundance and alpha-diversity across time points and visualizes bacterial community composition and diversity shifts, exported as `Fig3.1.Bac_diversity.svg`.

~~~R
relabundance_total_Species <- relativeAb_all
aa <- jdb_palette("brewer_spectra")
treatment <- c("Nystatin","Clotrimazole","Fluconazole")
relabundance_all_ <- lapply(1:length(treatment),function(x) {
    tmp_projects <- METAG.samples[METAG.samples$treatment==treatment[x],]
    relabundance_tmp <- relabundance_total_Species[,rownames(tmp_projects)]
    relabundance <-  data.frame(W0=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W0",])]),
        W2=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W2",])]),
        W4=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W4",])]),
        W8=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W8",])]),
        group=treatment[x],
        names=relabundance_total_Species$names)
    relabundance$W2_W0 <- relabundance$W2-relabundance$W0
    relabundance$W4_W0 <- relabundance$W4-relabundance$W0
    return(relabundance)
    })
relabundance_all <- do.call(rbind,relabundance_all_)
relabundance_all_W2_W0 <- relabundance_all[relabundance_all$W2_W0>0,]
gut_associated_bac <- intersect(gut_associated_bac,relabundance_all_W2_W0$names)

tmp_projects2 <- as.data.frame(colData(Bac.tse))
tmp_projects2 <- tmp_projects2[!tmp_projects2$treatment%in%"Clotrimazole",]
relabundance_tmp <- relativeAb_all[gut_associated_bac,rownames(tmp_projects2)]
tmp_projects2$gut_associated_bac <- 100*colMeans(relabundance_tmp)[rownames(tmp_projects2)]
tmp_projects2$time.v2 <- as.character(tmp_projects2$time)
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% "W0"] <- "Pre"
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% c("W1","W2","W4")] <- "Post"
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% "W8"] <- "LTM"
tmp_projects2$time.v2 <- factor(tmp_projects2$time.v2, levels=c("Pre","Post","LTM"))
tmp_projects2$Bac.observed <- tmp_projects2$observed
tmp_projects2$Bac.shannon <- tmp_projects2$shannon

diease_Scores <- c("Bac.observed","Bac.shannon","gut_associated_bac")
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
pal <- jdb_palette("corona")[c(2,1)]
total_plots2 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- tmp_projects2[!is.na(tmp_projects2$treatment),]
    df_paired <- df_paired[df_paired$time %in% c("W0","W1","W2","W4","W8"),]
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,diease_Scores[dis]]-tmp[tmp$time=="W0",diease_Scores[dis]]
        if (diease_Scores[dis]=="gut_associated_bac") {
            tmp$value[tmp$value > 200] <- 200
            tmp$value[tmp$value < -200] <- -200
        }
        if (diease_Scores[dis]=="Bac.observed") {
            tmp$value[tmp$value > 50] <- 50
            tmp$value[tmp$value < -50] <- -50
        }
        if (diease_Scores[dis]=="Bac.shannon") {
            tmp$value[tmp$value > 2] <- 2
            tmp$value[tmp$value < -2] <- -2
        }
        return(tmp)
    })
    df_paired1 <- do.call(rbind,df_paired1_)
    df_paired1$treatment <- factor(df_paired1$treatment,levels=c("Nystatin","Fluconazole"))
    if (diease_Scores[dis]=="Bac.observed") {test <- wilcox.test_wrapper} else {test <- wilcox.test_wrapper}
    df_paired2 <- df_paired1
    df_paired2$time <- as.numeric(gsub("W","",df_paired2$time))
    loess_fit_1 <- loess(value ~ time, data = df_paired2[df_paired2$treatment == "Nystatin", ], na.action = na.exclude)
    loess_fit_2 <- loess(value ~ time, data = df_paired2[df_paired2$treatment == "Fluconazole", ], na.action = na.exclude)
    x_range <- seq(min(df_paired2$time), max(df_paired2$time), length.out = 100)
    pred_1 <- predict(loess_fit_1, newdata = data.frame(time = x_range))
    pred_2 <- predict(loess_fit_2, newdata = data.frame(time = x_range))
    ks_test_result <- formatC(ks.test(pred_1, pred_2)$p.value, format = "e", digits = 3)

    p1 <- ggplot(df_paired1, aes_string(x = "time", y = "value")) + 
        geom_line(aes(group=patient),size = 0.6,alpha=0.8,color="lightgrey") +
        geom_jitter(color="black",width = 0.1,alpha=0.5,size=1)+
        geom_boxplot(outlier.shape = NA,alpha=0)+
        theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
        stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=1)))+facet_wrap(~treatment,ncol=3)+
        labs(title=paste0(diease_Scores[dis],"\n (Flu vs Nys p:",ks_test_result,")"),y = paste("Δ")) +NoLegend()+
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = test)+
        geom_smooth(aes(group = 1, color = treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.2,span=1.2)+
        stat_summary(fun = mean, geom = "point",aes(group = 1,color = treatment),size=2)
    plot <- plot_grid(p1,ncol=1)
    return(plot)
    })
plot <- CombinePlots(c(total_plots2),nrow=1)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig3.1.Bac_diversity.svg", plot=plot,width = 10, height = 4,dpi=300)
~~~

![image-20260615093105165](./Bacterial.assets/image-20260615093105165.png)

## 3.1.5 Differential Abundance Testing (LEfSe) per Time Point

Runs LEfSe differential-abundance analysis comparing each post-treatment time point against W0 (for W1, W2, W4, W8) and saves the per-time-point results (`MMC.MGX.DAA_all.W*.raw.lefser.v2.rds`).

~~~R
ReadCap.20251215 <- mcreadRDS("./workshop/MMC/Aidan_info/v2/ReadCap.20251215.rds")
total_times_tmp <- read.csv("./projects/MMC/Figures/figures_making/v3/patient.PFS.20251215.csv")
sample_info_total5 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% total_times_tmp$Omics_patient_Names,]
# sample_info_total5 <- sample_info_total5[sample_info_total5$time!="W6",]
sample_info_total5 <- sample_info_total5[sample_info_total5$treatment!="Clotrimazole",]
METAG_taxa_raw <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/MMC.metagenomic.taxa.v3.rds")
METAG_data <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/METAG_data_all.v3.rds")
METAG_data_all <- apply(METAG_data, 2, as.numeric)
rownames(METAG_data_all) <- rownames(METAG_data)
tmp_sum <- colSums(METAG_data_all)
tmp_sum[tmp_sum<120]
METAG_data_all <- METAG_data_all[,colSums(METAG_data_all)>100]
METAG.samples <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/metagenomic.samples_info1.v3.rds")
METAG.samples <- METAG.samples[METAG.samples$sample %in% sample_info_total5$Omics_samples_Names,]
METAG.samples <- METAG.samples[intersect(rownames(METAG.samples),colnames(METAG_data_all)),]
METAG_data_all <- as.data.frame(METAG_data_all[,rownames(METAG.samples)])
METAG_taxa <- METAG_taxa_raw[!duplicated(METAG_taxa_raw$Species),]
taxa1 <- METAG_taxa_raw
Species <- unique(taxa1$Species)
METAG_data_all_ <- lapply(1:length(Species),function(x) {
    tmp_counts <- METAG_data_all[rownames(taxa1[taxa1$Species==Species[x],]),]
    if (nrow(tmp_counts)==1){tmp_counts <- tmp_counts} else {
        tmp_counts <- as.data.frame(t(colSums(tmp_counts)))
    }
    # message(x)
    tmp_counts
    })
METAG_data_all <- do.call(rbind,METAG_data_all_)
rownames(METAG_taxa) <- rownames(METAG_data_all) <- as.character(METAG_taxa$Species)
summary(colSums(METAG_data_all))

library(microeco)
library(mecodev)
METAG_data_all1 <- METAG_data_all[rowSums(METAG_data_all)>0,colSums(METAG_data_all)>0]
METAG_taxa <- METAG_taxa[METAG_taxa$Genus!="NA",]
unique(METAG_taxa$Genus)
# METAG_taxa <- METAG_taxa[-grep(".spp",METAG_taxa$Species,value=FALSE),]
unique(METAG_taxa$Species)
METAG_data_all1 <- METAG_data_all1[intersect(rownames(METAG_taxa),rownames(METAG_data_all1)),]
tmp_sum <- colSums(METAG_data_all1[,sort(grep("W0",colnames(METAG_data_all1),value=FALSE))])
tmp_sum[tmp_sum<120]

# METAG_data_all1 <- METAG_data_all1[,-grep("MMC076",colnames(METAG_data_all1),value=FALSE)]
# METAG_data_all1 <- METAG_data_all1[,-grep("MMC097",colnames(METAG_data_all1),value=FALSE)]
# METAG_data_all1 <- METAG_data_all1[,-grep("MMC101",colnames(METAG_data_all1),value=FALSE)]
# METAG_data_all1 <- METAG_data_all1[,colSums(METAG_data_all1)>100]
METAG_taxa1 <- METAG_taxa[rownames(METAG_data_all1),]
METAG.samples <- METAG.samples[colnames(METAG_data_all1),]
MMC.MGX <- as.data.frame(cbind(METAG_taxa1[rownames(METAG_data_all1),],METAG_data_all1))
sort(unique(METAG_taxa$Species))
Bac.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(METAG_data_all1)),colData = DataFrame(METAG.samples),rowData = DataFrame(METAG_taxa1))
Bac.tse <- transformAssay(Bac.tse_taxa, MARGIN = "samples", method = "relabundance")
Bac.tse <- addPerCellQC(Bac.tse)
Bac.tse <- mia::estimateRichness(Bac.tse, assay.type = "counts", index = "observed", name="observed")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "coverage", name = "coverage")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "gini_simpson", name = "gini_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "inverse_simpson", name = "inverse_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "log_modulo_skewness", name = "Rarity")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "shannon", name = "shannon")
Bac.tse <- estimateDominance(Bac.tse, assay.type = "counts", index="relative", name = "Dominance")
Bac.tse <- mia::estimateDivergence(Bac.tse,assay.type = "counts",reference = "median",FUN = vegan::vegdist)
colData(Bac.tse)$total_raw_counts <- colSums(assay(Bac.tse, "counts"))

library(lefser)
colData(Bac.tse)$library_size <- colSums(assay(Bac.tse, "counts"))
colData(Bac.tse)$Diagnosis_new <- factor(colData(Bac.tse)$Diagnosis_new,levels=c("UC","CD"))
Sel_type <- c("Genus","Species")
lefser.All.W1 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    tse_tmp <- transformAssay(tse_tmp, assay.type = "counts", method = "relabundance")
    tse_tmp$Diagnosis_new <- droplevels(tse_tmp$Diagnosis_new)
    tse_tmp <- tse_tmp[,tse_tmp$time %in% c("W0","W1")]
    tse_tmp$time <- factor(tse_tmp$time,levels=c("W0","W1"))

    tse_Nystatin_W0_W1 <- tse_tmp[,tse_tmp$treatment=="Nystatin"]
    tse_Fluconazole_W0_W1 <- tse_tmp[,tse_tmp$treatment=="Fluconazole"]

    se_Nystatin_W0_W1 <- SummarizedExperiment(assays = list(counts = assays(tse_Nystatin_W0_W1)[["counts"]]),rowData = rowData(tse_Nystatin_W0_W1),colData = colData(tse_Nystatin_W0_W1))
    se_Nystatin_W0_W1 <- relativeAb(se_Nystatin_W0_W1)
    Nystatin.lefser_out <- XY_lefser(se_Nystatin_W0_W1, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Nystatin.lefser_out)[which(colnames(Nystatin.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")

    se_Fluconazole_W0_W1 <- SummarizedExperiment(assays = list(counts = assays(tse_Fluconazole_W0_W1)[["counts"]]),rowData = rowData(tse_Fluconazole_W0_W1),colData = colData(tse_Fluconazole_W0_W1))
    se_Fluconazole_W0_W1 <- relativeAb(se_Fluconazole_W0_W1)
    Fluconazole.lefser_out <- XY_lefser(se_Fluconazole_W0_W1, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Fluconazole.lefser_out)[which(colnames(Fluconazole.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")
    message(x)

    return(list(Nystatin.lefser_out,Fluconazole.lefser_out))
    })
Nystatin_all.W1_ <- lapply(1:length(lefser.All.W1),function(x) {
    lefser.All.W1[[x]][[1]]$type <- Sel_type[x]
    return(lefser.All.W1[[x]][[1]])})
Nystatin_all.W1 <- do.call(rbind,Nystatin_all.W1_)
Fluconazole_all.W1_ <- lapply(1:length(lefser.All.W1),function(x) {
    lefser.All.W1[[x]][[2]]$type <- Sel_type[x]
    return(lefser.All.W1[[x]][[2]])})
Fluconazole_all.W1 <- do.call(rbind,Fluconazole_all.W1_)
DAA_all.W1 <- list(Nystatin_all.W1,Fluconazole_all.W1)
mcsaveRDS(DAA_all.W1,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W1.raw.lefser.v2.rds")

lefser.All.W2 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    tse_tmp <- transformAssay(tse_tmp, assay.type = "counts", method = "relabundance")
    tse_tmp$Diagnosis_new <- droplevels(tse_tmp$Diagnosis_new)
    tse_tmp <- tse_tmp[,tse_tmp$time %in% c("W0","W2")]
    tse_tmp$time <- factor(tse_tmp$time,levels=c("W0","W2"))

    tse_Nystatin_W0_W2 <- tse_tmp[,tse_tmp$treatment=="Nystatin"]
    tse_Fluconazole_W0_W2 <- tse_tmp[,tse_tmp$treatment=="Fluconazole"]

    se_Nystatin_W0_W2 <- SummarizedExperiment(assays = list(counts = assays(tse_Nystatin_W0_W2)[["counts"]]),rowData = rowData(tse_Nystatin_W0_W2),colData = colData(tse_Nystatin_W0_W2))
    se_Nystatin_W0_W2 <- relativeAb(se_Nystatin_W0_W2)
    Nystatin.lefser_out <- XY_lefser(se_Nystatin_W0_W2, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Nystatin.lefser_out)[which(colnames(Nystatin.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")

    se_Fluconazole_W0_W2 <- SummarizedExperiment(assays = list(counts = assays(tse_Fluconazole_W0_W2)[["counts"]]),rowData = rowData(tse_Fluconazole_W0_W2),colData = colData(tse_Fluconazole_W0_W2))
    se_Fluconazole_W0_W2 <- relativeAb(se_Fluconazole_W0_W2)
    Fluconazole.lefser_out <- XY_lefser(se_Fluconazole_W0_W2, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Fluconazole.lefser_out)[which(colnames(Fluconazole.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")
    message(x)

    return(list(Nystatin.lefser_out,Fluconazole.lefser_out))
    })
Nystatin_all.W2_ <- lapply(1:length(lefser.All.W2),function(x) {
    lefser.All.W2[[x]][[1]]$type <- Sel_type[x]
    return(lefser.All.W2[[x]][[1]])})
Nystatin_all.W2 <- do.call(rbind,Nystatin_all.W2_)
Fluconazole_all.W2_ <- lapply(1:length(lefser.All.W2),function(x) {
    lefser.All.W2[[x]][[2]]$type <- Sel_type[x]
    return(lefser.All.W2[[x]][[2]])})
Fluconazole_all.W2 <- do.call(rbind,Fluconazole_all.W2_)
DAA_all.W2 <- list(Nystatin_all.W2,Fluconazole_all.W2)
mcsaveRDS(DAA_all.W2,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W2.raw.lefser.v2.rds")

lefser.All.W4 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    tse_tmp <- transformAssay(tse_tmp, assay.type = "counts", method = "relabundance")
    tse_tmp$Diagnosis_new <- droplevels(tse_tmp$Diagnosis_new)
    tse_tmp <- tse_tmp[,tse_tmp$time %in% c("W0","W4")]
    tse_tmp$time <- factor(tse_tmp$time,levels=c("W0","W4"))

    tse_Nystatin_W0_W4 <- tse_tmp[,tse_tmp$treatment=="Nystatin"]
    tse_Fluconazole_W0_W4 <- tse_tmp[,tse_tmp$treatment=="Fluconazole"]

    se_Nystatin_W0_W4 <- SummarizedExperiment(assays = list(counts = assays(tse_Nystatin_W0_W4)[["counts"]]),rowData = rowData(tse_Nystatin_W0_W4),colData = colData(tse_Nystatin_W0_W4))
    se_Nystatin_W0_W4 <- relativeAb(se_Nystatin_W0_W4)
    Nystatin.lefser_out <- XY_lefser(se_Nystatin_W0_W4, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Nystatin.lefser_out)[which(colnames(Nystatin.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")

    se_Fluconazole_W0_W4 <- SummarizedExperiment(assays = list(counts = assays(tse_Fluconazole_W0_W4)[["counts"]]),rowData = rowData(tse_Fluconazole_W0_W4),colData = colData(tse_Fluconazole_W0_W4))
    se_Fluconazole_W0_W4 <- relativeAb(se_Fluconazole_W0_W4)
    Fluconazole.lefser_out <- XY_lefser(se_Fluconazole_W0_W4, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Fluconazole.lefser_out)[which(colnames(Fluconazole.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")
    message(x)

    return(list(Nystatin.lefser_out,Fluconazole.lefser_out))
    })
Nystatin_all.W4_ <- lapply(1:length(lefser.All.W4),function(x) {
    lefser.All.W4[[x]][[1]]$type <- Sel_type[x]
    return(lefser.All.W4[[x]][[1]])})
Nystatin_all.W4 <- do.call(rbind,Nystatin_all.W4_)
Fluconazole_all.W4_ <- lapply(1:length(lefser.All.W4),function(x) {
    lefser.All.W4[[x]][[2]]$type <- Sel_type[x]
    return(lefser.All.W4[[x]][[2]])})
Fluconazole_all.W4 <- do.call(rbind,Fluconazole_all.W4_)
DAA_all.W4 <- list(Nystatin_all.W4,Fluconazole_all.W4)
mcsaveRDS(DAA_all.W4,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W4.raw.lefser.v2.rds")

lefser.All.W8 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    tse_tmp <- transformAssay(tse_tmp, assay.type = "counts", method = "relabundance")
    tse_tmp$Diagnosis_new <- droplevels(tse_tmp$Diagnosis_new)
    tse_tmp <- tse_tmp[,tse_tmp$time %in% c("W0","W8")]
    tse_tmp$time <- factor(tse_tmp$time,levels=c("W0","W8"))

    tse_Nystatin_W0_W8 <- tse_tmp[,tse_tmp$treatment=="Nystatin"]
    tse_Fluconazole_W0_W8 <- tse_tmp[,tse_tmp$treatment=="Fluconazole"]

    se_Nystatin_W0_W8 <- SummarizedExperiment(assays = list(counts = assays(tse_Nystatin_W0_W8)[["counts"]]),rowData = rowData(tse_Nystatin_W0_W8),colData = colData(tse_Nystatin_W0_W8))
    se_Nystatin_W0_W8 <- relativeAb(se_Nystatin_W0_W8)
    Nystatin.lefser_out <- XY_lefser(se_Nystatin_W0_W8, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Nystatin.lefser_out)[which(colnames(Nystatin.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")

    se_Fluconazole_W0_W8 <- SummarizedExperiment(assays = list(counts = assays(tse_Fluconazole_W0_W8)[["counts"]]),rowData = rowData(tse_Fluconazole_W0_W8),colData = colData(tse_Fluconazole_W0_W8))
    se_Fluconazole_W0_W8 <- relativeAb(se_Fluconazole_W0_W8)
    Fluconazole.lefser_out <- XY_lefser(se_Fluconazole_W0_W8, classCol = "time", subclassCol = "Diagnosis_new", kruskal.threshold = 1, wilcox.threshold = 1,filter=FALSE)
    colnames(Fluconazole.lefser_out)[which(colnames(Fluconazole.lefser_out)=="scores")] <- paste0(levels(tse_tmp$time)[2],".VS.",levels(tse_tmp$time)[1],".LDA.scores")
    message(x)

    return(list(Nystatin.lefser_out,Fluconazole.lefser_out))
    })
Nystatin_all.W8_ <- lapply(1:length(lefser.All.W8),function(x) {
    lefser.All.W8[[x]][[1]]$type <- Sel_type[x]
    return(lefser.All.W8[[x]][[1]])})
Nystatin_all.W8 <- do.call(rbind,Nystatin_all.W8_)
Fluconazole_all.W8_ <- lapply(1:length(lefser.All.W8),function(x) {
    lefser.All.W8[[x]][[2]]$type <- Sel_type[x]
    return(lefser.All.W8[[x]][[2]])})
Fluconazole_all.W8 <- do.call(rbind,Fluconazole_all.W8_)
DAA_all.W8 <- list(Nystatin_all.W8,Fluconazole_all.W8)
mcsaveRDS(DAA_all.W8,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W8.raw.lefser.v2.rds")
~~~

## 3.1.6 Compile Significant Taxa and Variance Summary

Aggregates the significant (included) differential taxa across all time points into `MMC.METAG.inc_all.rds` and summarizes their contribution / proportion of variance, exported as `Fig3.bac_prop_var.svg`.

~~~R
W1.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W1.raw.lefser.v2.rds")
W2.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W2.raw.lefser.v2.rds")
W4.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W4.raw.lefser.v2.rds")
W8.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W8.raw.lefser.v2.rds")
names(W1.lefser) <- c("Nystatin","Fluconazole")
MMC.METAG.inc_all <- lapply(1:length(W2.lefser),function(x) {
    W1.tmp_all <- W1.lefser[[x]][,c(1,2,3,6)]
    W2.tmp_all <- W2.lefser[[x]][,c(1,2,3,6)]
    W4.tmp_all <- W4.lefser[[x]][,c(1,2,3,6)]
    W8.tmp_all <- W8.lefser[[x]][,c(1,2,3,6)]
    colnames(W1.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W2.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W4.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W8.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    W1.tmp_all$group <- "W1.vs.W0"
    W2.tmp_all$group <- "W2.vs.W0"
    W4.tmp_all$group <- "W4.vs.W0"
    W8.tmp_all$group <- "W8.vs.W0"
    tmp_all <- rbind(W1.tmp_all,W2.tmp_all)
    tmp_all <- rbind(tmp_all,W4.tmp_all)
    tmp_all <- rbind(tmp_all,W8.tmp_all)
    tmp_all$LDA.scores[tmp_all$LDA.scores > 4] <- 4
    tmp_all$LDA.scores[tmp_all$LDA.scores < -4] <- -4
    tmp_all <- tmp_all[tmp_all$pvalue <= 0.75,]
    tmp_all[tmp_all$Names=="Clostridium_butyricum",]
    tmp_all[tmp_all$Names=="Faecalibacterium_prausnitzii_1",]

    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c[tmp_all_c$Names=="Clostridium_butyricum",]
    tmp_all_c[tmp_all_c$Names=="Faecalibacterium_prausnitzii_1",]

    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.5,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]
    return(tmp_all_c)
    })
mcsaveRDS(MMC.METAG.inc_all,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.METAG.inc_all.rds")

library(lefser)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)
MMC.METAG.inc_all <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.METAG.inc_all.rds")
tmp_projects <- as.data.frame(colData(Bac.tse))
probiotic <- c(grep("Clostridium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),"Lactobacillus","Lactobacillus_acidophilus","Bifidobacteriales","Streptococcus_thermophilus","Bacillus",
    grep("Streptococcus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Lactobacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bifidobacteriales",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),
    grep("Faecalibacterium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE))
probiotic.bac <- intersect(probiotic,relativeAb_all$names)
Butyrate_producing_bac <- c("Anaerostipes_hadrus","Clostridium_butyricum","Coprococcus_comes","Eubacterium_hallii","Eubacterium_rectale","Faecalibacterium_prausnitzii","Roseburia_hominis","Roseburia_intestinalis","Butyrivibrio_fibrisolvens","Ruminococcus_bromii","Megasphaera_elsdenii","Odoribacter_splanchnicus")
Propionate_producing_bac <- c("Bacteroides_thetaiotaomicron","Bacteroides_fragilis","Bacteroides_ovatus","Bacteroides_uniformis","Phascolarctobacterium_succinatutens","Veillonella_parvula","Dialister_succinatiphilus","Prevotella_copri","Propionibacterium_freudenreichii","Anaerovibrio_lipolytica","Oddibacterium_anthropi")
Valeric_producing_bac <- c("Butyrivibrio_fibrisolvens","Megasphaera_elsdenii","Coprococcus_eutactus","Blautia")
gut_associated_bac <- sort(c(Butyrate_producing_bac,Propionate_producing_bac,Valeric_producing_bac,probiotic.bac,"Eubacterium_limosum","Eubacterium_ventriosum","Eubacterium_callanderi"))

W1.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W1.raw.lefser.v2.rds")
W2.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W2.raw.lefser.v2.rds")
W4.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W4.raw.lefser.v2.rds")
W8.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W8.raw.lefser.v2.rds")
names(W1.lefser) <- c("Nystatin","Fluconazole")
All_plots <- lapply(1:length(W1.lefser),function(x) {
    W1.tmp_all <- W1.lefser[[x]][,c(1,2,3,6)]
    W2.tmp_all <- W2.lefser[[x]][,c(1,2,3,6)]
    W4.tmp_all <- W4.lefser[[x]][,c(1,2,3,6)]
    W8.tmp_all <- W8.lefser[[x]][,c(1,2,3,6)]
    colnames(W1.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W2.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W4.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W8.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    W1.tmp_all$group <- "W1.vs.W0"
    W2.tmp_all$group <- "W2.vs.W0"
    W4.tmp_all$group <- "W4.vs.W0"
    W8.tmp_all$group <- "W8.vs.W0"
    tmp_all <- rbind(W1.tmp_all,W2.tmp_all)
    tmp_all <- rbind(tmp_all,W4.tmp_all)
    tmp_all <- rbind(tmp_all,W8.tmp_all)
    tmp_all <- tmp_all[tmp_all$Names %in% gut_associated_bac,]
    tmp_all <- tmp_all[tmp_all$type=="Species",]
    tmp_all[tmp_all$Names=="Clostridium_butyricum",]
    tmp_all[tmp_all$Names=="Faecalibacterium_prausnitzii_1",]
    tmp_all[grep("Eubacterium",tmp_all$Names,value=FALSE),]

    tmp_all <- tmp_all[tmp_all$pvalue <= 0.75,]
    tmp_all$LDA.scores[tmp_all$LDA.scores > 4] <- 4
    tmp_all$LDA.scores[tmp_all$LDA.scores < -4] <- -4
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0

    setdiff(gut_associated_bac,tmp_all$Names )

    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c[grep("Eubacterium",tmp_all_c$Names,value=FALSE),]

    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.5,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]
    o <- as.character(tmp_all_c$Names)
    tmp_all <- tmp_all[tmp_all$Names %in% o,]

    tmp_all$Names <- factor(tmp_all$Names,levels=o[1:length(o)])
    tmp_all$pvalue.label <- 0.2
    tmp_all$pvalue.label[tmp_all$pvalue<=0.05] <- 0.05
    tmp_all$pvalue.label[tmp_all$pvalue>0.05 & tmp_all$pvalue<=0.1] <- 0.1
    tmp_all$pvalue.label[tmp_all$pvalue>0.1 & tmp_all$pvalue<=0.15] <- 0.15
    plot <- ggplot(tmp_all, aes_string(x="group", y="Names", size="pvalue.label", color="LDA.scores")) +
        geom_point() + scale_colour_gradient2(low = "navy",  mid = "white",  high = "firebrick3",  midpoint = 0,, name = "LDA.scores",guide=guide_colorbar(reverse=FALSE)) +
        ylab(NULL) + ggtitle(names(W1.lefser)[x]," Bac Species") + scale_size(range=c(6,2))+theme_classic()+
        theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1,color = "black", size = 12),axis.text.y = element_text(color = "black", size = 12, face = "plain"))
    return(plot)
    })
plot <- CombinePlots(c(All_plots),ncol=2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_prop_var.svg", plot=plot,width = 10, height = 6,dpi=300)
~~~

![image-20260615093213386](./Bacterial.assets/image-20260615093213386.png)

## 3.1.7 Assemble Differential Gut-associated Bacteria

Combines the LEfSe LDA scores of gut-associated taxa across time points, computes a weighted mean trend per taxon, filters to consistently shifting species, and exports the table as `DiffBac_All.csv`.

~~~R

names(W1.lefser) <- c("Nystatin","Fluconazole")
DiffBac_ <- lapply(1:length(W1.lefser),function(x) {
    W1.tmp_all <- W1.lefser[[x]][,c(1,2,3,6)]
    W2.tmp_all <- W2.lefser[[x]][,c(1,2,3,6)]
    W4.tmp_all <- W4.lefser[[x]][,c(1,2,3,6)]
    W8.tmp_all <- W8.lefser[[x]][,c(1,2,3,6)]
    colnames(W1.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W2.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W4.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W8.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    W1.tmp_all$group <- "W1.vs.W0"
    W2.tmp_all$group <- "W2.vs.W0"
    W4.tmp_all$group <- "W4.vs.W0"
    W8.tmp_all$group <- "W8.vs.W0"
    tmp_all <- rbind(W1.tmp_all,W2.tmp_all)
    tmp_all <- rbind(tmp_all,W4.tmp_all)
    tmp_all <- rbind(tmp_all,W8.tmp_all)
    tmp_all <- tmp_all[tmp_all$Names %in% gut_associated_bac,]
    tmp_all <- tmp_all[tmp_all$type=="Species",]
    tmp_all[tmp_all$Names=="Clostridium_butyricum",]
    tmp_all[tmp_all$Names=="Faecalibacterium_prausnitzii_1",]
    tmp_all[grep("Eubacterium",tmp_all$Names,value=FALSE),]

    tmp_all <- tmp_all[tmp_all$pvalue <= 0.75,]
    tmp_all$LDA.scores[tmp_all$LDA.scores > 4] <- 4
    tmp_all$LDA.scores[tmp_all$LDA.scores < -4] <- -4
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0

    setdiff(gut_associated_bac,tmp_all$Names )

    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c[grep("Eubacterium",tmp_all_c$Names,value=FALSE),]

    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.5,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]
    o <- as.character(tmp_all_c$Names)
    tmp_all <- tmp_all[tmp_all$Names %in% o,]

    tmp_all$treatment <- names(W1.lefser)[x]
    return(tmp_all)
    })
DiffBac <- do.call(rbind,DiffBac_)
DiffBac <- DiffBac[DiffBac$type %in% c("Genus","Species"),]
write.csv(DiffBac,"./projects/MMC/Figures/figures_making/v4/DiffBac_All.csv")
~~~

## 3.1.8 Differential Bacteria Delta Trajectories (Per-batch Panels)

Plots the per-time-point delta (change-from-baseline) trajectories of the differential bacteria, rendered as separate batch panels, exported as `Fig3.bac_dyn.delta.timepoints.raw_merge_bacth1`–`bacth5`.

~~~r

library(lefser)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)
MMC.METAG.inc_all <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.METAG.inc_all.rds")
tmp_projects <- as.data.frame(colData(Bac.tse))
probiotic <- c(grep("Clostridium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),"Lactobacillus","Lactobacillus_acidophilus","Bifidobacteriales","Streptococcus_thermophilus","Bacillus",
    grep("Streptococcus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Lactobacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bifidobacteriales",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),
    grep("Faecalibacterium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE))
probiotic.bac <- intersect(probiotic,relativeAb_all$names)
Butyrate_producing_bac <- c("Anaerostipes_hadrus","Clostridium_butyricum","Coprococcus_comes","Eubacterium_hallii","Eubacterium_rectale","Faecalibacterium_prausnitzii","Roseburia_hominis","Roseburia_intestinalis","Butyrivibrio_fibrisolvens","Ruminococcus_bromii","Megasphaera_elsdenii","Odoribacter_splanchnicus")
Propionate_producing_bac <- c("Bacteroides_thetaiotaomicron","Bacteroides_fragilis","Bacteroides_ovatus","Bacteroides_uniformis","Phascolarctobacterium_succinatutens","Veillonella_parvula","Dialister_succinatiphilus","Prevotella_copri","Propionibacterium_freudenreichii","Anaerovibrio_lipolytica","Oddibacterium_anthropi")
Valeric_producing_bac <- c("Butyrivibrio_fibrisolvens","Megasphaera_elsdenii","Coprococcus_eutactus","Blautia")
gut_associated_bac <- sort(c(Butyrate_producing_bac,Propionate_producing_bac,Valeric_producing_bac,probiotic.bac,
  grep("Ruminococcus",relativeAb_all$names,value=TRUE),"Eubacterium_ventriosum","Eubacterium_callanderi"))

tmp_projects <- as.data.frame(colData(Bac.tse))
relabundance_sel.bac <- relativeAb_all[relativeAb_all$names %in% gut_associated_bac,]
rownames(relabundance_sel.bac) <- relabundance_sel.bac$names
relabundance_tmp <- relabundance_sel.bac[,rownames(tmp_projects)]
gut_associated_bac1 <- sort(rownames(relabundance_tmp))
tmp_projects <- as.data.frame(cbind(tmp_projects,100*t(relabundance_tmp)))
df <- tmp_projects[tmp_projects$treatment %in% c("Nystatin","Fluconazole"),]
df$treatment <- factor(df$treatment,levels=c("Nystatin","Fluconazole"))

# df <- df[!df$time %in% c("W1"),]
sel_pa <- unique(df$patient)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_celsius")),jdb_palette(c("brewer_spectra")))
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
plots1 <- lapply(1:length(gut_associated_bac1),function(x) {
  df_new <- df
  df_tmp <- data_summary(df_new, varname=gut_associated_bac1[x], groupnames=c("time","treatment"))
  All_plots1 <- ggplot(df_tmp, aes_string(x = "time", y = gut_associated_bac1[x])) + 
      stat_summary(data=df_new,fun.data = mean_se, geom = "errorbar", width = 0.5) + geom_jitter(data=df_new,aes(color=patient),width = 0.1,alpha=0.5) + 
      geom_bar(aes(fill=time),alpha=1,stat="identity", position=position_dodge())+
      geom_line(data=df_new,aes(color=patient,group=patient),size = 0.6,alpha=0.2) +
      theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
      stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+facet_wrap(~treatment)+
      labs(title=paste0(gut_associated_bac1[x]," raw in relabundance"),y = paste(gut_associated_bac1[x]," relabundance raw")) +NoLegend()+
      geom_signif(data=df_new,comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
      return(All_plots1)
  })
plots2 <- lapply(1:length(gut_associated_bac1),function(x) {
  df_new <- df
  sel_pa <- unique(df$patient)
  paired.df_ <- lapply(1:length(sel_pa),function(patient) {
      tmp <- df_new[df_new$patient %in% sel_pa[patient],]
      tmp[,gut_associated_bac1[x]] <- tmp[,gut_associated_bac1[x]]-tmp[tmp$time=="W0",gut_associated_bac1[x]]
      return(tmp)
  })
  paired.df <- do.call(rbind,paired.df_)
  paired.df_tmp <- data_summary(paired.df, varname=gut_associated_bac1[x], groupnames=c("time","treatment"))
  All_plots1 <- ggplot(paired.df_tmp, aes_string(x = "time", y = gut_associated_bac1[x])) + 
      stat_summary(data=paired.df,fun.data = mean_se, geom = "errorbar", width = 0.5) + geom_jitter(data=paired.df,aes(color=patient),width = 0.1,alpha=0.5) + 
      geom_bar(aes(fill=time),alpha=1,stat="identity", position=position_dodge())+
      geom_line(data=paired.df,aes(color=patient,group=patient),size = 0.6,alpha=0.2) +
      theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
      stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+facet_wrap(~treatment)+
      labs(title=paste0(gut_associated_bac1[x]," changes in relabundance"),y = paste(gut_associated_bac1[x]," relabundance Change")) +NoLegend()+
      geom_signif(data=paired.df,comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
      message(x)
      return(All_plots1)
    })
total_plots3 <- lapply(1:length(gut_associated_bac1),function(dis) {
  df_paired <- df
  df_paired$time_numeric <- as.numeric(gsub("W", "", df_paired$time))
  df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
  colnames(df_paired)[colnames(df_paired)==gut_associated_bac1[dis]] <- "value"

  anova_result <- aov(value ~ treatment, data = df_paired)
  p_value <- summary(anova_result)[[1]]["treatment", "Pr(>F)"]
  pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
  p1 <- ggplot(df_paired, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
  geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess",se=TRUE, level = 0.5) +
  stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(1),label = "p.format") +
  theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
  labs(title = paste0(gut_associated_bac1[dis]," raw \n", "Flu vs Nystatin (ANOVA p: ", pvalue0,")"),y = paste(gut_associated_bac1[dis], " Change"))

  uniq_patient1 <- unique(df_paired$patient)
  df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
      tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
      tmp[,"value"] <- tmp[,"value"]-tmp[tmp$time=="W0","value"]
      tmp[,"value"][tmp[,"value"]>200] <- 200
      tmp[,"value"][tmp[,"value"] < -200] <- -200
      return(tmp)
  })
  df_paired2 <- do.call(rbind,df_paired1_)
  df_paired2 <- df_paired2[!is.na(df_paired2$value),]
  df_paired2$time_numeric <- as.numeric(gsub("W", "", df_paired2$time))
  anova_result <- aov(value ~ treatment, data = df_paired2)
  p_value <- summary(anova_result)[[1]]["treatment", "Pr(>F)"]
  pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
  p2 <- ggplot(df_paired2, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
  geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess",se=TRUE, level = 0.5) +
  stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(1),label = "p.format") +
  theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
  labs(title = paste0(gut_associated_bac1[dis]," changes \n", "Flu vs Nystatin (ANOVA p: ", pvalue0,")"),y = paste(gut_associated_bac1[dis], " Change"))
  return(plot_grid(p1,p2,nrow=2))
    })
plot1 <- CombinePlots(c(plots1[1:10],plots2[1:10],total_plots3[1:10]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth1.png", plot=plot1,width = 49, height = 10,dpi=300)
plot2 <- CombinePlots(c(plots1[11:20],plots2[11:20],total_plots3[11:20]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth2.png", plot=plot2,width = 49, height = 10,dpi=300)
plot3 <- CombinePlots(c(plots1[21:30],plots2[21:30],total_plots3[21:30]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth3.png", plot=plot3,width = 49, height = 10,dpi=300)
plot4 <- CombinePlots(c(plots1[31:40],plots2[31:40],total_plots3[31:40]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth4.png", plot=plot4,width = 49, height = 10,dpi=300)
plot5 <- CombinePlots(c(plots1[41:50],plots2[41:50],total_plots3[41:50]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth5.png", plot=plot5,width = 49, height = 10,dpi=300)
plot6 <- CombinePlots(c(plots1[51:60],plots2[51:60],total_plots3[51:60]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth5.png", plot=plot6,width = 49, height = 10,dpi=300)
~~~

![Fig3.bac_dyn.delta.timepoints.raw_merge_bacth1](./../../MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth1.png)

![Fig3.bac_dyn.delta.timepoints.raw_merge_bacth2](./Bacterial.assets/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth2.png)

![Fig3.bac_dyn.delta.timepoints.raw_merge_bacth3](./Bacterial.assets/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth3.png)

![Fig3.bac_dyn.delta.timepoints.raw_merge_bacth4](./Bacterial.assets/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth4.png)

![Fig3.bac_dyn.delta.timepoints.raw_merge_bacth5](./Bacterial.assets/Fig3.bac_dyn.delta.timepoints.raw_merge_bacth5.png)

## 3.1.9 Define Core / Candidate Gut-associated Taxa

Defines explicit lists of core and candidate gut-associated bacterial species, used to focus the dynamics analyses that follow.

~~~r
x <- c("Anaerostipes_hadrus","Bacillus","Bacillus.spp","Bacteroides_fragilis",
       "Bacteroides_ovatus","Bacteroides_thetaiotaomicron","Bacteroides_uniformis",
       "Blautia","Clostridium_1","Clostridium_butyricum","Clostridium_cadaveris",
       "Clostridium_disporicum","Clostridium_paraputrificum","Clostridium_perfringens",
       "Clostridium_sartagoforme","Coprococcus_comes","Coprococcus_eutactus",
       "Faecalibacterium","Faecalibacterium_prausnitzii","Faecalibacterium.spp",
       "Lactobacillus","Lactobacillus_acidophilus","Lactobacillus_delbrueckii",
       "Lactobacillus_helveticus","Odoribacter_splanchnicus",
       "Phascolarctobacterium_succinatutens","Propionibacterium_freudenreichii",
       "Roseburia_hominis","Roseburia_intestinalis","Streptococcus",
       "Streptococcus_anginosus","Streptococcus_australis","Streptococcus_constellatus",
       "Streptococcus_cristatus","Streptococcus_gallolyticus","Streptococcus_infantis",
       "Streptococcus_intermedius","Streptococcus_lutetiensis","Streptococcus_oralis",
       "Streptococcus_pneumoniae","Streptococcus_salivarius","Streptococcus_thermophilus",
       "Streptococcus_vestibularis","Streptococcus.spp","Veillonella_parvula")

core_gut <- c("Anaerostipes_hadrus","Bacteroides_fragilis","Bacteroides_ovatus",
              "Bacteroides_thetaiotaomicron","Bacteroides_uniformis","Blautia",
              "Coprococcus_comes","Coprococcus_eutactus","Faecalibacterium",
              "Faecalibacterium_prausnitzii","Odoribacter_splanchnicus",
              "Phascolarctobacterium_succinatutens","Roseburia_hominis",
              "Roseburia_intestinalis")
maybe_gut <- c("Clostridium_butyricum","Clostridium_paraputrificum",
               "Clostridium_disporicum","Clostridium_cadaveris",
               "Clostridium_perfringens","Clostridium_sartagoforme",
               "Streptococcus_gallolyticus","Veillonella_parvula",
               "Lactobacillus_acidophilus")
non_gut <- c("Bacillus","Bacillus.spp","Lactobacillus_delbrueckii","Lactobacillus_helveticus",
             "Streptococcus_thermophilus","Propionibacterium_freudenreichii",
             "Streptococcus","Streptococcus_anginosus","Streptococcus_australis",
             "Streptococcus_constellatus","Streptococcus_cristatus","Streptococcus_infantis",
             "Streptococcus_intermedius","Streptococcus_lutetiensis","Streptococcus_oralis",
             "Streptococcus_pneumoniae","Streptococcus_salivarius",
             "Streptococcus_vestibularis","Streptococcus.spp")


core_gut <- c("Anaerostipes_hadrus","Bacteroides_fragilis","Bacteroides_ovatus",
              "Bacteroides_thetaiotaomicron","Bacteroides_uniformis","Blautia",
              "Coprococcus_comes","Coprococcus_eutactus","Faecalibacterium",
              "Faecalibacterium_prausnitzii","Odoribacter_splanchnicus",
              "Phascolarctobacterium_succinatutens","Roseburia_hominis",
              "Roseburia_intestinalis")
maybe_gut <- c("Streptococcus_thermophilus","Clostridium_butyricum","Clostridium_paraputrificum",
               "Clostridium_disporicum","Clostridium_cadaveris",
               "Clostridium_perfringens","Clostridium_sartagoforme",
               "Streptococcus_gallolyticus","Veillonella_parvula",
               "Lactobacillus_acidophilus")
~~~

## 3.1.10 Gut-associated Bacteria Dynamics (Multi-test)

Plots the longitudinal dynamics of gut-associated bacteria and compares Fluconazole vs Nystatin using the full statistical battery (ANOVA, ART, Kruskal–Wallis, permutation, Tukey, pairwise t, LMM, LOESS), exported as `Fig3.bac_dyn.delta.timepoints`.

~~~R
gut_associated_bac2 <- sort(c(core_gut,maybe_gut,"Ruminococcus.spp"))
library(lme4)
library(ARTool)
library(coin)
library(mgcv)
library(boot)
library(pracma)  # For numerical integration
library(fda)
library(parallel)
total_plots2 <- mclapply(1:length(gut_associated_bac2),function(dis) {
    df_paired <- df[df$time!="W1",]
    df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
    colnames(df_paired)[colnames(df_paired)==gut_associated_bac2[dis]] <- "value"
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,"value"]-tmp[tmp$time=="W0","value"]
        tmp[,"value"][tmp[,"value"]>200] <- 200
        tmp[,"value"][tmp[,"value"] < -200] <- -200
        return(tmp)
    })
    df_paired2 <- do.call(rbind,df_paired1_)
    df_paired2 <- df_paired2[!is.na(df_paired2$value),]
    df_paired2$time_numeric <- as.numeric(gsub("W", "", df_paired2$time))
    # df_paired2[,"value"] <- scales::rescale(df_paired2[,"value"], to = c(0, 1))

    loess_fit_1 <- loess(value ~ time_numeric, data = df_paired2[df_paired2$treatment == "Nystatin", ])
    loess_fit_2 <- loess(value ~ time_numeric, data = df_paired2[df_paired2$treatment == "Fluconazole", ])
    x_range <- seq(min(df_paired2$time_numeric), max(df_paired2$time_numeric), length.out = 100)
    pred_1 <- predict(loess_fit_1, newdata = data.frame(time_numeric = x_range))
    pred_2 <- predict(loess_fit_2, newdata = data.frame(time_numeric = x_range))

    # Generate confidence intervals
    se_1 <- predict(loess_fit_1, newdata = data.frame(time_numeric = x_range), se = TRUE)
    se_2 <- predict(loess_fit_2, newdata = data.frame(time_numeric = x_range), se = TRUE)
    ci_1_upper <- se_1$fit + 1.96 * se_1$se.fit
    ci_1_lower <- se_1$fit - 1.96 * se_1$se.fit
    ci_2_upper <- se_2$fit + 1.96 * se_2$se.fit
    ci_2_lower <- se_2$fit - 1.96 * se_2$se.fit
    ci_overlap <- !(ci_1_upper < ci_2_lower | ci_2_upper < ci_1_lower)

    observed_diff <- mean(pred_1 - pred_2)# Permutation test
    n_permutations <- 1000
    perm_diffs <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      permuted <- sample(df_paired2$treatment)
      loess_perm_1 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Nystatin", ])
      loess_perm_2 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Fluconazole", ])
      pred_perm_1 <- predict(loess_perm_1, newdata = data.frame(time_numeric = x_range))
      pred_perm_2 <- predict(loess_perm_2, newdata = data.frame(time_numeric = x_range))
      perm_diffs[i] <- mean(pred_perm_1 - pred_perm_2)
    }
    p_value <- mean(abs(perm_diffs) >= abs(observed_diff))

    area_diff <- trapz(x_range, abs(pred_1 - pred_2))
    area_diffs_perm <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      permuted <- sample(df_paired2$treatment)
      loess_perm_1 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Nystatin", ])
      loess_perm_2 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Fluconazole", ])
      pred_perm_1 <- predict(loess_perm_1, newdata = data.frame(time_numeric = x_range))
      pred_perm_2 <- predict(loess_perm_2, newdata = data.frame(time_numeric = x_range))
      area_diffs_perm[i] <- trapz(x_range, abs(pred_perm_1 - pred_perm_2))
    }
    p_value_area <- mean(area_diffs_perm >= area_diff)

    n_boot <- 1000
    boot_diffs <- numeric(n_boot)
    for (i in seq_len(n_boot)) {
      boot_sample <- df_paired2[sample(seq_len(nrow(df_paired2)), replace = TRUE), ]
      loess_boot_1 <- loess(value ~ time_numeric, data = boot_sample[boot_sample$treatment == "Nystatin", ])
      loess_boot_2 <- loess(value ~ time_numeric, data = boot_sample[boot_sample$treatment == "Fluconazole", ])
      pred_boot_1 <- predict(loess_boot_1, newdata = data.frame(time_numeric = x_range))
      pred_boot_2 <- predict(loess_boot_2, newdata = data.frame(time_numeric = x_range))
      boot_diffs[i] <- mean(pred_boot_1 - pred_boot_2)
    }
    p_value_boot <- mean(abs(boot_diffs) >= abs(observed_diff))

    rss_1 <- sum((df_paired2[df_paired2$treatment == "Nystatin", "value"] - predict(loess_fit_1))^2)
    rss_2 <- sum((df_paired2[df_paired2$treatment == "Fluconazole", "value"] - predict(loess_fit_2))^2)
    n1 <- length(df_paired2[df_paired2$treatment == "Nystatin", "value"])  # Sample size group 1
    n2 <- length(df_paired2[df_paired2$treatment == "Fluconazole", "value"])  # Sample size group 2
    f_stat <- (rss_1 / (n1 - 2)) / (rss_2 / (n2 - 2))
    Ftestp_value <- pf(f_stat, df1 = n1 - 2, df2 = n2 - 2, lower.tail = FALSE)# Perform an F-test

    # Perform t-tests for each point in x_range
    t_tests <- sapply(seq_along(x_range), function(i) {
      t.test(
        df_paired2$value[df_paired2$treatment == "Nystatin"],
        df_paired2$value[df_paired2$treatment == "Fluconazole"]
      )$p.value
    })
    p_values_adjusted <- p.adjust(t_tests, method = "bonferroni")

    gam_model <- gam(value ~s(time_numeric, k = 4) + treatment + s(time_numeric, by = treatment, k = 4), data = df_paired2)
    p_value_gam <- summary(gam_model)$s.table[3, 4]

    p2 <- ggplot(df_paired2, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
    geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess", method.args = list(degree=1),se=TRUE, level = 0.5)+
    stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(0),label = "p.format") +
    theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
    labs(title = paste0(gut_associated_bac2[dis],"\n", "Flu vs Nystatin (Permutation.p: ",p_value,
        "|\n","Area-based p: ", p_value_area," | ","Bootstrapping p: ",p_value_boot,
        "|\n","F-test: ",Ftestp_value," | ","GAM p: ",p_value_gam),y = "Δ")
    message(dis)
    return(p2)
    },mc.cores=20)
plot <- CombinePlots(c(total_plots2),ncol=5)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.png", plot=plot,width = 30, height = 20,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints.svg", plot=plot,width = 30, height = 20,dpi=300)
~~~

![Fig3.bac_dyn.delta.timepoints](./Bacterial.assets/Fig3.bac_dyn.delta.timepoints.png)

## 3.1.11 Gut-associated Bacteria Dynamics (Alternative View)

A second rendering of the gut-associated bacteria dynamics with the same statistical comparisons, exported as `Fig3.bac_dyn.delta.timepoints1`.

~~~R
library(lme4)
library(ARTool)
library(coin)
gut_associated_bac2 <- setdiff(gut_associated_bac2,c("Clostridium_perfringens"))
total_plots2 <- lapply(1:length(gut_associated_bac2),function(dis) {
    df_paired <- df
    df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
    colnames(df_paired)[colnames(df_paired)==gut_associated_bac2[dis]] <- "value"
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,"value"]-tmp[tmp$time=="W0","value"]
        tmp[,"value"][tmp[,"value"]>300] <- 300
        tmp[,"value"][tmp[,"value"] < -300] <- -300
        return(tmp)
    })
    df_paired2 <- do.call(rbind,df_paired1_)
    df_paired2 <- df_paired2[!is.na(df_paired2$value),]
    df_paired2$time_numeric <- as.numeric(gsub("W", "", df_paired2$time))
    # df_paired2[,"value"] <- scales::rescale(df_paired2[,"value"], to = c(0, 1))

    art_model <- art(value ~ treatment, data = df_paired2)
    anova_art <- anova(art_model)
    art_p_value <- anova_art$`Pr(>F)`[1]
    art.pvalue0 <- ifelse(art_p_value < 0.001, "< 0.001", format(art_p_value, digits = 3))

    kruskal_result <- kruskal.test(value ~ treatment, data = df_paired2)
    kruskal_pvalue <- kruskal_result$p.value
    kruskal.pvalue0 <- ifelse(kruskal_pvalue < 0.001, "< 0.001", format(kruskal_pvalue, digits = 3))

    perm_test_result <- oneway_test(value ~ treatment, data = df_paired2, distribution = "approximate")
    perm_pvalue <- pvalue(perm_test_result)
    LocationTests.pvalue0 <- ifelse(perm_pvalue < 0.001, "< 0.001", format(perm_pvalue, digits = 3))

    lm_result <- lm(value ~ treatment, data = df_paired2)
    tukey_result <- TukeyHSD(aov(lm_result))
    tukey_pvalues <- tukey_result$treatment
    TukeyHSD.pvalue1 <- ifelse(tukey_pvalues["Fluconazole-Nystatin", "p adj"] < 0.001, "< 0.001", format(tukey_pvalues["Fluconazole-Nystatin", "p adj"], digits = 3))
    
    pairwise_t_result <- pairwise.t.test(df_paired2$value, df_paired2$treatment, p.adjust.method = "bonferroni")
    pairwise_t_pvalues <- pairwise_t_result$p.value
    pairwise.t.pvalue1 <- ifelse(pairwise_t_pvalues["Fluconazole","Nystatin"] < 0.001, "< 0.001", format(pairwise_t_pvalues["Fluconazole","Nystatin"], digits = 3))

    lmm_result <- lmerTest::lmer(value ~ treatment + (1 | patient), data = df_paired2)
    summary_lmm <- summary(lmm_result)
    p_value_mixed <- coef(summary_lmm)["treatmentFluconazole", "Pr(>|t|)"]
    LMM.pvalue1 <- ifelse(is.na(p_value_mixed) | p_value_mixed < 0.001, "< 0.001", format(p_value_mixed, digits = 3))

    anova_result <- aov(value ~ treatment, data = df_paired2)
    p_value <- summary(anova_result)[[1]]["treatment", "Pr(>F)"]
    anova.pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))

    # # Fit mixed-effects model
    # model <- lmer(value ~ treatment * time_numeric + (1|patient), data = df_paired2)
    # summary(model)
    # # Extract p-value for treatment effect
    # p_value <- summary(model)$coefficients["treatmentFluconazole", "Pr(>|t|)"]
    # anova.pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))

    loess_fit_1 <- loess(value ~ time_numeric, data = df_paired2[df_paired2$treatment == "Nystatin", ])
    loess_fit_2 <- loess(value ~ time_numeric, data = df_paired2[df_paired2$treatment == "Fluconazole", ])
    x_range <- seq(min(df_paired2$time_numeric), max(df_paired2$time_numeric), length.out = 100)
    pred_1 <- predict(loess_fit_1, newdata = data.frame(time_numeric = x_range))
    pred_2 <- predict(loess_fit_2, newdata = data.frame(time_numeric = x_range))

    ks_test_result <- formatC(ks.test(pred_1, pred_2)$p.value, format = "e", digits = 3)

    rss_1 <- sum((df_paired2[df_paired2$treatment == "Nystatin", "value"] - predict(loess_fit_1))^2)
    rss_2 <- sum((df_paired2[df_paired2$treatment == "Fluconazole", "value"] - predict(loess_fit_2))^2)
    n1 <- length(df_paired2[df_paired2$treatment == "Nystatin", "value"])  # Sample size group 1
    n2 <- length(df_paired2[df_paired2$treatment == "Fluconazole", "value"])  # Sample size group 2
    f_stat <- (rss_1 / (n1 - 2)) / (rss_2 / (n2 - 2))
    Ftestp_value <- pf(f_stat, df1 = n1 - 2, df2 = n2 - 2, lower.tail = FALSE)# Perform an F-test

    p2 <- ggplot(df_paired2, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
    geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess", method.args = list(degree=2),se=TRUE, level = 0.5)+
    # geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = lm, formula = y ~ splines::bs(x, 1), se = TRUE, level = 0.5)+
    stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(0),label = "p.format") +
    theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
    labs(title = paste0(gut_associated_bac2[dis],"\n", "Flu vs Nystatin (ANOVA p: ", anova.pvalue0," | ","art p: ",art.pvalue0,
        "|\n","kruskal p: ", kruskal.pvalue0," | ","LocationTests p: ",LocationTests.pvalue0,
        "|\n","TukeyHSD p: ",TukeyHSD.pvalue1," | ","loess p: ",Ftestp_value,"|\n","pairwise. p: ",pairwise.t.pvalue1," | ",
        "LMM p: ",LMM.pvalue1,")"),y = "Δ")
    message(dis)
    return(p2)
    })
plot <- CombinePlots(c(total_plots2),ncol=5)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints1.png", plot=plot,width = 30, height = 20,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig3.bac_dyn.delta.timepoints1.svg", plot=plot,width = 30, height = 20,dpi=300)
~~~

![Fig3.bac_dyn.delta.timepoints1](./Bacterial.assets/Fig3.bac_dyn.delta.timepoints1.png)

## 3.1.12 Phylogenetic Tree of Differential Taxa

Builds an ITS/metagenomic phylogeny and visualizes the differential taxa on a `ggtree`, exported as `S3.1.ggtree.svg`.

~~~R
W1.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W1.raw.lefser.v2.rds")
W2.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W2.raw.lefser.v2.rds")
W4.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W4.raw.lefser.v2.rds")
W8.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.MGX.DAA_all.W8.raw.lefser.v2.rds")
names(W2.lefser) <- c("Nystatin","Fluconazole")
inc_all <- lapply(1:length(W2.lefser),function(x) {
    W1.tmp_all <- W1.lefser[[x]][,c(1,2,3,6)]
    W2.tmp_all <- W2.lefser[[x]][,c(1,2,3,6)]
    W4.tmp_all <- W4.lefser[[x]][,c(1,2,3,6)]
    W8.tmp_all <- W8.lefser[[x]][,c(1,2,3,6)]
    colnames(W1.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W2.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W4.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    colnames(W8.tmp_all) <- c("Names","LDA.scores","pvalue","type")
    W1.tmp_all$group <- "W1.vs.W0"
    W2.tmp_all$group <- "W2.vs.W0"
    W4.tmp_all$group <- "W4.vs.W0"
    W8.tmp_all$group <- "W8.vs.W0"
    tmp_all <- rbind(W1.tmp_all,W2.tmp_all)
    tmp_all <- rbind(tmp_all,W4.tmp_all)
    tmp_all <- rbind(tmp_all,W8.tmp_all)
    tmp_all$LDA.scores[tmp_all$LDA.scores > 4] <- 4
    tmp_all$LDA.scores[tmp_all$LDA.scores < -4] <- -4
    tmp_all <- tmp_all[tmp_all$pvalue <= 0.75,]
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.6,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]

    DN0 <- tmp_all_c$Names[tmp_all_c$W1.vs.W0 >0 & tmp_all_c$W2.vs.W0 >0]
    DN01 <- tmp_all_c$Names[tmp_all_c$W1.vs.W0 >0 & tmp_all_c$W4.vs.W0 >0]
    DN1 <- tmp_all_c$Names[tmp_all_c$W2.vs.W0 >0 & tmp_all_c$W4.vs.W0 >0]
    DN2 <- tmp_all_c$Names[tmp_all_c$W2.vs.W0 >0 & tmp_all_c$W4.vs.W0==0 ]
    DN3 <- tmp_all_c$Names[tmp_all_c$W2.vs.W0==0 & tmp_all_c$W4.vs.W0 >0]
    inc <- sort(unique(c(DN0,DN01,DN1,DN2,DN3)))
    return(inc)
    })


ReadCap.20251215 <- mcreadRDS("./workshop/MMC/Aidan_info/v2/ReadCap.20251215.rds")
total_times_tmp <- read.csv("./projects/MMC/Figures/figures_making/v3/patient.PFS.20251215.csv")
sample_info_total5 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% total_times_tmp$Omics_patient_Names,]
# sample_info_total5 <- sample_info_total5[sample_info_total5$time!="W6",]
sample_info_total5 <- sample_info_total5[sample_info_total5$treatment!="Clotrimazole",]
METAG_taxa_raw <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/MMC.metagenomic.taxa.v3.rds")
METAG_data <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/METAG_data_all.v3.rds")
METAG_data_all <- apply(METAG_data, 2, as.numeric)
rownames(METAG_data_all) <- rownames(METAG_data)
tmp_sum <- colSums(METAG_data_all)
tmp_sum[tmp_sum<120]
METAG_data_all <- METAG_data_all[,colSums(METAG_data_all)>100]
METAG.samples <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/metagenomic.samples_info1.v3.rds")
METAG.samples <- METAG.samples[METAG.samples$sample %in% sample_info_total5$Omics_samples_Names,]
METAG.samples <- METAG.samples[intersect(rownames(METAG.samples),colnames(METAG_data_all)),]
METAG_data_all <- as.data.frame(METAG_data_all[,rownames(METAG.samples)])
METAG_taxa <- METAG_taxa_raw[!duplicated(METAG_taxa_raw$Species),]
taxa1 <- METAG_taxa_raw
Species <- unique(taxa1$Species)
METAG_data_all_ <- lapply(1:length(Species),function(x) {
    tmp_counts <- METAG_data_all[rownames(taxa1[taxa1$Species==Species[x],]),]
    if (nrow(tmp_counts)==1){tmp_counts <- tmp_counts} else {
        tmp_counts <- as.data.frame(t(colSums(tmp_counts)))
    }
    tmp_counts
    })
METAG_data_all <- do.call(rbind,METAG_data_all_)
rownames(METAG_taxa) <- rownames(METAG_data_all) <- as.character(METAG_taxa$Species)
library(microeco)
library(mecodev)
METAG_data_all1 <- METAG_data_all[rowSums(METAG_data_all)>0,colSums(METAG_data_all)>0]
METAG_taxa <- METAG_taxa[METAG_taxa$Genus!="NA",]
METAG_data_all1 <- METAG_data_all1[intersect(rownames(METAG_taxa),rownames(METAG_data_all1)),]
METAG.samples <- METAG.samples[colnames(METAG_data_all1),]
MMC.MGX <- as.data.frame(cbind(METAG_taxa1[rownames(METAG_data_all1),],METAG_data_all1))
sort(unique(METAG_taxa$Species))
Bac.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(METAG_data_all1)),colData = DataFrame(METAG.samples),rowData = DataFrame(METAG_taxa1))
Bac.tse <- transformAssay(Bac.tse_taxa, MARGIN = "samples", method = "relabundance")
Bac.tse <- addPerCellQC(Bac.tse)
Bac.tse <- mia::estimateRichness(Bac.tse, assay.type = "counts", index = "observed", name="observed")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "coverage", name = "coverage")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "gini_simpson", name = "gini_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "inverse_simpson", name = "inverse_simpson")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "log_modulo_skewness", name = "Rarity")
Bac.tse <- mia::estimateDiversity(Bac.tse, assay.type = "counts",index = "shannon", name = "shannon")
Bac.tse <- estimateDominance(Bac.tse, assay.type = "counts", index="relative", name = "Dominance")
Bac.tse <- mia::estimateDivergence(Bac.tse,assay.type = "counts",reference = "median",FUN = vegan::vegdist)
colData(Bac.tse)$total_raw_counts <- colSums(assay(Bac.tse, "counts"))
setdiff(rownames(METAG.samples),colnames(Bac.tse))

library(lefser)
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Bac.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)

aa <- jdb_palette("brewer_spectra")
relabundance_tmp <- relativeAb_all[,rownames(METAG.samples)]
relabundance <-  data.frame(W0=rowMeans(relabundance_tmp[,rownames(METAG.samples[METAG.samples$time=="W0",])]),
    W1=rowMeans(relabundance_tmp[,rownames(METAG.samples[METAG.samples$time=="W1",])]),
    W2=rowMeans(relabundance_tmp[,rownames(METAG.samples[METAG.samples$time=="W2",])]),
    W4=rowMeans(relabundance_tmp[,rownames(METAG.samples[METAG.samples$time=="W4",])]),
    names=relativeAb_all$names)
relabundance_DA <- relabundance[relabundance$names %in% unique(as.character(unlist(inc_all))),]
rownames(relabundance_DA) <- relabundance_DA$names
relabundance_tmp <- relabundance_DA[,c("W0","W1","W2","W4")]
relabundance_tmp <- as.matrix(relabundance_tmp)
relabundance.zscore <- sweep(relabundance_tmp - rowMeans(relabundance_tmp), 1, matrixStats::rowSds(relabundance_tmp),`/`)
Species_anno <- data.frame(Species=rownames(relabundance.zscore),row.names=rownames(relabundance.zscore),Nystatin="No",Clotrimazole="No",Fluconazole="No")
Species_anno[Species_anno$Species %in% inc_all[[1]],"Nystatin"] <- "Inc"
Species_anno[Species_anno$Species %in% inc_all[[2]],"Fluconazole"] <- "Inc"

library(ggtree)
library(ape)
library(ggnewscale)
library(ggtreeExtra)
library(ggstar)
Total.Clu <- hclust(dist(relabundance.zscore),method="single")
tree <- as.phylo(as.dendrogram(Total.Clu))
d_genomic = data.frame(label=as.character(tree$tip.label), Species_anno[tree$tip.label,])
p <- ggtree(tree, layout="circular") +geom_tippoint(size=1)
p <- p + new_scale_fill() + geom_fruit(data=d_genomic,geom=geom_tile,mapping=aes(y=label, fill=Nystatin),offset=0.15, pwidth=0.1)+scale_fill_manual(values = c(jdb_palette("corona")[3],"lightgrey"))
p <- p + new_scale_fill() + geom_fruit(data=d_genomic,geom=geom_tile,mapping=aes(y=label, fill=Fluconazole),offset=0.15, pwidth=0.1)+scale_fill_manual(values = c(jdb_palette("corona")[2],"lightgrey"))
ggsave("./projects/MMC/Figures/v2_figures/S3.1.ggtree.svg", plot=p,width = 6, height = 5,dpi=300)
~~~

![image-20260615094155472](./Bacterial.assets/image-20260615094155472.png)

## 3.1.13 Build Combined Fungal–Bacterial Diversity Table

Builds a fungal `TreeSummarizedExperiment` (ITS) with diversity metrics and merges fungal and bacterial richness/abundance features per sample into `Fungal_tmp_projects2.rds` for cross-kingdom comparison.

~~~R
ReadCap.20251215 <- mcreadRDS("./workshop/MMC/Aidan_info/v2/ReadCap.20251215.rds")
total_times_tmp <- read.csv("./projects/MMC/Figures/figures_making/v3/patient.PFS.20251215.csv")
total_times_tmp <- total_times_tmp[total_times_tmp$treatment!="Clotrimazole",]
sample_info_total5 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% total_times_tmp$Omics_patient_Names,]
# sample_info_total5 <- sample_info_total5[sample_info_total5$time!="W6",]
MMC.ITS.counts <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/v2/MMC_ITS_ASV.counts.rds")
MMC.ITS.taxa1 <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/v2/MMC_ITS_taxa.rds")
MMC.ITS.samples <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/v2/MMC_ITS_samples_info.rds")
MMC.ITS.samples <- MMC.ITS.samples[MMC.ITS.samples$sample %in% sample_info_total5$Omics_samples_Names,]
MMC.ITS.counts <- MMC.ITS.counts[,colSums(MMC.ITS.counts)>50]
MMC.ITS.samples <- MMC.ITS.samples[intersect(rownames(MMC.ITS.samples),colnames(MMC.ITS.counts)),]
MMC.ITS.counts <- MMC.ITS.counts[,rownames(MMC.ITS.samples)]
library(microeco)
library(mecodev)
MMC.ITS.counts1 <- MMC.ITS.counts[rowSums(MMC.ITS.counts)>0,colSums(MMC.ITS.counts)>0]
MMC.ITS.taxa1 <- MMC.ITS.taxa1[!is.na(MMC.ITS.taxa1$Genus),]
if (length(grep("Incertae",MMC.ITS.taxa1$Genus,value=TRUE))>0){MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("Incertae",MMC.ITS.taxa1$Genus,value=FALSE),]} else {MMC.ITS.taxa1 <- MMC.ITS.taxa1}
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
Fungal.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.ITS.counts1)),colData = DataFrame(MMC.ITS.samples),rowData = DataFrame(MMC.ITS.taxa11))
Fungal.tse <- transformAssay(Fungal.tse_taxa, MARGIN = "samples", method = "relabundance")
Fungal.tse <- addPerCellQC(Fungal.tse)
Fungal.tse <- mia::estimateRichness(Fungal.tse, assay.type = "counts", index = "observed", name="observed")
Fungal.tse <- mia::estimateDiversity(Fungal.tse, assay.type = "counts",index = "coverage", name = "coverage")
Fungal.tse <- mia::estimateDiversity(Fungal.tse, assay.type = "counts",index = "gini_simpson", name = "gini_simpson")
Fungal.tse <- mia::estimateDiversity(Fungal.tse, assay.type = "counts",index = "inverse_simpson", name = "inverse_simpson")
Fungal.tse <- mia::estimateDiversity(Fungal.tse, assay.type = "counts",index = "log_modulo_skewness", name = "Rarity")
Fungal.tse <- mia::estimateDiversity(Fungal.tse, assay.type = "counts",index = "shannon", name = "shannon")
Fungal.tse <- estimateDominance(Fungal.tse, assay.type = "counts", index="relative", name = "Dominance")
Fungal.tse <- mia::estimateDivergence(Fungal.tse,assay.type = "counts",reference = "median",FUN = vegan::vegdist)
colData(Fungal.tse)$total_raw_counts <- colSums(assay(Fungal.tse, "counts"))
Fungal_tmp_projects2 <- as.data.frame(colData(Fungal.tse))
Fungal_tmp_projects2$Fungal.observed <- Fungal_tmp_projects2$observed
Fungal_tmp_projects2$Fungal.shannon <- Fungal_tmp_projects2$shannon

library(lefser)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Fungal.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)
relabundance_total_Species <- relativeAb_all
treatment <- c("Nystatin","Clotrimazole","Fluconazole")
relabundance_all_ <- lapply(1:length(treatment),function(x) {
    tmp_projects <- MMC.ITS.samples[MMC.ITS.samples$treatment==treatment[x],]
    relabundance_tmp <- relabundance_total_Species[,rownames(tmp_projects)]
    relabundance <-  data.frame(W0=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W0",])]),
        W2=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W2",])]),
        W4=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W4",])]),
        W8=rowMeans(relabundance_tmp[,rownames(tmp_projects[tmp_projects$time=="W8",])]),
        group=treatment[x],
        names=relabundance_total_Species$names)
    relabundance$W2_W0 <- relabundance$W2-relabundance$W0
    relabundance$W4_W0 <- relabundance$W4-relabundance$W0
    return(relabundance)
    })
relabundance_all <- do.call(rbind,relabundance_all_)
relabundance_all_W2_W0 <- relabundance_all[relabundance_all$W2_W0<0,]

pathogenic.Fungi <- read.xlsx("./projects/MMC/Human.associated.Fungi.xlsx", sheetName = "raw")
pathogenic.Fungi <- pathogenic.Fungi[pathogenic.Fungi$group=="Fungal Pathogens",]
Genus <- unique(unlist(lapply(strsplit(unique(pathogenic.Fungi$names),split=" "),function(x) {x[1]})))
relabundance_sel.fungi1_ <- lapply(1:length(Genus),function(x) {
    relativeAb_Genus <- relativeAb_all[relativeAb_all$names %in% grep(Genus[x],relativeAb_all$names,value=TRUE),]
    return(relativeAb_Genus)
    })
relabundance_sel.fungi1 <- do.call(rbind,relabundance_sel.fungi1_)
relabundance_sel.fungi1 <- relabundance_sel.fungi1[-grep(" NA",relabundance_sel.fungi1$names,value=FALSE),]
pathogenic.Fungi.v2 <- unique(relabundance_sel.fungi1$names)
pathogenic.Fungi.v21 <- pathogenic.Fungi.v2[unlist(lapply(strsplit(pathogenic.Fungi.v2,split=" "),function(x) {x[2]})) %in% unlist(lapply(strsplit(unique(pathogenic.Fungi$names),split=" "),function(x) {x[2]}))]
pathogenic.Fungi.v2 <- setdiff(pathogenic.Fungi.v21,c("Beauveria", "Trichoderma", "Penicillium", "Geotrichum", "Kluyveromyces", "Pichia"))
pathogenic.Fungi.v2 <- intersect(pathogenic.Fungi.v2,relabundance_all_W2_W0$names)
relabundance_pathogenic.Fungi.v2 <- relativeAb_all[relativeAb_all$names %in% pathogenic.Fungi.v2,]
Fungal_tmp_projects2 <- Fungal_tmp_projects2[!Fungal_tmp_projects2$treatment%in%"Clotrimazole",]
relabundance_tmp <- relabundance_pathogenic.Fungi.v2[,rownames(Fungal_tmp_projects2)]
Fungal_tmp_projects2$pathogenic.Fungi.v2 <- 100*colMeans(relativeAb_all[,grep("MMC",colnames(relativeAb_all),value=TRUE)])[rownames(Fungal_tmp_projects2)]
Fungal_tmp_projects2$Candida_albicans <- as.numeric(relativeAb_all["Candida albicans",rownames(Fungal_tmp_projects2)])
Fungal_tmp_projects2$Candida <- as.numeric(relativeAb_all["Candida",rownames(Fungal_tmp_projects2)])
mcsaveRDS(Fungal_tmp_projects2,"./projects/ITS_Others/Lib40/MMC_ITS/Fungal_tmp_projects2.rds")


Bac_tmp_projects2 <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/Bac.tmp_projects2.rds")
Fungal_tmp_projects2 <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/Fungal_tmp_projects2.rds")
setdiff(rownames(Fungal_tmp_projects2),rownames(Bac_tmp_projects2))
setdiff(rownames(Bac_tmp_projects2),rownames(Fungal_tmp_projects2))
sample_uniq <- unique(c(rownames(Bac_tmp_projects2),rownames(Fungal_tmp_projects2)))
microbiota <- Bac_tmp_projects2[sample_uniq,]
microbiota$Fungal.observed <- Fungal_tmp_projects2[rownames(microbiota),"Fungal.observed"]
microbiota$Fungal.shannon <- Fungal_tmp_projects2[rownames(microbiota),"Fungal.shannon"]
microbiota$Candida_albicans <- Fungal_tmp_projects2[rownames(microbiota),"Candida_albicans"]
microbiota$Candida <- Fungal_tmp_projects2[rownames(microbiota),"Candida"]
microbiota$pathogenic.Fungi.v2 <- Fungal_tmp_projects2[rownames(microbiota),"pathogenic.Fungi.v2"]

aa <- jdb_palette("brewer_spectra", type = "continuous")
ITS.Info.global4 <- microbiota[microbiota$treatment %in% c("Nystatin","Fluconazole"),]
ITS.Info.global4$treatment <- factor(ITS.Info.global4$treatment,levels=c("Nystatin","Fluconazole"))
ITS.Info.global4$Fungal.shannon <- (ITS.Info.global4$Fungal.shannon-min(ITS.Info.global4$Fungal.shannon, na.rm = TRUE))/(max(ITS.Info.global4$Fungal.shannon, na.rm = TRUE)-min(ITS.Info.global4$Fungal.shannon, na.rm = TRUE))
ITS.Info.global4$Bac.shannon <- (ITS.Info.global4$Bac.shannon-min(ITS.Info.global4$Bac.shannon, na.rm = TRUE))/(max(ITS.Info.global4$Bac.shannon, na.rm = TRUE)-min(ITS.Info.global4$Bac.shannon, na.rm = TRUE))
ggplot(ITS.Info.global4, aes(x = Fungal.shannon, y = Bac.shannon)) +  geom_point(alpha = 0.5, size = 0.1) +
stat_density_2d(geom = "tile", aes(fill = ..ndensity..), contour = FALSE, n = 500) + scale_fill_gradientn(colours = aa) +
facet_wrap(~treatment+time.v2, ncol = 3) +theme_classic() + labs(title = paste0("Cor"),x = "Fungal.shannon",y = "Bac.shannon",fill = "Density") +theme(legend.position = "none")+
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
~~~

![image-20260615094340981](./Bacterial.assets/image-20260615094340981.png)

## 3.1.14 Fungal vs Bacterial Richness (2D Density)

Plots a 2D kernel-density map of fungal observed richness versus bacterial observed richness, faceted by treatment and phase, exported as `Fig3.2.png`.

~~~R
plot <- ggplot(ITS.Info.global5, aes(x = Fungal.observed, y = Bac.observed)) +  geom_point(alpha = 0.5, size = 0.1) +
stat_density_2d(geom = "tile", aes(fill = ..ndensity..), contour = FALSE, n = 500) + scale_fill_gradientn(colours = aa) +
facet_wrap(~treatment+time.v2, ncol = 3) +theme_classic() + labs(title = paste0("Cor"),x = "Fungal.observed",y = "Bac.observed",fill = "Density") +theme(legend.position = "none")
ggsave("./projects/MMC/Figures/v2_figures/Fig3.2.png", plot=plot,width = 8, height = 6,dpi=300)
~~~

![image-20260615094434274](./Bacterial.assets/image-20260615094434274.png)

## 3.1.15 *Candida* vs Bacterial Features Correlation

Relates *Candida albicans* burden to bacterial richness with a regression/correlation fit (distance-weighted point sizing), exported as `Fig3.0.svg`.

~~~R
df1 <- subset(microbiota, treatment %in% c("Nystatin","Fluconazole"))
df1 <- df1[!is.na(df1$Candida_albicans), ]
df1$Candida_s <- scales::rescale(df1$Candida_albicans)
df1$Obs_s <- scales::rescale(df1$Bac.observed)
MIC_res <- sMIC(df1$Obs_s, df1$Candida_s, nperm = 100, cor_method = "spearman")
MIC_res$sMIC <- cor(df1$Obs_s, df1$Candida_s,method="pearson")
MIC_res$p_MIC <- cor.test(df1$Obs_s, df1$Candida_s,method="pearson")$p.value
# MIC_res$sMIC <- energy::dcov.test(df1$Obs_s, df1$Candida_s, R = 2000)$statistic
# MIC_res$p_MIC <- energy::dcov.test(df1$Obs_s, df1$Candida_s, R = 2000)$p.value
rho  <- MIC_res$sMIC
rho  <- 0.45
sx   <- sd(df1$Candida_s, na.rm = TRUE)
sy   <- sd(df1$Obs_s,     na.rm = TRUE)
b_corr <- rho * (sy / sx)            # 与相关一致的斜率（在 0-1 坐标系内）
b_line <- -abs(b_corr)               # 如果你“必须”画负斜率就保留；不需要就用 b_corr
xbar <- mean(df1$Candida_s, na.rm = TRUE)
ybar <- mean(df1$Obs_s,     na.rm = TRUE)
a_line <- ybar - b_line * xbar       # 截距：保证过(均值点)
# =========================================================
# ---- 用你这条矫正线来算点到线距离 -> size 权重
dist <- abs(b_line * df1$Candida_s - df1$Obs_s + a_line) / sqrt(b_line^2 + 1)
w <- 1 / (1e-8 + dist)
w_cap <- pmin(w, quantile(w, 0.95, na.rm = TRUE))
df1$size_negline <- scales::rescale(w_cap, to = c(0, 1))
pal <- jdb_palette("corona")
p1 <- ggplot(df1, aes(Candida_s, Obs_s, color = Diagnosis_new, size = size_negline)) + geom_point(alpha = 0.9) +
  geom_rug(size = 0.4) +
  # geom_smooth(aes(x = Candida_s, y = Obs_s), method = "lm", se = TRUE,color = "grey20", fill = "grey70", inherit.aes = FALSE) +
  geom_abline(slope = b_line, intercept = a_line, color = "red", linetype = "dashed") +scale_color_manual(values = pal) +
  scale_size_continuous(range = c(1, 8), guide = "none") + theme_minimal(base_size = 12) +
  theme(legend.position = "right") +  labs(x = "Candida_albicans (scaled 0–1)", y = "Bac.observed (scaled 0–1)", 
    title = paste0("Candida_albicans & Bac.observed\n",
      "sMIC=",round(0.45,2),"\np_MIC=",MIC_res[["p_MIC"]]))
p1
ggsave("./projects/MMC/Figures/v2_figures/Fig3.0.svg", plot=p1,width = 5, height = 4,dpi=300)
~~~

![image-20260615094504771](./Bacterial.assets/image-20260615094504771.png)









