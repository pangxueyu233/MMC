# 3.1 Main Bacterial Analysis: Data Processing, Visualization, and Statistical Testing

## 3.1.1. Overview of Bacterial ITS and Metagenomic Data Processing

The bacterial analysis incorporates **ITS sequencing and metagenomic profiling** from the **MMC cohort**, focusing on **treatment responses and gut microbial shifts**.

Key steps:

- **Load and clean data**: Filtering out low-abundance taxa and samples with missing data.
- **Normalize abundance**: Convert raw counts into **relative abundances**.
- **Diversity estimation**: Calculate **alpha-diversity metrics** such as **Shannon index, Simpson index, and species richness**.
- **Dominance and divergence analysis**: Investigate shifts in microbial community structure.

~~~R
MMC.ITS.counts <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.ITS.counts.v3.rds")
MMC.ITS.taxa1 <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.ITS.taxa1.v2.rds")
MMC.ITS.samples <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/ITS.samples_info1.v2.rds")
library(microeco)
library(mecodev)
MMC.ITS.counts1 <- MMC.ITS.counts[rowSums(MMC.ITS.counts)>0,colSums(MMC.ITS.counts)>0]
MMC.ITS.taxa1 <- MMC.ITS.taxa1[!is.na(MMC.ITS.taxa1$Genus),]
if (length(grep("Incertae",MMC.ITS.taxa1$Genus,value=TRUE))>0){MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("Incertae",MMC.ITS.taxa1$Genus,value=FALSE),]} else {MMC.ITS.taxa1 <- MMC.ITS.taxa1}
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
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
ITS.Info.global3 <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.All.features.Info.gloibal.rds")
tmp_projects2 <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.ITS.All.Info.gloibal.groups.rds")
tmp_projects2 <- tmp_projects2[,c("DAI","UC.score.v2","CD.score.raw","time","time.v2","treatment","patient","Diagnosis_new")]
tmp_projects2 <- as.data.frame(cbind(tmp_projects2,ITS.Info.global3[rownames(tmp_projects2),c("Bac.observed","Bac.shannon","Candida_albicans","Candida","Butyrate_producing_bac")]))

relabundance_total_Species <- relativeAb_all
aa <- jdb_palette("brewer_spectra")
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
library(xlsx)
pathogenic.Fungi <- read.xlsx("/mnt/d/xiangyu.ubuntu/projects/MMC/Human.associated.Fungi.xlsx", sheetName = "raw")
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
relabundance_tmp <- relabundance_pathogenic.Fungi.v2[,rownames(tmp_projects2)]
tmp_projects2$pathogenic.Fungi.v2 <- 100*colMeans(relabundance_tmp)[rownames(tmp_projects2)]
diease_Scores <- c("Bac.observed","Bac.shannon","Candida_albicans","Candida","pathogenic.Fungi.v2")
df1 <-tmp_projects2
df1$Candida_albicans <- log(df1$Candida_albicans+1,2)
df1$Butyrate_producing_bac <- log(df1$Butyrate_producing_bac+1,2)
sel_pa <- unique(df1$patient)
paired.df_ <- lapply(1:length(sel_pa),function(x) {
    tmp <- df1[df1$patient %in% sel_pa[x],]
    tmp[,"Butyrate_producing_bac"] <- tmp[,"Butyrate_producing_bac"]-tmp[tmp$time=="W0","Butyrate_producing_bac"]
    tmp[,"Bac.observed"] <- tmp[,"Bac.observed"]-tmp[tmp$time=="W0","Bac.observed"]
    tmp[,"Bac.shannon"] <- tmp[,"Bac.shannon"]-tmp[tmp$time=="W0","Bac.shannon"]
    tmp[,"Candida_albicans"] <- tmp[,"Candida_albicans"]-tmp[tmp$time=="W0","Candida_albicans"]
    tmp[,"Candida"] <- tmp[,"Candida"]-tmp[tmp$time=="W0","Candida"]
    tmp[,"pathogenic.Fungi.v2"] <- tmp[,"pathogenic.Fungi.v2"]-tmp[tmp$time=="W0","pathogenic.Fungi.v2"]
    return(tmp)
})
paired.df <- do.call(rbind,paired.df_)
paired.df$abs.Candida <- abs(paired.df$Candida)
paired.df$abs.Ca <- abs(paired.df$Candida_albicans)
paired.df$abs.Fungi.observed <- abs(paired.df$pathogenic.Fungi.v2)
paired.df$abs.Bac.observed <- abs(paired.df$Bac.observed)
paired.df$abs.Bac.shannon <- abs(paired.df$Bac.shannon)
pal <- jdb_palette("corona")
total.lm <- lm(Candida_albicans ~ Bac.shannon, data = df1)
p1 <- ggplot(df1, aes_string(x = "Candida_albicans",y = "Bac.shannon",  color = "Diagnosis_new",size="Candida_albicans"))+ geom_rug(size=0.5)+ geom_point(alpha = 0.8) +  
geom_smooth(method = "lm", se = TRUE, color="red") + scale_color_manual(values = pal)+  theme_minimal() + theme(legend.position = "right")+
scale_size_continuous(range = c(1, 6), guide = "none") + labs(x = "Candida_albicans", y = "Bac.shannon", title = paste0("Cor Change-Bac.shannon and Change-Candida_albicans")) +  
annotate("text", x =3, y = 1, color = "#E51718", label = paste0("rho = ", round(cor(total.lm$model$Bac.shannon, total.lm$model$Candida_albicans, method = "pearson"),2),
    "; p = ", round(cor.test(total.lm$model$Bac.shannon, total.lm$model$Candida_albicans, method = "pearson")$p.value,4)))
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/Fig3.0.svg", plot=p1,width = 5, height = 4,dpi=300)
~~~

![image-20250205130257604](./Chapter3.assets/image-20250205130257604.png)

## 3.1.2. Alpha and Beta Diversity Analysis

- **Figure 3.0**: Correlation between **Candida albicans abundance and Shannon index of bacterial diversity**, stratified by disease status.
- **Figure 3.1**: **PCoA analysis** using **Bray-Curtis dissimilarity**, demonstrating temporal and treatment-based clustering of metagenomic samples.

~~~R
MMC.metagenomic.counts <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/metagenomic_data1.v2.rds")
MMC.metagenomic.taxa1 <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.metagenomic.taxa.v2.rds")
MMC.metagenomic.samples <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/metagenomic.samples_info1.v2.rds")
taxa1 <- MMC.metagenomic.taxa1
Species <- unique(taxa1$Species)
MMC.metagenomic.counts_ <- lapply(1:length(Species),function(x) {
    tmp_counts <- MMC.metagenomic.counts[rownames(taxa1[taxa1$Species==Species[x],]),]
    if (nrow(tmp_counts)==1){tmp_counts <- tmp_counts} else {
        tmp_counts <- as.data.frame(t(colSums(tmp_counts)))
    }
    tmp_counts
    })
MMC.metagenomic.counts <- do.call(rbind,MMC.metagenomic.counts_)
MMC.metagenomic.taxa1 <- MMC.metagenomic.taxa1[!duplicated(MMC.metagenomic.taxa1$Species),]
rownames(MMC.metagenomic.counts) <- rownames(MMC.metagenomic.taxa1) <- MMC.metagenomic.taxa1$Species
MMC.metagenomic.taxa1[grep("_1",MMC.metagenomic.taxa1$Species,value=FALSE),]
MMC.metagenomic.taxa1[grep("_NA",MMC.metagenomic.taxa1$Species,value=FALSE),]
MMC.metagenomic.taxa1[grep("Streptococcus_thermophilus",MMC.metagenomic.taxa1$Species,value=FALSE),]
MMC.metagenomic.taxa1[grep("Leuconostoc.spp",MMC.metagenomic.taxa1$Species,value=FALSE),]
library(microeco)
library(mecodev)
MMC.metagenomic.taxa1 <- MMC.metagenomic.taxa1[MMC.metagenomic.taxa1$Genus!="NA",]
MMC.metagenomic.counts1 <- MMC.metagenomic.counts[intersect(rownames(MMC.metagenomic.taxa1),rownames(MMC.metagenomic.counts)),]
MMC.metagenomic.counts1 <- MMC.metagenomic.counts1[rowSums(MMC.metagenomic.counts1)>0,colSums(MMC.metagenomic.counts1)>0]
MMC.metagenomic.taxa11 <- MMC.metagenomic.taxa1[rownames(MMC.metagenomic.counts1),]
MMC.metagenomic.samples <- MMC.metagenomic.samples[colnames(MMC.metagenomic.counts1),]
MMC.metagenomic.samples$time <- factor(MMC.metagenomic.samples$time,levels=c("W0","W1","W2","W4","W8"))
MMC.metagenomic.samples$group.v2 <- factor(MMC.metagenomic.samples$group.v2,levels=c("Resp","NR"))

tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.metagenomic.counts1)),colData = DataFrame(MMC.metagenomic.samples),rowData = DataFrame(MMC.metagenomic.taxa11))
METAG_tse <- transformAssay(tse_taxa, MARGIN = "samples", method = "relabundance")
METAG_tse <- addPerCellQC(METAG_tse)
METAG_tse <- mia::estimateRichness(METAG_tse, assay.type = "counts", index = "observed", name="observed")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "coverage", name = "coverage")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "gini_simpson", name = "gini_simpson")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "inverse_simpson", name = "inverse_simpson")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "log_modulo_skewness", name = "Rarity")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "shannon", name = "shannon")
METAG_tse <- estimateDominance(METAG_tse, assay.type = "counts", index="relative", name = "Dominance")
METAG_tse <- mia::estimateDivergence(METAG_tse,assay.type = "counts",reference = "median",FUN = vegan::vegdist)
colData(METAG_tse)$total_raw_counts <- colSums(assay(METAG_tse, "counts"))
METAG_tse$time.v2 <- as.character(METAG_tse$time)
METAG_tse$time.v2[METAG_tse$time.v2 %in% c("W0")] <- "Pre"
METAG_tse$time.v2[METAG_tse$time.v2 %in% c("W1","W2","W4")] <- "Post"
METAG_tse$time.v2[METAG_tse$time.v2 %in% c("W8")] <- "LTM"
METAG_tse$time.v2 <- factor(METAG_tse$time.v2,levels=c("Pre","Post","LTM"))
library(ggside)
generate_ellipse <- function(median_x, mad_x, median_y, mad_y, n = 500, factor) {
  t <- seq(0, 2*pi, length.out = n)
  a <- mad_x * factor  # Semi-major axis length
  b <- mad_y * factor  # Semi-minor axis length
  x <- median_x + a * cos(t)
  y <- median_y + b * sin(t)
  return(data.frame(x = x, y = y))
}
treatment <- c("Nystatin", "Clotrimazole", "Fluconazole")
pal <- jdb_palette("corona")
ALL_PLOTS <- lapply(treatment, function(treat) {
    tse_tmp <- METAG_tse[,which(METAG_tse$treatment %in% c(treat))]
    tse_tmp <- tse_tmp[,which(tse_tmp$time %in%  c("W0","W1","W2","W4","W8"))]
    tse_tmp <- subsetByPrevalentFeatures(tse_tmp, rank = "Species", detection = 0.01, prevalence = 0.01)
    tse_tmp <- transformAssay(tse_tmp, assay.type = "counts", method = "relabundance")
    tse_tmp <- runMDS(tse_tmp, FUN = vegan::vegdist, method = "bray", assay.type = "relabundance", name = "MDS_bray")
    p <- plotReducedDim(tse_tmp, "MDS_bray", colour_by = "time.v2", point_size = 5)
    PcoA <- ggplot_build(p)$data[[1]]
    PcoA$time.v2 <- unique(tse_tmp$time.v2)[PcoA$group]
    sbg <- PcoA %>% dplyr::group_by(time.v2) %>% dplyr::summarise(median_x = median(x),mad_x = mad(x),median_y = median(y),mad_y = mad(y),.groups = 'drop')
    ellipse_data <- do.call(rbind, lapply(1:nrow(sbg), function(i) {
        df <- generate_ellipse(sbg$median_x[i], sbg$mad_x[i], sbg$median_y[i], sbg$mad_y[i], factor = 0.8)
        df$time.v2 <- sbg$time.v2[i]
        df
    }))
    p <- ggplot(PcoA, aes(x = x, y = y, color = time.v2)) +
    geom_point(aes(shape = time.v2, fill = time.v2), size = 5, alpha = 0.5) +
    geom_polygon(data = ellipse_data, aes(x = x, y = y, fill = time.v2), alpha = 0.2) +
    geom_errorbar(data = sbg, aes(x = median_x, y = median_y, ymin = median_y - 0.5 * mad_y, ymax = median_y + 0.5 * mad_y), width = 0) +
    geom_errorbarh(data = sbg, aes(x = median_x, y = median_y, xmin = median_x - 0.5 * mad_x, xmax = median_x + 0.5 * mad_x), height = 0) +
    geom_point(data = sbg, aes(x = median_x, y = median_y, fill = time.v2), color = "black", shape = 22, size = 5, alpha = 0.7, show.legend = FALSE) +
    theme_bw() + scale_fill_manual(values = pal) + scale_color_manual(values = pal) + labs(x = "PCoA Axis 1", y = "PCoA Axis 2", title = paste0(treat, " in METAG.PCoA")) +
    geom_xsideboxplot(orientation = "y") + geom_ysideboxplot(orientation = "x")
    message(treat)
    return(p)
})
plot <- CombinePlots(ALL_PLOTS,nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/Fig3.1.svg", plot=plot,width = 15, height = 4,dpi=300)
~~~

![image-20250205130501684](./Chapter3.assets/image-20250205130501684.png)

## 3.1.3. Longitudinal Changes in Bacterial and Fungal Composition

- **Figure 3.2**: Density plots depicting **fungal vs. bacterial diversity across treatments and time**.
- **Figure 3.3**: Paired comparisons of **bacterial richness, diversity, and short-chain fatty acid (SCFA)-producing bacteria** over time.

~~~R
ITS.Info.global3 <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.All.features.Info.gloibal.rds")
tmp_projects2 <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.ITS.All.Info.gloibal.groups.rds")
tmp_projects2 <- tmp_projects2[,c("DAI","UC.score.v2","CD.score.raw","time","time.v2","treatment","patient")]
tmp_projects2 <- as.data.frame(cbind(tmp_projects2,ITS.Info.global3[rownames(tmp_projects2),c("Bac.observed","Fungi.observed","Bac.shannon","Propionate_producing_bac","Butyrate_producing_bac","Valeric_producing_bac","probiotic.bac")]))
aa <- jdb_palette("brewer_spectra", type = "continuous")
ITS.Info.global4 <- tmp_projects2[tmp_projects2$treatment %in% c("Nystatin","Fluconazole"),]
ITS.Info.global4$treatment <- factor(ITS.Info.global4$treatment,levels=c("Nystatin","Fluconazole"))
ITS.Info.global4$Fungi.observed <- (ITS.Info.global4$Fungi.observed-min(ITS.Info.global4$Fungi.observed, na.rm = TRUE))/(max(ITS.Info.global4$Fungi.observed, na.rm = TRUE)-min(ITS.Info.global4$Fungi.observed, na.rm = TRUE))
ITS.Info.global4$Bac.shannon <- (ITS.Info.global4$Bac.shannon-min(ITS.Info.global4$Bac.shannon, na.rm = TRUE))/(max(ITS.Info.global4$Bac.shannon, na.rm = TRUE)-min(ITS.Info.global4$Bac.shannon, na.rm = TRUE))
ggplot(ITS.Info.global4, aes(x = Fungi.observed, y = Bac.shannon)) +  geom_point(alpha = 0.5, size = 0.1) +
stat_density_2d(geom = "tile", aes(fill = ..ndensity..), contour = FALSE, n = 500) + scale_fill_gradientn(colours = aa) +
facet_wrap(~treatment+time.v2, ncol = 3) +theme_classic() + labs(title = paste0("Cor"),x = "Fungi.observed",y = "Bac.shannon",fill = "Density") +theme(legend.position = "none")+
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
plot <- ggplot(ITS.Info.global4, aes(x = Fungi.observed, y = Bac.shannon)) +  geom_point(alpha = 0.5, size = 0.1) +
stat_density_2d(geom = "tile", aes(fill = ..ndensity..), contour = FALSE, n = 500) + scale_fill_gradientn(colours = aa) +
facet_wrap(~treatment+time.v2, ncol = 3) +theme_classic() + labs(title = paste0("Cor"),x = "Fungi.observed",y = "Bac.shannon",fill = "Density") +theme(legend.position = "none")
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/Fig3.2.png", plot=plot,width = 8, height = 6,dpi=300)
~~~

![image-20250205130555420](./Chapter3.assets/image-20250205130555420.png)

~~~R
diease_Scores <- c("Bac.observed","Bac.shannon","Propionate_producing_bac","Butyrate_producing_bac","Valeric_producing_bac","probiotic.bac")
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
pal <- jdb_palette("corona")[c(1,3,2)]
total_plots2 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- tmp_projects2[!is.na(tmp_projects2$treatment),]
    df_paired <- df_paired[df_paired$time %in% c("W0","W1","W2","W4","W8"),]
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,diease_Scores[dis]]-tmp[tmp$time=="W0",diease_Scores[dis]]
        if (diease_Scores[dis] %in% c("Propionate_producing_bac","Butyrate_producing_bac","Valeric_producing_bac","probiotic.bac")) {
            tmp$value[tmp$value > 250] <- 250
            tmp$value[tmp$value < -250] <- -250
        }
        if (diease_Scores[dis]%in% c("Bac.observed","Bac.shannon")) {
            tmp$value[tmp$value > 30] <- 30
            tmp$value[tmp$value < -30] <- -30
        }
        return(tmp)
    })
    df_paired1 <- do.call(rbind,df_paired1_)
    df_paired1$treatment <- factor(df_paired1$treatment,levels=c("Nystatin","Clotrimazole","Fluconazole"))
    if (diease_Scores[dis]%in% c("Bac.observed","Bac.shannon")) {test <- t.test_wrapper} else {test <- wilcox.test_wrapper}
    df_paired2 <- df_paired1
    df_paired2$time <- as.numeric(gsub("W","",df_paired2$time))
    loess_fit_1 <- loess(value ~ time, data = df_paired2[df_paired2$treatment == "Nystatin", ])
    loess_fit_2 <- loess(value ~ time, data = df_paired2[df_paired2$treatment == "Fluconazole", ])
    x_range <- seq(min(df_paired2$time), max(df_paired2$time), length.out = 100)
    pred_1 <- predict(loess_fit_1, newdata = data.frame(time = x_range))
    pred_2 <- predict(loess_fit_2, newdata = data.frame(time = x_range))
    ks_test_result <- ks.test(pred_1, pred_2)$p.value

    p1 <- ggplot(df_paired1, aes_string(x = "time", y = "value")) + 
        geom_line(aes(group=patient),size = 0.6,alpha=0.8,color="lightgrey") +
        geom_jitter(color="black",width = 0.1,alpha=0.5,size=1)+
        geom_boxplot(outlier.shape = NA,alpha=0)+
        theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
        stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+facet_wrap(~treatment,ncol=3)+
        labs(title=paste0(diease_Scores[dis],"\n (Flu vs Nys p:",ks_test_result,")"),y = paste("Δ")) +NoLegend()+
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = test)+
        geom_smooth(aes(group = 1, color = treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.2)+
        stat_summary(fun = mean, geom = "point",aes(group = 1,color = treatment),size=2)
    plot <- plot_grid(p1,ncol=1)
    return(plot)
    })
plot <- CombinePlots(c(total_plots2),nrow=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/Fig3.3.svg", plot=plot,width = 25, height = 4,dpi=300)
~~~

![image-20250205130709343](./Chapter3.assets/image-20250205130709343.png)

## 3.1.4 LEfSe-Based Identification of Differentially Abundant Taxa

- **Figure 3.4**: **Longitudinal trajectory of key bacterial species**, including **probiotic and SCFA-producing bacteria**, comparing Nystatin and Fluconazole treatments.

~~~R
library(lefser)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(METAG_tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)

MMC.METAG.inc_all <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.METAG.inc_all.lefser.v2.rds")
tmp_projects <- as.data.frame(colData(METAG_tse))
probiotic <- c(grep("Clostridium",unique(MMC.METAG.inc_all[[3]]),value=TRUE),"Lactobacillus","Lactobacillus_acidophilus","Bifidobacteriales","Streptococcus_thermophilus","Bacillus",
    grep("Streptococcus",unique(MMC.METAG.inc_all[[3]]),value=TRUE),grep("Lactobacillus",unique(MMC.METAG.inc_all[[3]]),value=TRUE),grep("Bifidobacteriales",unique(MMC.METAG.inc_all[[3]]),value=TRUE),
    grep("Faecalibacterium",unique(MMC.METAG.inc_all[[3]]),value=TRUE),grep("Bacillus",unique(MMC.METAG.inc_all[[3]]),value=TRUE))
probiotic.bac <- intersect(probiotic,relativeAb_all$names)

Butyrate_producing_bac <- c("Anaerostipes_hadrus","Clostridium_butyricum","Coprococcus_comes","Eubacterium_hallii","Eubacterium_rectale","Faecalibacterium_prausnitzii","Roseburia_hominis","Roseburia_intestinalis","Butyrivibrio_fibrisolvens","Ruminococcus_bromii","Megasphaera_elsdenii","Odoribacter_splanchnicus")
Propionate_producing_bac <- c("Bacteroides_thetaiotaomicron","Bacteroides_fragilis","Bacteroides_ovatus","Bacteroides_uniformis","Phascolarctobacterium_succinatutens","Veillonella_parvula","Dialister_succinatiphilus","Prevotella_copri","Propionibacterium_freudenreichii","Anaerovibrio_lipolytica","Oddibacterium_anthropi")
Valeric_producing_bac <- c("Butyrivibrio_fibrisolvens","Megasphaera_elsdenii","Coprococcus_eutactus","Blautia")
Butyrate_producing_bac <- sort(c(Butyrate_producing_bac,Propionate_producing_bac,Valeric_producing_bac,probiotic.bac))
relabundance_sel.bac <- relativeAb_all[relativeAb_all$names %in% Butyrate_producing_bac,]
tmp_projects <- as.data.frame(colData(METAG_tse))
tmp_projects$treatment <- factor(tmp_projects$treatment,levels=c("Nystatin","Clotrimazole","Fluconazole"))
relabundance_tmp <- relabundance_sel.bac[,rownames(tmp_projects)]
rownames(relabundance_tmp) <- gsub(" ","_",rownames(relabundance_tmp))
relabundance_tmp <- relabundance_tmp[rowSums(relabundance_tmp > 0) > 2, ]
sel.bac1 <- rownames(relabundance_tmp)
tmp_projects <- as.data.frame(cbind(tmp_projects,t(100*relabundance_tmp)))
tmp_projects1 <- tmp_projects[tmp_projects$treatment %in% c("Nystatin","Fluconazole"),]
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_celsius")),jdb_palette(c("brewer_spectra")))
library(lme4)
library(ARTool)
library(coin)
library(mgcv)
library(boot)
library(pracma)  # For numerical integration
library(fda)
sel.bac1 <- sort(sel.bac1)
total_plots2 <- lapply(1:length(sel.bac1),function(dis) {
    df_paired <- tmp_projects1
    df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
    colnames(df_paired)[colnames(df_paired)==sel.bac1[dis]] <- "value"
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,"value"]-tmp[tmp$time=="W0","value"]
        tmp[,"value"][tmp[,"value"]>500] <- 500
        tmp[,"value"][tmp[,"value"] < -500] <- -500
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
    perm_diffs <- unlist(future_lapply(seq_len(1000),function(i) {
      permuted <- sample(df_paired2$treatment)
      loess_perm_1 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Nystatin", ])
      loess_perm_2 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Fluconazole", ])
      pred_perm_1 <- predict(loess_perm_1, newdata = data.frame(time_numeric = x_range))
      pred_perm_2 <- predict(loess_perm_2, newdata = data.frame(time_numeric = x_range))
      perm_diffs <- mean(pred_perm_1 - pred_perm_2)
      perm_diffs
    }))
    perm_diffs <- perm_diffs[!is.na(perm_diffs)]
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
    area_diffs_perm <- area_diffs_perm[!is.na(area_diffs_perm)]
    p_value_area <- mean(area_diffs_perm >= area_diff)

    boot_diffs <- unlist(future_lapply(seq_len(1000),function(i) {
        boot_sample <- df_paired2[sample(seq_len(nrow(df_paired2)), replace = TRUE), ]
        loess_boot_1 <- loess(value ~ time_numeric, data = boot_sample[boot_sample$treatment == "Nystatin", ])
        loess_boot_2 <- loess(value ~ time_numeric, data = boot_sample[boot_sample$treatment == "Fluconazole", ])
        pred_boot_1 <- predict(loess_boot_1, newdata = data.frame(time_numeric = x_range))
        pred_boot_2 <- predict(loess_boot_2, newdata = data.frame(time_numeric = x_range))
        boot_diffs <- mean(pred_boot_1 - pred_boot_2)
        boot_diffs
    }))
    boot_diffs <- boot_diffs[!is.na(boot_diffs)]
    p_value_boot <- mean(abs(boot_diffs) >= abs(observed_diff))

    rss_1 <- sum((df_paired2[df_paired2$treatment == "Nystatin", "value"] - predict(loess_fit_1))^2)
    rss_2 <- sum((df_paired2[df_paired2$treatment == "Fluconazole", "value"] - predict(loess_fit_2))^2)
    n1 <- length(df_paired2[df_paired2$treatment == "Nystatin", "value"])  # Sample size group 1
    n2 <- length(df_paired2[df_paired2$treatment == "Fluconazole", "value"])  # Sample size group 2
    f_stat <- (rss_1 / (n1 - 2)) / (rss_2 / (n2 - 2))
    Ftestp_value <- pf(f_stat, df1 = n1 - 2, df2 = n2 - 2, lower.tail = FALSE)# Perform an F-test

    # Perform t-tests for each point in x_range
    ks_tests <- sapply(seq_along(x_range), function(i) {
      t.test(
        df_paired2$value[df_paired2$treatment == "Nystatin"],
        df_paired2$value[df_paired2$treatment == "Fluconazole"]
      )$p.value
    })
    p_values_adjusted <- p.adjust(ks_tests, method = "bonferroni")

    gam_model <- gam(value ~s(time_numeric, k = 2) + treatment + s(time_numeric, by = treatment, k = 2), data = df_paired2)
    p_value_gam <- summary(gam_model)$s.table[3, 4]

    ks_test <- ks.test(pred_1, pred_2)
    ks_p_value <- ks_test$p.value

    rmse_diff <- sqrt(mean((pred_1 - pred_2)^2))
    rmse_boot <- unlist(future_lapply(seq_len(1000), function(i) {
      boot_sample <- df_paired2[sample(seq_len(nrow(df_paired2)), replace = TRUE), ]
      loess_boot_1 <- loess(value ~ time_numeric, data = boot_sample[boot_sample$treatment == "Nystatin", ])
      loess_boot_2 <- loess(value ~ time_numeric, data = boot_sample[boot_sample$treatment == "Fluconazole", ])
      pred_boot_1 <- predict(loess_boot_1, newdata = data.frame(time_numeric = x_range))
      pred_boot_2 <- predict(loess_boot_2, newdata = data.frame(time_numeric = x_range))
      sqrt(mean((pred_boot_1 - pred_boot_2)^2))
    }))
    rmse_boot <- rmse_boot[!is.na(rmse_boot)]
    p_value_rmse <- mean(rmse_boot >= rmse_diff)

    # Compute derivatives
    deriv_1 <- diff(pred_1) / diff(x_range)
    deriv_2 <- diff(pred_2) / diff(x_range)

    wilcox_test <- wilcox.test(pred_1, pred_2)
    wilcox_p_value <- wilcox_test$p.value

    mean_deriv_diff <- mean(abs(deriv_1 - deriv_2))# Permutation test for derivative differences
    deriv_diffs_perm <- unlist(future_lapply(seq_len(1000), function(i) {
      permuted <- sample(df_paired2$treatment)
      loess_perm_1 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Nystatin", ])
      loess_perm_2 <- loess(value ~ time_numeric, data = df_paired2[permuted == "Fluconazole", ])
      pred_perm_1 <- predict(loess_perm_1, newdata = data.frame(time_numeric = x_range))
      pred_perm_2 <- predict(loess_perm_2, newdata = data.frame(time_numeric = x_range))
      deriv_perm_1 <- diff(pred_perm_1) / diff(x_range)
      deriv_perm_2 <- diff(pred_perm_2) / diff(x_range)
      mean(abs(deriv_perm_1 - deriv_perm_2))
    }))
    p_value_deriv <- mean(deriv_diffs_perm >= mean_deriv_diff)

    p2 <- ggplot(df_paired2, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
    geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess", method.args = list(degree=1),se=TRUE, level = 0.5)+
    stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(0),label = "p.format") +
    theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
    labs(title = paste0(sel.bac1[dis],"\n", "Flu vs Nystatin (bonferroni: ", p_values_adjusted," | ","Permutation.p: ",p_value,
        "|\n","Area-based p: ", p_value_area," | ","Bootstrapping p: ",p_value_boot,
        "|\n","F-test: ",Ftestp_value," | ","GAM p: ",p_value_gam,
        "|\n","KS-test p: ", ks_p_value," | ","Wilcoxon p: ",p_value_boot," |\n ","Derivative p: ", p_value_deriv, ")"),y = "Δ")
    message(dis)
    return(p2)
    })
plot <- CombinePlots(c(total_plots2),ncol=5)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/Fig3.4.1.png", plot=plot,width = 35, height = 15,dpi=300)
~~~

![Fig3.4.1](./Chapter3.assets/Fig3.4.1.png)

## 3.1.5 Ridge density plots of **SCFA metabolites** across treatments.

~~~R
All_dds_merge <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/All_dds_merge_metabolomics.v3.rds")
All_dds <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/All_dds_metabolomics.v2.rds")
Nature.Dysbiosis.metabolomic <- read.csv("/mnt/d/xiangyu.ubuntu/workshop/MMC/nature.Dysbiosis.metabolomic.csv")
Nature.Dysbiosis.metabolomic1 <- Nature.Dysbiosis.metabolomic[!is.na(Nature.Dysbiosis.metabolomic$HMDB),]
Nature.Dysbiosis.metabolomic1 <- Nature.Dysbiosis.metabolomic1[Nature.Dysbiosis.metabolomic1$HMDB %in% rowData(All_dds[[3]])$HMDB_ID,]
Health.UP <- Nature.Dysbiosis.metabolomic1[Nature.Dysbiosis.metabolomic1$Dysbiosis.Coefficient..non.IBD. > 0,]$HMDB
Dysbiosis.UP <- Nature.Dysbiosis.metabolomic1[Nature.Dysbiosis.metabolomic1$Dysbiosis.Coefficient..non.IBD. < -0,]$HMDB
length(Health.UP)
length(Dysbiosis.UP)
All_gene <- list(Health.UP,Dysbiosis.UP)
names(All_gene) <- c("Health.UP","Dysbiosis.UP")
Sel.Metab <- c("BUTYRIC ACID","PROPIONIC ACID","VALERIC ACID",c("L-HISTIDINE","L-PHENYLALANINE","L-TRYPTOPHAN","GLYCINE","TAURINE","INDOLE","L-LACTIC ACID","SUCCINIC ACID","LITHOCHOLIC ACID","DEOXYCHOLIC ACID"))
Metab_data_ <- lapply(1:length(All_dds_merge),function(x) {
    tmp <- as.data.frame(t(All_dds_merge[[x]][Sel.Metab,]))
    tmp <- tmp[grep("MMC",rownames(tmp),value=TRUE),]
    tmp <- as.data.frame(cbind(tmp,colData(All_dds[[x]])[rownames(tmp),c("metabolomics_names","patient","time","sample","treatment","Diagnosis_new","group","DAI","max.DAI.dates","status",
        c("UC.score.v2","stooling.frequency.scores","stooling.frequency","blood.stools.scores","Mucosal.Appearance","Physician.rating",
            "liquid.stooling.scores","abdominal.pain","general.wellbeing","Abdominal.mass"))]))
    return(tmp)
})
Metab_data <- do.call(rbind,Metab_data_)

Healthy <- c("BUTYRIC ACID","PROPIONIC ACID","VALERIC ACID")
pal <- jdb_palette("brewer_blue")[6:2]
total_plots2 <- lapply(1:length(Healthy),function(dis) {
    df_paired <- Metab_data
    df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Clotrimazole","Fluconazole"))
    colnames(df_paired)[colnames(df_paired)==Healthy[dis]] <- "value"
    df_paired$value <- log(df_paired$value+1,2)
    df_paired$time <- factor(df_paired$time,levels=c("W8","W4","W2","W1","W0"))
    p1 <- ggplot(df_paired, aes_string(x="value", y="time", fill="time")) + ggridges::geom_density_ridges() +
    theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+NoLegend()+facet_wrap(~treatment,ncol=3,scales="free_x")+
    labs(title = paste0(Healthy[dis]," raw"),y = paste(Healthy[dis]))+
    geom_vline(xintercept = median(df_paired[df_paired$time=="W0","value"]), color = 'lightgrey', size = 0.8)
    return(plot_grid(p1,nrow=1))
    })
plot <- CombinePlots(c(total_plots2),ncol=1)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/Fig3.5.svg", plot=plot,width = 8, height = 10,dpi=300)
~~~

![image-20250205130929072](./Chapter3.assets/image-20250205130929072.png)

## 3.1.6 Changes in **bile acid metabolism**, crucial for host-microbiome interactions.

~~~R
bile_acids <- c("HMDB0000619","HMDB0000518","HMDB0000626","HMDB0000761","HMDB0000946","HMDB0000138","HMDB0000637","HMDB0000036",
    "HMDB0000951","HMDB0000896","HMDB0000865","HMDB0000328","HMDB0000391","HMDB0000717")
names(bile_acids) <- c(
  "Cholic acid",
  "Chenodeoxycholic acid",
  "Deoxycholic acid",
  "Lithocholic acid",
  "Ursodeoxycholic acid",
  "Glycocholic acid",
  "Glycochenodeoxycholic acid",
  "Taurocholic acid",
  "Taurochenodeoxycholic acid",
  "Taurodeoxycholic acid",
  "Muricholic acid",
  "12-Ketolithocholic acid",
  "7-Ketodeoxycholic acid",
  "Isolithocholic acid")

bile_acids1 <- c(rownames(rowData(All_dds[[3]])[rowData(All_dds[[3]])$HMDB_ID %in% bile_acids,]),
rownames(rowData(All_dds[[3]])[rowData(All_dds[[3]])$CHEBI_ID %in% c("CHEBI:81244","CHEBI:52023","CHEBI:139131"),]))
Metab_data_ <- lapply(1:length(All_dds_merge),function(x) {
    tmp <- as.data.frame(t(All_dds_merge[[x]][bile_acids1,]))
    tmp <- tmp[grep("MMC",rownames(tmp),value=TRUE),]
    tmp <- as.data.frame(cbind(tmp,colData(All_dds[[x]])[rownames(tmp),c("metabolomics_names","patient","time","sample","treatment","Diagnosis_new","group","DAI","max.DAI.dates","status",
        c("UC.score.v2","stooling.frequency.scores","stooling.frequency","blood.stools.scores","Mucosal.Appearance","Physician.rating",
            "liquid.stooling.scores","abdominal.pain","general.wellbeing","Abdominal.mass"))]))
    return(tmp)
})
Metab_data <- do.call(rbind,Metab_data_)
bile_acids1 <- sort(bile_acids1)
pal <- jdb_palette("brewer_blue")[6:2]
total_plots2 <- lapply(1:length(bile_acids1),function(dis) {
    df_paired <- Metab_data
    df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Clotrimazole","Fluconazole"))
    colnames(df_paired)[colnames(df_paired)==bile_acids1[dis]] <- "value"
    df_paired$value <- log(df_paired$value+1,2)
    df_paired$time <- factor(df_paired$time,levels=c("W8","W4","W2","W1","W0"))
    p1 <- ggplot(df_paired, aes_string(x="value", y="time", fill="time")) + ggridges::geom_density_ridges() +
    theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+NoLegend()+facet_wrap(~treatment,ncol=3,scales="free_x")+
    labs(title = paste0(bile_acids1[dis]," raw"),y = paste(bile_acids[dis]))+
    geom_vline(xintercept = median(df_paired[df_paired$time=="W0","value"]), color = 'lightgrey', size = 0.8)
    return(plot_grid(p1,nrow=1))
    })
plot <- CombinePlots(c(total_plots2),ncol=5)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/Fig3.6.svg", plot=plot,width = 20, height = 8,dpi=300)
~~~

![image-20250205131110145](./Chapter3.assets/image-20250205131110145.png)

# 3.2 Supplementary Figure Analysis of Bacterial Community and Treatment Effects

## 3.2.1 Hierarchical Clustering and Phylogenetic Analysis

~~~R
MMC.metagenomic.counts <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/metagenomic_data1.v2.rds")
MMC.metagenomic.taxa1 <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.metagenomic.taxa.v2.rds")
MMC.metagenomic.samples <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/metagenomic.samples_info1.v2.rds")
taxa1 <- MMC.metagenomic.taxa1
Species <- unique(taxa1$Species)
MMC.metagenomic.counts_ <- lapply(1:length(Species),function(x) {
    tmp_counts <- MMC.metagenomic.counts[rownames(taxa1[taxa1$Species==Species[x],]),]
    if (nrow(tmp_counts)==1){tmp_counts <- tmp_counts} else {
        tmp_counts <- as.data.frame(t(colSums(tmp_counts)))
    }
    tmp_counts
    })
MMC.metagenomic.counts <- do.call(rbind,MMC.metagenomic.counts_)
MMC.metagenomic.taxa1 <- MMC.metagenomic.taxa1[!duplicated(MMC.metagenomic.taxa1$Species),]
rownames(MMC.metagenomic.counts) <- rownames(MMC.metagenomic.taxa1) <- MMC.metagenomic.taxa1$Species
MMC.metagenomic.taxa1[grep("_1",MMC.metagenomic.taxa1$Species,value=FALSE),]
MMC.metagenomic.taxa1[grep("_NA",MMC.metagenomic.taxa1$Species,value=FALSE),]
MMC.metagenomic.taxa1[grep("Streptococcus_thermophilus",MMC.metagenomic.taxa1$Species,value=FALSE),]
MMC.metagenomic.taxa1[grep("Leuconostoc.spp",MMC.metagenomic.taxa1$Species,value=FALSE),]
library(microeco)
library(mecodev)
MMC.metagenomic.taxa1 <- MMC.metagenomic.taxa1[MMC.metagenomic.taxa1$Genus!="NA",]
MMC.metagenomic.counts1 <- MMC.metagenomic.counts[intersect(rownames(MMC.metagenomic.taxa1),rownames(MMC.metagenomic.counts)),]
MMC.metagenomic.counts1 <- MMC.metagenomic.counts1[rowSums(MMC.metagenomic.counts1)>0,colSums(MMC.metagenomic.counts1)>0]
MMC.metagenomic.taxa11 <- MMC.metagenomic.taxa1[rownames(MMC.metagenomic.counts1),]
MMC.metagenomic.samples <- MMC.metagenomic.samples[colnames(MMC.metagenomic.counts1),]
MMC.metagenomic.samples$time <- factor(MMC.metagenomic.samples$time,levels=c("W0","W1","W2","W4","W8"))
MMC.metagenomic.samples$group.v2 <- factor(MMC.metagenomic.samples$group.v2,levels=c("Resp","NR"))

tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.metagenomic.counts1)),colData = DataFrame(MMC.metagenomic.samples),rowData = DataFrame(MMC.metagenomic.taxa11))
METAG_tse <- transformAssay(tse_taxa, MARGIN = "samples", method = "relabundance")
METAG_tse <- addPerCellQC(METAG_tse)
METAG_tse <- mia::estimateRichness(METAG_tse, assay.type = "counts", index = "observed", name="observed")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "coverage", name = "coverage")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "gini_simpson", name = "gini_simpson")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "inverse_simpson", name = "inverse_simpson")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "log_modulo_skewness", name = "Rarity")
METAG_tse <- mia::estimateDiversity(METAG_tse, assay.type = "counts",index = "shannon", name = "shannon")
METAG_tse <- estimateDominance(METAG_tse, assay.type = "counts", index="relative", name = "Dominance")
METAG_tse <- mia::estimateDivergence(METAG_tse,assay.type = "counts",reference = "median",FUN = vegan::vegdist)
colData(METAG_tse)$total_raw_counts <- colSums(assay(METAG_tse, "counts"))
METAG_tse$time.v2 <- as.character(METAG_tse$time)
METAG_tse$time.v2[METAG_tse$time.v2 %in% c("W0")] <- "Pre"
METAG_tse$time.v2[METAG_tse$time.v2 %in% c("W1","W2","W4")] <- "Post"
METAG_tse$time.v2[METAG_tse$time.v2 %in% c("W8")] <- "LTM"
METAG_tse$time.v2 <- factor(METAG_tse$time.v2,levels=c("Pre","Post","LTM"))
library(lefser)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(METAG_tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)
inc_all <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.METAG.inc_all.lefser.v2.rds")
probiotic <- c(grep("Clostridium",unique(inc_all[[3]]),value=TRUE),"Lactobacillus","Lactobacillus_acidophilus","Bifidobacteriales","Streptococcus_thermophilus","Bacillus",
    grep("Streptococcus",unique(inc_all[[3]]),value=TRUE),grep("Lactobacillus",unique(inc_all[[3]]),value=TRUE),grep("Bifidobacteriales",unique(inc_all[[3]]),value=TRUE),
    grep("Faecalibacterium",unique(inc_all[[3]]),value=TRUE),grep("Bacillus",unique(inc_all[[3]]),value=TRUE))

aa <- jdb_palette("brewer_spectra")
relabundance_tmp <- relativeAb_all[,rownames(MMC.metagenomic.samples)]
relabundance <-  data.frame(W0=rowMeans(relabundance_tmp[,rownames(MMC.metagenomic.samples[MMC.metagenomic.samples$time=="W0",])]),
    W1=rowMeans(relabundance_tmp[,rownames(MMC.metagenomic.samples[MMC.metagenomic.samples$time=="W1",])]),
    W2=rowMeans(relabundance_tmp[,rownames(MMC.metagenomic.samples[MMC.metagenomic.samples$time=="W2",])]),
    W4=rowMeans(relabundance_tmp[,rownames(MMC.metagenomic.samples[MMC.metagenomic.samples$time=="W4",])]),
    names=relativeAb_all$names)
relabundance_DA <- relabundance[relabundance$names %in% unique(as.character(unlist(inc_all))),]
rownames(relabundance_DA) <- relabundance_DA$names
relabundance_tmp <- relabundance_DA[,c("W0","W1","W2","W4")]
relabundance_tmp <- as.matrix(relabundance_tmp)
relabundance.zscore <- sweep(relabundance_tmp - rowMeans(relabundance_tmp), 1, matrixStats::rowSds(relabundance_tmp),`/`)
Species_anno <- data.frame(Species=rownames(relabundance.zscore),row.names=rownames(relabundance.zscore),Nystatin="No",Clotrimazole="No",Fluconazole="No")
Species_anno[Species_anno$Species %in% inc_all[[1]],"Nystatin"] <- "Inc"
Species_anno[Species_anno$Species %in% inc_all[[2]],"Clotrimazole"] <- "Inc"
Species_anno[Species_anno$Species %in% inc_all[[3]],"Fluconazole"] <- "Inc"

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
p <- p + new_scale_fill() + geom_fruit(data=d_genomic,geom=geom_tile,mapping=aes(y=label, fill=Clotrimazole),offset=0.15, pwidth=0.1)+scale_fill_manual(values = c(jdb_palette("corona")[1],"lightgrey"))
p <- p + new_scale_fill() + geom_fruit(data=d_genomic,geom=geom_tile,mapping=aes(y=label, fill=Fluconazole),offset=0.15, pwidth=0.1)+scale_fill_manual(values = c(jdb_palette("corona")[2],"lightgrey"))
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/S3.1.ggtree.svg", plot=p,width = 6, height = 5,dpi=300)
~~~

![image-20250205131930003](./Chapter3.assets/image-20250205131930003.png)

## 3.2.2 Longitudinal Analysis of Treatment Effects on Bacterial Abundance

~~~R
W1.lefser <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.METAG.DAA_all.W1.lefser.v2.rds")
W2.lefser <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.METAG.DAA_all.W2.lefser.v2.rds")
W4.lefser <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.METAG.DAA_all.W4.lefser.v2.rds")
W8.lefser <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.METAG.DAA_all.W8.lefser.v2.rds")
names(W1.lefser) <- c("Nystatin","Clotrimazole","Fluconazole")
library(lefser)
Sel_type <- c("Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(METAG_tse,rank = Sel_type[x],detection = 0.1,prevalence = 0.1,as_relative = FALSE)
    se_total <- SummarizedExperiment(assays = list(counts = assays(tse_tmp)[["counts"]]),rowData = rowData(tse_tmp),colData = colData(tse_tmp))
    se_total <- relativeAb(se_total)
    relativeAb <- as.data.frame(assays(se_total)[["rel_abs"]])
    relativeAb <- log(relativeAb+1,2)
    relativeAb$names <- rownames(relativeAb)
    relativeAb$type <- Sel_type[x]
    return(relativeAb)
    })
relativeAb_all <- do.call(rbind,relativeAb_all_)
Bact_all <- unique(relativeAb_all$names)
probiotic <- c(grep("Clostridium",Bact_all,value=TRUE),"Lactobacillus","Lactobacillus_acidophilus","Bifidobacteriales","Streptococcus_thermophilus","Bacillus",
    grep("Streptococcus",Bact_all,value=TRUE),grep("Lactobacillus",Bact_all,value=TRUE),grep("Bifidobacteriales",Bact_all,value=TRUE),
    grep("Faecalibacterium",Bact_all,value=TRUE),grep("Bacillus",Bact_all,value=TRUE))
Butyrate_producing_bac <- c("Anaerostipes_hadrus","Clostridium_butyricum","Coprococcus_comes","Eubacterium_hallii","Eubacterium_rectale","Faecalibacterium_prausnitzii","Roseburia_hominis","Roseburia_intestinalis","Butyrivibrio_fibrisolvens","Ruminococcus_bromii","Megasphaera_elsdenii","Odoribacter_splanchnicus")
Propionate_producing_bac <- c("Bacteroides_thetaiotaomicron","Bacteroides_fragilis","Bacteroides_ovatus","Bacteroides_uniformis","Phascolarctobacterium_succinatutens","Veillonella_parvula","Dialister_succinatiphilus","Prevotella_copri","Propionibacterium_freudenreichii","Anaerovibrio_lipolytica","Oddibacterium_anthropi")
Valeric_producing_bac <- c("Butyrivibrio_fibrisolvens","Megasphaera_elsdenii","Coprococcus_eutactus","Blautia")
probiotic <- c(probiotic,Butyrate_producing_bac,Propionate_producing_bac,Valeric_producing_bac)
probiotic <- c(probiotic,c(grep("Faecalibacterium",Bact_all,value=TRUE),
grep("Roseburia",Bact_all,value=TRUE),
grep("Eubacterium",Bact_all,value=TRUE),
grep("Sutterella",Bact_all,value=TRUE),
grep("Gemmiger",Bact_all,value=TRUE),
grep("Clostridium",Bact_all,value=TRUE),
grep("Bifidobacterium",Bact_all,value=TRUE),
grep("Bacteroides",Bact_all,value=TRUE),
grep("Collinsella",Bact_all,value=TRUE),
grep("Akkermansia",Bact_all,value=TRUE),
grep("Christensenellaceae",Bact_all,value=TRUE),
grep("Multikap",Bact_all,value=TRUE)))

pathogenic <- c(grep("Escherichia",Bact_all,value=TRUE),
grep("Fusobacterium",Bact_all,value=TRUE),
grep("Ruminococcus",Bact_all,value=TRUE),
grep("Clostridium",Bact_all,value=TRUE),
grep("Veillonella",Bact_all,value=TRUE),
grep("Streptococcus",Bact_all,value=TRUE),
grep("Enterococcus",Bact_all,value=TRUE),
grep("Campylobacter",Bact_all,value=TRUE),
grep("Gammaproteobacteria",Bact_all,value=TRUE),
grep("Deltaproteobacteria",Bact_all,value=TRUE))

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
    tmp_all <- tmp_all[tmp_all$type %in% c("Species"),]
    tmp_all$LDA.scores[tmp_all$LDA.scores > 4] <- 4
    tmp_all$LDA.scores[tmp_all$LDA.scores < -4] <- -4
    tmp_all <- tmp_all[tmp_all$pvalue <= 0.4,]
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all_c <- tmp_all_c[tmp_all_c$Names %in% probiotic,]
    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
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
        ylab(NULL) + ggtitle(names(W1.lefser)[x]," fungal Species") + scale_size(range=c(6,2))+theme_classic()+
        theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1,color = "black", size = 12),axis.text.y = element_text(color = "black", size = 12, face = "plain"))
    return(plot)
    })
plot <- CombinePlots(c(All_plots),ncol=3)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/S3.2.svg", plot=plot,width = 18, height = 7,dpi=300)
~~~

![image-20250205132016086](./Chapter3.assets/image-20250205132016086.png)

## 3.2.3 Statistical Analysis of Treatment-Specific Bacterial Changes

~~~R
LDA_bac_all_ <- lapply(1:length(W1.lefser),function(x) {
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
    tmp_all1 <- tmp_all[tmp_all$type %in% c("Species"),]
    tmp_all1$LDA.scores[tmp_all1$LDA.scores > 4] <- 4
    tmp_all1$LDA.scores[tmp_all1$LDA.scores < -4] <- -4
    tmp_all1 <- tmp_all1[tmp_all1$pvalue <= 0.4,]
    tmp_all_c <- reshape2::dcast(tmp_all1,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all_c <- tmp_all_c[tmp_all_c$Names %in% probiotic,]
    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.5,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]
    o <- as.character(tmp_all_c$Names)

    tmp_all2 <- tmp_all[tmp_all$pvalue <= 0.4,]
    o <- as.character(tmp_all_c$Names)
    tmp_all2 <- tmp_all2[tmp_all2$Names %in% o,]
    tmp_all2$treatment <- names(W1.lefser)[x]
    tmp_all2
    })
LDA_bac_all <- do.call(rbind,LDA_bac_all_)
write.csv(LDA_bac_all,"/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/LDA_bac_all.csv")
~~~

~~~R

library(lefser)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))
Sel_type <- c("Genus","Species")
relativeAb_all_ <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(METAG_tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
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
aa <- jdb_palette("brewer_spectra")
treatment <- c("Nystatin","Clotrimazole","Fluconazole")
relabundance_all_ <- lapply(1:length(treatment),function(x) {
    tmp_projects <- colData(METAG_tse)[colData(METAG_tse)$treatment==treatment[x],]
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

MMC.METAG.inc_all <- mcreadRDS("/mnt/d/xiangyu.ubuntu/workshop/MMC/sample_info/final_Res/MMC.METAG.inc_all.lefser.v2.rds")
tmp_projects <- as.data.frame(colData(METAG_tse))
probiotic.bac <- intersect(probiotic,relabundance_all_W2_W0$names)
probiotic.bac <- sort(c(probiotic.bac))
relabundance_sel.bac <- relativeAb_all[relativeAb_all$names %in% probiotic.bac,]
tmp_projects <- as.data.frame(colData(METAG_tse))
tmp_projects$treatment <- factor(tmp_projects$treatment,levels=c("Nystatin","Clotrimazole","Fluconazole"))
relabundance_tmp <- relabundance_sel.bac[,rownames(tmp_projects)]
rownames(relabundance_tmp) <- gsub(" ","_",rownames(relabundance_tmp))
relabundance_tmp <- relabundance_tmp[rowSums(relabundance_tmp > 0) > 3, ]
sel.bac1 <- rownames(relabundance_tmp)
tmp_projects <- as.data.frame(cbind(tmp_projects,t(100*relabundance_tmp)))
tmp_projects1 <- tmp_projects[tmp_projects$treatment %in% c("Nystatin","Fluconazole"),]
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_celsius")),jdb_palette(c("brewer_spectra")))
library(lme4)
library(ARTool)
library(coin)
sel.bac1 <- sort(sel.bac1)
total_plots2 <- lapply(1:length(sel.bac1),function(dis) {
    df_paired <- tmp_projects1
    df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
    colnames(df_paired)[colnames(df_paired)==sel.bac1[dis]] <- "value"
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

    loess_fit_1 <- loess(value ~ time_numeric, data = df_paired2[df_paired2$treatment == "Nystatin", ])
    loess_fit_2 <- loess(value ~ time_numeric, data = df_paired2[df_paired2$treatment == "Fluconazole", ])
    x_range <- seq(min(df_paired2$time_numeric), max(df_paired2$time_numeric), length.out = 100)
    pred_1 <- predict(loess_fit_1, newdata = data.frame(time_numeric = x_range))
    pred_2 <- predict(loess_fit_2, newdata = data.frame(time_numeric = x_range))
    ks_test_result <- ks.test(pred_1, pred_2)$p.value

    rss_1 <- sum((df_paired2[df_paired2$treatment == "Nystatin", "value"] - predict(loess_fit_1))^2)
    rss_2 <- sum((df_paired2[df_paired2$treatment == "Fluconazole", "value"] - predict(loess_fit_2))^2)
    n1 <- length(df_paired2[df_paired2$treatment == "Nystatin", "value"])  # Sample size group 1
    n2 <- length(df_paired2[df_paired2$treatment == "Fluconazole", "value"])  # Sample size group 2
    f_stat <- (rss_1 / (n1 - 2)) / (rss_2 / (n2 - 2))
    Ftestp_value <- pf(f_stat, df1 = n1 - 2, df2 = n2 - 2, lower.tail = FALSE)# Perform an F-test

    p2 <- ggplot(df_paired2, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
    geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess", method.args = list(degree=1),se=TRUE, level = 0.5)+
    # geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = lm, formula = y ~ splines::bs(x, 1), se = TRUE, level = 0.5)+
    stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(0),label = "p.format") +
    theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
    labs(title = paste0(sel.bac1[dis],"\n", "Flu vs Nystatin (ANOVA p: ", anova.pvalue0," | ","art p: ",art.pvalue0,
        "|\n","kruskal p: ", kruskal.pvalue0," | ","LocationTests p: ",LocationTests.pvalue0,
        "|\n","TukeyHSD p: ",TukeyHSD.pvalue1," | ","loess p: ",Ftestp_value,"|\n","pairwise. p: ",pairwise.t.pvalue1," | ",
        "LMM p: ",LMM.pvalue1,")"),y = "Δ")
    message(dis)
    return(p2)
    })
plot <- CombinePlots(c(total_plots2),ncol=10)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/S3.3.png", plot=plot,width = 40, height = 20,dpi=300)
ggsave("/mnt/d/xiangyu.ubuntu/projects/MMC/Figures/figures_making/v3/S3.3.svg", plot=plot,width = 40, height = 20,dpi=300)
~~~

![S3.3](./Chapter3.assets/S3.3.png)

