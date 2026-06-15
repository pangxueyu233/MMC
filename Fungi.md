# 2.1 Main Figure Fungal Analysis: Data Processing, Visualization, and Statistical Testing

This chapter integrates fungal ITS sequencing with clinical metadata to assess treatment-induced shifts in the gut mycobiota. It covers environment setup, fungal community dynamics over time, LEfSe differential abundance, beta-diversity (PERMANOVA/PCoA), alpha-diversity and rarefaction QC, baseline (W0) composition, and *Candida* burden.

## 2.1.1 Environment Setup and Libraries

Loads all required R packages (mia/miaViz, vegan, ComplexHeatmap, DESeq2, lefser, etc.), sources the custom helper functions, and configures parallel processing.

~~~R
suppressPackageStartupMessages({
    library(miaViz)
    library(scater)
    library(ggplot2)
    library(cowplot)
    library(trqwe)
    library(BuenColors)
    library(dplyr)
    library(vegan)
    library(ComplexHeatmap)
    library(circlize)
    library(ggplotify)
    library(cowplot)
    library(bluster)
    library(mia)
    library(NbClust)
    library(dendextend)
    library(ALDEx2)
    library(Maaslin2)
    library(MicrobiomeStat)
    library(DESeq2)
    library(mikropml)
    library(MLeval)
    library(caret) 
    library(patchwork)
    library(ggpubr)
    library(Seurat)
    library(magrittr)
    library(reshape2)
})
source("./code/log-summery/MyBestFunction_scRNA.R.v4.R")
library(future)
library(future.apply)
options(future.globals.maxSize = 60 * 1024^3)
plan("multicore", workers = 15)
plan()
require(rJava)
require(xlsx)
require(trqwe)
~~~

## 2.1.2 Fungal Community Dynamics Over Time

Builds the fungal `TreeSummarizedExperiment` from ITS counts and plots the per-time-point dynamics (change from baseline) of fungal features, exported as `Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth`.



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
colSums(MMC.ITS.counts1[,sort(grep("W0",colnames(MMC.ITS.counts1),value=FALSE))])

sort(unique(MMC.ITS.taxa1$Species))
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

JCM_fungal <- c("Candida albicans", "Candida parapsilosis", "Candida tropicalis", "Candida glabrata", "Candida krusei", "Candida lusitaniae", "Candida dubliniensis", "Candida guilliermondii", "Candida rugosa", 
"Trichosporon asahii","Trichosporon mucoides","Trichosporon asteroides","Trichosporon cutaneum","Trichosporon ovoides","Trichosporon inkin",
"Rhodotorula glutinis", "Rhodotorula mucilaginosa", "Rhodotorula rubra", "Rhodotorula minuta",
"Geotrichum capitatum","Aspergillus fumigatus", "Aspergillus terreus","Zygomycetes","Fusarium","Acremonium","Scedosporium","Paecilomyces","Trichoderma longibrachiatum")
WHO_fungal <- c("Cryptococcus neoformans","Candida auris","Aspergillus fumigatus","Candida albicans",
"Candida glabrata","Histoplasma","Eumycetoma","Mucorales","Fusarium","Candida tropicalis","Candida parapsilosis",
"Scedosporium","Lomentospora prolificans","Coccidioides","Candida krusei","Cryptococcus gattii","Talaromyces mamneffei","Pneumocystis jirovecii","Paracoccidioides")
pathogenic.fungi <- gsub(" ","_",unique(c(JCM_fungal,WHO_fungal)))
relativeAb_all$names <- gsub(" ","_",relativeAb_all$names)

tmp_projects <- as.data.frame(colData(Fungal.tse))
relabundance_sel.fungi <- relativeAb_all[relativeAb_all$names %in% pathogenic.fungi,]
rownames(relabundance_sel.fungi) <- relabundance_sel.fungi$names
relabundance_tmp <- relabundance_sel.fungi[,rownames(tmp_projects)]
pathogenic.fungi <- rownames(relabundance_tmp)
tmp_projects <- as.data.frame(cbind(tmp_projects,t(100*relabundance_tmp)))
df <- tmp_projects[tmp_projects$treatment %in% c("Nystatin","Fluconazole"),]
df$treatment <- factor(df$treatment,levels=c("Nystatin","Fluconazole"))
# df <- df[!df$time %in% c("W1"),]
sel_pa <- unique(df$patient)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_celsius")),jdb_palette(c("brewer_spectra")))
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
# comb <- list(c("W0","W2"),c("W0","W4"),c("W0","W8"))
plots1 <- lapply(1:length(pathogenic.fungi),function(x) {
  df_new <- df
  df_tmp <- data_summary(df_new, varname=pathogenic.fungi[x], groupnames=c("time","treatment"))
  All_plots1 <- ggplot(df_tmp, aes_string(x = "time", y = pathogenic.fungi[x])) + 
      stat_summary(data=df_new,fun.data = mean_se, geom = "errorbar", width = 0.5) + geom_jitter(data=df_new,aes(color=patient),width = 0.1,alpha=0.5) + 
      geom_bar(aes(fill=time),alpha=1,stat="identity", position=position_dodge())+
      geom_line(data=df_new,aes(color=patient,group=patient),size = 0.6,alpha=0.2) +
      theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
      stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+facet_wrap(~treatment)+
      labs(title=paste0(pathogenic.fungi[x]," raw in relabundance"),y = paste(pathogenic.fungi[x]," relabundance raw")) +NoLegend()+
      geom_signif(data=df_new,comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
      return(All_plots1)
  })
plots2 <- lapply(1:length(pathogenic.fungi),function(x) {
  df_new <- df
  sel_pa <- unique(df$patient)
  paired.df_ <- lapply(1:length(sel_pa),function(patient) {
      tmp <- df_new[df_new$patient %in% sel_pa[patient],]
      tmp[,pathogenic.fungi[x]] <- tmp[,pathogenic.fungi[x]]-tmp[tmp$time=="W0",pathogenic.fungi[x]]
      return(tmp)
  })
  paired.df <- do.call(rbind,paired.df_)
  paired.df_tmp <- data_summary(paired.df, varname=pathogenic.fungi[x], groupnames=c("time","treatment"))
  All_plots1 <- ggplot(paired.df_tmp, aes_string(x = "time", y = pathogenic.fungi[x])) + 
      stat_summary(data=paired.df,fun.data = mean_se, geom = "errorbar", width = 0.5) + geom_jitter(data=paired.df,aes(color=patient),width = 0.1,alpha=0.5) + 
      geom_bar(aes(fill=time),alpha=1,stat="identity", position=position_dodge())+
      geom_line(data=paired.df,aes(color=patient,group=patient),size = 0.6,alpha=0.2) +
      theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
      stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+facet_wrap(~treatment)+
      labs(title=paste0(pathogenic.fungi[x]," changes in relabundance"),y = paste(pathogenic.fungi[x]," relabundance Change")) +NoLegend()+
      geom_signif(data=paired.df,comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
      message(x)
      return(All_plots1)
    })
total_plots3 <- lapply(1:length(pathogenic.fungi),function(dis) {
  df_paired <- df
  df_paired$time_numeric <- as.numeric(gsub("W", "", df_paired$time))
  df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
  colnames(df_paired)[colnames(df_paired)==pathogenic.fungi[dis]] <- "value"
  # df_paired[,"value"][df_paired[,"value"]>1000] <- 1000
  anova_result <- aov(value ~ treatment, data = df_paired)
  p_value <- summary(anova_result)[[1]]["treatment", "Pr(>F)"]
  pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
  p1 <- ggplot(df_paired, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
  geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess",se=TRUE, level = 0.5) +
  stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(1),label = "p.format") +
  theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
  labs(title = paste0(pathogenic.fungi[dis]," raw \n", "Flu vs Nystatin (ANOVA p: ", pvalue0,")"),y = paste(pathogenic.fungi[dis], " Change"))

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
  anova_result <- aov(value ~ treatment, data = df_paired2)
  p_value <- summary(anova_result)[[1]]["treatment", "Pr(>F)"]
  pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
  df_paired2[,c("ITS_names","treatment","value","time")]
  p2 <- ggplot(df_paired2, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
  geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess",se=TRUE, level = 0.5,span=1.2) +
  stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(1),label = "p.format") +
  theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
  labs(title = paste0(pathogenic.fungi[dis]," changes \n", "Flu vs Nystatin (ANOVA p: ", pvalue0,")"),y = paste(pathogenic.fungi[dis], " Change"))
  return(plot_grid(p1,p2,nrow=2))
    })
plot <- CombinePlots(c(plots1,plots2,total_plots3),nrow=3)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth.png", plot=plot,width = 40, height = 10,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth.svg", plot=plot,width = 40, height = 10,dpi=300)
~~~

![Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth](./Chapter2.assets/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth.png)

## 2.1.3 Fungal Dynamics — Per-batch Panels (Raw)

Renders the raw fungal dynamics as separate per-batch panels, exported as `Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth1`–`bacth5`.

~~~R

JCM_fungal <- c("Candida albicans", "Candida parapsilosis", "Candida tropicalis", "Candida glabrata", "Candida krusei", "Candida lusitaniae", "Candida dubliniensis", "Candida guilliermondii", "Candida rugosa", 
"Trichosporon asahii","Trichosporon mucoides","Trichosporon asteroides","Trichosporon cutaneum","Trichosporon ovoides","Trichosporon inkin",
"Rhodotorula glutinis", "Rhodotorula mucilaginosa", "Rhodotorula rubra", "Rhodotorula minuta",
"Geotrichum capitatum","Aspergillus fumigatus", "Aspergillus terreus","Zygomycetes","Fusarium","Acremonium","Scedosporium","Paecilomyces","Trichoderma longibrachiatum")
WHO_fungal <- c("Cryptococcus neoformans","Candida auris","Aspergillus fumigatus","Candida albicans",
"Candida glabrata","Histoplasma","Eumycetoma","Mucorales","Fusarium","Candida tropicalis","Candida parapsilosis",
"Scedosporium","Lomentospora prolificans","Coccidioides","Candida krusei","Cryptococcus gattii","Talaromyces mamneffei","Pneumocystis jirovecii","Paracoccidioides")
pathogenic.fungi <- unique(c(JCM_fungal,WHO_fungal))
gut_fungi <- c("Candida albicans","Candida tropicalis","Candida zeylanoides","Candida sake","Candida inconspicua",
  "Candida glaebosa","Candida pseudoglaebosa","Candida tetrigidarum","Candida saitoana","Kazachstania bulder",
  "Malassezia furfur","Trichosporon asahii","Yarrowia lipolytica","Mucor circinelloides","Schizophyllum commune",
  "Apiotrichum domesticum")
gut_associated_fungi <- c("Candida albicans","Candida tropicalis","Candida inconspicua","Candida glaebosa",
  "Candida pseudoglaebosa","Candida sake","Candida zeylanoides","Candida saitoana","Malassezia furfur",
  "Trichosporon asahii","Yarrowia lipolytica","Mucor circinelloides")
Sig_fungi <- c("Aspergillus","Aspergillus terreus","Candida","Candida albicans","Candida dubliniensis","Candida parapsilosis","Saccharomyces cerevisiae",
    "Cryptococcus","Fusarium","Talaromyces","Malassezia","Malassezia globosa","Clavispora","Clavispora lusitaniae","Meyerozyma",
    "Meyerozyma guilliermondii","Mucor","Acremonium","Exophiala","Rhodotorula","Trichosporon","Trichosporon asahii","Debaryomyces",
    "Exserohilum","Saccharomyces cerevisiae","Candida sake","Candida tropicalis","Candida albicans","Candida parapsilosis",
    "Meyerozyma guilliermondii","Clavispora lusitaniae","Candida albicans","Aspergillus","Fusarium","Pichia mandshurica","Aspergillus penicillioides","Penicillium","Aspergillus reticulatus","Fusarium equiseti","Malassezia","Mucor circinelloides","Mucor","Malassezia restricta","Candida argentea","Candida sake","Penicillium citrinum","Rhodotorula","Aspergillus piperis","Aspergillus ruber","Geotrichum silvicola")
gut_associated_fungi <- unique(c(gut_associated_fungi,gut_fungi,pathogenic.fungi,Sig_fungi))

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
pathogenic.fungi <- gsub(" ","_",gut_associated_fungi)
relativeAb_all$names <- gsub(" ","_",relativeAb_all$names)

tmp_projects <- as.data.frame(colData(Fungal.tse))
relabundance_sel.fungi <- relativeAb_all[relativeAb_all$names %in% pathogenic.fungi,]
rownames(relabundance_sel.fungi) <- relabundance_sel.fungi$names
relabundance_tmp <- relabundance_sel.fungi[,rownames(tmp_projects)]
pathogenic.fungi <- sort(rownames(relabundance_tmp))
tmp_projects <- as.data.frame(cbind(tmp_projects,100*t(relabundance_tmp)))
df <- tmp_projects[tmp_projects$treatment %in% c("Nystatin","Fluconazole"),]
df$treatment <- factor(df$treatment,levels=c("Nystatin","Fluconazole"))
# df <- df[!df$time %in% c("W1"),]
sel_pa <- unique(df$patient)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_celsius")),jdb_palette(c("brewer_spectra")))
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
plots1 <- lapply(1:length(pathogenic.fungi),function(x) {
  df_new <- df
  df_tmp <- data_summary(df_new, varname=pathogenic.fungi[x], groupnames=c("time","treatment"))
  All_plots1 <- ggplot(df_tmp, aes_string(x = "time", y = pathogenic.fungi[x])) + 
      stat_summary(data=df_new,fun.data = mean_se, geom = "errorbar", width = 0.5) + geom_jitter(data=df_new,aes(color=patient),width = 0.1,alpha=0.5) + 
      geom_bar(aes(fill=time),alpha=1,stat="identity", position=position_dodge())+
      geom_line(data=df_new,aes(color=patient,group=patient),size = 0.6,alpha=0.2) +
      theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
      stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+facet_wrap(~treatment)+
      labs(title=paste0(pathogenic.fungi[x]," raw in relabundance"),y = paste(pathogenic.fungi[x]," relabundance raw")) +NoLegend()+
      geom_signif(data=df_new,comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
      return(All_plots1)
  })
plots2 <- lapply(1:length(pathogenic.fungi),function(x) {
  df_new <- df
  sel_pa <- unique(df$patient)
  paired.df_ <- lapply(1:length(sel_pa),function(patient) {
      tmp <- df_new[df_new$patient %in% sel_pa[patient],]
      tmp[,pathogenic.fungi[x]] <- tmp[,pathogenic.fungi[x]]-tmp[tmp$time=="W0",pathogenic.fungi[x]]
      return(tmp)
  })
  paired.df <- do.call(rbind,paired.df_)
  paired.df_tmp <- data_summary(paired.df, varname=pathogenic.fungi[x], groupnames=c("time","treatment"))
  All_plots1 <- ggplot(paired.df_tmp, aes_string(x = "time", y = pathogenic.fungi[x])) + 
      stat_summary(data=paired.df,fun.data = mean_se, geom = "errorbar", width = 0.5) + geom_jitter(data=paired.df,aes(color=patient),width = 0.1,alpha=0.5) + 
      geom_bar(aes(fill=time),alpha=1,stat="identity", position=position_dodge())+
      geom_line(data=paired.df,aes(color=patient,group=patient),size = 0.6,alpha=0.2) +
      theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
      stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+facet_wrap(~treatment)+
      labs(title=paste0(pathogenic.fungi[x]," changes in relabundance"),y = paste(pathogenic.fungi[x]," relabundance Change")) +NoLegend()+
      geom_signif(data=paired.df,comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
      message(x)
      return(All_plots1)
    })
total_plots3 <- lapply(1:length(pathogenic.fungi),function(dis) {
  df_paired <- df
  df_paired$time_numeric <- as.numeric(gsub("W", "", df_paired$time))
  df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
  colnames(df_paired)[colnames(df_paired)==pathogenic.fungi[dis]] <- "value"
  # df_paired[,"value"][df_paired[,"value"]>1000] <- 1000
  anova_result <- aov(value ~ treatment, data = df_paired)
  p_value <- summary(anova_result)[[1]]["treatment", "Pr(>F)"]
  pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
  p1 <- ggplot(df_paired, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
  geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess",se=TRUE, level = 0.5) +
  stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(1),label = "p.format") +
  theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
  labs(title = paste0(pathogenic.fungi[dis]," raw \n", "Flu vs Nystatin (ANOVA p: ", pvalue0,")"),y = paste(pathogenic.fungi[dis], " Change"))

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
  anova_result <- aov(value ~ treatment, data = df_paired2)
  p_value <- summary(anova_result)[[1]]["treatment", "Pr(>F)"]
  pvalue0 <- ifelse(p_value < 0.001, "< 0.001", format(p_value, digits = 3))
  df_paired2[,c("ITS_names","treatment","value","time")]
  p2 <- ggplot(df_paired2, aes_string(x = "time_numeric", y = "value")) + geom_jitter(aes(color=treatment),size = 1)+ 
  geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess",se=TRUE, level = 0.5,span=1.2) +
  stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(1),label = "p.format") +
  theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
  labs(title = paste0(pathogenic.fungi[dis]," changes \n", "Flu vs Nystatin (ANOVA p: ", pvalue0,")"),y = paste(pathogenic.fungi[dis], " Change"))
  return(plot_grid(p1,p2,nrow=2))
    })
plot1 <- CombinePlots(c(plots1[1:10],plots2[1:10],total_plots3[1:10]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth1.png", plot=plot1,width = 49, height = 10,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth1.svg", plot=plot1,width = 49, height = 10,dpi=300)
plot2 <- CombinePlots(c(plots1[11:20],plots2[11:20],total_plots3[11:20]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth2.png", plot=plot2,width = 49, height = 10,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth2.svg", plot=plot2,width = 49, height = 10,dpi=300)
plot3 <- CombinePlots(c(plots1[21:30],plots2[21:30],total_plots3[21:30]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth3.png", plot=plot3,width = 49, height = 10,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth3.svg", plot=plot3,width = 49, height = 10,dpi=300)
plot4 <- CombinePlots(c(plots1[31:40],plots2[31:40],total_plots3[31:40]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth4.png", plot=plot4,width = 49, height = 10,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth4.svg", plot=plot4,width = 49, height = 10,dpi=300)
plot5 <- CombinePlots(c(plots1[41:50],plots2[41:50],total_plots3[41:50]),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth5.png", plot=plot5,width = 49, height = 10,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth5.svg", plot=plot5,width = 49, height = 10,dpi=300)
~~~

![Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth1](./Chapter2.assets/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth1.png)

![Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth2](./Chapter2.assets/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth2.png)

![Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth3](./Chapter2.assets/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth3.png)

![Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth4](./Chapter2.assets/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth4.png)

![Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth5](./Chapter2.assets/Fig2.fungi_dyn.delta.timepoints.raw_merge_bacth5.png)

## 2.1.4 Fungal Dynamics — Merged Panel

Combines the fungal dynamics into a single merged panel, exported as `Fig2.fungi_dyn.delta.timepoints_merge_bacth1`.

~~~R

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
pathogenic.fungi <- gsub(" ","_",c(gut_associated_fungi,"Saccharomyces","Saccharomyces cerevisiae"))
relativeAb_all$names <- gsub(" ","_",relativeAb_all$names)

tmp_projects <- as.data.frame(colData(Fungal.tse))
relabundance_sel.fungi <- relativeAb_all[relativeAb_all$names %in% pathogenic.fungi,]
rownames(relabundance_sel.fungi) <- relabundance_sel.fungi$names
relabundance_tmp <- relabundance_sel.fungi[,rownames(tmp_projects)]
pathogenic.fungi <- sort(c(rownames(relabundance_tmp)))
tmp_projects <- as.data.frame(cbind(tmp_projects,100*t(relabundance_tmp)))
df <- tmp_projects[tmp_projects$treatment %in% c("Nystatin","Fluconazole"),]
df$treatment <- factor(df$treatment,levels=c("Nystatin","Fluconazole"))
# df <- df[!df$time %in% c("W1"),]
sel_pa <- unique(df$patient)
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_celsius")),jdb_palette(c("brewer_spectra")))
library(lme4)
library(ARTool)
library(coin)
total_plots2 <- lapply(1:length(pathogenic.fungi),function(dis) {
    df_paired <- df
    df_paired$treatment <- factor(df_paired$treatment,levels=c("Nystatin","Fluconazole"))
    colnames(df_paired)[colnames(df_paired)==pathogenic.fungi[dis]] <- "value"
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,"value"]-tmp[tmp$time=="W0","value"]
        if (pathogenic.fungi[dis] %in% 
          c("Candida_albicans","Candida_parapsilosis","Saccharomyces","Candida_sake","Candida_tropicalis")) {num <- 500} else {num <- 200}
        tmp[,"value"][tmp[,"value"]> num] <- num
        tmp[,"value"][tmp[,"value"] < (-1)*num] <- (-1)*num
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
    anova.pvalue0 <- sprintf("%.2e", p_value)

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
    geom_smooth(aes(color=treatment,fill=treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.1,span=1.2)+
    # geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess", method.args = list(degree=0.5),se=TRUE, level = 0.5)+
    # geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = lm, formula = y ~ splines::bs(x, 1), se = TRUE, level = 0.5)+
    stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(0),label = "p.format") +
    theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
    labs(title = paste0(pathogenic.fungi[dis],"\n", "Flu vs Nystatin (ANOVA p: ", anova.pvalue0," | ","art p: ",art.pvalue0,
        "|\n","kruskal p: ", kruskal.pvalue0," | ","LocationTests p: ",LocationTests.pvalue0,
        "|\n","TukeyHSD p: ",TukeyHSD.pvalue1," | ","loess p: ",Ftestp_value,"|\n","pairwise. p: ",pairwise.t.pvalue1," | ",
        "LMM p: ",LMM.pvalue1,")"),y = "Δ")+NoLegend()
    message(dis)
    return(p2)
    })
plot <- CombinePlots(c(total_plots2),ncol=10)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints_merge_bacth1.png", plot=plot,width = 49, height = 20,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_dyn.delta.timepoints_merge_bacth1.svg", plot=plot,width = 49, height = 20,dpi=300)
~~~

![Fig2.fungi_dyn.delta.timepoints_merge_bacth1](./Chapter2.assets/Fig2.fungi_dyn.delta.timepoints_merge_bacth1.png)

## 2.1.5 Differential Abundance Testing (LEfSe) per Time Point

Runs LEfSe differential-abundance analysis comparing each post-treatment time point (W1, W2, W4, W8) against W0 and saves the results (`MMC.ITS.DAA_all.W*.raw.lefser.v2.rds`).

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
colSums(MMC.ITS.counts1[,sort(grep("W0",colnames(MMC.ITS.counts1),value=FALSE))])
table(MMC.ITS.samples[MMC.ITS.samples$time=="W0","treatment"])
sort(unique(MMC.ITS.taxa1$Species))
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
colData(Fungal.tse)$library_size <- colSums(assay(Fungal.tse, "counts"))
colData(Fungal.tse)$Diagnosis_new <- factor(colData(Fungal.tse)$Diagnosis_new,levels=c("UC","CD"))
Sel_type <- c("Class","Order","Family","Genus","Species")
lefser.All.W1 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Fungal.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
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
mcsaveRDS(DAA_all.W1,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W1.raw.lefser.v2.rds")

lefser.All.W2 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Fungal.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
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
mcsaveRDS(DAA_all.W2,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W2.raw.lefser.v2.rds")

lefser.All.W4 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Fungal.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
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
mcsaveRDS(DAA_all.W4,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W4.raw.lefser.v2.rds")

lefser.All.W8 <- future_lapply(1:length(Sel_type),function(x) {
    tse_tmp <- subsetByPrevalentFeatures(Fungal.tse,rank = Sel_type[x],detection = 0,prevalence = 0,as_relative = FALSE)
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
mcsaveRDS(DAA_all.W8,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W8.raw.lefser.v2.rds")
~~~

## 2.1.6 Compile Significant Fungi and Variance Summary

Aggregates the significant differential fungi across time points and summarizes their contribution / proportion of variance, exported as `Fig2.fungi_prop_var.svg`.

~~~R
W1.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W1.raw.lefser.v2.rds")
W2.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W2.raw.lefser.v2.rds")
W4.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W4.raw.lefser.v2.rds")
W8.lefser <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.DAA_all.W8.raw.lefser.v2.rds")
names(W1.lefser) <- c("Nystatin","Fluconazole")
library(lefser)
JCM_fungal <- c("Candida albicans", "Candida parapsilosis", "Candida tropicalis", "Candida glabrata", "Candida krusei", "Candida lusitaniae", "Candida dubliniensis", "Candida guilliermondii", "Candida rugosa", 
"Trichosporon asahii","Trichosporon mucoides","Trichosporon asteroides","Trichosporon cutaneum","Trichosporon ovoides","Trichosporon inkin",
"Rhodotorula glutinis", "Rhodotorula mucilaginosa", "Rhodotorula rubra", "Rhodotorula minuta",
"Geotrichum capitatum","Aspergillus fumigatus", "Aspergillus terreus","Zygomycetes","Fusarium","Acremonium","Scedosporium","Paecilomyces","Trichoderma longibrachiatum")
WHO_fungal <- c("Cryptococcus neoformans","Candida auris","Aspergillus fumigatus","Candida albicans",
"Candida glabrata","Histoplasma","Eumycetoma","Mucorales","Fusarium","Candida tropicalis","Candida parapsilosis",
"Scedosporium","Lomentospora prolificans","Coccidioides","Candida krusei","Cryptococcus gattii","Talaromyces mamneffei","Pneumocystis jirovecii","Paracoccidioides")
pathogenic.fungi <- unique(c(JCM_fungal,WHO_fungal))
gut_fungi <- c("Candida albicans","Candida tropicalis","Candida zeylanoides","Candida sake","Candida inconspicua",
  "Candida glaebosa","Candida pseudoglaebosa","Candida tetrigidarum","Candida saitoana","Kazachstania bulder",
  "Malassezia furfur","Trichosporon asahii","Yarrowia lipolytica","Mucor circinelloides","Schizophyllum commune",
  "Apiotrichum domesticum")
gut_associated_fungi <- c("Candida albicans","Candida tropicalis","Candida inconspicua","Candida glaebosa",
  "Candida pseudoglaebosa","Candida sake","Candida zeylanoides","Candida saitoana","Malassezia furfur",
  "Trichosporon asahii","Yarrowia lipolytica","Mucor circinelloides")
gut_associated_fungi <- unique(c(gut_associated_fungi,gut_fungi,pathogenic.fungi))
gut_associated_fungi <- setdiff(gut_associated_fungi,c("Candida saitoana"))
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
    tmp_all <- tmp_all[tmp_all$Names %in% gut_associated_fungi,]
    tmp_all$LDA.scores[tmp_all$LDA.scores > 4] <- 4
    tmp_all$LDA.scores[tmp_all$LDA.scores < -4] <- -4
    tmp_all <- tmp_all[tmp_all$pvalue <= 1,]
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all[tmp_all$Names=="Candida albicans",]
    tmp_all[tmp_all$Names=="Candida parapsilosis",]
    tmp_all[tmp_all$Names=="Malassezia furfur",]

    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.18,]
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.4,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]
    o <- as.character(tmp_all_c$Names)
    tmp_all <- tmp_all[tmp_all$Names %in% o,]

    tmp_all_c[tmp_all_c$Names=="Candida albicans",]
    tmp_all_c[tmp_all_c$Names=="Candida parapsilosis",]

    tmp_all$Names <- factor(tmp_all$Names,levels=o[length(o):1])
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
plot <- CombinePlots(c(All_plots),ncol=2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig2.fungi_prop_var.svg", plot=plot,width = 10, height = 6,dpi=300)
~~~

![image-20260612181831356](./Chapter2.assets/image-20260612181831356.png)

## 2.1.7 High-variance Fungi and Prevalence Analysis

Assembles the differential fungi across time points, identifies high-variance fungal taxa (`HVFungi`) and decreasing taxa (`dec_all`), and visualizes their prevalence, exported as `Fig2.prevelance.svg`.

~~~R

names(W1.lefser) <- c("Nystatin","Fluconazole")
DiffFungi_ <- lapply(1:length(W1.lefser),function(x) {
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
    tmp_all <- tmp_all[tmp_all$Names %in% gut_associated_fungi,]
    tmp_all$LDA.scores[tmp_all$LDA.scores > 4] <- 4
    tmp_all$LDA.scores[tmp_all$LDA.scores < -4] <- -4
    tmp_all <- tmp_all[tmp_all$pvalue <= 1,]
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all[tmp_all$Names=="Candida albicans",]
    tmp_all[tmp_all$Names=="Candida parapsilosis",]
    tmp_all[tmp_all$Names=="Malassezia furfur",]

    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.18,]
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.4,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]
    o <- as.character(tmp_all_c$Names)
    tmp_all <- tmp_all[tmp_all$Names %in% o,]
    tmp_all$treatment <- names(W1.lefser)[x]
    return(tmp_all)
    })
DiffFungi <- do.call(rbind,DiffFungi_)
DiffFungi <- DiffFungi[DiffFungi$type %in% c("Genus","Species"),]
write.csv(DiffFungi,"./projects/MMC/Figures/figures_making/v4/DiffFungi_All.csv")



HVFungi <- lapply(1:length(W2.lefser),function(x) {
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
    tmp_all <- tmp_all[tmp_all$pvalue <= 1,]
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>=0.4,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]
    return(tmp_all_c)
    })
mcsaveRDS(HVFungi,"./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.HVFungi.lefser.Genus.v2.rds")



names(W2.lefser) <- c("Nystatin","Fluconazole")
dec_all <- lapply(1:length(W2.lefser),function(x) {
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
    tmp_all <- tmp_all[tmp_all$pvalue <= 1,]
    tmp_all_c <- reshape2::dcast(tmp_all,Names~group,value.var = "LDA.scores")
    tmp_all_c[is.na(tmp_all_c)] <- 0
    tmp_all_c$weight.mean <- (tmp_all_c$W1.vs.W0+tmp_all_c$W2.vs.W0+tmp_all_c$W4.vs.W0+tmp_all_c$W8.vs.W0)/4
    sd <- rowSds(as.matrix(tmp_all_c[,c("W1.vs.W0","W2.vs.W0","W4.vs.W0","W8.vs.W0")]))
    sd[sd==0] <- 0.1
    tmp_all_c$weight <- tmp_all_c$weight.mean/sd
    tmp_all_c <- tmp_all_c[abs(tmp_all_c$weight)>0,]
    tmp_all_c <- tmp_all_c[order(tmp_all_c$weight.mean,tmp_all_c$W1.vs.W0,tmp_all_c$W2.vs.W0,tmp_all_c$W4.vs.W0,tmp_all_c$W8.vs.W0),]

    DN0 <- tmp_all_c$Names[tmp_all_c$W1.vs.W0 <0 & tmp_all_c$W2.vs.W0 <0]
    DN01 <- tmp_all_c$Names[tmp_all_c$W1.vs.W0 <0 & tmp_all_c$W4.vs.W0 <0]
    DN1 <- tmp_all_c$Names[tmp_all_c$W2.vs.W0 <0 & tmp_all_c$W4.vs.W0 <0]
    DN2 <- tmp_all_c$Names[tmp_all_c$W2.vs.W0 <0 & tmp_all_c$W4.vs.W0==0 ]
    DN3 <- tmp_all_c$Names[tmp_all_c$W2.vs.W0==0 & tmp_all_c$W4.vs.W0 <0]
    dec <- sort(unique(c(DN0,DN01,DN1,DN2,DN3)))
    return(dec)
    })
mcsaveRDS(dec_all[[2]],"./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.dec_all.lefser.All.v2.rds")

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
require(rJava)
require(xlsx)
pathogenic.Fungi_all <- read.xlsx("./projects/MMC/Human.associated.Fungi.xlsx", sheetName = "raw")
pathogenic.Fungi <- pathogenic.Fungi_all[pathogenic.Fungi_all$group=="Fungal Pathogens",]
unique(pathogenic.Fungi$names)
Genus <- unique(unlist(lapply(strsplit(unique(pathogenic.Fungi$names),split=" "),function(x) {x[1]})))
relabundance_sel.fungi1_ <- lapply(1:length(Genus),function(x) {
    relativeAb_Genus <- relativeAb_all[relativeAb_all$names %in% grep(Genus[x],relativeAb_all$names,value=TRUE),]
    return(relativeAb_Genus)
    })
relabundance_sel.fungi1 <- do.call(rbind,relabundance_sel.fungi1_)
relabundance_sel.fungi1 <- relabundance_sel.fungi1[-grep(" NA",relabundance_sel.fungi1$names,value=FALSE),]
pathogenic.Fungi2 <- pathogenic.Fungi_all[pathogenic.Fungi_all$group=="human gut associated",]
relabundance_sel.fungi2 <- relativeAb_all[relativeAb_all$names %in% setdiff(pathogenic.Fungi2$names,relabundance_sel.fungi1$names),]
relabundance_sel.fungi <- rbind(relabundance_sel.fungi1,relabundance_sel.fungi2)
unique(relabundance_sel.fungi$names)
Species <- relabundance_sel.fungi[relabundance_sel.fungi$type=="Species","names"]
dec_all <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.dec_all.lefser.All.v2.rds")
Genus <- intersect(Genus,dec_all)
Species <- intersect(Species,dec_all)

library(lefser)
JCM_fungal <- c("Candida albicans", "Candida parapsilosis", "Candida tropicalis", "Candida glabrata", "Candida krusei", "Candida lusitaniae", "Candida dubliniensis", "Candida guilliermondii", "Candida rugosa", 
"Trichosporon asahii","Trichosporon mucoides","Trichosporon asteroides","Trichosporon cutaneum","Trichosporon ovoides","Trichosporon inkin",
"Rhodotorula glutinis", "Rhodotorula mucilaginosa", "Rhodotorula rubra", "Rhodotorula minuta",
"Geotrichum capitatum","Aspergillus fumigatus", "Aspergillus terreus","Zygomycetes","Fusarium","Acremonium","Scedosporium","Paecilomyces","Trichoderma longibrachiatum")
WHO_fungal <- c("Cryptococcus neoformans","Candida auris","Aspergillus fumigatus","Candida albicans",
"Candida glabrata","Histoplasma","Eumycetoma","Mucorales","Fusarium","Candida tropicalis","Candida parapsilosis",
"Scedosporium","Lomentospora prolificans","Coccidioides","Candida krusei","Cryptococcus gattii","Talaromyces mamneffei","Pneumocystis jirovecii","Paracoccidioides")
pathogenic.fungi <- unique(c(JCM_fungal,WHO_fungal))
gut_fungi <- c("Candida albicans","Candida tropicalis","Candida zeylanoides","Candida sake","Candida inconspicua",
  "Candida glaebosa","Candida pseudoglaebosa","Candida tetrigidarum","Candida saitoana","Kazachstania bulder",
  "Malassezia furfur","Trichosporon asahii","Yarrowia lipolytica","Mucor circinelloides","Schizophyllum commune",
  "Apiotrichum domesticum")
gut_associated_fungi <- c("Candida albicans","Candida tropicalis","Candida inconspicua","Candida glaebosa",
  "Candida pseudoglaebosa","Candida sake","Candida zeylanoides","Candida saitoana","Malassezia furfur",
  "Trichosporon asahii","Yarrowia lipolytica","Mucor circinelloides")
gut_associated_fungi <- unique(c(gut_associated_fungi,gut_fungi,pathogenic.fungi))
gut_associated_fungi_genus <- unlist(lapply(strsplit(gut_associated_fungi,split=" "),function(x) {return(x[1])}))
gut_associated_fungi <- unique(c(gut_associated_fungi,gut_associated_fungi_genus))
HVFungi <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.ITS.HVFungi.lefser.Genus.v2.rds")
gut_associated_fungi1 <- intersect(HVFungi[[2]][HVFungi[[2]]$weight < 0,]$Names,gut_associated_fungi)
gut_associated_fungi1 <- intersect(dec_all[[2]],gut_associated_fungi)

Genus1 <- intersect(HVFungi[[2]][HVFungi[[2]]$weight < 0,]$Names,Genus)

sel.time <- c("W0","W1","W2","W4","W8")
sel.treat <- c("Nystatin","Fluconazole")
Sel.type <- c("Phylum","Class","Order","Family","Genus","Species")
Prevalence_all_time_ <- future_lapply(1:length(sel.time),function(tmp) {
    tse_tmp <- Fungal.tse[,Fungal.tse$time==sel.time[tmp]]
    Prevalence_all_tmp_ <- future_lapply(1:length(Sel.type),function(i) {
        Prevalence_all_tmp_tmp_ <- lapply(1:length(sel.treat),function(x) {
            tse.tmp <- tse_tmp[,tse_tmp$treatment %in% sel.treat[x]]
            altExp(tse.tmp,Sel.type[i]) <- mergeFeaturesByRank(tse.tmp, Sel.type[i])
            Prevalence <- as.data.frame(getPrevalence(altExp(tse.tmp,Sel.type[i]), detection = 0.28, sort = FALSE,assay.type = "counts", as_relative = TRUE))
            colnames(Prevalence) <- "Prevalence"
            Prevalence$Fungi <- rownames(Prevalence)
            Prevalence$treatment <- sel.treat[x]
            return(Prevalence)
            })
        Prevalence_all_tmp_tmp <- do.call(rbind,Prevalence_all_tmp_tmp_)
        Prevalence_tmp <- reshape2::dcast(Prevalence_all_tmp_tmp,Fungi~treatment,value.var = "Prevalence")
        Prevalence_tmp$type <- Sel.type[i]
        return(Prevalence_tmp)
        })
    Prevalence_all_tmp <- do.call(rbind,Prevalence_all_tmp_)
    Prevalence_all_tmp$time <- sel.time[tmp]
    return(Prevalence_all_tmp)
    })
Prevalence_all_time_high_pre <- do.call(rbind,Prevalence_all_time_)
Prevalence_all_time_high_pre[Prevalence_all_time_high_pre$Fungi=="Candida",]
Prevalence_all_time_high_pre[Prevalence_all_time_high_pre$Fungi=="Fusarium",]

Prevalence_all_time_ <- future_lapply(1:length(sel.time),function(tmp) {
    tse_tmp <- Fungal.tse[,Fungal.tse$time==sel.time[tmp]]
    Prevalence_all_tmp_ <- future_lapply(1:length(Sel.type),function(i) {
        Prevalence_all_tmp_tmp_ <- lapply(1:length(sel.treat),function(x) {
            tse.tmp <- tse_tmp[,tse_tmp$treatment %in% sel.treat[x]]
            altExp(tse.tmp,Sel.type[i]) <- mergeFeaturesByRank(tse.tmp, Sel.type[i])
            Prevalence <- as.data.frame(getPrevalence(altExp(tse.tmp,Sel.type[i]), detection = 1/100, sort = FALSE,assay.type = "counts", as_relative = TRUE))
            colnames(Prevalence) <- "Prevalence"
            Prevalence$Fungi <- rownames(Prevalence)
            Prevalence$treatment <- sel.treat[x]
            return(Prevalence)
            })
        Prevalence_all_tmp_tmp <- do.call(rbind,Prevalence_all_tmp_tmp_)
        Prevalence_tmp <- reshape2::dcast(Prevalence_all_tmp_tmp,Fungi~treatment,value.var = "Prevalence")
        Prevalence_tmp$type <- Sel.type[i]
        return(Prevalence_tmp)
        })
    Prevalence_all_tmp <- do.call(rbind,Prevalence_all_tmp_)
    Prevalence_all_tmp$time <- sel.time[tmp]
    return(Prevalence_all_tmp)
    })
Prevalence_all_time_low_pre <- do.call(rbind,Prevalence_all_time_)
Prevalence_all_time_low_pre_W0 <- Prevalence_all_time_low_pre[Prevalence_all_time_low_pre$time=="W0",]
Prevalence_all_time_low_pre_W0 <- Prevalence_all_time_low_pre_W0[Prevalence_all_time_low_pre_W0$Fluconazole<0.5 | Prevalence_all_time_low_pre_W0$Nystatin<0.5,]
Prevalence_all_time_low_pre1 <- Prevalence_all_time_low_pre[Prevalence_all_time_low_pre$Fungi %in% Prevalence_all_time_low_pre_W0$Fungi,]
Prevalence_all_time_high_pre1 <- Prevalence_all_time_high_pre[Prevalence_all_time_high_pre$Fungi %in% setdiff(Prevalence_all_time_high_pre$Fungi,Prevalence_all_time_low_pre1$Fungi),]
Prevalence_all_time <- rbind(Prevalence_all_time_high_pre1,Prevalence_all_time_low_pre1)


Prevalence_all_time_high_pre[Prevalence_all_time_high_pre$Fungi=="Rhodotorula",]
Prevalence_all_time_low_pre[Prevalence_all_time_low_pre$Fungi=="Rhodotorula",]


index = c("Nystatin","Fluconazole")
label_genus <- lapply(1:length(index),function(x) {
    Prevalence_tmp <- Prevalence_all_time
    Prevalence_tmp <- Prevalence_tmp[,c("Fungi","time",index[x])]
    Prevalence <- reshape2::dcast(Prevalence_tmp,Fungi~time,value.var = index[x])
    Prevalence[Prevalence$Fungi=="Candida",]
    Prevalence$W2_W0 <- Prevalence$W2-Prevalence$W0
    Prevalence$W0_W2 <- Prevalence$W0-Prevalence$W2
    Prevalence$W2_W0[Prevalence$W2_W0>0.2] <- 0.2
    Prevalence$W2_W0[Prevalence$W2_W0< -0.2] <- -0.2
    Prevalence$dynamic <- "No_Changed"
    Prevalence$dynamic[Prevalence$W2_W0 < 0] <- "W0_enriched"
    Prevalence$dynamic[Prevalence$W2_W0 >0] <- "W2_enriched"
    Prevalence <- Prevalence[order(Prevalence$dynamic,decreasing=FALSE),]
    Prevalence1 <- Prevalence[Prevalence$dynamic!="No_Changed",]
    Prevalence1 <- Prevalence1[Prevalence1$Fungi %in% Genus,]
    return(Prevalence1$Fungi)
})
label_genus <- intersect(label_genus[[1]],label_genus[[2]])

color <- c("#ff4757", "lightgrey","#3742fa")
names(color) <- c("W2_enriched","No_Changed","W0_enriched")
index = c("Nystatin","Fluconazole")
All.plot0 <- lapply(1:length(index),function(x) {
    Prevalence_tmp <- Prevalence_all_time
    Prevalence_tmp <- Prevalence_tmp[,c("Fungi","time",index[x])]
    Prevalence <- reshape2::dcast(Prevalence_tmp,Fungi~time,value.var = index[x])
    Prevalence[Prevalence$Fungi=="Candida",]
    Prevalence$W2_W0 <- Prevalence$W2-Prevalence$W0
    Prevalence$W0_W2 <- Prevalence$W0-Prevalence$W2
    Prevalence$W2_W0[Prevalence$W2_W0>0.3] <- 0.3
    Prevalence$W2_W0[Prevalence$W2_W0< -0.3] <- -0.3
    Prevalence$dynamic <- "No_Changed"
    Prevalence$dynamic[Prevalence$W2_W0 <  0] <- "W0_enriched"
    Prevalence$dynamic[Prevalence$W2_W0 > 0] <- "W2_enriched"
    Prevalence <- Prevalence[order(Prevalence$dynamic,decreasing=FALSE),]
    Prevalence1 <- Prevalence[Prevalence$dynamic!="No_Changed",]
    aa <- jdb_palette("brewer_spectra")[1:length(jdb_palette("brewer_spectra"))]
    p1 <- ggplot(Prevalence, aes(x=W0, y=W2, color=W2_W0,size=abs(W0_W2))) + geom_point(alpha=1)+ scale_size(range = c(2, 6))+
    scale_colour_gradientn(colours = colorRampPalette(aa)(100))+
    theme_classic() +  ylab("Prevalence (W2)") + xlab("Prevalence (W0)")+ 
    labs(title=paste0("Prevalence of fungi genus ",index[x]))+
    geom_abline(slope=1, intercept = 0,linetype=2,color="lightgrey",size=1.5)+
    ggnewscale::new_scale_color()+ggrepel::geom_label_repel(data = Prevalence1[Prevalence1$Fungi %in% label_genus,], 
        aes(label = Fungi,color=dynamic), segment.color = 'black', show.legend = FALSE,size=3,max.overlaps=50,
        min.segment.length = unit(0.1, 'lines'))+ scale_color_manual(values=color)+xlim(0,1)+ylim(0,1)
    return(p1)
})
plot <- CombinePlots(c(All.plot0),ncol=2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig2.prevelance.svg", plot=plot,width = 10, height = 4,dpi=300)
~~~

![image-20260615084429630](./Chapter2.assets/image-20260615084429630.png)

## 2.1.8 Rebuild Fungal TSE for Beta-diversity

Reconstructs the fungal `TreeSummarizedExperiment` (relative abundance + QC) as input for the beta-diversity ordination that follows.

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
# MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("NA",MMC.ITS.taxa1$Species,value=FALSE),]
if (length(grep("Incertae",MMC.ITS.taxa1$Genus,value=TRUE))>0){MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("Incertae",MMC.ITS.taxa1$Genus,value=FALSE),]} else {MMC.ITS.taxa1 <- MMC.ITS.taxa1}
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
colSums(MMC.ITS.counts1[,sort(grep("W0",colnames(MMC.ITS.counts1),value=FALSE))])

sort(unique(MMC.ITS.taxa1$Species))
Fungal.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.ITS.counts1)),colData = DataFrame(MMC.ITS.samples),rowData = DataFrame(MMC.ITS.taxa11))
Fungal.tse <- transformAssay(Fungal.tse_taxa, MARGIN = "samples", method = "relabundance")
Fungal.tse <- addPerCellQC(Fungal.tse)
colData(Fungal.tse)$total_raw_counts <- colSums(assay(Fungal.tse, "counts"))
~~~

## 2.1.9 Beta-diversity: PERMANOVA and PCoA

Loads per-treatment PERMANOVA results and draws the PCoA ordination of fungal community composition, exported as `Fig2.1.Pcoa.svg`.

~~~R

All_sig.Flu <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC_ITS_PERMANOVA_res.Fluconazole.rds")
All_sig.Nys <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC_ITS_PERMANOVA_res.Nystatin.rds")
All_sig.Flu[All_sig.Flu$pvl_jac <= 0.05 | All_sig.Flu$pvl_bray <= 0.05,]
All_sig.Nys[All_sig.Nys$pvl_jac <= 0.05 | All_sig.Nys$pvl_bray <= 0.05,]
All_sig.Nys[All_sig.Nys$prevalence==0.5,]
All_sig.Flu[All_sig.Flu$prevalence==0.5,]

library(ggside)
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggside)
library(patchwork)
generate_ellipse <- function(median_x, mad_x, median_y, mad_y, n = 500, factor) {
  t <- seq(0, 2*pi, length.out = n)
  a <- mad_x * factor
  b <- mad_y * factor
  x <- median_x + a * cos(t)
  y <- median_y + b * sin(t)
  return(data.frame(x = x, y = y))
}
treatment <- c("Nystatin", "Fluconazole")
All_sig <- list(All_sig.Nys,All_sig.Flu)
pal <- jdb_palette("corona")
Sel_type <- c("Genus","Species")
ALL_PLOTS <- lapply(1:length(treatment), function(treat) {
  tse_tmp <- Fungal.tse[, Fungal.tse$treatment %in% treatment[treat]]
  tse_tmp <- tse_tmp[, tse_tmp$time %in% c("W0","W2","W4","W1","W8")]
  tse_tmp$time.v2 <- as.character(tse_tmp$time)
  tse_tmp$time.v2[tse_tmp$time.v2 %in% "W0"] <- "Pre"
  tse_tmp$time.v2[tse_tmp$time.v2 %in% c("W2","W4","W1")] <- "Post"
  tse_tmp$time.v2[tse_tmp$time.v2 %in% "W8"] <- "LTM"
  tse_tmp$time.v2 <- factor(tse_tmp$time.v2, levels=c("Pre","Post","LTM"))
  ## 转换为相对丰度
  mat_rel_ <- lapply(1:length(Sel_type),function(t) {
    tse_tmp1 <- subsetByPrevalentFeatures(tse_tmp, rank = Sel_type[t], detection = 0.25, prevalence = 0.25)
    ## 相对丰度 + Jaccard（presence/absence）
    tse_tmp1 <- transformAssay(tse_tmp1, assay.type = "counts", method = "relabundance")
    mat_rel1 <- t(assay(tse_tmp1, "relabundance"))
    return(mat_rel1)
    })
  mat_rel <- as.data.frame(do.call(cbind,mat_rel_))

  dist_jac <- vegdist(mat_rel > 0, method = "jaccard")
  ## PCoA
  pcoa_res <- pcoa(dist_jac)
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
  All_sig_tmp <- All_sig[[treat]][48,]
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
ggsave("./projects/MMC/Figures/v2_figures/Fig2.1.Pcoa.svg", plot=plot,width = 10, height = 4,dpi=300)
~~~

![image-20260615091654521](./Chapter2.assets/image-20260615091654521.png)

## 2.1.10 Fungal Alpha-diversity Over Time

Computes fungal alpha-diversity metrics and plots their changes across treatment time points, exported as `Fig2.1.fungal_diversity.svg`.

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
# MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("NA",MMC.ITS.taxa1$Species,value=FALSE),]
if (length(grep("Incertae",MMC.ITS.taxa1$Genus,value=TRUE))>0){MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("Incertae",MMC.ITS.taxa1$Genus,value=FALSE),]} else {MMC.ITS.taxa1 <- MMC.ITS.taxa1}
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
colSums(MMC.ITS.counts1[,sort(grep("W0",colnames(MMC.ITS.counts1),value=FALSE))])

sort(unique(MMC.ITS.taxa1$Species))
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

tmp_projects2 <- as.data.frame(colData(Fungal.tse))
relabundance_tmp <- relabundance_pathogenic.Fungi.v2[,rownames(tmp_projects2)]
tmp_projects2$pathogenic.Fungi.v2 <- 100*colMeans(relabundance_tmp)[rownames(tmp_projects2)]
tmp_projects2$time.v2 <- as.character(tmp_projects2$time)
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% "W0"] <- "Pre"
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% c("W1","W2","W4")] <- "Post"
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% "W8"] <- "LTM"
tmp_projects2$time.v2 <- factor(tmp_projects2$time.v2, levels=c("Pre","Post","LTM"))
tmp_projects2$Fungi.observed <- tmp_projects2$observed
tmp_projects2$Fungi.shannon <- tmp_projects2$shannon

diease_Scores <- c("Fungi.observed","Fungi.shannon","pathogenic.Fungi.v2")
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
pal <- jdb_palette("corona")[c(2,1)]
total_plots2 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- tmp_projects2[!is.na(tmp_projects2$treatment),]
    df_paired <- df_paired[df_paired$time %in% c("W0","W1","W2","W4","W8"),]
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,diease_Scores[dis]]-tmp[tmp$time=="W0",diease_Scores[dis]]
        if (diease_Scores[dis]=="pathogenic.Fungi.v2") {
            tmp$value[tmp$value > 500] <- 500
            tmp$value[tmp$value < -500] <- -500
        }
        if (diease_Scores[dis]=="Fungi.observed") {
            tmp$value[tmp$value > 40] <- 40
            tmp$value[tmp$value < -40] <- -40
        }
        if (diease_Scores[dis]=="Fungi.shannon") {
            tmp$value[tmp$value > 1.5] <- 1.5
            tmp$value[tmp$value < -1.5] <- -1.5
        }
        return(tmp)
    })
    df_paired1 <- do.call(rbind,df_paired1_)
    df_paired1$treatment <- factor(df_paired1$treatment,levels=c("Nystatin","Fluconazole"))
    if (diease_Scores[dis]=="Fungi.observed") {test <- wilcox.test_wrapper} else {test <- wilcox.test_wrapper}
    df_paired2 <- df_paired1
    df_paired2$time <- as.numeric(gsub("W","",df_paired2$time))
    loess_fit_1 <- loess(value ~ time, data = df_paired2[df_paired2$treatment == "Nystatin", ])
    loess_fit_2 <- loess(value ~ time, data = df_paired2[df_paired2$treatment == "Fluconazole", ])
    x_range <- seq(min(df_paired2$time), max(df_paired2$time), length.out = 100)
    pred_1 <- predict(loess_fit_1, newdata = data.frame(time = x_range))
    pred_2 <- predict(loess_fit_2, newdata = data.frame(time = x_range))
    ks_test_result <- formatC(ks.test(pred_1, pred_2)$p.value, format = "e", digits = 3)

    rss_1 <- sum((df_paired2[df_paired2$treatment == "Nystatin", "value"] - predict(loess_fit_1))^2)
    rss_2 <- sum((df_paired2[df_paired2$treatment == "Fluconazole", "value"] - predict(loess_fit_2))^2)
    n1 <- length(df_paired2[df_paired2$treatment == "Nystatin", "value"])  # Sample size group 1
    n2 <- length(df_paired2[df_paired2$treatment == "Fluconazole", "value"])  # Sample size group 2
    f_stat <- (rss_1 / (n1 - 2)) / (rss_2 / (n2 - 2))
    Ftestp_value <- pf(f_stat, df1 = n1 - 2, df2 = n2 - 2, lower.tail = FALSE)# Perform an F-test

    p1 <- ggplot(df_paired1, aes_string(x = "time", y = "value")) + 
        geom_line(aes(group=patient),size = 0.6,alpha=0.8,color="lightgrey") +
        geom_jitter(color="black",width = 0.1,alpha=0.5,size=1)+
        geom_boxplot(outlier.shape = NA,alpha=0)+
        theme_bw()+ scale_fill_manual(values = pal,guide="none")+scale_color_manual(values = pal,guide="none")+theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
        stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=1)))+facet_wrap(~treatment,ncol=3)+
        labs(title=paste0(diease_Scores[dis],"\n (Flu vs Nys p:",Ftestp_value,")"),y = paste("Δ")) +NoLegend()+
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = test)+
        geom_smooth(aes(group = 1, color = treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.2,span=1.2)+
        stat_summary(fun = mean, geom = "point",aes(group = 1,color = treatment),size=2)
        # stat_summary(fun = median, geom = "line",aes(group = 1,color = treatment),size=1.5)
    plot <- plot_grid(p1,ncol=1)
    return(plot)
    })
plot <- CombinePlots(c(total_plots2),nrow=1)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig2.1.fungal_diversity.svg", plot=plot,width = 10, height = 4,dpi=300)
~~~

![image-20260615091743105](./Chapter2.assets/image-20260615091743105.png)

## 2.1.11 Normalize Pathogenic Fungi and Richness

Min–max normalizes pathogenic-fungi load and observed fungal richness within each treatment arm, preparing the indices for the density comparisons below.

~~~R
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% "W0"] <- "Pre"
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% c("W1","W2","W4")] <- "Post"
tmp_projects2$time.v2[tmp_projects2$time.v2 %in% "W8"] <- "LTM"
tmp_projects2$time.v2 <- factor(tmp_projects2$time.v2, levels=c("Pre","Post","LTM"))
aa <- jdb_palette("brewer_spectra", type = "continuous")
ITS.Info.global4 <- tmp_projects2[tmp_projects2$treatment %in% c("Nystatin","Fluconazole"),]
ITS.Info.global4$treatment <- factor(ITS.Info.global4$treatment,levels=c("Nystatin","Fluconazole"))
# ITS.Info.global4$pathogenic.Fungi.v2 <- (ITS.Info.global4$pathogenic.Fungi.v2-min(ITS.Info.global4$pathogenic.Fungi.v2, na.rm = TRUE))/(max(ITS.Info.global4$pathogenic.Fungi.v2, na.rm = TRUE)-min(ITS.Info.global4$pathogenic.Fungi.v2, na.rm = TRUE))
# ITS.Info.global4$Fungi.observed <- (ITS.Info.global4$Fungi.observed-min(ITS.Info.global4$Fungi.observed, na.rm = TRUE))/(max(ITS.Info.global4$Fungi.observed, na.rm = TRUE)-min(ITS.Info.global4$Fungi.observed, na.rm = TRUE))
sel_treat <- c("Nystatin","Fluconazole")
ITS.Info.global5_ <- lapply(1:length(sel_treat),function(x) {
  tmp <- ITS.Info.global4[ITS.Info.global4$treatment==sel_treat[x],]
  tmp$pathogenic.Fungi.v2 <- (tmp$pathogenic.Fungi.v2-min(tmp$pathogenic.Fungi.v2, na.rm = TRUE))/(max(tmp$pathogenic.Fungi.v2, na.rm = TRUE)-min(tmp$pathogenic.Fungi.v2, na.rm = TRUE))
  tmp$Fungi.observed <- (tmp$Fungi.observed-min(tmp$Fungi.observed, na.rm = TRUE))/(max(tmp$Fungi.observed, na.rm = TRUE)-min(tmp$Fungi.observed, na.rm = TRUE))
  tmp
  })
ITS.Info.global5 <- do.call(rbind,ITS.Info.global5_)

# ITS.Info.global5 <- ITS.Info.global4
# ITS.Info.global5$pathogenic.Fungi.v2 <- (ITS.Info.global5$pathogenic.Fungi.v2-min(ITS.Info.global5$pathogenic.Fungi.v2, na.rm = TRUE))/(max(ITS.Info.global5$pathogenic.Fungi.v2, na.rm = TRUE)-min(ITS.Info.global5$pathogenic.Fungi.v2, na.rm = TRUE))
# ITS.Info.global5$Fungi.observed <- (ITS.Info.global5$Fungi.observed-min(ITS.Info.global5$Fungi.observed, na.rm = TRUE))/(max(ITS.Info.global5$Fungi.observed, na.rm = TRUE)-min(ITS.Info.global5$Fungi.observed, na.rm = TRUE))
~~~

## 2.1.12 Pathogenic Fungi vs Richness Density (Exploratory)

Plots a 2D density map of normalized pathogenic-fungi load versus observed richness, faceted by treatment and phase, with a 1:1 reference line (exploratory view).

~~~R
ggplot(ITS.Info.global5, aes(x = pathogenic.Fungi.v2, y = Fungi.observed)) +  geom_point(alpha = 0.5, size = 0.1) +
stat_density_2d(geom = "tile", aes(fill = ..ndensity..), contour = FALSE, n = 500) + scale_fill_gradientn(colours = aa) +
facet_wrap(~treatment+time.v2, ncol = 3) +theme_classic() + labs(title = paste0("Cor"),x = "pathogenic.Fungi.v2",y = "Fungi.observed",fill = "Density") +theme(legend.position = "none")+
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+xlim(0,0.75)+ylim(0,0.75)
~~~

![image-20260615091902676](./Chapter2.assets/image-20260615091902676.png)

## 2.1.13 Pathogenic Fungi vs Richness Density (Final Figure)

Renders the final pathogenic-fungi vs richness density figure, exported as `Fig2.1.denisty.png`.

~~~R
plot <- ggplot(ITS.Info.global6, aes(x = pathogenic.Fungi.v2, y = Fungi.observed)) +  geom_point(alpha = 0.5, size = 0.1) +
stat_density_2d(geom = "tile", aes(fill = ..ndensity..), contour = FALSE, n = 500) + scale_fill_gradientn(colours = aa) +
facet_wrap(~treatment+time.v2, ncol = 3) +theme_classic() +xlim(0,0.75)+ylim(0,0.75)+ 
labs(title = paste0("Cor"),x = "pathogenic.Fungi.v2",y = "Fungi.observed",fill = "Density") +
theme(legend.position = "none")
ggsave("./projects/MMC/Figures/v2_figures/Fig2.1.denisty.png", plot=plot,width = 8, height = 6,dpi=300)
~~~

![image-20260615092037203](./Chapter2.assets/image-20260615092037203.png)

## 2.1.14 Diversity Index (Supplementary)

Computes and plots fungal diversity indices for the supplementary figure, exported as `S2.2.index.svg`.

~~~R
ReadCap.20251215 <- mcreadRDS("./workshop/MMC/Aidan_info/v2/ReadCap.20251215.rds")
total_times_tmp <- read.csv("./projects/MMC/Figures/figures_making/v3/patient.PFS.20251215.csv")
total_times_tmp <- total_times_tmp[total_times_tmp$treatment!="Clotrimazole",]
sample_info_total5 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% total_times_tmp$Omics_patient_Names,]
# sample_info_total5 <- sample_info_total5[sample_info_total5$time!="W6",]
sample_info_total5[sample_info_total5$Omics_patient_Names=="MMC008",c("Omics_patient_Names","Omics_samples_Names","DAI","treatment")]
sample_info_total5[sample_info_total5$Omics_patient_Names=="MMC107",c("Omics_patient_Names","Omics_samples_Names","DAI","treatment")]
sample_info_total5[sample_info_total5$Omics_patient_Names=="MMC063",c("Omics_patient_Names","Omics_samples_Names","DAI","treatment")]
sample_info_total5[sample_info_total5$Omics_patient_Names=="MMC099",c("Omics_patient_Names","Omics_samples_Names","DAI","treatment")]
sample_info_total5[sample_info_total5$Omics_patient_Names=="MMC110",c("Omics_patient_Names","Omics_samples_Names","DAI","treatment")]

MMC.ITS.counts <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/v2/MMC_ITS_ASV.counts.rds")
MMC.ITS.taxa1 <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/v2/MMC_ITS_taxa.rds")
MMC.ITS.samples <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/v2/MMC_ITS_samples_info.rds")
MMC.ITS.samples <- MMC.ITS.samples[MMC.ITS.samples$sample %in% sample_info_total5$Omics_samples_Names,]
# MMC.ITS.counts <- MMC.ITS.counts[,colSums(MMC.ITS.counts)>50]
MMC.ITS.samples <- MMC.ITS.samples[intersect(rownames(MMC.ITS.samples),colnames(MMC.ITS.counts)),]
MMC.ITS.counts <- MMC.ITS.counts[,rownames(MMC.ITS.samples)]
library(microeco)
library(mecodev)
MMC.ITS.counts1 <- MMC.ITS.counts
# MMC.ITS.taxa1 <- MMC.ITS.taxa1[!is.na(MMC.ITS.taxa1$Genus),]
# if (length(grep("Incertae",MMC.ITS.taxa1$Genus,value=TRUE))>0){MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("Incertae",MMC.ITS.taxa1$Genus,value=FALSE),]} else {MMC.ITS.taxa1 <- MMC.ITS.taxa1}
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
colSums(MMC.ITS.counts1[,sort(grep("W0",colnames(MMC.ITS.counts1),value=FALSE))])
sort(unique(MMC.ITS.taxa1$Species))
Fungal.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.ITS.counts1)),colData = DataFrame(MMC.ITS.samples),rowData = DataFrame(MMC.ITS.taxa11))
Fungal.tse <- transformAssay(Fungal.tse_taxa, MARGIN = "samples", method = "relabundance")
Fungal.tse <- addPerCellQC(Fungal.tse)
colData(Fungal.tse)$total_raw_counts <- colSums(assay(Fungal.tse, "counts"))
Fungal.tse_raw <- Fungal.tse
# 添加分组信息
cutoff <- 50
df <- as.data.frame(colData(Fungal.tse)) %>% arrange(sum) %>%
  mutate(index = 1:n(),group = ifelse(sum <= cutoff, "Filtered (<50)", "Kept (>=50)"))
n_filtered <- sum(df$sum <= cutoff)
n_total <- nrow(df)
df[df$sum <= cutoff,]
p1 <- ggplot(df, aes(x = sum, fill = group)) + geom_histogram(position = "stack", bins = 50, color = "black", alpha = 0.7) +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1) + scale_x_log10() +
  labs(x = "Library size (log10 scale)",y = "Count",
    title = paste0("ITS (", n_filtered, "/", n_total, " samples < ", cutoff, ")"),fill = "Group") + theme_bw()
p2 <- ggplot(df, aes(y = index, x = sum, color = group)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1) +
  scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = "Library size (log10 scale)", y = "Sample index",
       title = paste0("ITS (cutoff=", cutoff, ")"), color = "Group") +
  theme_bw()
plot <- plot_grid(p1, p2, ncol = 2)
plot
ggsave("./projects/MMC/Figures/v2_figures/S2.2.index.svg", plot=plot,width = 8, height = 5,dpi=300)
~~~

![image-20260615092112532](./Chapter2.assets/image-20260615092112532.png)

## 2.1.15 Rebuild Fungal TSE for Rarefaction QC

Reconstructs the fungal `TreeSummarizedExperiment` and inspects per-sample sequencing depth as preparation for the rarefaction analysis.

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
dim(MMC.ITS.counts)
MMC.ITS.counts <- MMC.ITS.counts[,colSums(MMC.ITS.counts)>50]
dim(MMC.ITS.counts)
MMC.ITS.samples <- MMC.ITS.samples[intersect(rownames(MMC.ITS.samples),colnames(MMC.ITS.counts)),]
MMC.ITS.counts <- MMC.ITS.counts[,rownames(MMC.ITS.samples)]
MMC.ITS.taxa1 <- MMC.ITS.taxa1[!is.na(MMC.ITS.taxa1$Genus),]
library(microeco)
library(mecodev)
MMC.ITS.counts1 <- MMC.ITS.counts[rowSums(MMC.ITS.counts)>0,colSums(MMC.ITS.counts)>0]
summary(colSums(MMC.ITS.counts1))
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
summary(colSums(MMC.ITS.counts1))
dim(MMC.ITS.counts1)
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
summary(colSums(MMC.ITS.counts1))
Fungal.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.ITS.counts1)),colData = DataFrame(MMC.ITS.samples),rowData = DataFrame(MMC.ITS.taxa11))
Fungal.tse <- transformAssay(Fungal.tse_taxa, MARGIN = "samples", method = "relabundance")
Fungal.tse <- addPerCellQC(Fungal.tse)
colData(Fungal.tse)$total_raw_counts <- colSums(assay(Fungal.tse, "counts"))
summary(colData(Fungal.tse)$total_raw_counts)
setdiff(colnames(Fungal.tse),colnames(MMC.ITS.counts1))
sort(colSums(MMC.ITS.counts1))
~~~

## 2.1.16 Rarefaction Curves

Builds a `microtable` and computes rarefaction curves of observed richness across sequencing depths, exported as `S2.1.rarefy.svg`.

~~~R
library(microeco)
library(mecodev)
dataset <- microtable$new(sample_table = MMC.ITS.samples, otu_table = MMC.ITS.counts1, tax_table = MMC.ITS.taxa11)
dataset$tidy_dataset()
t1 <- trans_rarefy$new(dataset, alphadiv = "Observed", depth = c(0, 10, 50, 500, 2000, 4000, 6000, 8000))
plot <- t1$plot_rarefy(color = "patient", show_point = FALSE, add_fitting = TRUE)+labs(title="ITS")
ggsave("./projects/MMC/Figures/v2_figures/S2.1.rarefy.svg", plot=plot,width = 8, height = 5,dpi=300)
~~~

![image-20260615092205922](./Chapter2.assets/image-20260615092205922.png)

## 2.1.17 Rarefied Richness Summary

Summarizes observed richness at a sequencing depth of 8,000 reads (mean ± SD / SE).

~~~R
df <- t1$res_rarefy
summary(df[df$seqnum==8000,"Observed"])
obs <- df[df$seqnum == 8000, "Observed"]
mean_val <- mean(obs, na.rm = TRUE)
sd_val <- sd(obs, na.rm = TRUE)
n <- sum(!is.na(obs))
se_val <- sd_val / sqrt(n)
cat("Mean ± SD:", mean_val, "±", sd_val, "\n")
cat("Mean ± SE:", mean_val, "±", se_val, "\n")
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   3.00   22.75   35.50   41.76   59.25  118.00
Mean ± SD: 41.75833 ± 24.48048
Mean ± SE: 41.75833 ± 2.234752
~~~

## 2.1.18 Library Size (Supplementary)

Computes and plots per-sample library sizes for the supplementary figure, exported as `S2.2.Library.svg`.

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
dim(MMC.ITS.counts)
MMC.ITS.counts <- MMC.ITS.counts[,colSums(MMC.ITS.counts)>50]
dim(MMC.ITS.counts)
MMC.ITS.samples <- MMC.ITS.samples[intersect(rownames(MMC.ITS.samples),colnames(MMC.ITS.counts)),]
MMC.ITS.counts <- MMC.ITS.counts[,rownames(MMC.ITS.samples)]
library(microeco)
library(mecodev)
MMC.ITS.counts1 <- MMC.ITS.counts[rowSums(MMC.ITS.counts)>0,colSums(MMC.ITS.counts)>0]
summary(colSums(MMC.ITS.counts1))
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
summary(colSums(MMC.ITS.counts1))
dim(MMC.ITS.counts1)
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
summary(colSums(MMC.ITS.counts1))
Fungal.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.ITS.counts1)),colData = DataFrame(MMC.ITS.samples),rowData = DataFrame(MMC.ITS.taxa11))
Fungal.tse <- transformAssay(Fungal.tse_taxa, MARGIN = "samples", method = "relabundance")
Fungal.tse <- addPerCellQC(Fungal.tse)
colData(Fungal.tse)$total_raw_counts <- colSums(assay(Fungal.tse, "counts"))
summary(colData(Fungal.tse)$total_raw_counts)

obs <- as.data.frame(colData(Fungal.tse))$sum
mean_val <- mean(obs, na.rm = TRUE)
sd_val <- sd(obs, na.rm = TRUE)
n <- sum(!is.na(obs))
se_val <- sd_val / sqrt(n)
cat("Mean ± SE:", mean_val, "±", se_val, "\n")
Mean ± SE: 82612.98 ± 8257.325

obs <- as.data.frame(colData(Fungal.tse)[colData(Fungal.tse)$treatment=="Nystatin",])$sum
mean_val <- mean(obs, na.rm = TRUE)
sd_val <- sd(obs, na.rm = TRUE)
n <- sum(!is.na(obs))
se_val <- sd_val / sqrt(n)
cat("Mean ± SE:", mean_val, "±", se_val, "\n")
Mean ± SE: 96066.19 ± 14669.21

obs <- as.data.frame(colData(Fungal.tse)[colData(Fungal.tse)$treatment=="Fluconazole",])$sum
mean_val <- mean(obs, na.rm = TRUE)
sd_val <- sd(obs, na.rm = TRUE)
n <- sum(!is.na(obs))
se_val <- sd_val / sqrt(n)
cat("Mean ± SE:", mean_val, "±", se_val, "\n")
Mean ± SE: 76053.16 ± 9974.55

Fungal.tse$log10sum <- log10(Fungal.tse$sum)
All_plots1 <- plotColData(Fungal.tse,"log10sum","time", colour_by = "time") + theme(axis.text.x = element_text(angle = 45, hjust=1))+labs(y = "Library size (N)", x = "time",title="ITS")+geom_hline(yintercept = log10(cutoff), color = "red", linetype = "dashed", size = 1)
All_plots2 <- plotColData(Fungal.tse,"log10sum","treatment", colour_by = "treatment") + theme(axis.text.x = element_text(angle = 45, hjust=1))+labs(y = "Library size (N)", x = "treatment",title="ITS")+geom_hline(yintercept = log10(cutoff), color = "red", linetype = "dashed", size = 1)
plot <- CombinePlots(list(All_plots1,All_plots2),nrow=1)
ggsave("./projects/MMC/Figures/v2_figures//S2.2.Library.svg", plot=plot,width = 5, height = 5,dpi=300
~~~

![image-20260615092335104](./Chapter2.assets/image-20260615092335104.png)

## 2.1.19 Baseline (W0) Fungal Abundance

Summarizes baseline (W0) fungal taxonomic abundance across patients, exported as `Fig2.W0abudan.svg`.

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
dim(MMC.ITS.counts)
MMC.ITS.counts <- MMC.ITS.counts[,colSums(MMC.ITS.counts)>50]
dim(MMC.ITS.counts)
MMC.ITS.samples <- MMC.ITS.samples[intersect(rownames(MMC.ITS.samples),colnames(MMC.ITS.counts)),]
MMC.ITS.counts <- MMC.ITS.counts[,rownames(MMC.ITS.samples)]
library(microeco)
library(mecodev)
MMC.ITS.counts1 <- MMC.ITS.counts[rowSums(MMC.ITS.counts)>0,colSums(MMC.ITS.counts)>0]
MMC.ITS.taxa1 <- MMC.ITS.taxa1[!is.na(MMC.ITS.taxa1$Genus),]
if (length(grep("Incertae",MMC.ITS.taxa1$Genus,value=TRUE))>0){MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("Incertae",MMC.ITS.taxa1$Genus,value=FALSE),]} else {MMC.ITS.taxa1 <- MMC.ITS.taxa1}
summary(colSums(MMC.ITS.counts1))
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
summary(colSums(MMC.ITS.counts1))
dim(MMC.ITS.counts1)
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
summary(colSums(MMC.ITS.counts1))
Fungal.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.ITS.counts1)),colData = DataFrame(MMC.ITS.samples),rowData = DataFrame(MMC.ITS.taxa11))
Fungal.tse <- transformAssay(Fungal.tse_taxa, MARGIN = "samples", method = "relabundance")
Fungal.tse <- addPerCellQC(Fungal.tse)
colData(Fungal.tse)$total_raw_counts <- colSums(assay(Fungal.tse, "counts"))
summary(colData(Fungal.tse)$total_raw_counts)

Fungal.tse_W0 <- Fungal.tse[,Fungal.tse$time=="W0"]
Fungal.tse_W0_Genus <- mergeFeaturesByRank(Fungal.tse_W0, rank ="Genus", onRankOnly=TRUE)
filtered_counts <- assay(Fungal.tse_W0_Genus, "relabundance")
filtered_counts[filtered_counts>0.01] <- 1
filtered_counts[filtered_counts<0.01] <- 0
filtered_counts["Kazachstania",]
sum(filtered_counts["Kazachstania",])
filtered_counts["Wickerhamomyces",]
sum(filtered_counts["Wickerhamomyces",])
keep_taxa <- rowSums(filtered_counts) > 3
Fungal.tse_W0_Genus_filtered <- Fungal.tse_W0_Genus[keep_taxa, ]
counts1 <- assay(Fungal.tse_W0_Genus_filtered, "counts")
counts1[c("Candida"),] <- counts1[c("Issatchenkia"),]+counts1[c("Candida"),]
counts2 <- counts1[setdiff(rownames(counts1),"Issatchenkia"),]
Fungal.tse_W0_Genus_filtered1 <- Fungal.tse_W0_Genus_filtered[-which(rownames(Fungal.tse_W0_Genus_filtered)=="Issatchenkia"),]
Fungal.tse_W0_Genus_filtered1@assays@data$counts <- counts2
Prevalence <- as.data.frame(getPrevalence(Fungal.tse_W0_Genus_filtered1, detection = 0,prevalence = 0, sort = TRUE,
    assay.type = "counts", as_relative = TRUE))
colnames(Prevalence) <- "Prevalence"
Prevalence$Genus <- rownames(Prevalence)
Prevalence$Rank <- nrow(Prevalence):1
Prevalence$Genus[1:14]
top_taxa.tmp_Genus <- Prevalence$Genus[c(1:13,15)]
top_taxa.tmp_Genus <- c("Candida","Saccharomyces","Cladosporium","Alternaria","Aspergillus",
  "Pichia","Papulaspora","Debaryomyces","Vishniacozyma","Wallemia","Clavispora","Fusarium","Dekkera","Penicillium")
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))[-c(8,15)]
pal <- pal[c(2,1,3:length(top_taxa.tmp_Genus))]
names(pal) <- top_taxa.tmp_Genus
pal <- c(pal,"#c9caca")
names(pal)[length(names(pal))] <- "others"


Fungal.tse_W0_Nystatin <- Fungal.tse[,Fungal.tse$time=="W0" & Fungal.tse$treatment=="Nystatin"]
Fungal.tse_W0_Nystatin_Genus <- mergeFeaturesByRank(Fungal.tse_W0_Nystatin, rank ="Genus", onRankOnly=TRUE)
counts1_Nystatin <- assay(Fungal.tse_W0_Nystatin_Genus, "relabundance")
counts1_Nystatin[c("Candida"),] <- counts1_Nystatin[c("Issatchenkia"),]+counts1_Nystatin[c("Candida"),]
counts2_Nystatin <- counts1_Nystatin[setdiff(rownames(counts1_Nystatin),"Issatchenkia"),]
Fungal.tse_W0_Nystatin_Genus1 <- Fungal.tse_W0_Nystatin_Genus[-which(rownames(Fungal.tse_W0_Nystatin_Genus)=="Issatchenkia"),]
Fungal.tse_W0_Nystatin_Genus1@assays@data$relabundance <- counts2_Nystatin
relabu_Nystatin <- assay(Fungal.tse_W0_Nystatin_Genus1, "relabundance")
relabu_Nystatin1 <- relabu_Nystatin[top_taxa.tmp_Genus,]
relabu_Nystatin1_tmp <- as.data.frame(t(colSums(relabu_Nystatin[setdiff(rownames(relabu_Nystatin),top_taxa.tmp_Genus),])))
relabu_Nystatin1 <- rbind(relabu_Nystatin1,relabu_Nystatin1_tmp)
rownames(relabu_Nystatin1) <- c(top_taxa.tmp_Genus,"others")
relabu_Nystatin1 <- sweep(as.matrix(relabu_Nystatin1), 2, colSums(as.matrix(relabu_Nystatin1)), "/")
relabu_Nystatin2 <- reshape2::melt(relabu_Nystatin1)
relabu_Nystatin2 <- as.data.frame(cbind(relabu_Nystatin2,colData(Fungal.tse_W0_Nystatin_Genus1)[as.character(relabu_Nystatin2$Var2),]))
o <- relabu_Nystatin2[relabu_Nystatin2$time=="W0",]
o1 <- as.data.frame(o[o$Var1=="Candida","value"])
o1$patient <- o[o$Var1=="Candida","patient"]
o1 <- o1[order(o1[,1],decreasing=TRUE),]
relabu_Nystatin2$patient <- factor(relabu_Nystatin2$patient,levels=c(o1$patient))
top_taxa.tmp_Genus1 <- unique(c(setdiff(c(top_taxa.tmp_Genus,"others"), "Candida"), "Candida"))
relabu_Nystatin2$Var1 <- factor(relabu_Nystatin2$Var1, levels = rev(top_taxa.tmp_Genus1))
relabu_Nystatin2$time <- factor(relabu_Nystatin2$time,levels=c("W0","W1","W2","W4","W8"))
p1 <- ggplot(relabu_Nystatin2, aes(x = patient, y = value, fill = Var1)) +geom_bar(stat = "identity") +
  labs(title = paste0(" Top 14 Fungal Genus (Relative Abundance) Nystatin"),x = "Sample",y = "Relative Abundance") +
  theme_minimal() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_color_manual(values = pal)+ 
  scale_fill_manual(values = pal)
Fungal.tse_W0_Fluconazole <- Fungal.tse[,Fungal.tse$time=="W0" & Fungal.tse$treatment=="Fluconazole"]
Fungal.tse_W0_Fluconazole_Genus <- mergeFeaturesByRank(Fungal.tse_W0_Fluconazole, rank ="Genus", onRankOnly=TRUE)
counts1_Fluconazole <- assay(Fungal.tse_W0_Fluconazole_Genus, "relabundance")
counts1_Fluconazole[c("Candida"),] <- counts1_Fluconazole[c("Issatchenkia"),]+counts1_Fluconazole[c("Candida"),]
counts2_Fluconazole <- counts1_Fluconazole[setdiff(rownames(counts1_Fluconazole),"Issatchenkia"),]
Fungal.tse_W0_Fluconazole_Genus1 <- Fungal.tse_W0_Fluconazole_Genus[-which(rownames(Fungal.tse_W0_Fluconazole_Genus)=="Issatchenkia"),]
Fungal.tse_W0_Fluconazole_Genus1@assays@data$relabundance <- counts2_Fluconazole
relabu_Fluconazole <- assay(Fungal.tse_W0_Fluconazole_Genus1, "relabundance")
relabu_Fluconazole1 <- relabu_Fluconazole[top_taxa.tmp_Genus,]
relabu_Fluconazole1_tmp <- as.data.frame(t(colSums(relabu_Fluconazole1[setdiff(rownames(relabu_Fluconazole1),top_taxa.tmp_Genus),])))
relabu_Fluconazole1 <- rbind(relabu_Fluconazole1,relabu_Fluconazole1_tmp)
rownames(relabu_Fluconazole1) <- c(top_taxa.tmp_Genus,"others")
relabu_Fluconazole1 <- sweep(as.matrix(relabu_Fluconazole1), 2, colSums(as.matrix(relabu_Fluconazole1)), "/")
relabu_Fluconazole2 <- reshape2::melt(relabu_Fluconazole1)
relabu_Fluconazole2 <- as.data.frame(cbind(relabu_Fluconazole2,colData(Fungal.tse_W0_Fluconazole_Genus1)[as.character(relabu_Fluconazole2$Var2),]))
o <- relabu_Fluconazole2[relabu_Fluconazole2$time=="W0",]
o1 <- as.data.frame(o[o$Var1=="Candida","value"])
o1$patient <- o[o$Var1=="Candida","patient"]
o1 <- o1[order(o1[,1],decreasing=TRUE),]
relabu_Fluconazole2$patient <- factor(relabu_Fluconazole2$patient,levels=c(o1$patient))
top_taxa.tmp_Genus1 <- unique(c(setdiff(c(top_taxa.tmp_Genus,"others"), "Candida"), "Candida"))
relabu_Fluconazole2$Var1 <- factor(relabu_Fluconazole2$Var1, levels = rev(top_taxa.tmp_Genus1))
relabu_Fluconazole2$time <- factor(relabu_Fluconazole2$time,levels=c("W0","W1","W2","W4","W8"))
p2 <- ggplot(relabu_Fluconazole2, aes(x = patient, y = value, fill = Var1)) +geom_bar(stat = "identity") +
  labs(title = paste0(" Top 14 Fungal Genus (Relative Abundance) Fluconazole"),x = "Sample",y = "Relative Abundance") +
  theme_minimal() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_color_manual(values = pal)+ 
  scale_fill_manual(values = pal)
plot <- plot_grid(p1,p2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig2.W0abudan.svg", plot=plot,width = 8, height = 5,dpi=300)
~~~

![image-20260615092423518](./Chapter2.assets/image-20260615092423518.png)

## 2.1.20 Global Fungal Prevalence

Computes the global prevalence of fungal taxa across all samples, exported as `Fig2.prevee_global.svg`.

~~~R
Fungal.tse_W0 <- Fungal.tse
Fungal.tse_W0_Genus <- mergeFeaturesByRank(Fungal.tse_W0, rank ="Genus", onRankOnly=TRUE)
filtered_counts <- assay(Fungal.tse_W0_Genus, "relabundance")
filtered_counts[filtered_counts>0.01] <- 1
filtered_counts[filtered_counts<0.01] <- 0
filtered_counts["Kazachstania",]
sum(filtered_counts["Kazachstania",])
filtered_counts["Wickerhamomyces",]
sum(filtered_counts["Wickerhamomyces",])
keep_taxa <- rowSums(filtered_counts) > 4
Fungal.tse_W0_Genus_filtered <- Fungal.tse_W0_Genus[keep_taxa, ]
counts1 <- assay(Fungal.tse_W0_Genus_filtered, "counts")
counts1[c("Candida"),] <- counts1[c("Issatchenkia"),]+counts1[c("Candida"),]
counts2 <- counts1[setdiff(rownames(counts1),"Issatchenkia"),]
Fungal.tse_W0_Genus_filtered1 <- Fungal.tse_W0_Genus_filtered[-which(rownames(Fungal.tse_W0_Genus_filtered)=="Issatchenkia"),]
Fungal.tse_W0_Genus_filtered1@assays@data$counts <- counts2
Prevalence <- as.data.frame(getPrevalence(Fungal.tse_W0_Genus_filtered1, detection = 0,prevalence = 0, sort = TRUE,
    assay.type = "counts", as_relative = TRUE))
colnames(Prevalence) <- "Prevalence"
Prevalence$Genus <- rownames(Prevalence)
Prevalence$Rank <- nrow(Prevalence):1
Prevalence$Genus[1:14]
top_taxa.tmp_Genus <- Prevalence$Genus[c(1:14)]
top_taxa.tmp_Genus <- c("Candida","Saccharomyces","Cladosporium","Alternaria","Aspergillus",
  "Pichia","Papulaspora","Debaryomyces","Vishniacozyma","Wallemia","Clavispora","Fusarium","Dekkera","Penicillium")
pal <- c(jdb_palette("corona"),jdb_palette(c("lawhoops")),jdb_palette(c("brewer_spectra")))[-c(8,15)]
pal <- pal[c(2,1,3:length(top_taxa.tmp_Genus))]
names(pal) <- top_taxa.tmp_Genus
pal <- c(pal,"#c9caca")
names(pal)[length(names(pal))] <- "others"

plot <- ggplot(Prevalence, aes(x = Prevalence, y = Rank)) +
geom_point(size = 1, color = "#1f78b4") + ggrepel::geom_label_repel(data = Prevalence[c(1:14),], 
aes(label = Genus), segment.color = 'black', show.legend = FALSE,size=3,max.overlaps=50,min.segment.length = unit(0.1, 'lines'))+
theme_classic()+labs(title = "Top 14 Fungal Genus by Prevalence (Ranked)",x = "Prevalence",y = "Genus Rank")
ggsave("./projects/MMC/Figures/v2_figures/Fig2.prevee_global.svg", plot=plot,width = 5, height = 5,dpi=300)
~~~

![image-20260615092456169](./Chapter2.assets/image-20260615092456169.png)

## 2.1.21 Fungal Abundance Across All Time Points

Summarizes fungal taxonomic abundance across all time points, exported as `Fig2.Allabudan.svg`.

~~~R
Fungal.tse_W0_Nystatin <- Fungal.tse[, Fungal.tse$treatment=="Nystatin"]
Fungal.tse_W0_Nystatin_Genus <- mergeFeaturesByRank(Fungal.tse_W0_Nystatin, rank ="Genus", onRankOnly=TRUE)
counts1_Nystatin <- assay(Fungal.tse_W0_Nystatin_Genus, "relabundance")
counts1_Nystatin[c("Candida"),] <- counts1_Nystatin[c("Issatchenkia"),]+counts1_Nystatin[c("Candida"),]
counts2_Nystatin <- counts1_Nystatin[setdiff(rownames(counts1_Nystatin),"Issatchenkia"),]
Fungal.tse_W0_Nystatin_Genus1 <- Fungal.tse_W0_Nystatin_Genus[-which(rownames(Fungal.tse_W0_Nystatin_Genus)=="Issatchenkia"),]
Fungal.tse_W0_Nystatin_Genus1@assays@data$relabundance <- counts2_Nystatin
relabu_Nystatin <- assay(Fungal.tse_W0_Nystatin_Genus1, "relabundance")
relabu_Nystatin1 <- relabu_Nystatin[top_taxa.tmp_Genus,]
relabu_Nystatin1_tmp <- as.data.frame(t(colSums(relabu_Nystatin[setdiff(rownames(relabu_Nystatin),top_taxa.tmp_Genus),])))
relabu_Nystatin1 <- rbind(relabu_Nystatin1,relabu_Nystatin1_tmp)
rownames(relabu_Nystatin1) <- c(top_taxa.tmp_Genus,"others")
relabu_Nystatin1 <- sweep(as.matrix(relabu_Nystatin1), 2, colSums(as.matrix(relabu_Nystatin1)), "/")
relabu_Nystatin2 <- reshape2::melt(relabu_Nystatin1)
relabu_Nystatin2 <- as.data.frame(cbind(relabu_Nystatin2,colData(Fungal.tse_W0_Nystatin_Genus1)[as.character(relabu_Nystatin2$Var2),]))
o <- relabu_Nystatin2
o1 <- as.data.frame(o[o$Var1=="Candida","value"])
o1$ITS_names <- o[o$Var1=="Candida","ITS_names"]
o1 <- o1[order(o1[,1],decreasing=TRUE),]
relabu_Nystatin2$ITS_names <- factor(relabu_Nystatin2$ITS_names,levels=c(o1$ITS_names))
top_taxa.tmp_Genus1 <- unique(c(setdiff(c(top_taxa.tmp_Genus,"others"), "Candida"), "Candida"))
relabu_Nystatin2$Var1 <- factor(relabu_Nystatin2$Var1, levels = rev(top_taxa.tmp_Genus1))
relabu_Nystatin2$time <- factor(relabu_Nystatin2$time,levels=c("W0","W1","W2","W4","W8"))
p1 <- ggplot(relabu_Nystatin2, aes(x = ITS_names, y = value, fill = Var1)) +geom_bar(stat = "identity") +
  labs(title = paste0(" Top 14 Fungal Genus (Relative Abundance) Nystatin"),x = "Sample",y = "Relative Abundance") +
  theme_minimal() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_color_manual(values = pal)+ 
  scale_fill_manual(values = pal)
Fungal.tse_W0_Fluconazole <- Fungal.tse[, Fungal.tse$treatment=="Fluconazole"]
Fungal.tse_W0_Fluconazole_Genus <- mergeFeaturesByRank(Fungal.tse_W0_Fluconazole, rank ="Genus", onRankOnly=TRUE)
counts1_Fluconazole <- assay(Fungal.tse_W0_Fluconazole_Genus, "relabundance")
counts1_Fluconazole[c("Candida"),] <- counts1_Fluconazole[c("Issatchenkia"),]+counts1_Fluconazole[c("Candida"),]
counts2_Fluconazole <- counts1_Fluconazole[setdiff(rownames(counts1_Fluconazole),"Issatchenkia"),]
Fungal.tse_W0_Fluconazole_Genus1 <- Fungal.tse_W0_Fluconazole_Genus[-which(rownames(Fungal.tse_W0_Fluconazole_Genus)=="Issatchenkia"),]
Fungal.tse_W0_Fluconazole_Genus1@assays@data$relabundance <- counts2_Fluconazole
relabu_Fluconazole <- assay(Fungal.tse_W0_Fluconazole_Genus1, "relabundance")
relabu_Fluconazole1 <- relabu_Fluconazole[top_taxa.tmp_Genus,-c(which(colnames(relabu_Fluconazole)=="MMC003W0F"),which(colnames(relabu_Fluconazole)=="MMC003W1F"))]
relabu_Fluconazole1_tmp <- as.data.frame(t(colSums(relabu_Fluconazole1[setdiff(rownames(relabu_Fluconazole1),top_taxa.tmp_Genus),])))
relabu_Fluconazole1 <- rbind(relabu_Fluconazole1,relabu_Fluconazole1_tmp)
rownames(relabu_Fluconazole1) <- c(top_taxa.tmp_Genus,"others")
relabu_Fluconazole1 <- sweep(as.matrix(relabu_Fluconazole1), 2, colSums(as.matrix(relabu_Fluconazole1)), "/")
relabu_Fluconazole2 <- reshape2::melt(relabu_Fluconazole1)
relabu_Fluconazole2 <- as.data.frame(cbind(relabu_Fluconazole2,colData(Fungal.tse_W0_Fluconazole_Genus1)[as.character(relabu_Fluconazole2$Var2),]))
o <- relabu_Fluconazole2
o1 <- as.data.frame(o[o$Var1=="Candida","value"])
o1$ITS_names <- o[o$Var1=="Candida","ITS_names"]
o1 <- o1[order(o1[,1],decreasing=TRUE),]
relabu_Fluconazole2$ITS_names <- factor(relabu_Fluconazole2$ITS_names,levels=c(o1$ITS_names))
top_taxa.tmp_Genus1 <- unique(c(setdiff(c(top_taxa.tmp_Genus,"others"), "Candida"), "Candida"))
relabu_Fluconazole2$Var1 <- factor(relabu_Fluconazole2$Var1, levels = rev(top_taxa.tmp_Genus1))
relabu_Fluconazole2$time <- factor(relabu_Fluconazole2$time,levels=c("W0","W1","W2","W4","W8"))
p2 <- ggplot(relabu_Fluconazole2, aes(x = ITS_names, y = value, fill = Var1)) +geom_bar(stat = "identity") +
  labs(title = paste0(" Top 14 Fungal Genus (Relative Abundance) Fluconazole"),x = "Sample",y = "Relative Abundance") +
  theme_minimal() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +scale_color_manual(values = pal)+ 
  scale_fill_manual(values = pal)
plot <- plot_grid(p1,p2)
ggsave("./projects/MMC/Figures/v2_figures/Fig2.Allabudan.svg", plot=plot,width = 8, height = 5,dpi=300)
~~~

![image-20260615092547535](./Chapter2.assets/image-20260615092547535.png)

## 2.1.22 *Candida* / *Saccharomyces* Relative Abundance Table

Rebuilds the fungal TSE on the extended cohort (including additional CD patients) and exports relative abundances of *Candida albicans*, *Saccharomyces cerevisiae*, and their genera as `relativeAb_forMD.csv`.

~~~R
ReadCap.20251215 <- mcreadRDS("./workshop/MMC/Aidan_info/v2/ReadCap.20251215.rds")
total_times_tmp <- read.csv("./projects/MMC/Figures/figures_making/v3/patient.PFS.20251215.csv")
Nys_CD <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% 
    setdiff(ReadCap.20251215$Omics_patient_Names,total_times_tmp$Omics_patient_Names) &
    ReadCap.20251215$Diagnosis_new=="CD" & ReadCap.20251215$treatment=="Nystatin" ,]
Flu_CD <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% 
    setdiff(ReadCap.20251215$Omics_patient_Names,total_times_tmp$Omics_patient_Names) &
    ReadCap.20251215$Diagnosis_new=="CD" & ReadCap.20251215$treatment=="Fluconazole" &
    ReadCap.20251215$CD.score.raw<=5 ,]
Flu_CD <- Flu_CD[!is.na(Flu_CD$treatment),]
c("MMC114",Flu_CD$Omics_patient_Names,Nys_CD$Omics_patient_Names)
sample_info_total2 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% c("MMC114",Flu_CD$Omics_patient_Names,Nys_CD$Omics_patient_Names,total_times_tmp$Omics_patient_Names),]
sample_info_total2 <- sample_info_total2[sample_info_total2$time!="W6",]
sample_info_total2 <- sample_info_total2[sample_info_total2$treatment!="Clotrimazole",]
sample_info_total5 <- ReadCap.20251215[ReadCap.20251215$Omics_patient_Names %in% sample_info_total2$Omics_patient_Names,]
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
# MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("NA",MMC.ITS.taxa1$Species,value=FALSE),]
if (length(grep("Incertae",MMC.ITS.taxa1$Genus,value=TRUE))>0){MMC.ITS.taxa1 <- MMC.ITS.taxa1[-grep("Incertae",MMC.ITS.taxa1$Genus,value=FALSE),]} else {MMC.ITS.taxa1 <- MMC.ITS.taxa1}
MMC.ITS.counts1 <- MMC.ITS.counts1[intersect(rownames(MMC.ITS.taxa1),rownames(MMC.ITS.counts1)),]
MMC.ITS.taxa11 <- MMC.ITS.taxa1[rownames(MMC.ITS.counts1),]
MMC.ITS.samples <- MMC.ITS.samples[colnames(MMC.ITS.counts1),]
MMC.ITS <- as.data.frame(cbind(MMC.ITS.taxa11[rownames(MMC.ITS.counts1),],MMC.ITS.counts1))
colSums(MMC.ITS.counts1[,sort(grep("W0",colnames(MMC.ITS.counts1),value=FALSE))])

sort(unique(MMC.ITS.taxa1$Species))
Fungal.tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = as.matrix(MMC.ITS.counts1)),colData = DataFrame(MMC.ITS.samples),rowData = DataFrame(MMC.ITS.taxa11))
Fungal.tse <- transformAssay(Fungal.tse_taxa, MARGIN = "samples", method = "relabundance")
Fungal.tse <- addPerCellQC(Fungal.tse)
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
relativeAb_forMD <- relativeAb_all[relativeAb_all$names%in% c("Candida albicans","Saccharomyces cerevisiae","Candida","Saccharomyces"),]
write.csv(relativeAb_forMD[,sort(colnames(relativeAb_forMD))],"./projects/MMC/Figures/v2_figures/relativeAb_forMD.csv")
~~~

## 2.1.23 Baseline *Candida albicans* Burden (Density)

Plots the W0 baseline distribution (density) of *Candida albicans* relative abundance, exported as `Candida_albicans.W0_onlyabudan_density.svg`.

~~~R
relativeAb_all_Ca <- relativeAb_all[relativeAb_all$names=="Candida albicans",]
relativeAb_all_Ca <- as.data.frame(t(relativeAb_all_Ca[,colnames(MMC.ITS.counts1)]))
relativeAb_all_Ca1 <- as.data.frame(cbind(relativeAb_all_Ca,MMC.ITS.samples[rownames(relativeAb_all_Ca),]))
colnames(relativeAb_all_Ca1) <- gsub(" ","_",colnames(relativeAb_all_Ca1))
relativeAb_all_Ca1 <- relativeAb_all_Ca1[relativeAb_all_Ca1$treatment %in% c("Nystatin","Fluconazole"),]
relativeAb_all_Ca1$treatment <- factor(relativeAb_all_Ca1$treatment,levels=c("Nystatin","Fluconazole"))
relativeAb_all_Ca2 <- relativeAb_all_Ca1[relativeAb_all_Ca1$time=="W0",]
comb <- list(c("Nystatin","Fluconazole"))
plot <- ggplot(relativeAb_all_Ca2, aes_string(x="Candida_albicans", y="treatment", fill="treatment")) + ggridges::geom_density_ridges() +
    theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+NoLegend()+
    labs(title = paste0("Candida_albicans"),y = paste("Candida_albicans"))+
      geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
plot
ggsave("./projects/MMC/Figures/v2_figures/Candida_albicans.W0_onlyabudan_density.svg", plot=plot,width = 5, height = 5,dpi=300)
~~~

![image-20260615092711566](./Chapter2.assets/image-20260615092711566.png)

## 2.1.24 Baseline *Candida albicans* Burden (Abundance)

Plots the W0 baseline *Candida albicans* abundance per group, exported as `Candida_albicans.W0_onlyabudan.svg`.

~~~R
relativeAb_all_Ca2 <- relativeAb_all_Ca1[relativeAb_all_Ca1$time=="W0",]
p1 <- ggbarplot(relativeAb_all_Ca2, x = "treatment", y = "Candida_albicans",add = c("mean_se","jitter"),legend = "none",rotate=FALSE,
  color = "treatment", fill = "treatment", alpha = 1,title=paste0("Candida_albicans"))+
    stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) + 
    stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+
    stat_summary(fun = mean, geom = "line",aes(group = 1),col = "red",size=1)+NoLegend()+
    theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
    geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
relativeAb_all_Ca <- relativeAb_all[relativeAb_all$names=="Candida",]
relativeAb_all_Ca <- as.data.frame(t(relativeAb_all_Ca[,colnames(MMC.ITS.counts1)]))
relativeAb_all_Ca1 <- as.data.frame(cbind(relativeAb_all_Ca,MMC.ITS.samples[rownames(relativeAb_all_Ca),]))
colnames(relativeAb_all_Ca1) <- gsub(" ","_",colnames(relativeAb_all_Ca1))
relativeAb_all_Ca1 <- relativeAb_all_Ca1[relativeAb_all_Ca1$treatment %in% c("Nystatin","Fluconazole"),]
relativeAb_all_Ca1$treatment <- factor(relativeAb_all_Ca1$treatment,levels=c("Nystatin","Fluconazole"))
relativeAb_all_Ca2 <- relativeAb_all_Ca1[relativeAb_all_Ca1$time=="W0",]
p2 <- ggbarplot(relativeAb_all_Ca2, x = "treatment", y = "Candida",add = c("mean_se","jitter"),legend = "none",rotate=FALSE,
  color = "treatment", fill = "treatment", alpha = 1,title=paste0("Candida"))+
    stat_summary(fun.y = mean, geom="point",colour="darkred", size=3) + 
    stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))+
    stat_summary(fun = mean, geom = "line",aes(group = 1),col = "red",size=1)+NoLegend()+
    theme(axis.text.x  = element_text(angle=45, vjust=1,hjust = 1))+
    geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = ks.test_wrapper)
plot_grid(p1,p2)
ggsave("./projects/MMC/Figures/v2_figures/Candida_albicans.W0_onlyabudan.svg", plot=p1,width = 5, height = 5,dpi=300)
~~~

![image-20260615092734826](./Chapter2.assets/image-20260615092734826.png)

## 2.1.25 Baseline (W0) PCoA by Treatment

Loads the W0 baseline PERMANOVA results and draws the baseline PCoA ordination colored by treatment group, exported as `Fig2.W0.Pcoa.svg`.

~~~R
All_sig <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC_ITS_PERMANOVA_res.W0.treatmente.rds")
library(ggside)
library(ape)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggside)
library(patchwork)
generate_ellipse <- function(median_x, mad_x, median_y, mad_y, n = 500, factor) {
  t <- seq(0, 2*pi, length.out = n)
  a <- mad_x * factor
  b <- mad_y * factor
  x <- median_x + a * cos(t)
  y <- median_y + b * sin(t)
  return(data.frame(x = x, y = y))
}

Fungal.tse_W0 <- Fungal.tse[,Fungal.tse$time=="W0"]
Fungal.tse_W0 <- Fungal.tse_W0[,Fungal.tse_W0$treatment %in% c("Nystatin","Fluconazole")]
Fungal.tse_W0$treatment <- factor(Fungal.tse_W0$treatment, levels = c("Nystatin","Fluconazole"))
tse_tmp <- Fungal.tse_W0
Sel_type <- c("Species")
mat_rel_ <- lapply(1:length(Sel_type),function(t) {
  tse_tmp1 <- subsetByPrevalentFeatures(tse_tmp, rank = Sel_type[t], detection = 0.1, prevalence = 0.1)
  ## 相对丰度 + Jaccard（presence/absence）
  tse_tmp1 <- transformAssay(tse_tmp1, assay.type = "counts", method = "relabundance")
  mat_rel1 <- t(assay(tse_tmp1, "relabundance"))
  return(mat_rel1)
  })
mat_rel <- as.data.frame(do.call(cbind,mat_rel_))

dist_jac <- vegdist(mat_rel > 0, method = "jaccard")
## PCoA
pcoa_res <- pcoa(dist_jac)
var_explained <- round(100 * pcoa_res$values$Relative_eig[1:2], 1)
PcoA <- data.frame(x = pcoa_res$vectors[,1],y = pcoa_res$vectors[,2],treatment = tse_tmp$treatment)
## 分组求 median 和 MAD
sbg <- PcoA %>% group_by(treatment) %>%
  summarise(median_x = median(x), mad_x = mad(x),median_y = median(y), mad_y = mad(y), .groups='drop')
ellipse_data <- do.call(rbind, lapply(1:nrow(sbg), function(i) {
  df <- generate_ellipse(sbg$median_x[i], sbg$mad_x[i],sbg$median_y[i], sbg$mad_y[i], factor=0.8)
  df$treatment <- sbg$treatment[i]
  df
}))
All_sig_tmp <- All_sig[All_sig$prevalence==0.1,]
p <- ggplot(PcoA, aes(x=x, y=y, color=treatment)) + geom_point(aes(shape=treatment, fill=treatment), size=5, alpha=0.5) +
  geom_polygon(data=ellipse_data, aes(x=x, y=y, fill=treatment), alpha=0.2) +
  geom_errorbar(data = sbg, aes(x = median_x, y = median_y, ymin = median_y - 0.5 * mad_y, ymax = median_y + 0.5 * mad_y), width = 0) +
  geom_errorbarh(data = sbg, aes(x = median_x, y = median_y, xmin = median_x - 0.5 * mad_x, xmax = median_x + 0.5 * mad_x), height = 0) +
  geom_point(data=sbg, aes(x=median_x, y=median_y, fill=treatment),color="black", shape=22, size=5, alpha=0.7, show.legend=FALSE) +
  theme_bw() + scale_fill_manual(values=pal[2:1]) + scale_color_manual(values=pal[2:1]) +
  labs(x = paste0("PCoA1 (", var_explained[1], "%)"),y = paste0("PCoA2 (", var_explained[2], "%)"),
      title = All_sig_tmp$title_lab_bray) +
  geom_xsideboxplot(orientation="y") + geom_ysideboxplot(orientation="x")
p
ggsave("./projects/MMC/Figures/v2_figures/Fig2.W0.Pcoa.svg", plot=p,width = 5, height = 4,dpi=300)
~~~

![image-20260615092804928](./Chapter2.assets/image-20260615092804928.png)





