# 6.3 Correlation: Clinical Indices vs Microbiome/Metabolome Features

This script relates microbiome and metabolome features to clinical disease indices (DAI, UC and CD scores). It computes SCFA-producer abundances, derives per-indicator weekly change tables, correlates features with clinical scores (raw and W0-delta, via permutation MIC), and produces longitudinal bar/LOESS plots comparing Nystatin (ORNT) vs Fluconazole (GIFT) in UC and CD subsets.

## 6.3.1 SCFA-producer / Probiotic Abundances

Defines butyrate-, propionate-, and valerate-producing (and probiotic) taxa sets, sums their relative abundance per sample, and appends them to the merged multi-omics table.

~~~R
relativeAb_all <-mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/MMC.metagenomic.relativeAb_all.rds")
MMC.METAG.inc_all <- mcreadRDS("./projects/ITS_Others/Lib40/MMC_ITS/MMC.METAG.inc_all.rds")
probiotic <- c(grep("Clostridium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),"Lactobacillus","Lactobacillus_acidophilus","Bifidobacteriales","Streptococcus_thermophilus","Bacillus",
    grep("Streptococcus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Lactobacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bifidobacteriales",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),
    grep("Faecalibacterium",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE),grep("Bacillus",unique(MMC.METAG.inc_all[[2]]$Names),value=TRUE))
Butyrate_producing_bac <- c("Anaerostipes_hadrus","Clostridium_butyricum","Coprococcus_comes","Eubacterium_hallii","Eubacterium_rectale","Faecalibacterium_prausnitzii","Roseburia_hominis","Roseburia_intestinalis","Butyrivibrio_fibrisolvens","Ruminococcus_bromii","Megasphaera_elsdenii","Odoribacter_splanchnicus")
Propionate_producing_bac <- c("Bacteroides_thetaiotaomicron","Bacteroides_fragilis","Bacteroides_ovatus","Bacteroides_uniformis","Phascolarctobacterium_succinatutens","Veillonella_parvula","Dialister_succinatiphilus","Prevotella_copri","Propionibacterium_freudenreichii","Anaerovibrio_lipolytica","Oddibacterium_anthropi")
Valeric_producing_bac <- c("Butyrivibrio_fibrisolvens","Megasphaera_elsdenii","Coprococcus_eutactus","Blautia")
intersect(unique(MMC.METAG.inc_all[[2]]$Names),Butyrate_producing_bac)
intersect(unique(MMC.METAG.inc_all[[2]]$Names),Propionate_producing_bac)
intersect(unique(MMC.METAG.inc_all[[2]]$Names),Valeric_producing_bac)

microbiota3 <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/All_omics_values.rds")
relabundance_tmp <- relativeAb_all[intersect(Butyrate_producing_bac,rownames(relativeAb_all)),intersect(colnames(relativeAb_all),rownames(microbiota3))]
microbiota3$Butyrate_producing_bac <- as.numeric(100*colMeans(relabundance_tmp)[rownames(microbiota3)])
relabundance_tmp <- relativeAb_all[intersect(Propionate_producing_bac,rownames(relativeAb_all)),intersect(colnames(relativeAb_all),rownames(microbiota3))]
microbiota3$Propionate_producing_bac <- as.numeric(100*colMeans(relabundance_tmp)[rownames(microbiota3)])
relabundance_tmp <- relativeAb_all[intersect(Valeric_producing_bac,rownames(relativeAb_all)),intersect(colnames(relativeAb_all),rownames(microbiota3))]
microbiota3$Valeric_producing_bac <- as.numeric(100*colMeans(relabundance_tmp)[rownames(microbiota3)])

mcsaveRDS(microbiota3,"./workshop/MMC/sample_info/final_Res/v2/All_omics_values.rds")
~~~

## 6.3.2 Per-indicator Weekly Change Tables

Adds fecal SCFA metabolite levels (butyric/propionic/valeric acid) and, for each indicator, builds a per-patient table of log2 values at W0–W8 plus the W2-vs-W0 change (with W4/W8 fallback) and percent-of-baseline change, saved as `All_microbiome_diff.rds`.

~~~R
microbiota3 <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/All_omics_values.rds")
corrected <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/metabolomics.counts.corrected.rds")
Healthy <- c("BUTYRIC ACID","PROPIONIC ACID","VALERIC ACID")
corrected_sel <- as.data.frame(t(corrected[Healthy,]))
colnames(corrected_sel) <- gsub(" ","_",colnames(corrected_sel))
setdiff(rownames(microbiota3),rownames(corrected_sel))
setdiff(rownames(corrected_sel),rownames(microbiota3))
microbiota3 <- as.data.frame(cbind(microbiota3,corrected_sel[rownames(microbiota3),]))
microbiota3[grep("MMC096",microbiota3$patient,value=FALSE),]
indicators <- c("pathogenic.Fungi.v2","Fungal.observed","Fungal.shannon","Bac.shannon","Bac.observed","gut_associated_bac",
        "Candida","Candida_albicans","Dysbiosis.UP","Butyrate_producing_bac","BUTYRIC_ACID","PROPIONIC_ACID","VALERIC_ACID")
outdir <- "./projects/MMC/Figures/figures_making/v3/"
process_indicator <- function(ReadCap, indicator, outdir){
  df <- microbiota3[!is.na(microbiota3[[indicator]]), ]
  df <- df[df$time %in% c("W0","W1","W2","W4","W8"), ]
  rownames(df) <- df$MGX_names
  patient_info <- df[!duplicated(df$patient), ]
  uniq_patient <- unique(df$patient)
  
  All_weeks_ <- lapply(uniq_patient, function(pid) {
    tmp <- df[df$patient == pid, ]
    W0 <- tmp[tmp$time=="W0", indicator]
    W1 <- tmp[tmp$time=="W1", indicator]
    W2 <- tmp[tmp$time=="W2", indicator]
    W4 <- tmp[tmp$time=="W4", indicator]
    W8 <- tmp[tmp$time=="W8", indicator]
    data.frame(W0=ifelse(length(W0)==0,NA,W0),
               W1=ifelse(length(W1)==0,NA,W1),
               W2=ifelse(length(W2)==0,NA,W2),
               W4=ifelse(length(W4)==0,NA,W4),
               W8=ifelse(length(W8)==0,NA,W8))
  })
  All_weeks <- do.call(rbind, All_weeks_)
  All_weeks <- log(All_weeks+1,2)
  All_weeks$patient <- uniq_patient
  total_times_tmp <- All_weeks
  
  # 各个差值
  total_times_tmp$W1_W0 <- total_times_tmp$W1 - total_times_tmp$W0
  total_times_tmp$W2_W0 <- total_times_tmp$W2 - total_times_tmp$W0
  total_times_tmp$W4_W0 <- total_times_tmp$W4 - total_times_tmp$W0
  total_times_tmp$W8_W0 <- total_times_tmp$W8 - total_times_tmp$W0
  
  # 主 diff 列
  diff_col <- paste0(indicator, "_diff")
  total_times_tmp[[diff_col]] <- total_times_tmp$W2 - total_times_tmp$W0
  total_times_tmp[is.na(total_times_tmp$W2_W0) & is.na(total_times_tmp[[diff_col]]), diff_col] <- 
    total_times_tmp[is.na(total_times_tmp$W2_W0) & is.na(total_times_tmp[[diff_col]]), "W4_W0"]
  total_times_tmp[is.na(total_times_tmp$W4_W0) & is.na(total_times_tmp[[diff_col]]), diff_col] <- 
    total_times_tmp[is.na(total_times_tmp$W4_W0) & is.na(total_times_tmp[[diff_col]]), "W8_W0"]
  
  # ⚡ 新增 per baseline 百分比变化
  per_col <- paste0(indicator, "_per_baseline")
  total_times_tmp[[per_col]] <- 100 * abs(total_times_tmp[[diff_col]]) / total_times_tmp$W0
  
  # 删除 NA
  total_times_tmp <- total_times_tmp[!is.na(total_times_tmp[[diff_col]]), ]
  
  # 加诊断和治疗信息
  total_times_tmp$Diagnosis_new <- patient_info[as.character(total_times_tmp$patient), "Diagnosis_new"]
  total_times_tmp$treatment <- patient_info[as.character(total_times_tmp$patient), "treatment"]
  total_times_tmp$treatment <- factor(total_times_tmp$treatment, levels=c("Nystatin","Clotrimazole","Fluconazole"))
  total_times_tmp$treatment.v2 <- as.character(total_times_tmp$treatment)
  total_times_tmp$treatment.v2[total_times_tmp$treatment.v2 %in% c("Clotrimazole","Fluconazole")] <- "Clo_Flu"
  total_times_tmp$treatment.v2 <- factor(total_times_tmp$treatment.v2, levels=c("Nystatin","Clo_Flu"))
  
  if("group" %in% colnames(total_times_tmp)){
    total_times_tmp$group.v2 <- as.character(total_times_tmp$group)
    total_times_tmp$group.v2[total_times_tmp$group.v2 %in% c("PD","SD")] <- "NR"
  }
  
  total_times_tmp <- total_times_tmp[order(total_times_tmp$patient), ]
  
  outfile <- file.path(outdir, paste0("patient.", indicator, ".csv"))
  write.csv(total_times_tmp, outfile, row.names=FALSE)
  message(indicator)
  return(total_times_tmp)
}

# 循环跑
results_list <- lapply(indicators, function(ind){
  process_indicator(microbiota3, ind, outdir)
})
names(results_list) <- indicators
mcsaveRDS(results_list,"./workshop/MMC/Aidan_info/v2/All_microbiome_diff.rds")
~~~

## 6.3.3 Features vs DAI (Raw)

Defines the signed permutation-MIC function (`sMIC_new`) and the scatter helper (`plot_module_cor_xASm`), then correlates each microbiome/metabolome feature against the raw DAI score, exported as `Fig5.1.clincal_raw.cor.all`.

~~~R
sMIC_new <- function(x, y, nperm = 2e4, seed = NULL, method = "spearman",
                 return_null = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  rho  <- suppressWarnings(cor(x, y, method = method))
  mic0 <- mine(x, y)$MIC
  if(length(nperm)>1){nperm=sample(nperm,1)}
  null_mic <- vapply(seq_len(nperm), function(i) mine(x, sample(y))$MIC, 0.0)
  k <- sum(null_mic >= mic0)
  p_mic <- (k + 1) / (nperm + 1)   # 永不为 0 的置换 p

  sign_rho <- if (is.na(rho) || rho == 0) 0 else sign(rho)
  out <- c(sMIC = sign_rho * mic0, MIC = mic0, rho = unname(rho),
           p_MIC = p_mic, k = k, B = nperm)
  if (return_null) attr(out, "null_mic") <- null_mic
  out
}
microbiota3 <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/All_omics_values.rds")
microbiota3$time <- gsub("F","",microbiota3$time)
microbiota3$time.v2 <- as.character(microbiota3$time)
microbiota3$time.v2[microbiota3$time.v2 %in% "W0"] <- "Pre"
microbiota3$time.v2[microbiota3$time.v2 %in% c("W1","W2","W4")] <- "Post"
microbiota3$time.v2[microbiota3$time.v2 %in% "W8"] <- "LTM"
microbiota3$time.v2 <- factor(microbiota3$time.v2, levels=c("Pre","Post","LTM"))
plot_module_cor_xASm <- function(df, module, yvar = "Candida_albicans",MIC=NULL, size = "abs.Ca",color_by="Diagnosis_new") {
  pal <- jdb_palette("corona")
  f <- as.formula(paste(yvar, "~", module))
  total.lm <- lm(f, data = df)
  if (is.null(MIC)) {
    rho <- cor(total.lm$model[[module]], total.lm$model[[yvar]], method = "pearson")
    pval <- cor.test(total.lm$model[[module]], total.lm$model[[yvar]], method = "pearson")$p.value
  } else {
    rho <- MIC[[module]][["sMIC"]]
    pval <- MIC[[module]][["p_MIC"]]
  }
  if(is.list(size)) {
    size1 <- size[[module]]
  } else {size1=size}
  ypos <- max(df[[yvar]], na.rm = TRUE)
  xpos <- median(df[[module]], na.rm = TRUE)
  p <- ggplot(df, aes_string(module, yvar, color = color_by, size = size1)) + geom_point(alpha = 0.9) +
    geom_rug(size = 0.4) + geom_smooth(aes_string(x = module, y = yvar), method = "lm", se = TRUE,color = "grey20", fill = "grey70", inherit.aes = FALSE) +
    scale_color_manual(values = pal) +scale_size_continuous(range = c(1, 8), guide = "none") + theme_minimal(base_size = 12) +
    labs(x = module, y = yvar,title = paste0("Cor: ", module, " vs ", yvar)) +
    annotate("text", x = xpos, y = ypos, color = "#E51718",
      label = paste0("rho = ", round(rho, 2), "; p = ", signif(pval, 3)))
  p
  return(p)
}
diease_Scores <- c("pathogenic.Fungi.v2","Fungal.observed","Fungal.shannon","Bac.shannon","Bac.observed","gut_associated_bac",
        "Candida","Candida_albicans","Bile_Acid_Biosynthesis","Fatty_Acid_Biosynthesis","Health.UP","Dysbiosis.UP","Module1","Module2","Module3","Module4","Butyrate_producing_bac","Valeric_producing_bac")
DAI_MIC <- future_lapply(1:length(diease_Scores),function(x) {
  tmp <- microbiota3[!is.na(microbiota3[[diease_Scores[x]]]),]
  res <- sMIC_new(tmp[[diease_Scores[x]]], tmp$DAI,method = "spearman")
  res
  })
names(DAI_MIC) <- diease_Scores

size_name <- as.list(diease_Scores)
names(size_name) <- diease_Scores

plots <- lapply(diease_Scores, function(m) 
  plot_module_cor_xASm(microbiota3, m, yvar = "DAI",MIC=DAI_MIC,size=size_name,color_by="time.v2"))
plot <- CombinePlots(c(plots),nrow=4)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_raw.cor.all.svg", plot=plot,width = 25, height = 16,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_raw.cor.all.png", plot=plot,width = 25, height = 16,dpi=300)
~~~

![Fig5.1.clincal_raw.cor.all](./Cor_Clincal_micro.assets/Fig5.1.clincal_raw.cor.all.png)

## 6.3.4 Features vs UC Score (Raw)

Restricts to UC patients and correlates each feature with the raw UC score, exported as `Fig5.1.UC.score.v2_raw.cor.all`.

~~~R
microbiota3_UC <- microbiota3[microbiota3$Diagnosis_new=="UC",]
UC.score.v2_MIC <- future_lapply(1:length(diease_Scores),function(x) {
  tmp <- microbiota3_UC[!is.na(microbiota3_UC[[diease_Scores[x]]]),]
  res <- sMIC_new(tmp[[diease_Scores[x]]], tmp$UC.score.v2,method = "spearman")
  res
  })
names(UC.score.v2_MIC) <- diease_Scores

plots <- lapply(diease_Scores, function(m) 
  plot_module_cor_xASm(microbiota3_UC, m, yvar = "UC.score.v2",MIC=UC.score.v2_MIC,size=size_name))
plot <- CombinePlots(c(plots),nrow=4)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.UC.score.v2_raw.cor.all.svg", plot=plot,width = 30, height = 16,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.UC.score.v2_raw.cor.all.png", plot=plot,width = 30, height = 16,dpi=300)
~~~

![Fig5.1.UC.score.v2_raw.cor.all](./Cor_Clincal_micro.assets/Fig5.1.UC.score.v2_raw.cor.all.png)

## 6.3.5 Features vs CD Score (Raw)

Restricts to CD patients and correlates each feature with the raw CD score, exported as `Fig5.1.CD.score.raw_raw.cor.all`.

~~~R
microbiota3_CD <- microbiota3[microbiota3$Diagnosis_new=="CD",]
CD.score.raw_MIC <- future_lapply(1:length(diease_Scores),function(x) {
  tmp <- microbiota3_CD[!is.na(microbiota3_CD[[diease_Scores[x]]]),]
  res <- sMIC_new(tmp[[diease_Scores[x]]], tmp$CD.score.raw,method = "spearman")
  res
  })
names(CD.score.raw_MIC) <- diease_Scores

plots <- lapply(diease_Scores, function(m) 
  plot_module_cor_xASm(microbiota3_CD, m, yvar = "CD.score.raw",MIC=CD.score.raw_MIC,size=size_name))
plot <- CombinePlots(c(plots),nrow=4)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.CD.score.raw_raw.cor.all.svg", plot=plot,width = 30, height = 16,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.CD.score.raw_raw.cor.all.png", plot=plot,width = 30, height = 16,dpi=300)
~~~

![Fig5.1.CD.score.raw_raw.cor.all](./Cor_Clincal_micro.assets/Fig5.1.CD.score.raw_raw.cor.all.png)

## 6.3.6 Features vs DAI (W0-delta)

Builds a per-patient baseline-subtracted (delta) feature table (excluding W0), recomputes the MIC function to also return raw correlation rho/p, and correlates each delta feature with DAI change, exported as `Fig5.1.clincal_delta.cor.all`.

~~~R
Sel_features <- c("DAI","UC.score.v2","CD.score.raw","pathogenic.Fungi.v2","Fungal.observed","Fungal.shannon","Bac.shannon","Bac.observed","gut_associated_bac",
        "Candida","Candida_albicans","Bile_Acid_Biosynthesis","Fatty_Acid_Biosynthesis","Health.UP","Dysbiosis.UP","Module1","Module2","Module3","Module4","Butyrate_producing_bac","Valeric_producing_bac")
df1 <- microbiota3
df1$time <- gsub("F","",df1$time)
df1$time.v2 <- as.character(df1$time)
df1$time.v2[df1$time.v2 %in% "W0"] <- "Pre"
df1$time.v2[df1$time.v2 %in% c("W1","W2","W4")] <- "Post"
df1$time.v2[df1$time.v2 %in% "W8"] <- "LTM"
df1$time.v2 <- factor(df1$time.v2, levels=c("Pre","Post","LTM"))
sel_pa <- names(table(df1$patient))[table(df1$patient)>1]
paired.df_ <- lapply(1:length(sel_pa),function(x) {
    tmp <- df1[df1$patient %in% sel_pa[x],]
    All_tmp_ <- lapply(1:length(Sel_features),function(i) {
      tmp_col <- tmp[,Sel_features[[i]]]-tmp[tmp$time=="W0",Sel_features[[i]]]
      tmp_col
      })
    All_tmp <- do.call(cbind,All_tmp_)
    colnames(All_tmp) <- Sel_features
    rownames(All_tmp) <- rownames(tmp)
    abs_d <- abs(All_tmp)
    colnames(abs_d) <- paste0("abs.",colnames(abs_d))
    All_tmp <- as.data.frame(cbind(All_tmp,abs_d))
    message(x)
    return(All_tmp)
})
paired.df <- do.call(rbind,paired.df_)
paired.df1 <- as.data.frame(cbind(paired.df,microbiota3[rownames(paired.df),setdiff(colnames(microbiota3),colnames(paired.df))]))
paired.df1$time <- gsub("F","",paired.df1$time)
paired.df1$time.v2 <- as.character(paired.df1$time)
paired.df1$time.v2[paired.df1$time.v2 %in% "W0"] <- "Pre"
paired.df1$time.v2[paired.df1$time.v2 %in% c("W1","W2","W4")] <- "Post"
paired.df1$time.v2[paired.df1$time.v2 %in% "W8"] <- "LTM"
paired.df1$time.v2 <- factor(paired.df1$time.v2, levels=c("Pre","Post","LTM"))
paired.df1 <- paired.df1[paired.df1$time!="W0",]

sMIC_new <- function(x, y, nperm = 2e4, seed = NULL, method = "spearman",
                 return_null = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]

  rho  <- suppressWarnings(cor(x, y, method = method))
  mic0 <- mine(x, y)$MIC
  if(length(nperm)>1){nperm=sample(nperm,1)}
  null_mic <- vapply(seq_len(nperm), function(i) mine(x, sample(y))$MIC, 0.0)
  k <- sum(null_mic >= mic0)
  p_mic <- (k + 1) / (nperm + 1)   # 永不为 0 的置换 p

  sign_rho <- if (is.na(rho) || rho == 0) 0 else sign(rho)
  cor.rho <- cor(x, y, method = method)
  cor.pval <- cor.test(x, y, method = method)$p.value
  out <- c(sMIC = sign_rho * mic0, MIC = mic0, rho = unname(rho),
           p_MIC = p_mic, k = k, B = nperm, cor.rho = cor.rho, cor.pval = cor.pval)
  if (return_null) attr(out, "null_mic") <- null_mic
  out
}

DAI_MIC_delta <- future_lapply(1:length(diease_Scores),function(x) {
  tmp <- paired.df1[!is.na(paired.df1[[diease_Scores[x]]]),]
  res <- sMIC_new(tmp[[diease_Scores[x]]], tmp$DAI,method = "spearman", nperm = 2e4)
  res
  })
names(DAI_MIC_delta) <- diease_Scores
mcsaveRDS(DAI_MIC_delta, "./projects/MMC/Figures/v2_figures/DAI_MIC_delta.rds", mc.cores=20)

DAI_MIC_delta <- mcreadRDS("./projects/MMC/Figures/v2_figures/DAI_MIC_delta.rds", mc.cores=20)
size_name <- as.list(paste0("abs.",diease_Scores))
names(size_name) <- diease_Scores

plots <- lapply(diease_Scores, function(m) 
  plot_module_cor_xASm(paired.df1, m, yvar = "DAI",MIC=DAI_MIC_delta,size=size_name,color_by="treatment"))
plot <- CombinePlots(c(plots),nrow=4)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.svg", plot=plot,width = 25, height = 16,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.png", plot=plot,width = 25, height = 16,dpi=300)
~~~

![Fig5.1.clincal_delta.cor.all](./Cor_Clincal_micro.assets/Fig5.1.clincal_delta.cor.all.png)

## 6.3.7 Features vs DAI (W0-delta, lm/MIC Line Option)

Recomputes the DAI delta MIC with a smaller permutation budget and redefines `plot_module_cor_xASm` to support either an `lm` fit or a MIC-derived fixed slope, exported as `Fig5.1.clincal_delta.cor.all.keep`.

~~~R

DAI_MIC_delta <- future_lapply(1:length(diease_Scores),function(x) {
  tmp <- paired.df[!is.na(paired.df[[diease_Scores[x]]]),]
  res <- sMIC_new(tmp[[diease_Scores[x]]], tmp$DAI,method = "spearman", nperm = seq(1e2,3e2,by=10))
  res
  })
names(DAI_MIC_delta) <- diease_Scores
mcsaveRDS(DAI_MIC_delta, "./projects/MMC/Figures/v2_figures/DAI_MIC_delta.all.rds", mc.cores=20)

DAI_MIC_delta <- mcreadRDS("./projects/MMC/Figures/v2_figures/DAI_MIC_delta.all.rds", mc.cores=20)
size_name <- as.list(paste0("abs.",diease_Scores))
names(size_name) <- diease_Scores

plot_module_cor_xASm <- function(df, module, yvar = "Candida_albicans", 
                                 MIC = NULL, size = "abs.Ca", 
                                 color_by = "Diagnosis_new", 
                                 line_mode = c("lm", "mic")) {
  line_mode <- match.arg(line_mode)  # 默认参数检查
  pal <- jdb_palette("corona")
  f <- as.formula(paste(yvar, "~", module))
  total.lm <- lm(f, data = df)
  
  # 相关系数和 p 值
  if (is.null(MIC)) {
    rho <- cor(total.lm$model[[module]], total.lm$model[[yvar]], method = "pearson")
    pval <- cor.test(total.lm$model[[module]], total.lm$model[[yvar]], method = "pearson")$p.value
  } else {
    rho <- MIC[[module]][["sMIC"]]
    pval <- MIC[[module]][["p_MIC"]]
  }
  
  # 点大小设置
  if (is.list(size)) {
    size1 <- size[[module]]
  } else {
    size1 <- size
  }
  # 文本位置
  ypos <- max(df[[yvar]], na.rm = TRUE)
  xpos <- median(df[[module]], na.rm = TRUE)
  # 平均点坐标
  xbar <- mean(df[[module]], na.rm = TRUE)
  ybar <- mean(df[[yvar]], na.rm = TRUE)
  # lm 斜率和截距
  b_hat <- coef(total.lm)[2]
  a_hat <- coef(total.lm)[1]
  # MIC 映射斜率（人为设置：|rho| * k，保证方向和 rho 一致）
  slope_mic <- rho * 1   # ⚠️ 这里的 “2” 是缩放因子，可以调
  intercept_mic <- ybar - slope_mic * xbar
    p <- ggplot(df, aes_string(module, yvar, color = color_by, size = size1)) + 
    geom_point(alpha = 0.9) +
    geom_rug(size = 0.4) + 
    scale_color_manual(values = pal) +
    scale_size_continuous(range = c(1, 8), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(x = module, y = yvar, title = paste0("Cor: ", module, " vs ", yvar)) +
    annotate("text", x = xpos, y = ypos, color = "#E51718",
             label = paste0("rho = ", round(rho, 2), 
                            "; p = ", signif(pval, 3)))
  # 添加直线
  if (line_mode == "lm") {
    p <- p + geom_smooth(aes_string(x = module, y = yvar), 
                         method = "lm", se = TRUE,
                         color = "grey20", fill = "grey70", inherit.aes = FALSE)
  } else if (line_mode == "mic") {
    p <- p + geom_abline(slope = slope_mic, intercept = intercept_mic, 
                         color = "red", linetype = "dashed")
  }
  
  return(p)
}

paired.df_tmp <- as.data.frame(cbind(paired.df,microbiota3[rownames(paired.df),setdiff(colnames(microbiota3),colnames(paired.df))]))
plots <- lapply(diease_Scores, function(m) 
  plot_module_cor_xASm(paired.df_tmp, m, yvar = "DAI",MIC=DAI_MIC_delta,size=size_name,color_by="treatment"))
plot <- CombinePlots(c(plots),nrow=4)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.keep.svg", plot=plot,width = 25, height = 16,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.keep.png", plot=plot,width = 25, height = 16,dpi=300)
~~~

![Fig5.1.clincal_delta.cor.all.keep](./../../MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.keep.png)

## 6.3.8 Features vs DAI (W0-delta, Colored by Diagnosis)

Same DAI delta correlations as above, but colored by diagnosis (UC vs CD), exported as `Fig5.1.clincal_delta.cor.all.dieasse`.

~~~R
plots <- lapply(diease_Scores, function(m) 
  plot_module_cor_xASm(paired.df1, m, yvar = "DAI",MIC=DAI_MIC_delta,size=size_name,color_by="Diagnosis_new"))
plot <- CombinePlots(c(plots),nrow=4)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.dieasse.svg", plot=plot,width = 30, height = 16,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.dieasse.png", plot=plot,width = 30, height = 16,dpi=300)
~~~

![Fig5.1.clincal_delta.cor.all.dieasse](./../../MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.dieasse.png)

## 6.3.9 Features vs DAI (W0-delta, MIC Slope with 95% CI)

Defines `plot_module_cor_xASm_new` (MIC-based slope plus a 95% confidence band) and re-plots the DAI delta correlations with this regression band, exported as `Fig5.1.clincal_delta.cor.all.keep`.

~~~R
plot_module_cor_xASm_new <- function(df, module, yvar = "Candida_albicans",
                                 MIC = NULL, size = "abs.Ca",
                                 color_by = "Diagnosis_new",
                                 line_mode = c("lm", "mic"),
                                 mic_ci = TRUE,
                                 scale_factor = 1) {   # 新增参数
  line_mode <- match.arg(line_mode)
  pal <- jdb_palette("corona")
  f <- as.formula(paste(yvar, "~", module))
  total.lm <- lm(f, data = df)

  # Pearson or MIC
  if (is.null(MIC)) {
    rho <- cor(total.lm$model[[module]], total.lm$model[[yvar]], method = "pearson")
    pval <- cor.test(total.lm$model[[module]], total.lm$model[[yvar]], method = "pearson")$p.value
  } else {
    rho <- MIC[[module]][["sMIC"]]
    pval <- MIC[[module]][["p_MIC"]]
  }

  # 文本位置
  ypos <- max(df[[yvar]], na.rm = TRUE)
  xpos <- median(df[[module]], na.rm = TRUE)

  # MIC-based 斜率: rho × sd(y)/sd(x)，再乘以 scale_factor
  slope_mic <- rho * (sd(df[[yvar]], na.rm = TRUE) / sd(df[[module]], na.rm = TRUE)) * scale_factor
  intercept_mic <- mean(df[[yvar]], na.rm = TRUE) - slope_mic * mean(df[[module]], na.rm = TRUE)

  # 点大小
  if (is.list(size)) {
    size1 <- size[[module]]
  } else {
    size1 <- size
  }

  # 基础图
  p <- ggplot(df, aes_string(module, yvar, color = color_by, size = size1)) +
    geom_point(alpha = 0.9) +
    geom_rug(size = 0.4) +
    scale_color_manual(values = pal) +
    scale_size_continuous(range = c(1, 8), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(x = module, y = yvar, title = paste0("Cor: ", module, " vs ", yvar)) +
    annotate("text", x = xpos, y = ypos, color = "#E51718",
             label = paste0("rho = ", round(rho, 2),
                            "; p = ", signif(pval, 3)))

  if (line_mode == "lm") {
    p <- p + geom_smooth(aes_string(x = module, y = yvar),method = "lm", se = TRUE,color = "grey20", fill = "grey70", inherit.aes = FALSE)
  } else if (line_mode == "mic") {
    # 固定斜率直线
    p <- p + geom_abline(slope = slope_mic, intercept = intercept_mic,color = "red", linetype = "dashed")
    # 95% CI 阴影
    if (mic_ci) {
      x_seq <- seq(min(df[[module]], na.rm = TRUE),
                   max(df[[module]], na.rm = TRUE), length.out = 200)
      y_pred <- intercept_mic + slope_mic * x_seq
      
      # 用残差标准差生成随 x 变化的置信区间（更像真实回归带）
      resid_sd <- sd(df[[yvar]] - (intercept_mic + slope_mic * df[[module]]), na.rm = TRUE)
      n <- nrow(df)
      se_fit <- resid_sd * sqrt(1/n + (x_seq - mean(df[[module]], na.rm = TRUE))^2 /
                                   sum((df[[module]] - mean(df[[module]], na.rm = TRUE))^2, na.rm = TRUE))
      
      df_ribbon <- data.frame(
        x = x_seq,
        y = y_pred,
        ymin = y_pred - 1.96 * se_fit,
        ymax = y_pred + 1.96 * se_fit
      )
      
      p <- p + geom_ribbon(data = df_ribbon,
                           aes(x = x, ymin = ymin, ymax = ymax),
                           inherit.aes = FALSE, alpha = 0.2, fill = "grey70")
    }
  }
  return(p)
}
plots <- lapply(diease_Scores, function(m) 
  plot_module_cor_xASm_new(paired.df1, m, yvar = "DAI",MIC=DAI_MIC_delta,size=size_name,color_by="treatment", 
    line_mode="mic", scale_factor=0.8))
plot <- CombinePlots(c(plots),nrow=4)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.keep.svg", plot=plot,width = 25, height = 16,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.clincal_delta.cor.all.keep.png", plot=plot,width = 25, height = 16,dpi=300)
~~~

![Fig5.1.clincal_delta.cor.all.keep](./Cor_Clincal_micro.assets/Fig5.1.clincal_delta.cor.all.keep-1781535736984-50.png)

## 6.3.10 Longitudinal Microbiome Features in UC (Weekly, Multi-test)

For UC patients, plots each feature over weeks (W0–W8) as raw bars plus per-patient delta LOESS trajectories, annotating Fluconazole-vs-Nystatin differences with a battery of tests (ANOVA, ART, Kruskal–Wallis, permutation, Tukey, pairwise t, LMM, LOESS F-test), exported as `Fig5.1.micrbiome_Features.UC`.

~~~R
diease_Scores <- c("Candida_albicans","pathogenic.Fungi.v2","Fungal.observed","Bac.observed","gut_associated_bac",
        "Health.UP","Dysbiosis.UP","Butyrate_producing_bac")
microbiota3 <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/All_omics_values.rds")
microbiota3$Candida_albicans <- 100*microbiota3$Candida_albicans
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
pal <- jdb_palette("corona")
microbiota3_UC <- microbiota3[microbiota3$treatment!="Clotrimazole" & microbiota3$Diagnosis_new=="UC",]
microbiota3_UC <- microbiota3_UC[microbiota3_UC$time!="W6",]
microbiota3_UC$treatment <- factor(microbiota3_UC$treatment,levels=c("Nystatin","Fluconazole"))
microbiota3_UC$UC.score.v2[microbiota3_UC$UC.score.v2 > 9 & microbiota3_UC$Diagnosis_new=="UC"] <- 9
microbiota3_UC$time <- gsub("F","",microbiota3_UC$time)
sel_pa <- names(table(microbiota3_UC$patient))[table(microbiota3_UC$patient)>1]
total_plots1 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_UC
    colnames(df_paired)[which(colnames(df_paired)==diease_Scores[dis])] <- "value"
    bar_data <- df_paired %>%dplyr::group_by(time, treatment) %>%  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    p1 <- ggplot(df_paired, aes_string(x = "time", y = "value")) +
        theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+ 
        geom_bar(data=bar_data,aes(fill=time,color=time),alpha=0.5,stat="identity", position=position_dodge())+
        NoLegend()+labs(title = paste0(diease_Scores[dis]),y ="raw")+
        geom_jitter(color="black",width = 0.1,alpha=0.5)+ stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) + 
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = wilcox.test_wrapper)+facet_wrap(~treatment,ncol=3)+
        stat_summary(fun.y=mean, colour="black", geom="text", size=6,show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=1)))
    return(p1)
    })
library(lme4)
library(ARTool)
library(coin)
total_plots3 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_UC[microbiota3_UC$patient %in% sel_pa,]
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,diease_Scores[dis]]-tmp[tmp$time=="W0",diease_Scores[dis]]
        tmp$value[tmp$value > 500] <- 500
        tmp$value[tmp$value < -500] <- -500
        return(tmp)
    })
    df_paired2 <- do.call(rbind,df_paired1_)
    df_paired2 <- df_paired2[!is.na(df_paired2$value),]
    df_paired2$time_numeric <- as.numeric(gsub("W", "", df_paired2$time))

    art_model <- art(value ~ treatment, data = df_paired2)
    anova_art <- anova(art_model)
    art_p_value <- anova_art$`Pr(>F)`[1]
    art.pvalue0 <- sprintf("%.2e", art_p_value)

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
    # geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess", method.args = list(degree=1),se=TRUE, level = 0.5)+
    geom_smooth(aes(color=treatment,fill=treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.2,span=1.2)+
    stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(0),label = "p.format") +
    theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
    labs(title = paste0(diease_Scores[dis],"\n", "Flu vs Nystatin (ANOVA p: ", anova.pvalue0," | ","art p: ",art.pvalue0,
        "|\n","kruskal p: ", kruskal.pvalue0," | ","LocationTests p: ",LocationTests.pvalue0,
        "|\n","TukeyHSD p: ",TukeyHSD.pvalue1," | ","loess p: ",Ftestp_value,"|\n","pairwise. p: ",pairwise.t.pvalue1," | ",
        "LMM p: ",LMM.pvalue1,")"),y = "Δ")+  NoLegend()
    message(dis)
    return(p2)
    })
plot <- CombinePlots(c(total_plots1,total_plots3),nrow=2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.micrbiome_Features.UC.svg", plot=plot,width = 45, height = 8,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.micrbiome_Features.UC.png", plot=plot,width = 45, height = 8,dpi=300)
~~~

![Fig5.1.micrbiome_Features.UC](./Cor_Clincal_micro.assets/Fig5.1.micrbiome_Features.UC.png)

## 6.3.11 Longitudinal Microbiome Features in CD (Weekly, Multi-test)

Same weekly raw-bar plus delta-LOESS analysis and Fluconazole-vs-Nystatin test battery as above, applied to the CD subset, exported as `Fig5.1.micrbiome_Features.CD`.

~~~r
diease_Scores <- c("Candida_albicans","pathogenic.Fungi.v2","Fungal.observed","Bac.observed","gut_associated_bac",
        "Health.UP","Dysbiosis.UP","Butyrate_producing_bac")
microbiota3 <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/All_omics_values.rds")
microbiota3$Candida_albicans <- 100*microbiota3$Candida_albicans
comb <- list(c("W0","W1"),c("W0","W2"),c("W0","W4"),c("W0","W8"))
pal <- jdb_palette("corona")
microbiota3_CD <- microbiota3[microbiota3$treatment!="Clotrimazole" & microbiota3$Diagnosis_new=="CD",]
microbiota3_CD <- microbiota3_CD[microbiota3_CD$time!="W6",]
microbiota3_CD$treatment <- factor(microbiota3_CD$treatment,levels=c("Nystatin","Fluconazole"))
microbiota3_CD$time <- gsub("F","",microbiota3_CD$time)
sel_pa <- names(table(microbiota3_CD$patient))[table(microbiota3_CD$patient)>1]
total_plots1 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_CD
    colnames(df_paired)[which(colnames(df_paired)==diease_Scores[dis])] <- "value"
    bar_data <- df_paired %>%dplyr::group_by(time, treatment) %>%  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    p1 <- ggplot(df_paired, aes_string(x = "time", y = "value")) +
        theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+ 
        geom_bar(data=bar_data,aes(fill=time,color=time),alpha=0.5,stat="identity", position=position_dodge())+
        NoLegend()+labs(title = paste0(diease_Scores[dis]),y ="raw")+
        geom_jitter(color="black",width = 0.1,alpha=0.5)+ stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) + 
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = wilcox.test_wrapper)+facet_wrap(~treatment,ncol=3)+
        stat_summary(fun.y=mean, colour="black", geom="text", size=6,show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=1)))
    return(p1)
    })
library(lme4)
library(ARTool)
library(coin)
total_plots3 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_CD[microbiota3_CD$patient %in% sel_pa,]
    uniq_patient1 <- unique(df_paired$patient)
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,diease_Scores[dis]]-tmp[tmp$time=="W0",diease_Scores[dis]]
        tmp$value[tmp$value > 500] <- 500
        tmp$value[tmp$value < -500] <- -500
        return(tmp)
    })
    df_paired2 <- do.call(rbind,df_paired1_)
    df_paired2 <- df_paired2[!is.na(df_paired2$value),]
    df_paired2$time_numeric <- as.numeric(gsub("W", "", df_paired2$time))

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
    # geom_smooth(aes(color=treatment,fill=treatment),size = 2,alpha=0.2,method = "loess", method.args = list(degree=1),se=TRUE, level = 0.5)+
    geom_smooth(aes(color=treatment,fill=treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.1,span=1.2)+
    stat_compare_means(aes(group = treatment), method = "wilcox.test", label.y = c(0),label = "p.format") +
    theme_bw()+ scale_color_manual(values = pal[c(2,1)])+scale_fill_manual(values = pal[c(2,1)])+
    labs(title = paste0(diease_Scores[dis],"\n", "Flu vs Nystatin (ANOVA p: ", anova.pvalue0," | ","art p: ",art.pvalue0,
        "|\n","kruskal p: ", kruskal.pvalue0," | ","LocationTests p: ",LocationTests.pvalue0,
        "|\n","TukeyHSD p: ",TukeyHSD.pvalue1," | ","loess p: ",Ftestp_value,"|\n","pairwise. p: ",pairwise.t.pvalue1," | ",
        "LMM p: ",LMM.pvalue1,")"),y = "Δ")+  NoLegend()
    message(dis)
    return(p2)
    })
plot <- CombinePlots(c(total_plots1,total_plots3),nrow=2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.micrbiome_Features.CD.svg", plot=plot,width = 45, height = 8,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.micrbiome_Features.CD.png", plot=plot,width = 45, height = 8,dpi=300)
~~~

![Fig5.1.micrbiome_Features.CD](./Cor_Clincal_micro.assets/Fig5.1.micrbiome_Features.CD.png)

## 6.3.12 UC Features by Treatment Phase (Pre/Post/LTM)

Collapses time into Pre (W0), Post (W1–W4), and LTM (W8) phases and plots raw and delta feature levels per phase for UC, exported as `Fig5.1.microbiomics.time.v2.UC`.

~~~R

diease_Scores <- c("Candida_albicans","pathogenic.Fungi.v2","Fungal.observed","Bac.observed","gut_associated_bac",
        "Health.UP","Dysbiosis.UP","Butyrate_producing_bac")
microbiota3 <- mcreadRDS("./workshop/MMC/sample_info/final_Res/v2/All_omics_values.rds")
pal <- jdb_palette("corona")
microbiota3_UC <- microbiota3[microbiota3$treatment!="Clotrimazole" & microbiota3$Diagnosis_new=="UC",]
microbiota3_UC <- microbiota3_UC[microbiota3_UC$time!="W6",]
microbiota3_UC$treatment <- factor(microbiota3_UC$treatment,levels=c("Nystatin","Fluconazole"))
microbiota3_UC$UC.score.v2[microbiota3_UC$UC.score.v2 > 9 & microbiota3_UC$Diagnosis_new=="UC"] <- 9
microbiota3_UC$time <- gsub("F","",microbiota3_UC$time)
microbiota3_UC$time.v2 <- as.character(microbiota3_UC$time)
microbiota3_UC$time.v2[microbiota3_UC$time.v2 %in% "W0"] <- "Pre"
microbiota3_UC$time.v2[microbiota3_UC$time.v2 %in% c("W1","W2","W4")] <- "Post"
microbiota3_UC$time.v2[microbiota3_UC$time.v2 %in% "W8"] <- "LTM"
microbiota3_UC$time.v2 <- factor(microbiota3_UC$time.v2, levels=c("Pre","Post","LTM"))
sel_pa <- names(table(microbiota3_UC$patient))[table(microbiota3_UC$patient)>1]
pal <- jdb_palette("corona")
comb <- list(c("Pre","Post"),c("Pre","LTM"),c("Pre","LTM"))
total_plots1 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_UC
    colnames(df_paired)[which(colnames(df_paired)==diease_Scores[dis])] <- "value"
    bar_data <- df_paired %>%dplyr::group_by(time.v2, treatment) %>%  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    p1 <- ggplot(df_paired, aes_string(x = "time.v2", y = "value")) +
        theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+ 
        geom_bar(data=bar_data,aes(fill=time.v2,color=time.v2),alpha=0.5,stat="identity", position=position_dodge())+
        NoLegend()+labs(title = paste0(diease_Scores[dis]),y ="raw")+
        geom_jitter(color="black",width = 0.1,alpha=0.5)+ stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) + 
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = wilcox.test_wrapper)+facet_wrap(~treatment,ncol=3)+
        stat_summary(fun.y=mean, colour="black", geom="text", size=6,show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=1)))
    return(p1)
    })
total_plots2 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_UC[microbiota3_UC$patient %in% sel_pa,]
    uniq_patient1 <- unique(df_paired$patient)
    if (diease_Scores[dis]=="Fungal.observed") {thresh=25} else {thresh=500}
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,diease_Scores[dis]]-tmp[tmp$time.v2=="Pre",diease_Scores[dis]]
        tmp$value[tmp$value > thresh] <- thresh
        tmp$value[tmp$value < -thresh] <- -thresh
        return(tmp)
    })
    df_paired1 <- do.call(rbind,df_paired1_)
    bar_data <- df_paired1 %>%dplyr::group_by(time.v2, treatment) %>%  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    p1 <- ggplot(df_paired1, aes_string(x = "time.v2", y = "value")) +
        theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+ 
        geom_bar(data=bar_data,alpha=0.5,stat="identity", position=position_dodge())+
        NoLegend()+labs(title = paste0(diease_Scores[dis]," changes"),y ="Δ")+
        geom_smooth(aes(group = 1, color = treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.2,span=1.2)+
        geom_jitter(color="black",width = 0.1,alpha=0.5)+ stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) + 
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = wilcox.test_wrapper)+facet_wrap(~treatment,ncol=3)+
        stat_summary(fun.y=mean, colour="black", geom="text", size=6,show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))
    return(p1)
    })
plot <- CombinePlots(c(total_plots1,total_plots2),nrow=2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.microbiomics.time.v2.UC.svg", plot=plot,width = 30, height = 8,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.microbiomics.time.v2.UC.png", plot=plot,width = 30, height = 8,dpi=300)
~~~

![Fig5.1.microbiomics.time.v2.UC](./Cor_Clincal_micro.assets/Fig5.1.microbiomics.time.v2.UC.png)

## 6.3.13 CD Features by Treatment Phase (Pre/Post/LTM)

Same Pre/Post/LTM phase comparison of raw and delta feature levels applied to the CD subset, exported as `Fig5.1.microbiomics.time.v2.CD`.

~~~R
microbiota3_CD <- microbiota3[microbiota3$treatment!="Clotrimazole" & microbiota3$Diagnosis_new=="CD",]
microbiota3_CD <- microbiota3_CD[microbiota3_CD$time!="W6",]
microbiota3_CD$treatment <- factor(microbiota3_CD$treatment,levels=c("Nystatin","Fluconazole"))
microbiota3_CD$time <- gsub("F","",microbiota3_CD$time)
microbiota3_CD$time.v2 <- as.character(microbiota3_CD$time)
microbiota3_CD$time.v2[microbiota3_CD$time.v2 %in% "W0"] <- "Pre"
microbiota3_CD$time.v2[microbiota3_CD$time.v2 %in% c("W1","W2","W4")] <- "Post"
microbiota3_CD$time.v2[microbiota3_CD$time.v2 %in% "W8"] <- "LTM"
microbiota3_CD$time.v2 <- factor(microbiota3_CD$time.v2, levels=c("Pre","Post","LTM"))
sel_pa <- names(table(microbiota3_CD$patient))[table(microbiota3_CD$patient)>1]
pal <- jdb_palette("corona")
comb <- list(c("Pre","Post"),c("Pre","LTM"),c("Pre","LTM"))
total_plots1 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_CD
    colnames(df_paired)[which(colnames(df_paired)==diease_Scores[dis])] <- "value"
    bar_data <- df_paired %>%dplyr::group_by(time.v2, treatment) %>%  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    p1 <- ggplot(df_paired, aes_string(x = "time.v2", y = "value")) +
        theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+ 
        geom_bar(data=bar_data,aes(fill=time.v2,color=time.v2),alpha=0.5,stat="identity", position=position_dodge())+
        NoLegend()+labs(title = paste0(diease_Scores[dis]),y ="raw")+
        geom_jitter(color="black",width = 0.1,alpha=0.5)+ stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) + 
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = wilcox.test_wrapper)+facet_wrap(~treatment,ncol=3)+
        stat_summary(fun.y=mean, colour="black", geom="text", size=6,show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=1)))
    return(p1)
    })
total_plots2 <- lapply(1:length(diease_Scores),function(dis) {
    df_paired <- microbiota3_CD[microbiota3_CD$patient %in% sel_pa,]
    uniq_patient1 <- unique(df_paired$patient)
    if (diease_Scores[dis]=="Fungal.observed") {thresh=25} else {thresh=500}
    df_paired1_ <- lapply(1:length(uniq_patient1),function(i) {
        tmp <- df_paired[df_paired$patient %in% uniq_patient1[i],]
        tmp[,"value"] <- tmp[,diease_Scores[dis]]-tmp[tmp$time.v2=="Pre",diease_Scores[dis]]
        tmp$value[tmp$value > thresh] <- thresh
        tmp$value[tmp$value < -thresh] <- -thresh
        return(tmp)
    })
    df_paired1 <- do.call(rbind,df_paired1_)
    bar_data <- df_paired1 %>%dplyr::group_by(time.v2, treatment) %>%  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    p1 <- ggplot(df_paired1, aes_string(x = "time.v2", y = "value")) +
        theme_bw()+ scale_color_manual(values = pal)+scale_fill_manual(values = pal)+ 
        geom_bar(data=bar_data,alpha=0.5,stat="identity", position=position_dodge())+
        NoLegend()+labs(title = paste0(diease_Scores[dis]," changes"),y ="Δ")+
        geom_smooth(aes(group = 1, color = treatment), method = "loess", size = 1.5, se = TRUE,alpha=0.2,span=1.2)+
        geom_jitter(color="black",width = 0.1,alpha=0.5)+ stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) + 
        geom_signif(comparisons = comb,step_increase = 0.1,map_signif_level = FALSE,test = wilcox.test_wrapper)+facet_wrap(~treatment,ncol=3)+
        stat_summary(fun.y=mean, colour="black", geom="text", size=6,show_guide = FALSE,  vjust=-0.7, aes( label=round(..y.., digits=2)))
    return(p1)
    })
plot <- CombinePlots(c(total_plots1,total_plots2),nrow=2)
plot
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.microbiomics.time.v2.CD.svg", plot=plot,width = 30, height = 8,dpi=300)
ggsave("./projects/MMC/Figures/v2_figures/Fig5.1.microbiomics.time.v2.CD.png", plot=plot,width = 30, height = 8,dpi=300)
~~~

![Fig5.1.microbiomics.time.v2.CD](./Cor_Clincal_micro.assets/Fig5.1.microbiomics.time.v2.CD.png)



