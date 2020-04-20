setwd("/Users/ltozzi/Dropbox (PanLab)/Meta-MDD_metaanalysis_project")

library(meta)
library(tidyverse)
library(grid)
library(plyr)

# Create merged dataframe of clinical and nws
clinical_data=read_csv('clinical_CG_LT.csv') # clinical data only for people included in Yan (2019)
nw_data=read_csv('FC_nws_LT.csv') # average network data
merged=merge(clinical_data, nw_data)

# Labeling networks and sites
nws=c('AUD', 'CER', 'COP', 'DMN', 'DAN', 'FPN', 'MEM', 'SAL', 'SMH', 'SMM', 'SC', 'UNC', 'VAN', 'VIS')
sites=sort(unique(merged$Site))

###### DEFINE FUNCTIONS

# Function to create group mean and std dataframe for a network
create_ES_df_binintvar<-function(df, sites, interestvar, nw){
  
  vars=c('g1_n', 'g1_mean', 'g1_std', 'g2_n', 'g2_mean', 'g2_std')
  nw_df=data.frame(matrix(nrow = length(sites), ncol = length(vars)+1))
  colnames(nw_df) <- c('Site', vars)
  nw_df$Site=sites
  interest1=unique(df[,interestvar])[1]
  interest2=unique(df[,interestvar])[2]
  
  for (site in sites)
  {
    nw_df[nw_df$Site==site,vars[1]]<-length(df[df[,interestvar]==interest1 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[2]]<-mean(df[df[,interestvar]==interest1& df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[3]]<-sd(df[df[,interestvar]==interest1 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[4]]<-length(df[df[,interestvar]==interest2 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[5]]<-mean(df[df[,interestvar]==interest2 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[6]]<-sd(df[df[,interestvar]==interest2 & df$Site==site,nw])
  }
  return(nw_df)
}

# Function to create correlation dataframe for a network
create_ES_df_corr<-function(df, sites, nw, interestvar){
  
  vars=c('n', 'corr')
  nw_df=data.frame(matrix(nrow = length(sites), ncol = length(vars)+1))
  colnames(nw_df) <- c('Site', vars)
  nw_df$Site=sites

  for (site in sites)
  {
    nw_df[nw_df$Site==site,vars[1]]<-length(df[df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[2]]<-cor(df[df$Site==site,nw], df[df$Site==site,interestvar],method = 'spearman', use='complete.obs')
  }
  
  return(nw_df)
}

# Function to run between groups meta-analysis for a network
run_btw_group_meta <- function(dataframe) {
      m.hksj <- metacont(
      g1_n,
      g1_mean,
      g1_std,
      g2_n,
      g2_mean,
      g2_std,
      studlab = dataframe$Site,
      data = dataframe,
      comb.fixed = FALSE,
      comb.random = TRUE,
      method.tau = "REML",
      #choosing REML according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950030/
      hakn = TRUE,
      prediction = TRUE,
      sm = "SMD"
    )
  return(m.hksj)
}

# Function to run between correlation meta-analysis for a network
run_corr_meta <- function(dataframe){
  m.cor<-metacor(cor = corr,
          n = n,
          studlab = dataframe$Site,
          data=dataframe,
          method.tau = "REML",
          hakn = TRUE,
          comb.fixed = FALSE,
          comb.random = TRUE,
          prediction = TRUE,
          sm = "ZCOR")
  return(m.cor)
}

# Function to generate forest plots and add title
drawforplot<-function(metaandata, title)
{
  postscript(file=paste('plots/', title,'.ps',sep=''))
  forest(metaandata, main=title,layout = "RevMan5",
         digits.sd = 2)
  grid.text(nw, .5, .80, gp=gpar(cex=2))
  dev.off()
}

# Function to export meta-analysis details for groups
export_meta_group <-function(metaandata, nw){
  df=data.frame(network = nw)
  df$ngroup1<-sum(metaandata$n.e)
  df$ngroup2<-sum(metaandata$n.c)
  df$effect<-metaandata$TE.random
  df$effect_l<-metaandata$lower.random
  df$effect_u<-metaandata$upper.random
  df$effect_p<-metaandata$pval.random
  df$predict_l<-metaandata$lower.predict
  df$predict_u<-metaandata$upper.predict
  df$tau2<-metaandata$tau2
  df$q<-metaandata$Q
  df$qp<-metaandata$pval.Q
  df$I2<-metaandata$I2
  return(df)
}

# Function to export meta-analysis details for correlations
export_meta_corr <-function(metaandata, nw){
  df=data.frame(network = nw)
  df$ngroup<-sum(metaandata$n)
  df$effect<-metaandata$TE.random
  df$effect_l<-metaandata$lower.random
  df$effect_u<-metaandata$upper.random
  df$effect_p<-metaandata$pval.random
  df$predict_l<-metaandata$lower.predict
  df$predict_u<-metaandata$upper.predict
  df$tau2<-metaandata$tau2
  df$q<-metaandata$Q
  df$qp<-metaandata$pval.Q
  df$I2<-metaandata$I2
  return(df)
}

### AGE AND SEX REGRESSION

# Regress age and sex from data at each site for each network
merged_agesexreg<-merged
for (site in sites){
  sitedata=merged[merged$Site==site,]
  for (nw in nws){
  agesex.mdl<-lm(paste(nw,'~Age+Sex'),data=sitedata)
  merged_agesexreg[merged_agesexreg$Site==site, nw]<-agesex.mdl$residuals
    }
}

# Remove sites with less than 20 people per group
factor_count<-by(merged_agesexreg$Site, merged_agesexreg$Group, count)
sites_keep1<-factor_count[['-1']][factor_count[["-1"]]$freq>=20, 'x']
sites_keep2<-factor_count[['1']][factor_count[["1"]]$freq>=20, 'x']
sites_keep<-intersect(sites_keep1,sites_keep2)

# Run meta-analysis for all networks for comparison between groups
MA_list_hcmdd_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_binintvar(df=merged_agesexreg, sites = sites_keep, interestvar = 'Group', nw=nw)
  MA_list_hcmdd_agesexreg[[nw]]<-run_btw_group_meta(nwdf)
}

# Remove sites with less than 20 values of HAMD
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$HAMD), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with HAMD
MA_list_hamdcorr_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'HAMD', nw=nw)
  MA_list_hamdcorr_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of HAMA
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$HAMA), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with HAMA
MA_list_hamacorr_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'HAMA', nw=nw)
  MA_list_hamacorr_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

###### EXPORT RESULTS

# Draw and save forest plots
for (nw in nws){
  drawforplot(MA_list_hcmdd_agesexreg[[nw]], paste(nw, '_hcmdd_agesexreg', sep=''))
  drawforplot(MA_list_hamdcorr_agesexreg[[nw]], paste(nw, '_hamdcorr_agesexreg', sep=''))
  drawforplot(MA_list_hamacorr_agesexreg[[nw]], paste(nw, '_hamacorr_agesexreg', sep=''))  
}

# Export group meta-analyses
meta_export_hcmdd_agesexreg <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(meta_export_hcmdd_agesexreg)<-colnames(export_meta_group(meta_export_hcmdd_agesexreg[[nw]], nw))

for (nw in nws){
  meta_export_hcmdd_agesexreg<-rbind(meta_export_hcmdd_agesexreg, export_meta_group(MA_list_hcmdd_agesexreg[[nw]], nw))
  }

# Export correlation meta-analyses
meta_export_hamd_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_hamd_agesexreg)<-colnames(export_meta_corr(MA_list_hamdcorr_agesexreg[[nw]], nw))
meta_export_hama_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_hama_agesexreg)<-colnames(export_meta_corr(MA_list_hamacorr_agesexreg[[nw]], nw))

for (nw in nws){
  meta_export_hamd_agesexreg<-rbind(meta_export_hamd_agesexreg, export_meta_corr(MA_list_hamdcorr_agesexreg[[nw]], nw))
  meta_export_hama_agesexreg<-rbind(meta_export_hama_agesexreg, export_meta_corr(MA_list_hamacorr_agesexreg[[nw]], nw))
}

# Save CSVs
write_csv(x = meta_export_hcmdd_agesexreg, path ='tables/summary_hcmdd_agesexreg.csv')
write_csv(x = meta_export_hamd_agesexreg, path ='tables/summary_hamd_agesexreg.csv')
write_csv(x = meta_export_hama_agesexreg, path ='tables/summary_hama_agesexreg.csv')



