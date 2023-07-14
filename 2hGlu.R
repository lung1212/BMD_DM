library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)
library(dplyr)
set.seed(123)

#read exposure and outcome datasets
exposure<-read_exposure_data(filename = "eBMD_clumped.txt",
                             sep = "\t",
                             snp_col = "RSID",
                             beta_col = "BETA",
                             se_col = "SE",
                             effect_allele_col = "EA",
                             other_allele_col = "NEA",
                             eaf_col = "EAF",
                             pval_col = "P.NI",
                             samplesize_col = "N",
)
outcome<-read_outcome_data(filename = "2hglu.txt",
                           sep = "\t",
                           snp_col = "variant",
                           beta_col = "beta",
                           se_col = "standard_error",
                           effect_allele_col = "effect_allele",
                           other_allele_col = "other_allele",
                           eaf_col = "effect_allele_frequency",
                           pval_col = "p_value",
                           samplesize_col = "sample_size"
                           
)

#harmonise
dat<-harmonise_data(exposure, outcome, action = 2)

#find SNPs to proxy
if(nrow(dat)>=1){dat<-subset(dat, dat$mr_keep==T)}
exposure_DF<-as.data.frame(exposure)
prox_tot<-subset(exposure_DF, !exposure_DF$SNP %in% dat$SNP)

#read all significant exposures as dataset to identify proxies from
exposure2<-read_exposure_data(filename = "eBMD_sig.txt",
                              sep = "\t",
                              snp_col = "RSID",
                              beta_col = "BETA",
                              se_col = "SE",
                              effect_allele_col = "EA",
                              other_allele_col = "NEA",
                              eaf_col = "EAF",
                              pval_col = "P.NI",
                              samplesize_col = "N",
)


#6 SNPs duplicated, make sure the required A1/A2 combo appears first, as read_exposure_data only keeps the first instance. In this case, only one SNP is not in the outcome.
#rs11580218 TG
#rs491616   GG 
#rs2077111  AG 
#rs11773399 AC
#rs10858956 CG 
#rs2696531  NA

#find proxies, filter by R2 and Pvalue, then retain only those which are not palindromic.
proxies<-data.frame()
for(j in 1:nrow(prox_tot)){
  x<-system(paste0("plink --bfile ./1KGP/EUR --r2 yes-really --ld-snp ",prox_tot$SNP[j]," --ld-window 99999 --ld-window-r2 0.8"))
  y<-try(read.table("plink.ld"))
  if(is(y, 'try-error')) next
  system("rm plink.ld")
  colnames(y)<-y[1,]
  y<-y[-1,]
  y$R2<-as.numeric(y$R2)
  if (ncol(y)>1){
    y<-subset(y, y$SNP_B %in% exposure2$SNP)
    y<-subset(y, !y$SNP_B %in% exposure$SNP)
    if(nrow(y)==0){next}
    prox<-merge(exposure2, y,by.x="SNP", by.y = "SNP_B")
    dat_proxy<-harmonise_data(prox, outcome, action = 2)
    if(nrow(dat_proxy)>=1){
      dat_proxy<-subset(dat_proxy, dat_proxy$mr_keep==T)
      dat_proxy$R2<-as.numeric(dat_proxy$R2)
      dat_proxy<-dat_proxy[order(-dat_proxy$R2, dat_proxy$pval.exposure),]
      proxies<-rbind(proxies, dat_proxy[1,])
    }
  }
}
dat<-rbind(dat, proxies[,c(intersect(colnames(dat), colnames(proxies)))])
dat$id.exposure<-"eBMD"

#remove outliers using radial MR
library(RadialMR)
test_f<-format_radial(BXG=dat$beta.exposure, BYG=dat$beta.outcome, seBXG=dat$se.exposure, seBYG=dat$se.outcome, RSID=dat$SNP)
test_IVWradial<-ivw_radial(test_f, alpha=0.05, weights=3)
dat_radial<-dat[(!dat$SNP %in% test_IVWradial$outliers$SNP),]

#identify SNPs associated with waist circumference, waist to hip ratio, Coronary disease, BMI or LDL-C
library(phenoscanner)
res <- phenoscanner(snpquery=dat_radial$SNP[1:100], pvalue = 5E-8, proxies = "EUR", r2 = 0.8, build = 37)
res1 <- phenoscanner(snpquery=dat_radial$SNP[101:200], pvalue = 5E-8, proxies = "EUR", r2 = 0.8, build = 37)
res2 <- phenoscanner(snpquery=dat_radial$SNP[201:300], pvalue = 5E-8, proxies = "EUR", r2 = 0.8, build = 37)
res3 <- phenoscanner(snpquery=dat_radial$SNP[301:400], pvalue = 5E-8, proxies = "EUR", r2 = 0.8, build = 37)
res4 <- phenoscanner(snpquery=dat_radial$SNP[401:500], pvalue = 5E-8, proxies = "EUR", r2 = 0.8, build = 37)
res5 <- phenoscanner(snpquery=dat_radial$SNP[501:nrow(dat_radial)], pvalue = 5E-8, proxies = "EUR", r2 = 0.8, build = 37)
pheno<-rbind(res$results, res1$results, res2$results, res3$result, res4$results, res5$result)
  
#check for variations in trait names, identify terms that encompass only relevant traits

pheno<-subset(pheno, pheno$ancestry == "European" & as.numeric(pheno$p)<5e-8)
pheno<-subset(pheno, grepl("ldl", trait, ignore.case=T) |
                grepl("waist", trait, ignore.case=T) |
                grepl("coronary", trait, ignore.case=T) | 
                grepl("Body mass index", trait, ignore.case=T)|
                grepl("Low density lipoprotein", trait, ignore.case=T))

#remove SNPs from phenoscanner to create final dataset
dat_final<-subset(dat_radial, !dat_radial$SNP %in% pheno$snp) 

#conduct analyses
results_main<-mr(dat_final)
egger_intercept<-mr_pleiotropy_test(dat_final)
het<-mr_heterogeneity(dat_final)
library(MendelianRandomization)
results_conmix<-mr_conmix(mr_input(bx=dat_final$beta.exposure, by=dat_final$beta.outcome, bxse=dat_final$se.exposure, byse=dat_final$se.outcome), psi=0, CIMin=-5, CIMax=5, CIStep=0.01)

results_main$UCI<-results_main$b+(1.96*results_main$se)
results_main$LCI<-results_main$b-(1.96*results_main$se)

#calculate R2
dat_final2<-add_rsq(dat_final)
r2<-sum(dat_final2$rsq.exposure)

#output IV information
write.table(dat_final, "2hGLU_IVs.txt", sep="\t", quote = F, col.names = T, row.names = F)