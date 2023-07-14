library(TwoSampleMR)
library(ieugwasr)
library(genetics.binaRies)

  dat<-read.csv("eBMD_sig.txt", sep="\t", header=T)
  dat$id<-"id"
  
  clump<-ld_clump(
    dplyr::tibble(rsid=dat$RSID, pval=dat$P.NI, id=dat$id),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "./1KGP/EUR",
    clump_kb=10000,
    clump_r2=0.001
  )
  dat2<-subset(dat, dat$RSID %in% clump$rsid)
  write.table(dat2, "eBMD_clumped.txt", row.names = F, col.names = T, sep="\t", quote=F)

