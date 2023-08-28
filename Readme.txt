Download the respective GWAS files
eBMD: http://www.gefos.org/?q=content/data-release-2018
T2DM: https://t2d.hugeamp.org/dinspector.html?dataset=Mahajan2022_T2D_EU
2hGlu: https://magicinvestigators.org/downloads/
LD file: https://mrcieu.github.io/ieugwasr/articles/local_ld.html
 


Get significant and independent SNPs from Exposure, and retain only SNPs that we would use for DM and 2hGlu. (Linux)

#eBMD

->get significant SNPs with P.NI<5e-8
awk 'NR==1 {print $0} NR>1 {if ($13 < 5e-8) print $0}' <(zcat eBMD.txt.gz ) > eBMD_sig.txt

->get independent SNPs, 10000Kb 0.001r2
Run clump.R

#DM

->get a smaller DM file containing all SNPs in eBMD_sig.txt
awk 'NR==FNR {a[$2]; next} FNR==1 {print $0" ncase ncontrol"} FNR > 1 && $4 in a {print $0" 74124 824006"}' eBMD_sig.txt T2DM.txt > DM_short.txt

#2hGlu

->get a smaller 2hGlu file containing all SNPs in eBMD_sig.txt
awk 'NR==FNR {a[$2]; next} FNR==1 || FNR>1 && $1 in a' eBMD_sig.txt <(zcat 2hGlu.tsv.gz ) > 2hglu.txt


Run Rscripts 2hGlu.R and T2DM.R 

