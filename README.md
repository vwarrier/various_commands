# Various commands

## This is a repository for various small commands that have been used in analyses

### Command to loop through files in a folder and create files for PRSice

```{R}


file.names <- dir(path = " ", pattern =".TBL")

for(i in 1:length(file.names)){
  data1 = fread(file.names[i], fill = TRUE)
  data2 = data1[,c("MarkerName", "CHROM", "POS", "Allele1", "Allele2", "Effect", "StdErr", "P-value")]
  #Change the command above according to the actual column names
  setnames(data2, 1, "SNP")
  setnames(data2, 2, "CHR")
  setnames(data2, 3, "BP")
  setnames(data2, 4, "A1")
  setnames(data2, 5, "A2")
  setnames(data2, 6, "BETA")
  setnames(data2, 7, "SE")
  setnames(data2, 8, "P")
  data2$A1 = toupper(data2$A1)
  data2$A2 = toupper(data2$A2)
  write.table(data2, file = paste0(i, "forprsice.txt"), row.names = F, col.names = T, quote = F)}

```


### Command to loop through rows in a fam file and create seperate fam file for each individual in the pedigree

```{R}
family = fread("/directory", fill = TRUE)

for (row in 1:nrow(family)){
  i = as.character(family[row,2])
  data = family[row,]
  write.table(data, file = paste0("directory", as.character(i), ".snps.common.vcf.fam"), row.names = F, col.names = F, quote = F)
  write.table(data, file = paste0("directory", as.character(i), ".indels.common.vcf.fam"), row.names = F, col.names = F, quote = F)}
  
  ```
  
  
  ### Binomial sign test usign GWAS
  
  ```bash
  # Binomial sign test script

#./plink --bfile ./files_imputed/g1000_eur --clump ~/ALSPAC/pgssumstats/SQprsice.txt --clump-field P --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 1000 --out SQclumped
 
```

```{R}
library(data.table)

sqclumped = fread("SQclumped.clumped")
sqclumped = sqclumped[,c(1:5)]

sq = fread("~/ALSPAC/pgssumstats/SQprsice.txt")
sqclumped = sq[sq$SNP %in% sqclumped$SNP,]

rm(sq)

empathy = fread("~/ALSPAC/pgssumstats/empathyprsice.txt")

merged = merge(sqclumped, empathy, by = "SNP")

merged$Z.x = merged$BETA.x/merged$SE.x
merged$Z.y = merged$BETA.y/merged$SE.y

merged$Z.y2 = ifelse(merged$A1.x == merged$A1.y, merged$Z.y, "NA")

merged$Z.y2 = ifelse(merged$A1.x == merged$A2.y, (merged$Z.y *-1), merged$Z.y2)

merged2 = merged[!is.na(merged$Z.y2),]

#Binomial sign test
merged$concordant = ifelse(merged$Z.x > 0 & merged$Z.y2 > 0, "C", "D")
merged$concordant = ifelse(merged$Z.x < 0 & merged$Z.y2 < 0, "C", merged$concordant)

concordant = subset(merged, concordant == "C")
discordant = subset(merged, concordant == "D")

C = nrow(concordant)
T = nrow(merged)

binom.test(C,T)
```

