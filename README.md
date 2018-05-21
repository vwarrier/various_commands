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
  

