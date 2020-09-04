### Example for comparing two runs in the same project
### Oligo level count data as example

### read in count tables

duo_counts_2 <- read.delim("path_to_count_table_1", row.names = 1)
duo_counts_3 <- read.delim("path_to_count_table_2", row.names = 1)

### Add in column names if needed
colnames(duo_counts_2) <- c("library", c(replicate_ids), "plmean", "PlasmidsBCsum")
colnames(duo_counts_3) <- c("library", c(replicate_ids), "plmean", "PlasmidsBCsum")

### Generate condition matrix for the data - needs cell type
duo_condition <- as.data.frame(c(rep("LCL", 4), rep("DNA", 4)), stringsAsFactors=F)
colnames(duo_condition) <- "condition"
rownames(duo_condition) <- colnames(duo_counts_2)[2:9]

### Run the steps of the pipeline
names_list <- duoNames(duo_counts_2, duo_counts_3, libExcl1 = "library to exclude 1", libExcl2 = "library to exclude 2")

neg_list_e <- c("negative element 1", ... , "negative element n")
neg_list_s <- names_list[[element]][grepl("negative pattern", names_list[[element]])]

duo_2_DE <- duoSeq(duo_counts_2, duo_condition, 2, "project_name", namesList = names_list, libExcl = "library to exclude 1", negListE = neg_list_e, negListS = neg_list_s)
duo_3_DE <- duoSeq(duo_counts_3, duo_condition, 3, "project_name", namesList = names_list, libExcl = "library to exclude 2", negListE = neg_list_e, negListS = neg_list_s)

duo_norm_counts <- duoNorm(duo_2_DE, duo_3_DE, duo_condition, "original_barcodes", pairsList = c("S_2", "S_3", "S_2", "S_3", "SE_2", "ES_3"))


### Example for analyzing a single run within the project
### Barcode level count data as input
### Duo Only

# Read in count table and create condition table

count_init <- read.delim("path_to_count_table",stringsAsFactors=F)
attr_init <- read.delim("path_to_attributes_table", stringsAsFactors = F)

condition_table <- as.data.frame(c(rep("DNA",4), rep("GM12878",4), rep("K562",4), rep("HepG2",4), rep("SK.N.SH", 4)), stringsAsFactors=F)
colnames(condition_table) <- "condition"
rownames(condition_table) <- colnames(RESTscreen_14082020)[4:23]

## Convert to oligo level count data
oligo_count_table <- duoStats(count_init,condition_table)

## Get common names to include and filter 
oligo_names <- duoNames(oligo_count_table,oligo_count_table,libExcl1 = c("E","S"), libExcl2 = c("E","S"),duoOnly = T)

## Run DESeq analysis
ENC <- "En02"
all_names_split <- colsplit(oligo_names$ES, pattern = "\\^", names = c("enhancer","silencer"))
SNC <- attr_init$ID[which(attr_init$project=="NegCtrl")]

enh_all <- c("En02", "En09", "En11", "En19", "En21")

oligo_attr <- duoAttr(attr_init, enhList = enh_all)

oligo_seq <- duoSeq(oligo_attr, oligo_count_table, condition_table, 1, "example",oligo_names, negListE = ENC, negListS = SNC, libExcl = c("E","S"))

## Count Table correlation (cell types)
duoCor(dataCount=oligo_count_table,dataCond=condition_table,namesList = oligo_names,filePrefix = "example",run=1,libExcl = c("E","S"), libIncl = "ES_1")

## log2FoldChange Correlation
E02_randoms <- rownames(results(oligo_seq$ES_1))[grepl("En02\\^random", rownames(results(oligo_seq$ES_1)))]
all_randoms <- rownames(results(oligo_seq$ES_1))[grepl("random", rownames(results(oligo_seq$ES_1)))]
E02_all <- rownames(results(oligo_seq$ES_1))[grepl("En02", rownames(results(oligo_seq$ES_1)))]
E02_positive <- E02_all[grepl("Pos", E02_all)]
E02_negative <- E02_all[grepl("Neg", E02_all)]

duoLogCor(oligo_seq, condition_table, "En02^negative", E02_negative)
duoLogCor(oligo_seq, condition_table, "En02^positives", E02_positive)
duoLogCor(oligo_seq, condition_table, "all_randoms", all_randoms)
duoLogCor(oligo_seq, condition_table, "En02^random", E02_randoms)
