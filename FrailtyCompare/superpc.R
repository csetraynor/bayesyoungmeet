#Superpca
library(superpc)
library(tidyverse)

#obtain gene expression data
gene_expression <- read_table2("brca_metabric/data_expression.txt", col_names = TRUE)

#do you need to sample genes? 
genes <- gene_expression[sample(nrow(gene_expression), 100), ]

#imput (with knn) 
require(impute)

genes_impute <- genes %>%
  select(-Entrez_Gene_Id, -Hugo_Symbol) %>%
  as.matrix()
genes_impute <- impute.knn(genes_impute, k = 10)$data

#or drop NA values

#obtain id
id_gene_expression <- gsub("-", ".", colnames(genes)[3:1906])

#from clinical data keep observations with gene expression

#rownames to a column sample
metabric_clinicaldata <- rownames_to_column(metabric_clinicaldata, var = "sample")

#select outcome variables and id
metabric_clinicaldata <- metabric_clinicaldata %>%
  filter(!is.na(OS_STATUS)) %>%   
  mutate(IsDeceased = (OS_STATUS == "DECEASED")) %>% # convert to boolean
  select(OS_MONTHS, IsDeceased, sample) %>% 
  rename(time = OS_MONTHS, status = IsDeceased) 

#obtain sample with gene expression measurements and match 
Y <- metabric_clinicaldata[metabric_clinicaldata$sample %in% id_gene_expression,] 
Y <- arrange(Y,  match(Y$sample, id_gene_expression))

surv_time <- Y$time
surv_status <- Y$status
surv_id <- Y$sample

all(id_gene_expression == surv_id)

#prepare for spca
data <- list( x=genes_impute, y = surv_time, censoring.status = surv_status, 
              featurenames = genes$Hugo_Symbol)

#train superpc
fitsuperpc <- superpc.train(data, type="survival")
