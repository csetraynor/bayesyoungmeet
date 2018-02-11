#get data from cbioportal 

require(cgdsr)

mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)

listofcancerstudies = getCancerStudies(mycgds)
#View(listofcancerstudies)

#get  study id
glioblastome_2013_id_sutdy = getCancerStudies(mycgds)[55,1]
glioblastome_2008_id_sutdy = getCancerStudies(mycgds)[56,1]
glioblastome_2013_case_list = getCaseLists(mycgds, glioblastome_2013_id_sutdy)[2,1]

#get avialable genetic profiles
glioblastome_2013_geneticprofile = getGeneticProfiles(mycgds, glioblastome_2013_id_sutdy)[4,1]

#get data slices for a specified list of genes, genetic profile and case list
gene_data_slices = getProfileData(mycgds, c('PIK3CA', 'TP53'), glioblastome_2013_geneticprofile, glioblastome_2013_case_list)

#Get Clinical Data for the case list
glioblastome_2013_clinical_data <-  getClinicalData(mycgds, glioblastome_2013_case_list)
#documentation
#help('cgdsr')

#write.csv(glioblastome_2013_clinical_data, "glioblastome2013.csv")

#transform data to tibble
library(tidyverse)
g_clin_tib <- as_tibble(glioblastome_2013_clinical_data)
g_clin_tib[g_clin_tib == "" | g_clin_tib == " "] <- NA
g_clin_tib <- rownames_to_column(g_clin_tib, var = 'sample')


  

