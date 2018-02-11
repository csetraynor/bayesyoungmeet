#Download Colorectal Metastasic Cancer Dataset

#get data from cbioportal 

require(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)

listofcancerstudies = getCancerStudies(mycgds)
#View(listofcancerstudies)

#get study id
id_sutdy = getCancerStudies(mycgds)[135,1] #targeted 
#id_sutdy_history = getCancerStudies(mycgds)[56,1] Get another study for power prior?

#get case list
case_list = getCaseLists(mycgds, id_sutdy)[2,1] #All tumours include also NGS, and CNA
#case_list_ngs = getCaseLists(mycgds, id_sutdy)[2,2] #substract only NGS?

#get avialable genetic profiles
geneticprofile = getGeneticProfiles(mycgds, id_sutdy)[2,1]

#get data slices for a specified list of genes, genetic profile and case list
gene_mutated_genes = getMutationData(x = mycgds, genes =  c("TP53", "HIST1H3C", "HIST1H3B", "H3F3C", "KMT2D", "CDKN2A", "NOTCH1", "TERT", "SOX2", "MYCN", "MYC", "PLK2", "FGFR3", "NOTCH3", "PHOX2B", "MYOD1", "MEN1", "MCL1", "CEBPA", "SOX17", "CTLA4", "NFKBIA", "MSH6", "CD79B", "FAT1", "NKX2-1", "MED12", "NKX3-1", "NFE2L2", "RPS6KA4", "PIK3CG", "PIK3CA", "PNRC1", "BAP1", "AURKA", "TBX3", "JAK3", "IRS2", "ERCC5", "ARID1A", "PIM1", "SF3B1", "FAM46C", "IFNGR1", "ERBB3", "HNF1A", "FGF19", "MEF2B", "NRAS", "ASXL1", "RECQL4", "MYCL", "KDR", "CIC", "GATA2", "FLT4", "SPEN", "LATS1", "STK11", "CASP8", "GATA1", "PARP1", "ERCC4", "ROS1", "PTCH1", "AXIN2", "PDGFRB", "POLD1", "MPL", "PRDM1", "MAPK3", "ATR", "IDH1", "IRF4", "HGF", "TSC2", "RET", "NOTCH2", "EWSR1", "SMARCA4", "FGFR4", "IL7R", "MUTYH", "DICER1", "BRCA2", "RARA", "GLI1", "BCL6", "STAT5A", "CCNE1", "MDM2", "EP300", "TGFBR1", "PMS2", "POLE", "NBN", "NTRK1", "CCND1", "CTNNB1", "RBM10", "AXL", "IKZF1", "ATM", "RAC1", "DNAJB1", "KMT2C", "CDK12", "CD276", "RAD54L", "SMARCD1", "CREBBP", "TET1", "BRD4", "VEGFA", "VTCN1", "AXIN1", "DOT1L", "MST1R", "CSF3R", "ARID2", "EPCAM", "BCOR", "KMT2A", "PPP2R1A", "KDM6A", "NSD1", "PALB2", "FGFR1", "FLT3", "CSF1R", "BRCA1", "IGF2", "ERCC2", "SETD2", "TCF3", "DDR2", "KDM5C", "BIRC3", "EPHA5", "FH", "MTOR", "TET2", "TRAF7", "IRS1", "RICTOR", "CARD11", "STK40", "DNMT3B", "KRAS", "FGFR2", "MAX", "SYK", "CDH1", "SMO", "RASA1", "RB1", "CHEK1", "FANCA", "DIS3", "SH2D1A", "MAP3K1", "XIAP", "SMAD4", "KIT", "APC", "NF1", "PIK3R1", "ATF1", "MLH1", "CUL3", "ASXL2", "TGFBR2", "RAD50", "MAP3K13", "GATA3", "EPHA7", "FLCN", "INSR", "SRC", "PTPN11", "AR", "NCOR1", "MET", "NF2", "SUZ12", "FLT1", "RUNX1", "RNF43", "ETV1", "PTPRS", "ZFHX3", "GRIN2A", "CRKL", "KDM5A", "ERCC3", "INPP4A", "EPHA3", "RAD51", "EGFR", "STAT5B", "EZH2", "RAD52", "ERBB2", "BARD1", "PAK5", "ABL1", "TMPRSS2", "NTRK3", "LATS2", "LMO1", "PIK3CB", "BRIP1", "SDHA", "PMS1", "PBRM1", "TOP1", "ARID5B", "MGA", "PIK3CD", "TP63", "PAK1", "HRAS", "RHOA", "FBXW7", "TSC1", "MAPK1", "CHEK2", "DNMT3A", "AKT2", "BMPR1A", "JAK1", "XPO1", "DNMT1", "PIK3C3", "PDPK1", "SUFU", "CDC73", "PDGFRA", "E2F3", "ATRX", "STAG2", "JAK2", "ARID1B", "STAT3", "MITF", "NCOA3", "RAF1", "ETV6", "RPTOR", "FIP1L1", "IGF1", "COP1", "YES1", "TSHR", "EPHB1", "MRE11", "SMAD2", "BLM", "BRAF", "ALK", "FOXP1", "ERBB4", "PIK3C2G", "NTRK2", "MAP2K4", "YAP1", "PTPRT", "INPP4B", "IGF1R", "PTPRD", "PAX5", "FYN", "TCF7L2", "ESR1", "MSH2", "GSK3B", "ERG", "PRKN", "NEGR1", "RAD51B", "MDC1", "NOTCH4", "DAXX", "PTEN", "GREM1", "B2M", "AKT3", "PDCD1"), geneticProfile = geneticprofile, caseList = case_list)

#Get Clinical Data for the case list
clinical_data <-  getClinicalData(mycgds, case_list)
#documentation
#help('cgdsr')

