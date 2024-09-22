## comparing ES genes in mouse TB model and cattle bovis model

## import excel file from Bellerose
br_data <- read_csv(here("data/Bellarose_attenuating mutations.csv"), skip=9) 
# 574 but last 12 are NA
br_list <- br_data[1:562,] %>% select(gene) %>% pull()
length(br_list)
#562
#change formatting to match our data ('RVBD_0007' to 'Rv0007')
br_list <- sub("RVBD_", "Rv", br_list)

## import our supp file to get list of ess bovis genes, identical/variable
ess_data <- read_csv(here("Output/supp_table_08_02.csv"), skip=1)
nrow(ess_data)
#141

## import all_genes table to get info for genes that were NOT attenuated in bovis
all_data <- read_csv(here("Output/supp_table_all_genes_08_02.csv"), skip=1, name_repair = "unique")
all_data <- all_data %>% select(-c("...9", "...12"))
colnames(all_data) <- c("M. bovis locus","H37Rv Ortholog","Name","Number of TA sites", 
                        "Essentiality in input library",	"Functional category", 
                        "Average log2 Fold Change in lung",
                        "No of cattle with significant change in lung", 
                        "Average log2 Fold Change	in nodes", 
                        "No of cattle with significant change in nodes", 
                        "Average log2 Fold Change (both lungs and nodes)") 
#ess in our data and in bellarose
ess_data %>% dplyr::filter(`H37Rv Ortholog` %in% br_list)
#86
# this doesn't match Sharon's number of 109?

# in bellarose and in our data
in_supp <- br_list[br_list %in% ess_data$`H37Rv Ortholog`]
length(in_supp)
#86

# are there any duplicates in essential data?
ess_data %>% group_by(`H37Rv Ortholog`) %>% filter(n()>1)
# no duplicate Rv numbers

# in bellarose but not in our data
not_in_supp <- br_list[!br_list %in% ess_data$`H37Rv Ortholog`]
length(not_in_supp)
#476 not in our list of essential genes

# complete info for genes in mouse list but not in ours
# read in malone ortholog table for identical/variable
orthologs <- read_csv(here("data/mb_mtb_orthologs_1.txt"))
colnames(orthologs) <- c("bovis", "mtb", "id_var")
# add i/v to list of mouse essential not in bovis essential
not_in_df <- as_tibble(not_in_supp)
not_in_df <- not_in_df %>% mutate(snps = ifelse(value %in% orthologs$mtb, orthologs$id_var, NA))
colnames(not_in_df) <- c("Mtb", "id_var")
nrow(not_in_df)
#476


# do with complete list
not_in_df <- as_tibble(not_in_supp)
orthologs1 <- read_csv(here("data/Malone_orthologs.csv"), skip=2)
not_in_df2 <- not_in_df %>% mutate(snps = ifelse(value %in% orthologs1$`Rv tag`, orthologs1$Orthology, NA))
not_in_df2 %>% dplyr::filter(is.na(snps))
#12 missing, tb locus (Rv) not listed in Malone spreadsheet, these don't seem to have bovis orthologs 
colnames(not_in_df2) <- c("Mtb", "id_var")
nrow(not_in_df)
#476

# add in our logfold change
final_df <- left_join(not_in_df2, all_data, by=c("Mtb" = "H37Rv Ortholog"), na_matches="never")
nrow(final_df)
#481 
#there are 5 duplicate Rv numbers in this list, because associated with two Mbovis loci
final_df %>% group_by(Mtb) %>% 
  dplyr::filter(n()>1)

View(final_df)

# several don't have id_var values, these are in original spreadsheet, so not sure why?
final_df %>% dplyr::filter(is.na(id_var))
#12 missing, no bovis gene associated with it?
#one not in any tb annotation I can find, RVBD_0454Ac. typo? or is it 0454?

#clean up final column (too many ,,,,,,,)

#how many don't have mbovis loci
final_df %>% dplyr::filter(is.na(`M. bovis locus`) & !is.na(id_var))
#6 not in our list of tested genes
#Rv0590A could be Rv0590 instead in ours
#Rv0616 could be Rv0616A in ours or Rv0616c
#Rv0924c/MB0948c missing, we have Rv0925c associated with this gene (both in malone)
#Rv2098c/MB2125c missing, we have Rv2099c assoc with gene (both in malone)
#Rv2160A/MB2184c missing, we have Rv2160c (also associated with MB2184c, both in Malone)
#Rv3117/MB3144  missing, we have Rv3121 (both in Malone)
#all non-attenuating

final_df$`Average log2 Fold Change (both lungs and nodes)`<- sub(",,,,,,,,,,,,,,,", "", final_df$`Average log2 Fold Change (both lungs and nodes)`)

length(unique(final_df$`M. bovis locus`))
#464

#write csv
write_excel_csv(final_df, here("lung_tnseq/Output/bellarose_compare.csv"))
