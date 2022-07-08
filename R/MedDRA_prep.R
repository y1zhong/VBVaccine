# library(dplyr)
# 
# MedDRA <- read_csv("~/Dropbox (University of Michigan)/BayesianVaccine/Rcode/yuan/data/final_MedDRA.csv")
# MovePTterms <- read_xlsx("~/Dropbox (University of Michigan)/BayesianVaccine/Rcode/yuan/data/MovePTterms.xlsx", sheet = "Other") %>%
#   rename(PT = AE_NAME) %>%
#   mutate(HLGT = "Others",
#          soc = "Others")
# 
# MedDRA <- MedDRA %>%
#   select(soc, HLGT, PT) %>%
#   add_row(soc="Neoplasms benign, malignant and unspecified (incl cysts and polyps)", HLGT="Cancer", PT="Cancer")%>%
#   distinct()%>%
#   rbind(MovePTterms)
# save(MedDRA, file = "data/MedDRA.RData")
# 
# 
# neg_control<- read_xlsx("~/Dropbox (University of Michigan)/BayesianVaccine/Rcode/yuan/data/Negative controls Final.xlsx",
#                         sheet = "Final35") %>%
#   filter(PT != "Hypoxia") %>%
#   mutate(PT= replace(PT, PT == "Magnetic resonance imaging", "Nuclear magnetic resonance imaging"))
# save(neg_control, file = "data/neg_control.RData")
