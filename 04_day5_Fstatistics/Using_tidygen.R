## Examples with tidypopgen

# call libraries 
library(tidypopgen)
library(admixtools)
library(ggplot2)

# define working directory 
setwd("~/embo_popgen_2025_figt/Andrea_Manica/")

# look at files in data directory
dir("./data")

# convert files into a gen_tibble
modern_gt <- tidypopgen::gen_tibble(
  "./data/modern_samples.bed",
  valid_alleles = c("A","T","C","G"),
  missing_alleles = c("X")
  )

# inspect data
modern_gt

modern_gt %>% group_by(population) %>% tally()

# do QC on loci
loci_report <- modern_gt %>%
  group_by(population) %>%
  qc_report_loci()

autoplot(loci_report, type = "missing")

# remove loci that have missingness above 4%
modern_gt <- modern_gt %>% 
  select_loci_if(loci_missingness(genotypes)<0.04)

# read in ancient data
ancient_gt <- tidypopgen::gen_tibble(
  "./data/ancient_samples.vcf",
   valid_alleles = c("A","T","C","G","X"), # where X is unknown 
   quiet = TRUE
  )

# representing as pseudohaploid
ancient_gt <- gt_pseudohaploid(ancient_gt)

# look at IDs
ancient_gt$id

# rename to be more informative
ancient_gt$id[ancient_gt$id == "GB20"] <- "Mota"

# set populations
ancient_gt$population <- ancient_gt$id

# we need to check if it is compatible to merge our datasets: to a test run
merged_dry <- rbind_dry_run(
  modern_gt, ancient_gt, 
  flip_strand = TRUE
  )

# we still have a large overlap, so merge data
merged_gt <- rbind(modern_gt, ancient_gt, 
                     flip_strand = TRUE,
                     backingfile = "./data/merged_samples")

merged_gt <- merged_gt %>% group_by(population)
# data is now cleaned and merged! 

f2_dir <- "./data/f2_tidypopgen"

# compute f2
f2_tidypopgen <- gt_extract_f2(
  merged_gt, 
  outdir = "./data/f2_tidypopgen", 
  overwrite = TRUE
  )

f2_blocks = f2_from_precomp("./data/f2_tidypopgen")

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
  "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")

lbk_f3out <- f3(data = f2_blocks, 
                pop1 = "Mbuti", # the outgroup
                pop2 = "LBK",   # the population of interest
                pop3 = lbk_modern_panel) # the sister population
lbk_f3out

lbk_f3out %>% arrange(desc(est)) # order by the amount of shared drift


ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# order factor levels 
lbk_f3out$pop3<-factor(lbk_f3out$pop3, levels = lbk_f3out$pop3[order(lbk_f3out$est)])

# plot this with the ordered levels 
ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# looking at MA1 - run f3
ma1_f3out <- f3(
  data = f2_blocks, 
  pop1 = "Mbuti", 
  pop2 = "MA1", 
  pop3 = lbk_modern_panel
)

ma1_f3out %>% arrange(desc(est)) # arrange in order of shared drift (descending)

# order levels
ma1_f3out$pop3<-factor(ma1_f3out$pop3, levels = ma1_f3out$pop3[order(ma1_f3out$est)])

# plot this with the ordered levels 
ggplot(ma1_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# Now moving to the admixed population ----------------------------------------
aa_f3admix <- f3(data = f2_blocks,
                 pop1 = "AA",
                 pop2 = "Yoruba",
                 pop3 = "French")
aa_f3admix

# define eurasian samples
eurasian_sources <- c("French","Spanish","Sardinian","LBK")

somali_f3admix <- f3(data = f2_blocks,
                     pop1 = "Somali",
                     pop2 = eurasian_sources, 
                     pop3 = "Mota")
somali_f3admix

somali_f3admix %>% arrange(est)
somali_f3admix$pop3<-factor(somali_f3admix$pop3, levels = somali_f3admix$pop3[order(somali_f3admix$est)])

# plot this with the ordered levels 
ggplot(somali_f3admix, aes(pop2, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with Mota", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# repeat for Dinka 
dinka_f3admix <- f3(data = f2_blocks,
                     pop1 = "Dinka",
                     pop2 = eurasian_sources, 
                     pop3 = "Mota")

dinka_f3admix$pop3<-factor(dinka_f3admix$pop3, levels = dinka_f3admix$pop3[order(dinka_f3admix$est)])

# plot this with the ordered levels 
ggplot(dinka_f3admix, aes(pop2, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with Mota", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# -------------------------------------

# Looking at neandertals 
neand_french_f4 <- f4(data = f2_blocks,
                     pop1 = "pan_troglodytes",
                     pop2 = "AltaiNea",
                     pop3 = "Mbuti",
                     pop4 = "French")
neand_french_f4 # we can see an excess of neanderthal ancestry in the French

neand_han_f4 <- f4(data = f2_blocks,
                      pop1 = "pan_troglodytes",
                      pop2 = "AltaiNea",
                      pop3 = "Mbuti",
                      pop4 = "Han")
neand_han_f4 

# Looking at the extent/presence of Yamnaya gene flow into Europeans

eur_pops <- c("French", "Basque", "Spanish")

yamnaya_f4 <- f4(data = f2_blocks,
                 pop1 = "pan_troglodytes",
                 pop2 = "Yamnaya",
                 pop3 = "LBK",
                 pop4 = eur_pops)
yamnaya_f4


# looking at F4 ratios
pops <- c("Han", "pan_troglodytes", "AA","Yoruba","French")

qpf4ratio(data = f2_blocks, 
          pops = pops)

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
  "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")

modern_panel_gt <- merged_gt %>% filter (population %in% lbk_modern_panel)

# remove monomorphic sites (ungrouping the tibble, as we want to use global frequencies
# rather than population frequencies)
modern_panel_gt <- modern_panel_gt %>% ungroup() %>% select_loci_if(loci_maf(genotypes)>0)

modern_panel_gt <- gt_update_backingfile(modern_panel_gt)


# reset the ploidy for this tibble, as it is now all diploids
attr(modern_panel_gt$genotypes, "ploidy") <- 2
modern_panel_gt <- modern_panel_gt %>% gt_impute_simple(method="mode")

modern_panel_gt <- modern_panel_gt %>% select_loci_if(loci_ld_clump(genotypes))

modern_pca <- modern_panel_gt %>% gt_pca_randomSVD()

autoplot(modern_pca, type = "scores") 

library(ggplot2)
autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")


lbk_gt <- merged_gt %>% filter(id == "LBK")
lbk_pca_scores <- predict(modern_pca, new_data = lbk_gt, project_method = "least_square")
lbk_pca_scores

autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")+
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black")


autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population") +
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black") +
  lims(x=c(30, 70), y = c(-10, 15))

