# Quantifying relationships between many populations 

# load libraries
library(admixtools)
library(tidypopgen)

# load f2 blocks
f2_blocks = f2_from_precomp("./data/f2_tidypopgen", verbose = FALSE)

# consider hybridisation between archaic and modern humans
neand_euras_wave <- qpwave(
  data = f2_blocks,
  left = c("French","Spanish","Tuscan", "Han", "Onge"),
  right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave
# so we fail to reject the null (for data fit)
# adding in Han and Onge we see that rank 1 best fits the data

# looking at contributions into the french
french_adm <- qpadm(data = f2_blocks,
      left = c("Loschbour", "LBK", "Yamnaya"), # populations we think have contributed to our target
      right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"), # probable outgroups
      target= "French") # population of interest

french_adm$popdrop

# looking at the Basque 
basque_adm <- qpadm(
  data = f2_blocks,
  left = c("Loschbour", "LBK", "Yamnaya"), # populations we think have contributed to our target
  right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"), # probable outgroups
  target= "Basque" # population of interest
  ) 

basque_adm$popdrop
# hence yamnaya component is larger in the French (0.729>0.648)

# Making graph -------------------------------------------------------------
# define graph
base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

# define edges
base_edges %>% edges_to_igraph() %>% plot_graph()

# redefine tree
base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges

# convert to igraph object 
base_igraph <- base_edges %>% edges_to_igraph()
# check the tree is valid
is_valid(base_igraph)

base_igraph %>% plot_graph()
# to make interactive:
base_igraph %>% plotly_graph()

# add in f2_block data
base_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)
base_qpgraph$f3

# check that graph is compatible with the data 
base_qpgraph$f3 %>% filter(abs(z)>2)

base_qpgraph$edges %>% plot_graph()

# make fake graph 
base_swapped_qpgraph <- matrix(
  c("R",	"Dinka",
    "R", "eAfr",
    "eAfr",	"Mbuti",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,byrow = TRUE,
  dimnames=list(NULL, c("from","to"))
  )

# convert to igraph object 
base_swapped_qpgraph <- base_swapped_qpgraph %>% edges_to_igraph()
# check the tree is valid
is_valid(base_swapped_qpgraph)

base_swapped_qpgraph %>% plot_graph()
# to make interactive:
base_swapped_qpgraph %>% plotly_graph()

# add in f2_block data
base_swapped_qpgraph <- qpgraph(data = f2_blocks, graph = base_swapped_qpgraph)

fits = qpgraph_resample_multi(
  f2_blocks, 
  graphlist = list(base_qpgraph[[1]], base_swapped_qpgraph[[1]]),
  nboot = 100
  )
# look at the different fits 
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)

base_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

## Looking at Yamnaya edges ------------------------------
yamnaya_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "Dinka",
    "eAfr", "outAfrica",
    "outAfrica", "Han",
    "outAfrica", "wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "Loschbour"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
yamnaya_igraph <- yamnaya_edges %>% edges_to_igraph()
yamnaya_igraph %>% plot_graph()

# adding HGs
lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "outAfrica",
    "eAfr", "Dinka",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "Loschbour"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to"))
)
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

# add in f2_block data
lbk_extra_igraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)
lbk_extra_igraph$f3

# check that graph is compatible with the data 
lbk_extra_igraph$f3 %>% filter(abs(z)>2)

# looking at admixture
admixture_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "BasalEurasian",
    "eAfr", "Dinka",
    "BasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "wEurasian", "WHG",
    "pLoschbour", "pLBK",
    "BasalEurasian","pLBK",
    "pLBK", "LBK"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to"))
)
admixture_edges <- admixture_edges %>% edges_to_igraph()
admixture_edges %>% plot_graph()

is_valid(admixture_edges)



# add in more edges
lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to"))
  )
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)

lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)
lbk_extra_qpgraph$edges %>% plot_graph()

