Y <-read.csv("final return 2016.csv")
Y<-as.matrix(Y[1:818,2:72])
library(BigVAR)
library(frequencyConnectedness)
big_var_est<- function(Y) {Model1 = constructModel(as.matrix(Y),1,struct="Basic",gran=c(50,10),
                                                   h = 20,cv = "Rolling",verbose = FALSE,IC = TRUE,VARX = list(), 
                                                   T1 = floor(nrow(Y)/3),T2 = floor(2 * nrow(Y)/3),ONESE = FALSE,
                                                   ownlambdas = FALSE,recursive = FALSE,dates = as.character(NULL), 
                                                   window.size = 500,separate_lambdas = FALSE)
MYResults= cv.BigVAR(Model1)}
oo <- big_var_est(Y)
spilloverDY12(oo, n.ahead = 100, no.corr = F)
staticnet=spilloverDY12(oo, n.ahead = 100, no.corr = F)
W1<-as.matrix(staticnet$tables[[1]])
write.csv(W1,"subsamplematrix1.csv")
############################################
Y <-read.csv("final return 2016.csv")
Y<-as.matrix(Y[1177:1749,2:72])
library(BigVAR)
library(frequencyConnectedness)
big_var_est<- function(Y) {Model1 = constructModel(as.matrix(Y),1,struct="Basic",gran=c(50,10),
                                                   h = 20,cv = "Rolling",verbose = FALSE,IC = TRUE,VARX = list(), 
                                                   T1 = floor(nrow(Y)/3),T2 = floor(2 * nrow(Y)/3),ONESE = FALSE,
                                                   ownlambdas = FALSE,recursive = FALSE,dates = as.character(NULL), 
                                                   window.size = 500,separate_lambdas = FALSE)
MYResults= cv.BigVAR(Model1)}
oo <- big_var_est(Y)
spilloverDY12(oo, n.ahead = 100, no.corr = F)
staticnet=spilloverDY12(oo, n.ahead = 100, no.corr = F)
W2<-as.matrix(staticnet$tables[[1]])
write.csv(W2,"subsamplematrix2.csv")
##############################################
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

mat <- as.matrix(W1)

clean_names <- colnames(mat)
clean_names <- gsub("\\.\\.", ", ", clean_names)
clean_names <- gsub("\\.", " ", clean_names)     
colnames(mat) <- clean_names
rownames(mat) <- clean_names

commodities <- c(
  "Corn, China","Wheat, China","Japonica Rice, China","Pork, China","Beef, China",
  "Mutton, China","Cotton, China","Soybean, China","Soybean Meal, China",
  "Soybean Oil, China","Rapeseed Meal, China","Rapeseed Oil, China","Palm Oil, China",
  "Pure Benzene, China","Toluene, China","O Xylene, China","Styrene, China",
  "Butadiene, China","P Xylene, China","PTA, China","Epoxy Ethane, China","Propylene, China",
  "Polystyrene, China","Natural Rubber, China","PP, China","Causticsoda, China",
  "Potassic Fertilizers, China","Xylene, China","Glycol, China","Acrylamide, China",
  "Acrylic Acid, China","Acrylonitrile, China","Butyl Acrylate, China","PET Chips, China",
  "Methyl Methacrylate, China","PVC, China","Urea, China","Methanol, China",
  "Metallurgical Coke, China","Hard Coking Coal, China","Coking Coal, China",
  "Hot rolled Sheet Coil, China","Platinum, China","Silver, China","Palladium, China",
  "Gold, China","Copper, China","Nickel, China","Aluminum, China","Zinc, China",
  "Lead, China","Tin, China","Chrome Metal, China","Silicon Metal, China","Cobalt, China",
  "Ingot, China","Electrolytic Manganese, China","Nickel Sulfate, China","Antimony, China",
  "Soybean Oil, US","Diesel Oil, Singapore","Gasoline, Japan","Brent","WTI",
  "Henry Hub Natural Gas","Crude Oil, UAE Dubai","Naphtha, Korea",
  "Gold, UK","Silver, UK","Platinum, UK","Palladium, UK"
)

get_category <-  function(name) {
  n <- tolower(name) 
  if (str_detect(n, "brent|wti|dubai|crude|natural gas|diesel|gasoline|naphtha|coke|coal")) {
    return("Energy")
  }
  else if (str_detect(n, "gold|silver|platinum|palladium|copper|nickel|aluminum|zinc|lead|tin|chrome|silicon|cobalt|ingot|manganese|antimony|sheet|coil")) {
    return("Metals")
  }
  else if (str_detect(n, "corn|wheat|rice|pork|beef|mutton|cotton|soybean|rapeseed|palm|rubber|meal")) {
    return("Agriculture")
  }
  else {
    return("Chemicals") 
  }
}

node_data <- tibble(name = commodities) %>%
  mutate(Category = map_chr(name, get_category))

g <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE, diag = FALSE)

g_tidy <- as_tbl_graph(g) %>%
  left_join(node_data, by = "name") %>%
  mutate(
    Influence = centrality_degree(mode = "out", weights = weight)
  )


threshold <- quantile(E(g)$weight, 0.95)

g_pruned <- g_tidy %>%
  activate(edges) %>%
  filter(weight > threshold) %>%
  activate(nodes) %>%
  filter(!node_is_isolated()) 

set.seed(42)

my_colors <- c(
  "Agriculture" = "#2ca02c", 
  "Energy"      = "#d62728", 
  "Metals"      = "#ff7f0e", 
  "Chemicals"   = "#1f77b4"  
)

plot <- ggraph(g_pruned, layout = "fr") + 
  
  
  geom_edge_fan(aes(width = weight, alpha = weight), 
                arrow = arrow(length = unit(2, 'mm')), 
                end_cap = circle(3, 'mm'), 
                color = "grey20", 
                show.legend = FALSE) +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_edge_alpha(range = c(0.2, 0.6)) +
  
  
  geom_node_point(aes(size = Influence, color = Category), alpha = 0.9) +
  
  
  scale_color_manual(values = my_colors) +
  scale_size_continuous(range = c(3, 10)) + 
  
  
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 15) +
  
  
  labs(
    color = "Commodity Sector",
    size = "Market Influence"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey50")
  )

print(plot)
#########################################
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

mat <- as.matrix(W2)

clean_names <- colnames(mat)
clean_names <- gsub("\\.\\.", ", ", clean_names)
clean_names <- gsub("\\.", " ", clean_names)     
colnames(mat) <- clean_names
rownames(mat) <- clean_names

commodities <- c(
  "Corn, China","Wheat, China","Japonica Rice, China","Pork, China","Beef, China",
  "Mutton, China","Cotton, China","Soybean, China","Soybean Meal, China",
  "Soybean Oil, China","Rapeseed Meal, China","Rapeseed Oil, China","Palm Oil, China",
  "Pure Benzene, China","Toluene, China","O Xylene, China","Styrene, China",
  "Butadiene, China","P Xylene, China","PTA, China","Epoxy Ethane, China","Propylene, China",
  "Polystyrene, China","Natural Rubber, China","PP, China","Causticsoda, China",
  "Potassic Fertilizers, China","Xylene, China","Glycol, China","Acrylamide, China",
  "Acrylic Acid, China","Acrylonitrile, China","Butyl Acrylate, China","PET Chips, China",
  "Methyl Methacrylate, China","PVC, China","Urea, China","Methanol, China",
  "Metallurgical Coke, China","Hard Coking Coal, China","Coking Coal, China",
  "Hot rolled Sheet Coil, China","Platinum, China","Silver, China","Palladium, China",
  "Gold, China","Copper, China","Nickel, China","Aluminum, China","Zinc, China",
  "Lead, China","Tin, China","Chrome Metal, China","Silicon Metal, China","Cobalt, China",
  "Ingot, China","Electrolytic Manganese, China","Nickel Sulfate, China","Antimony, China",
  "Soybean Oil, US","Diesel Oil, Singapore","Gasoline, Japan","Brent","WTI",
  "Henry Hub Natural Gas","Crude Oil, UAE Dubai","Naphtha, Korea",
  "Gold, UK","Silver, UK","Platinum, UK","Palladium, UK"
)

get_category <-  function(name) {
  n <- tolower(name) 
  if (str_detect(n, "brent|wti|dubai|crude|natural gas|diesel|gasoline|naphtha|coke|coal")) {
    return("Energy")
  }
  else if (str_detect(n, "gold|silver|platinum|palladium|copper|nickel|aluminum|zinc|lead|tin|chrome|silicon|cobalt|ingot|manganese|antimony|sheet|coil")) {
    return("Metals")
  }
  else if (str_detect(n, "corn|wheat|rice|pork|beef|mutton|cotton|soybean|rapeseed|palm|rubber|meal")) {
    return("Agriculture")
  }
  else {
    return("Chemicals") 
  }
}

node_data <- tibble(name = commodities) %>%
  mutate(Category = map_chr(name, get_category))

g <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE, diag = FALSE)

g_tidy <- as_tbl_graph(g) %>%
  left_join(node_data, by = "name") %>%
  mutate(
    Influence = centrality_degree(mode = "out", weights = weight)
  )


threshold <- quantile(E(g)$weight, 0.95)

g_pruned <- g_tidy %>%
  activate(edges) %>%
  filter(weight > threshold) %>%
  activate(nodes) %>%
  filter(!node_is_isolated()) 

set.seed(42)

my_colors <- c(
  "Agriculture" = "#2ca02c", 
  "Energy"      = "#d62728", 
  "Metals"      = "#ff7f0e", 
  "Chemicals"   = "#1f77b4"  
)

plot <- ggraph(g_pruned, layout = "fr") + 
  
  
  geom_edge_fan(aes(width = weight, alpha = weight), 
                arrow = arrow(length = unit(2, 'mm')), 
                end_cap = circle(3, 'mm'), 
                color = "grey20", 
                show.legend = FALSE) +
  scale_edge_width(range = c(0.2, 1.2)) +
  scale_edge_alpha(range = c(0.2, 0.6)) +
  
  
  geom_node_point(aes(size = Influence, color = Category), alpha = 0.9) +
  
  
  scale_color_manual(values = my_colors) +
  scale_size_continuous(range = c(3, 10)) + 
  
  
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 15) +
  
  
  labs(
    color = "Commodity Sector",
    size = "Market Influence"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey50")
  )

print(plot)