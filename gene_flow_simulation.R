library(tidyverse)

# gene flow simulation

popnames <- c("A","B","C","D","E")

pfreq <- c(1,0.35, 0.76, 0.15, 0.45)

gene_f_df <- data.frame(cbind(popnames,pfreq))
gene_f_df$pfreq <- as.numeric(gene_f_df$pfreq)



# create a time vector
time <- c(0:100)

# create a function for calculating allele frequency
p_freq <- function(pAVG,pIN,m,t){
  pt = pAVG + ((pIN - pAVG)*(1-m)^t)
  return(pt)
}

#empty vector to store results
population <- NULL
years <- NULL
freq <- NULL
flow_rate <- NULL
# calculate p AVG
AVG <- mean(gene_f_df$pfreq)

# make simulation
m_list <- seq(0,0.1, by =0.01)
for (m in m_list) {
  for (i in time)
  {
    for (row in 1:nrow(gene_f_df)) 
    {
      pop = gene_f_df[row, "popnames"]
      p_initial = gene_f_df[row,"pfreq"]
      population = append(population, pop)
      years = append(years, i)
      flow_rate = append(flow_rate,m)
      freq = append(freq, p_freq(AVG,p_initial,m,i))
    }
  }
}

# save result into data frame
results <- data.frame(population,flow_rate,years,freq)

# plot results
results %>% 
  filter(flow_rate != 0) %>% 
  ggplot(aes(x = years, y =freq, colour = population)) +
  geom_line() + facet_wrap(flow_rate~.) + theme_light() +
  labs(y = "Allele Frequency")
