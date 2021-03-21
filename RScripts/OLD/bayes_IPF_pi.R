bayes_IPF_pi <- function(prior, contingency_table, runs){
  IPF_pi_chain <- matrix(nrow = nrow(contingency_table), ncol = runs)
  IPF_pi_chain[, 1] <- rep(1/nrow(contingency_table), nrow(contingency_table))
  
  contingency_table %<>% mutate("Freq" = Freq + prior)
  
  var1_levels <- unique(contingency_table$Var1)
  
  for(run in 2:runs){
    for(var in var1_levels){
      index <- which(contingency_table$Var1 == var)
      alpha_post <- contingency_table %>%  filter(Var1 == var) %>% 
        summarise_at("Freq", sum) %>% as.numeric()
      
      gamma_draws <- rgamma(n = length(index), shape = alpha_post)
      
      IPF_pi_chain[index, run] <- 
        IPF_pi_chain[index, (run - 1)]*
        ((gamma_draws/sum(gamma_draws))/sum(IPF_pi_chain[index, (run - 1)]))
    }
  }
  
  return(IPF_pi_chain)
}

#testing
prior <- alpha_0
contingency_table <- cross_class_label

chain <- bayes_IPF_pi(prior, contingency_table, runs = 100)

#look at plot
bayes_IPF_pi_chain_data <- chain %>% as.data.frame() %>% 
  set_colnames(seq(1, ncol(chain))) %>%
  mutate("Cell" = contingency_table$`Cell Label`) %>% 
  pivot_longer(-c("Cell"), names_to = "Run", values_to = "probability") %>%
  mutate_at("Run", as.numeric) %>%
  mutate_if(is.character, as.factor) 

extended_pallette6 <- colorRampPalette(wes_palette("Darjeeling1"))(6)

bayes_IPF_pi_chain_plot <- ggplot(data = bayes_IPF_pi_chain_data, 
                        aes(x = Run, y = probability, colour = Cell)) +       
  geom_line(aes(group = Cell)) + 
  theme_minimal() + xlab("Run") + ylab("Probability of cell membership") +  
  scale_color_manual(values = rev(extended_pallette6)) + 
  scale_x_continuous(breaks = seq(1, ncol(chain))) 

