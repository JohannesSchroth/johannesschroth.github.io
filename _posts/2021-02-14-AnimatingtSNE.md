---
layout: post
title: Animating tSNE iterations
subheading: A stepwise implementation of tSNE
author: Johannes
categories: Data-Science
banner: assets/images/banners/rnaseqinr_banner.png
---

### Environment Setup

```
library(dslabs)
library(rsvd)
library(dplyr)
library(tidyverse)
library(gganimate)

source('~/development/FIt-SNE-master/fast_tsne.R', chdir=T)

set.seed(1234)

```

### Data Import and Preprocessing

```
mnist <- read_mnist()
```
### Stepwise tSNE Iteration Function
```

tsne_iterate <- function(data, iteration_sequence) {
  
  start <- Sys.time()
  
  tsne_df_tmp <- data.frame()
  run <- 0
  
  for (i in iteration_sequence) {
    
    run <- run + 1
    
    message(paste('Run', run, 'of', length(iteration_sequence)))
    
    tsne_tmp <- data %>%
      fftRtsne(perplexity = nrow(.)/100,
               theta = 0.5,
               max_iter = i,
               learning_rate = nrow(.)/12,
               rand_seed = 1234) %>%
      as.data.frame %>%
      mutate('iter' = i)
    
    tsne_df_tmp <- rbind(tsne_df_tmp, tsne_tmp)
  }
  
  message(paste('Time to complete:', round(Sys.time() - start, 2), 'min'))
  
  return(tsne_df_tmp)
}
```

### Running the function
```
down_sample <- sample(nrow(mnist$train$images), 10000)
tsne_df <- tsne_iterate(mnist$train$images[down_sample,], seq(50,2000,50))

colnames(tsne_df) <- c('tSNE1', 'tSNE2', 'iteration')
```

### Animating with gganimate

```
tsne_df %>%
  mutate('labels' = rep(as.factor(mnist$train$labels[down_sample]), length(unique(iteration)))) %>%
  ggplot(aes(tSNE1, tSNE2, colour = labels)) +
  geom_point() +
  theme_void() +
  transition_states(iteration,
                    transition_length = 1,
                    state_length = 1) +
  labs(title = 'Max Iteration: {closest_state}')
```

<img src="/assets/gif/mnist.gif" width="500" height="500" style="display: block; margin-left: auto; margin-right: auto;"/>




