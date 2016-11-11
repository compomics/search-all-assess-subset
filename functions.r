library(tidyverse)



pi0_likelihood = function(n_targets, n_decoys, pi0) {
  dbeta(pi0 / (1 + pi0), n_decoys + 1, n_targets + 1) *
    1 / (1 + pi0) ^ 2 /
    pbeta(.5, n_decoys + 1, n_targets + 1)
}

#' Title
#'
#' @param df 
#' @param score_higher 
#'
#' @return
#' @export
#'
#' @examples
calculate_fdr = function(df,score_higher = TRUE) {
  mutate(
    df,
    index = row_number(),
    ## pi_0 of subset PSMs
    pi_0 = (sum(decoy[subset == 1] == 1)) /
      (sum(decoy[subset == 1] != 1)),
    ## conservative pi_0 of subset PSMs with upperbound of 1
    pi_0_cons = (sum(decoy[subset == 1] == 1) + 1) /
      (sum(decoy[subset == 1] != 1) + 1),
    pi_0_cons = ifelse(pi_0_cons > 1, 1, pi_0_cons)
  ) %>%
    # Sort the score depending on if higher scores are better or not
    arrange(desc(if (score_higher) score else -score)) %>%
    # Calculate classical FDR on subset
    mutate(FDR = ifelse((subset == 1) & (decoy == 0),
                        cummax(cumsum(decoy == 1) / 
                                 cumsum(decoy != 1)), NA)) %>%
    # calculate BH FDR on subset
    mutate(FDR_BH = ifelse((subset == 1) & (decoy == 0),
                           cummax((cumsum(decoy == 1) / max(decoy == 1)) / 
                                    (cumsum(decoy != 1) / decoy != 1)),
                           NA)) %>%
    # calculate stable FDR on subset
    mutate(FDR_stable = FDR_BH * pi_0_cons) %>%
    ## Does not allow any FDR to be above 1
    mutate_at(vars(FDR, FDR_BH, FDR_stable),
              funs(ifelse((. > 1) & !is.na(.), 1, .))) %>%
    # Put dataframe back in original order
    ungroup %>% arrange(index) %>%
    select(-index)
}

#' Title
#'
#' @param n_targets
#' @param n_decoys
#'
#' @return
#' @export
#'
#' @examples
pi0plot = function(n_targets, n_decoys) {
  grid = seq(0, 1, .001)
  dens = pi0_likelihood(n_targets, n_decoys, grid)
  df = data_frame(grid, dens)
  ggplot(df, aes(x = grid, y = dens)) + geom_line(col = 'dark grey') +
    geom_vline(xintercept = (n_decoys + 1) / (n_targets + 1),
               col = 'black') +
    labs(
      x = expression(pi[0]),
      y = 'Density',
      title = expression('Posterior probability of' ~ pi[0])
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.2)),
      axis.text = element_text(size = rel(1.2)),
      axis.title.y = element_text(angle = 0)
    )
}

PPplot = function(score, label, pi0 = 0,title = 'PP plot of target PSMs' ,
                  xlab = 'Decoy percentile' ,ylab = 'Target\npercentile'){
  Ft = ecdf(score[!label])
  Fd = ecdf(score[label])
  x = score[!label]
  df = data_frame(Fdp = Fd(x), Ftp = Ft(x))
  ggplot(df,aes(Fdp,Ftp)) + geom_point(,color = 'dark grey') +
    geom_abline(slope = pi0,color = 'black') +
    labs(x = xlab, y = ylab ,title = title) +
    xlim(0,1) + ylim(0,1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.2)),
      axis.text = element_text(size = rel(1.2)),
      axis.title.y = element_text(angle = 0))
}

plot_diag = function(df){
  df = mutate(df,subset = subset == 1,decoy = decoy ==1)
  dfsub = filter(df,subset)
  dfdec = filter(df,decoy)
  
  p1 =  pi0plot(sum(!dfsub$decoy),sum(dfsub$decoy))
  
  p2 = PPplot(df$score,df$decoy,sum(dfsub$decoy)/sum(!dfsub$decoy),
              xlab = 'All decoy percentile')
  
  d = filter(df,decoy)
  d2 = filter(d,subset ) %>%
    mutate(subset = FALSE) ## add subset decoys to rest because it's also part of all decoys'
  d = bind_rows(d,d2)
  p3 = PPplot(d$score,!d$subset,1,title = 'PP plot of subset decoy PSMs',
              xlab = 'All decoy percentile',ylab = 'Subset decoy\npercentile')
  
  d = filter(df,subset)
  p4 = PPplot(d$score,d$decoy, sum(dfsub$decoy)/sum(!dfsub$decoy),
              xlab = 'Subset decoy percentile')
  
  p_all = cowplot::plot_grid(p1,p2,p3,p4, align = 'v',labels = 'auto',hjust = -4, label_size = 16)
  #print(p_all)
  return(list(pi0plot = p1, decoyall_targetsubset = p2,
              decoyall_decoy_subset = p3, decoysubset_target_subset = p4, all = p_all))
}


# pi_0D = 2*par$pi_0/(1+par$pi_0)
# par$n_decoy = rbinom(1,par$n_all_min,.5)
# par$ns_H0 = rbinom(1,par$n,par$pi_0D)
# par$ns_decoy =rbinom(1,par$ns_H0,.5)
# par$ns_target = par$n-par$ns_decoy


sample_dataset = function(H1_n = 160,H0_n = 40, decoy_n = H0_n ,decoy_large_n = 2000,
                          H0_mean=2.75, H1_mean=3.31,H0_sd=.13,H1_sd=.28,
                          decoy_mean = H0_mean, decoy_sd = H0_sd,
                          decoy_large_mean = H0_mean, decoy_large_sd = H0_sd){
  bind_rows(
    data_frame(score = rnorm(decoy_large_n,decoy_large_mean,decoy_large_sd),
               decoy = TRUE,
               H0 = FALSE,
               subset = FALSE)
    ,
    data_frame(score = rnorm(decoy_n,decoy_mean,decoy_sd),
               decoy = TRUE,
               H0 = FALSE,
               subset = TRUE)
    ,
    data_frame(score = rnorm(H0_n,H0_mean,H0_sd),
               decoy = FALSE,
               H0 = TRUE,
               subset = TRUE)
    ,
    data_frame(score = rnorm(H1_n,H1_mean,H1_sd),
               decoy = FALSE,
               H0 = FALSE,
               subset = TRUE)
  )
}


plot_theo_dist = function(H1_n = 160,H0_n = 40, decoy_n = H0_n ,decoy_large_n = 2000,
                                 H0_mean=2.75, H1_mean=3.31,H0_sd=.13,H1_sd=.28,
                          decoy_mean = H0_mean, decoy_sd = H0_sd,
                          decoy_large_mean = H0_mean, decoy_large_sd = H0_sd){

    ## make grid of scores
  d = data_frame(score = seq(1,5,.01))
  ## calculate theoretical density for eacht dist

  p_H0 = H0_n/(H0_n + H1_n)
  p_H1 = H1_n/(H0_n + H1_n)
  p_decoy_large = p_H0 * decoy_large_n/(decoy_n + decoy_large_n)
  p_decoy = p_H0 * decoy_n/(decoy_n + decoy_large_n)

 d = bind_rows(
    mutate(d,dens = p_decoy_large * dnorm(score,decoy_large_mean,decoy_large_sd),
               dist = 'decoy_large')
    ,
    mutate(d,dens = p_decoy * dnorm(score,decoy_mean,decoy_sd),
           dist = 'decoy')
    ,
    mutate(d,dens = p_H0 * dnorm(score,H0_mean,H0_sd),
           dist = 'H0')
    ,
    mutate(d,dens = p_H1 * dnorm(score,H1_mean,H1_sd),
           dist = 'H1')
    ,
    mutate(d,dens = p_decoy * dnorm(score,decoy_mean,decoy_sd) + 
             p_decoy_large * dnorm(score,decoy_large_mean,decoy_large_sd),
           dist = 'decoy_all')
    ,
    mutate(d,dens = p_H0 * dnorm(score,H0_mean,H0_sd) + 
             p_H1 * dnorm(score,H1_mean,H1_sd),
           dist = 'target')
  )
  
  p1 = ggplot(d,aes(score,dens,col = dist)) + geom_line() +
    labs(x = 'Score',y = 'Density',colour= '') +
    theme(axis.title = element_text(size = rel(1.2)), axis.title.y = element_text(angle = 0),
          legend.position = 'top')
  print(p1)
  list(data = d,plot = p1)
  }


# d = read.csv('/home/adriaan/Dropbox/paper_subset/revision/Rpackage/cytoplasm.csv',sep = ',')
# write.csv(d,'/home/adriaan/Dropbox/paper_subset/revision/Rpackage/cytoplasm_back.csv',row.names = FALSE)
# nrow(d)
# d2 = filter(d, !(!decoy & !subset))
# write.csv(d2,'/home/adriaan/Dropbox/paper_subset/revision/Rpackage/cytoplasm.csv',row.names = FALSE)
# d = read_delim('/home/adriaan/Dropbox/paper_subset/revision/Rpackage/cytoplasm.csv',delim = ',')
# head(d)
# 

