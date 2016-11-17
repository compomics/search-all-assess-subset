#' Calculate density from pi0 distribution
#'
#' @param pi0
#' @param n_targets
#' @param n_decoys
#' @return
#' @export
#'
#' @examples
dpi0= function(pi0, n_targets, n_decoys) {
  dbeta(pi0 / (1 + pi0), n_decoys + 1, n_targets + 1) *
    1 / (1 + pi0) ^ 2 /
    pbeta(.5, n_decoys + 1, n_targets + 1)
}

##' sample from pi0 distribution
##'
##' .. content for \details{} ..
##' @param n
##' @param n_targets
##' @param n_decoys
##' @return
##' @export
rpi0 = function(n, n_targets, n_decoys) {
  x = as.numeric()
  ## sample until n samples has been reached
  while(n > 0){
    xt = rbeta(n,n_decoys + 1,n_targets + 1)
    ## remove samples outside 0-.5 boundary
    xt = xt[xt<.5]
    xt = xt/(1-xt)
    x = c(x,xt)
    n = n - length(xt)
  }
  return(x)
}

#' Title
#'
#' @param n
#' @param pi0
#' @param sims
#'
#' @return
#' @export
#' @import tibble
#' @examples
simulate_subset = function(n,pi0,sims = 1){
  pi0D = 2*pi0/(1+pi0)
  min_n = rbinom(sims ,n , pi0D)
  data_frame(
    n,
    pi0,
    decoy_n =rbinom(sims, min_n, .5),
    target_n = n - decoy_n,
    H0_n = min_n - decoy_n,
    H1_n = target_n - H0_n)
}

#' Calculate qvalues of the subset PSMs.
#'
#' @param df A dataframe with at least the columns scores from decoy and subset target PSMs (see details)
#' @param score_higher TRUE if higher scores mean a better peptide sequence match.
#'
#' @return A data frame containing the collumns from the original input,
#' A column with the conservative pi_0 estimation,
#' the estimated qvalues from the stable FDR estimation strategy using the conservative pi_0 estimator and a large decoy set.
#' Lastly, a column containing qvalues estimated with the Benjamini Hochbergh procedure.
#' @export
#' @import dplyr
#'
#' @examples
calculate_fdr = function(df,score_higher = TRUE) {
  mutate(
    df,
    index = row_number(),
    ## pi_0 of subset PSMs
    # pi_0 = (sum(decoy[subset == 1] == 1)) /
    #   (sum(decoy[subset == 1] != 1)),
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

#' Creates density plot of the pi0 distribution
#'
#' @param n_targets
#' @param n_decoys
#'
#' @return
#' @export
#' @import ggplot2
#' @examples
pi0plot = function(n_targets, n_decoys) {
  grid = seq(0, 1, .001)
  dens = dpi0(grid, n_targets, n_decoys)
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

#' Creates PP plot of two empirical distributions
#'
#' @param score
#' @param label
#' @param pi0
#' @param title
#' @param xlab
#' @param ylab
#'
#' @return
#' @export
#' @import ggplot2
#' @import tibble
#' @examples
PPplot = function(score, label, pi0 = 0,title = 'PP plot of target PSMs' ,
                  xlab = 'Decoy percentile' ,ylab = 'Target\npercentile'){
  p = ggplot()  +
    geom_abline(slope = pi0,color = 'black') +
    labs(x = xlab, y = ylab ,title = title) +
    xlim(0,1) + ylim(0,1) +
    theme_bw() +
    theme(
      plot.title = element_text(size = rel(1.5)),
      axis.title = element_text(size = rel(1.2)),
      axis.text = element_text(size = rel(1.2)),
      axis.title.y = element_text(angle = 0))

  if ((length(score[!label]) == 0) | (length(score[label]) == 0))
    return(p + annotate('text',label = 'NOT ENOUGH DATA TO PLOT',x = .5,y = .5))

  Ft = ecdf(score[!label])
  Fd = ecdf(score[label])
  x = score[!label]
  df = data_frame(Fdp = Fd(x), Ftp = Ft(x))

  p + geom_point(data = df,aes(Fdp,Ftp),color = 'dark grey')
}


#' Title
#'
#' @param df
#'
#' @return
#' @export
#' @import ggplot2
#' @import dplyr
#' @import cowplot
#' @examples
plot_diag = function(df){
  df = mutate(df,subset = subset == 1,decoy = decoy ==1) %>%
    ## remove non subset targets
    filter(decoy|subset)

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


#' Title
#'
#' @param H1_n
#' @param H0_n
#' @param decoy_n
#' @param decoy_large_n
#' @param H0_mean
#' @param H1_mean
#' @param H0_sd
#' @param H1_sd
#' @param decoy_mean
#' @param decoy_sd
#' @param decoy_large_mean
#' @param decoy_large_sd
#'
#' @return
#' @export
#' @import dplyr
#'
#' @examples
sample_dataset = function(H1_n = 160,H0_n = 40, decoy_n = H0_n ,decoy_large_n = 2000,
                          H0_mean=2.75, H1_mean=3.31,H0_sd=.13,H1_sd=.28,
                          decoy_mean = H0_mean, decoy_sd = H0_sd,
                          decoy_large_mean = H0_mean, decoy_large_sd = H0_sd){

  d1 = d2 = d3 = d4 = NULL
  if (decoy_large_n){
   d1 = data_frame(score = rnorm(decoy_large_n,decoy_large_mean,decoy_large_sd),
               decoy = TRUE,
               H0 = FALSE,
               subset = FALSE)
  }
  if (decoy_n){
  d2 = data_frame(score = rnorm(decoy_n,decoy_mean,decoy_sd),
               decoy = TRUE,
               H0 = FALSE,
               subset = TRUE)
  }
  if (H0_n){
  d3 = data_frame(score = rnorm(H0_n,H0_mean,H0_sd),
               decoy = FALSE,
               H0 = TRUE,
               subset = TRUE)
  }
  if (H1_n){
  d4 =   data_frame(score = rnorm(H1_n,H1_mean,H1_sd),
               decoy = FALSE,
               H0 = FALSE,
               subset = TRUE)
  bind_rows(d1,d2,d3,d4)
  }
  }

#' plots the theoretical distribution of all components in the PSM distribution
#' @param H0_mean
#' @param H1_mean
#' @param H0_sd
#' @param H1_sd
#' @param decoy_mean
#' @param decoy_sd
#' @param decoy_large_mean
#' @param decoy_large_sd
#' @return
#' @export
#' @import dplyr
#' @import ggplot2
#' @examples
plot_theo_dist2 = function(H0_mean=2.75, H1_mean=3.31,H0_sd=.13,H1_sd=.28,
                          decoy_mean = H0_mean, decoy_sd = H0_sd,
                          decoy_extra_mean = H0_mean, decoy_extra_sd = H0_sd){

    ## make grid of scores
  d = data_frame(score = seq(1,5,.01))
  ## calculate theoretical density for eacht dist

 d = bind_rows(
    mutate(d,dens = dnorm(score,H1_mean,H1_sd),
           PSMs = 'correct subset')
    ,
    mutate(d,dens = dnorm(score,H0_mean,H0_sd),
           PSMs = 'incorrect subset')
    ,
    mutate(d,dens = dnorm(score,decoy_mean,decoy_sd),
           PSMs = 'subset decoy')
    ,
    mutate(d,dens = dnorm(score,decoy_extra_mean,decoy_extra_sd),
           PSMs = 'extra decoy')
  )

 d = mutate(d,PSMs = factor(PSMs,levels = unique(d$PSMs)))
  p1 = ggplot(d,aes(score,dens,col = PSMs,size = PSMs,linetype = PSMs)) + geom_line() +
    labs(x = 'Score',y = 'Density') +
    theme(axis.title = element_text(size = rel(1.2)), axis.title.y = element_text(angle = 0),
          legend.position = 'top') +
    scale_size_manual(values=c(4,4,1,1))+
    scale_linetype_manual(values=c(1,1,1,2)) +
    scale_color_manual(values = c('green','red','orange','blue'))
  print(p1)
  list(data = d,plot = p1)
  }


#' Title
#'
#' @param H1_n
#' @param H0_n
#' @param decoy_n
#' @param decoy_large_n
#' @param H0_mean
#' @param H1_mean
#' @param H0_sd
#' @param H1_sd
#' @param decoy_mean
#' @param decoy_sd
#' @param decoy_large_mean
#' @param decoy_large_sd
#'
#' @return
#' @export
#' @import dplyr
#' @import ggplot2
#' @examples
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
               dist = 'extra_decoys')
    ,
    mutate(d,dens = p_decoy * dnorm(score,decoy_mean,decoy_sd),
           dist = 'subset_decoys')
    ,
    mutate(d,dens = p_H0 * dnorm(score,H0_mean,H0_sd),
           dist = 'H0')
    ,
    mutate(d,dens = p_H1 * dnorm(score,H1_mean,H1_sd),
           dist = 'H1')
    ,
    mutate(d,dens = p_decoy * dnorm(score,decoy_mean,decoy_sd) +
             p_decoy_large * dnorm(score,decoy_large_mean,decoy_large_sd),
           dist = 'all_decoys')
    ,
    mutate(d,dens = p_H0 * dnorm(score,H0_mean,H0_sd) +
             p_H1 * dnorm(score,H1_mean,H1_sd),
           dist = 'targets')
  )

 d = mutate(d,dist = factor(dist,levels = unique(d$dist)))
  p1 = ggplot(d,aes(score,dens,col = dist,size = dist,linetype = dist)) + geom_line() +
    labs(x = 'Score',y = 'Density', colour = '',linetype = '',size = '') +
    theme(axis.title = element_text(size = rel(1.2)), axis.title.y = element_text(angle = 0),
          legend.position = 'top') +
    scale_size_manual(values=c(1,1,1,1,2,2))+
    scale_linetype_manual(values=c(2,3,2,3,1,1)) +
    scale_color_manual(values = c('red','red','blue','blue','red','blue'))
  print(p1)
  list(data = d,plot = p1)
  }
