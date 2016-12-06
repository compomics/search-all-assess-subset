#' Density function for the pi0 distribution
#'
#' @param pi0 vector of quantiles.
#' @param n_targets vector of observed target PSMs.
#' @param n_decoys vector of observed decoy PSMs.
#' @return vector of densities.
#' The length is the maximum length of the numerical arguments.
#' Returns 'NaN' for 'pi0 < 0' and 'pi > 1'.
#' @export
#'
#' @examples
dpi0= function(pi0, n_targets, n_decoys) {
  pi0 = ifelse((pi0 > 1) | (pi0 < 0),NaN,pi0)
  dbeta(pi0 / (1 + pi0), n_decoys + 1, n_targets + 1) *
    1 / (1 + pi0) ^ 2 /
    pbeta(.5, n_decoys + 1, n_targets + 1)
}

##' Random generation for the pi0 distribution
##'
##' @param n number of observations.
##' @param n_targets number of observed target PSMs.
##' @param n_decoys number of observed decoy PSMs.
##' @return vector of random deviates.
##' The length equals `n'.
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

#' Random generation of a dataset after TDA.
#'
#' Random generation of number of decoy, correct target and incorrect target PSMs target PSMs after a competitive target-decoy search.
#'
#' @param n number of total PSMs.
#' @param pi0 theoretical \eqn{\pi_0}.
#' @param sims number of observations.
#'
#' @return A data frame with ``sims'' rows and 6 rows:
#' \describe{
#'   \item{n}{number of PSMs.}
#'   \item{pi0}{theoretical \eqn{\pi_0}.}
#'   \item{decoy_n}{number of decoy PSMs.}
#'   \item{target_n}{number of target PSMs.}
#'   \item{H0_n}{number of incorrect target PSMs.}
#'   \item{H1_n}{number of correct target PSMs.}
#' }
#'
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

#' Calculate qvalues on the subset PSMs.
#'
#' @param df dataframe with at least 3 columns:
#'\describe{
#' \item{score}{score assigned to the peptide to spectrum match (PSM).}
#' \item{subset}{TRUE or 1 if PSM belongs to the subset in interest, FALSE or 0 otherwise.}
#' \item{decoy}{TRUE or 1 if decoy PSM, FALSE or 0 otherwise.}
#' }
#'Additional columns are allowed but ignored.
#' Target and decoy PSMs are assumbed to be from a competitive target decoy database search.
#' @param score_higher TRUE if a higher score means a better PSM.
#'
#' @return A data frame containing all columns in ``df''.
#' Following columns are added:
#' \describe{
#'\item{pi_0_cons}{conservative estimation of \eqn{\pi_0}.}
#'\item{FDR}{estimated subset PSM qvalues calculated according the competitive target decoy approach.}
#'\item{FDR_BH}{estimated subset PSM qvalues calculated according the Benjamini Hochbergh procedure.
#' When provided, non-subset decoy PSMs are used to stabilize estimates in small subsets}
#'\item{FDR_stable}{estimated subset PSM qvalues calculated with ``pi_0_cons''.
#' When provided, non-subset decoy PSMs are used to stabilize estimates in small subsets}
#' }
#' @export
#' @import dplyr
#'
#' @examples
calculate_fdr = function(df,score_higher = TRUE) {
  ## TODO Look in to dplyr::cume_dist to calculate FDR while avoiding arranging data frames
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
    mutate(FDR = cummax(cumsum((decoy == 1) & (subset == 1)) /
                        cumsum((decoy != 1) & (subset == 1)))) %>%
    # calculate BH FDR on subset
    mutate(FDR_BH = cummax((cumsum(decoy == 1) / sum(decoy == 1)) /
                           (cumsum((decoy != 1) & (subset == 1)) / sum((decoy != 1) & (subset == 1))))) %>%
    # calculate stable FDR on subset
    mutate(FDR_stable = FDR_BH * pi_0_cons) %>%
    ## Does not allow any FDR to be above 1 and set fdr of decoys and non subset PSMs to NA
    mutate_at(vars(FDR, FDR_BH, FDR_stable),
              funs(ifelse(!((decoy == 0) & (subset == 1)), NA,
                   ifelse(. > 1 , 1, .)))) %>%
    # Put dataframe back in original order
    ungroup %>% arrange(index) %>%
    select(-index)
}

#' Creates density plot of the pi0 distribution
#'
#' @param n_targets vector of observed target PSMs.
#' @param n_decoys vector of observed decoy PSMs.
#'
#' @return ggplot object.
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
#' @param score vector of quantiles of distribution 1 and 2
#' @param label vector of logical values. TRUE if score belongs to distribution 1
#' @param pi0 mixture coefficient of distribution 1 in distribution 2
#' @param title main title.
#' @param xlab label on x-axis.
#' @param ylab label on y-axis.
#'
#' @return ggplot object
#' @export
#' @import ggplot2
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


#' Plot diagnostic plots to evaluate assumptions from the search all, search subset strategy.
#'
#' @param df dataframe with at least 3 columns:
#'\describe{
#' \item{score}{score assigned to the peptide to spectrum match (PSM).}
#' \item{subset}{TRUE or 1 if PSM belongs to the subset in interest, FALSE or 0 otherwise.}
#' \item{decoy}{TRUE or 1 if decoy PSM, FALSE or 0 otherwise.}
#' }
#'Additional columns are allowed but ignored.
#' Target and decoy PSMs are assumbed to be from a competitive target decoy database search.
#'
#' @return
#' @export
#' @import dplyr
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

  return(list(pi0plot = p1, decoyall_targetsubset = p2,
              decoyall_decoy_subset = p3, decoysubset_target_subset = p4, all = p_all))
}


#' Sample scores for a random dataset
#'
#' The scores are sampled from two-component mixture distribution of Gaussians.
#'
#' @param H0_n number of incorrect subset target PSMs.
#' @param H1_n number of correct subset target PSMs.
#' @param decoy_n number of subset decoy PSMs.
#' @param decoy_large_n number of non subset decoy PSMs.
#' @param H0_mean mean of the incorrect subset target PSM distribution.
#' @param H1_mean mean of the correct subset target PSM distribution.
#' @param H0_sd sd of the incorrect subset target PSM distribution.
#' @param H1_sd sd of the correct subset target PSM distribution.
#' @param decoy_mean mean of the subset decoy PSM distribution.
#' @param decoy_sd sd of the subset decoy PSM distribution.
#' @param decoy_large_mean mean of the non subset decoy PSM distribution.
#' @param decoy_large_sd sd of the non subset decoy PSM distribution.
#'
#' @return
#' @export
#' @import dplyr
#' @keywords internal
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
  }
  bind_rows(d1,d2,d3,d4)
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
#' @import dplyr
#' @import ggplot2
#' @examples
plot_theo_dist = function(H0_mean=2.75, H1_mean=3.31,H0_sd=.13,H1_sd=.28,
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


##' Parses a mzID file generated by MS-GF+.
##'
##' See \url{https://omics.pnl.gov/software/ms-gf}
##' for more info on how to perform a database search on MSMS dataset with MS-GF+ and how to generate a mzID file.
##' Note that most functions in these package require data from a competitive target decoy search.
##'
##' @param mzid_path Location of the mzID file.
##'
#' @return A data frame containing the following 7 columns:
#' \describe{
#'\item{spec_id}{Id of the spectrum from the searched dataset file.}
#'\item{sequence}{Amino acid sequence matching the spectra.}
#'\item{protein_id}{Id of the sequence from the database file.}
#' \item{score}{score assigned to the peptide to spectrum match (PSM).}
#' \item{database}{Name of the database file used to search the spectra.}
#' \item{decoy}{TRUE if decoy PSM, FALSE otherwise.}
#' \item{database_size}{Number of sequences in the database file.}
#' }
##' @import dplyr
##' @import stringr
##' @import mzR
##' @export
##' @author
parse_msgf_mzid = function(mzid_path){
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Package mzR needed for this function to work. Please install it from bioconductor.",
      call. = FALSE)
  }

  d = openIDfile(mzid_path)
  if (!any(grepl('MS-GF+',mzidInfo(d)$software)))
    stop('Only mzID files generated by MS-GF+ can be parsed')
  t =  cbind(psms(d),score(d)[,-1]) %>% as_data_frame %>%
    filter(rank == 1) %>%
    transmute(spec_id = acquisitionNum,
              sequence = as.character(sequence),
              protein_id = as.character(DatabaseAccess),
              protein_id = str_replace(protein_id, '(XXX_)',''), ##remove XXX_ in decoy protein
              score= MS.GF.SpecEValue,
              database = gsub(".*/+([^/]+).fasta","\\1",levels(database(d)[['location']])),
              decoy = isDecoy,
              database_size = database(d)$numDatabaseSequences) %>%
    distinct
}

##' Checks if protein id appears in the headers of a fasta file.
##'
##' @param protein_id Vector of protein ids.
##' @param fastapath Location of the fasta file.
##' @return Logical vector, TRUE if protein id is present in provided fasta file, FALSE otherwise.
##' @export
##' @author
id_is_present = function(protein_id,fastapath){
  fas = stringr::read_lines(fastapath)
  headers = fas[grepl('>',fas,fixed = TRUE)]
  ## reducing search time by only searching the unique
  protein_id_unique = unique(protein_id)
  is_subset = sapply(protein_id_unique,function(p){any(stringr::str_detect(headers,fixed(p)))})
  protein_id %in% protein_id_unique[is_subset]
}


##' Preprocess data from a MS-GF mzID file.
##'
##' The parsed data frame from saas::parse_msgf_mzid function contains sometimes multiple entries for a spectrum. (eg. if sequence can be assigned to multiple protein ids). This function takes care of this by default.
##'
##' @param dat Data frame generated by the saas::parse_msgf_mzid function.
##' @param remove_target_decoy_PSM TRUE to remove PSMs that match both a target and decoy sequence.
##' @param remove_multiple_proteins_PSM TRUE to remove PSMs that can be assigned to multiple protein ids.
##'
##' @return Data frame with the same columns as ``dat'' and 1 column is added:
##' \describe{\item{non_subset}{TRUE if sequence can be assigned to a non subset protein id}}
##'
##' Every spectrum haves only 1 row in the data frame.
##'
##' @import dplyr
##' @export
##' @author
preprocess = function(dat,remove_target_decoy_PSM = TRUE, remove_multiple_proteins_PSM = FALSE){
  if (remove_target_decoy_PSM == TRUE){
    dat = group_by(dat,spec_id) %>% filter(!(any(decoy) & any(!decoy)))
  }
  if (remove_multiple_proteins_PSM == TRUE){
    dat = group_by(dat,spec_id) %>% filter(n() > 1)
  }
  dat = group_by(dat,spec_id,sequence,decoy,database,database_size,score) %>%
    summarise(subset = any(subset),non_subset = any(!subset),
              protein_id = paste0(protein_id,collapse = ';')) %>%
    ungroup
  dat
}
