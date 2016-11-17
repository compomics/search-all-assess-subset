# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @param rawinput
##' @param sep
##' @return
##' @author
check_input = function(rawinput,sep = ','){
  df = try({read.csv(text = rawinput,sep = sep)})
  if (class(df) == 'try-error'){return('cannot read input')}
  if(any(!(c('id','score','decoy','subset') %in% colnames(df))))
    {return('missing columns or wrong column names')}
  if(!is.numeric(df$score)) {return('Score should be numeric')}
  if(any(!(df$decoy ) %in% c(0,1))) {return('Decoy variable should be 0 or 1')}
  if(any(!(df$subset ) %in% c(0,1))) {return('Subset variable should be 0 or 1')}
  if(!any((df$decoy  == 1))) {return('No decoys available')}
  if(!any((df$decoy  == 0))) {return('No targets available')}
  if(!any((df$subset  == 1))) {return('No subset PSMs indicated')}
  if(!any((df$decoy[df$subset == 1]  == 1))) {return('No subset decoys available')}
##  return(df)
  return('Data format correct')
}
##' server
##'
##' .. content for \details{} ..
##' @param input
##' @param output
##' @param session
##' @return
##' @author
server = function(input, output, session) {
  observe({cat(input$tabs, '\n')})
  ## logic data input tab
  #######################
  ## read data from given path
  rawfiledata = reactive({readLines(req(input$file1$datapath))})
  # ## or read example data
  rawexamplefiledata = reactive({
    req(input$action)
    readLines(system.file("extdata", "cytoplasm.csv", package = "saas"))})

  observe({updateTextInput(session,inputId = 'rawinput',value = paste0(rawfiledata(),collapse = '\n'))})
  observe({updateTextInput(session,inputId = 'rawinput',value = paste0(rawexamplefiledata(),collapse = '\n'))})

  rawdata = reactive({req(input$rawinput)})

  generate_errormessage = reactive({check_input(rawdata())})

  output$errormessage = reactive({generate_errormessage()})

  df = reactive({req(generate_errormessage() == 'Data format correct')
    read.csv(text = rawdata(),sep = ',')})

  ## logic calculation tab
  ########################
  df_calc = reactive({calculate_fdr(df(), score_higher = input$scorehigher)})

  output$contents <- renderDataTable({
    mutate_at(req(df_calc()),
              vars(score,pi_0_cons, FDR, FDR_BH, FDR_stable),
              funs(round(.*1000)/1000))
  }, options = list(
    lengthMenu = c(10, 20, 100, 1000),
    pageLength = 20
  ))

  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$file1$datapath, '_fdr.csv', sep = '')
    },
    content = function(file) {
      write.csv(df_calc(), file)
    })

  ## logic diagnistic tab
  #######################
  output$plot1 = renderPlot({
    plot_diag(df())$all
  }, width = 940,height = 410)

  ## logic simulation tab
  #######################
  observeEvent((input$make_subset_action),{
    par = simulate_subset(input$subset_n,input$subset_pi0)
    updateTextInput(session,inputId = 'H1_n',value = par$H1_n)
    updateTextInput(session,inputId = 'H0_n',value = par$H0_n)
    updateTextInput(session,inputId = 'decoy_n',value = par$decoy_n)
  })

  ## check if decoy dist are the same and change value of respective boxes
observe({
    if(input$check_subset_decoys_same){
      updateTextInput(session,inputId = 'decoy_mean',value = input$H0_mean)
      updateTextInput(session,inputId = 'decoy_sd',value = input$H0_sd)
    }
    if(input$check_subset_decoys_same){
      updateTextInput(session,inputId = 'decoy_extra_mean',value = input$decoy_mean)
      updateTextInput(session,inputId = 'decoy_extra_sd',value = input$decoy_sd)
    }
  })

## check if there are enough PSMs for diagnostics
#enough_targets_decoys = reactive({input$H1_n + input$H0_n) > 0) & (input$decoy_n > 0)})
#output$not_enough_decoys_targets = reactive({ifelse(enough_targets_decoys,'',
#                                                    'You need at least 1 decoy or target to visualize the diagnostic plots!')})

output$plot_theo = renderPlot({
  plot_theo_dist2(
                 H0_mean = input$H0_mean,
                 H1_mean = input$H1_mean,
                 H0_sd = input$H0_sd,
                 H1_sd = input$H1_sd,
                 decoy_mean = input$decoy_mean,
                 decoy_sd = input$decoy_sd,
                 decoy_extra_mean = input$decoy_extra_mean,
                 decoy_extra_sd = input$decoy_extra_sd)$plot
}, width = 400, height = 300)

  output$plot_sim_diag = renderPlot({
    req(input$action_simulate)
    isolate({
    decoy_n = input$decoy_n
    if (input$check_also_decoys){
      decoy_n =simulate_subset(input$subset_n,input$pi0)$decoy_n
      updateTextInput(session,inputId = 'decoy_n',value = decoy_n)
    }
    target_n = input$subset_n - decoy_n
    pi0 = rpi0(1, target_n, decoy_n)
    H0_n = round(pi0 * target_n)
    H1_n = target_n - H0_n
      plotdata = sample_dataset(
        H1_n = H1_n,
        decoy_n = decoy_n,
        decoy_large_n = input$decoy_extra_n,
        H0_n = H0_n,
        H0_mean = input$H0_mean,
        H1_mean = input$H1_mean,
        H0_sd = input$H0_sd,
        H1_sd = input$H1_sd,
        decoy_mean = input$decoy_mean,
        decoy_sd = input$decoy_sd,
        decoy_large_mean = input$decoy_extra_mean,
        decoy_large_sd = input$decoy_extra_sd)
       plot_diag(plotdata)$all
    })
  }, width = 940,height = 410)
}


##' UI shiny app
##'
##' .. content for \details{} ..
##' @return
##' @author
ui = function() fluidPage(navbarPage(
  "Search All, Assess Subset",
  tabPanel('Data input',
           sidebarLayout(
             sidebarPanel(width = 7,
                          HTML(markdownToHTML(text =
                                                '
This tool demonstrates the search-all-asses-subset method to control
the False Discovery Rate (FDR) when you are only interested in a subset
of the identified PSMs from a shotgun proteomics experiment.

The workflow is as follows:

1. search the spectra for all expected peptides
2. Removing irrelevant peptides
3. Calculate FDR

You can load the data from step **1.** into this webtool to perform
step **2.** and **3.**'
                          )), hr(),
                          fileInput ('file1',
                                     'Choose Text File',
                                     accept = c('text/csv',
                                                'text/comma-separated-values,text/plain',
                                                '.csv')
                          ),
                          checkboxInput("scorehigher", "Are higher scores better?", TRUE),
                          HTML(markdownToHTML(text =
                                                '
Use the **Choose File** button above to load a CSV file from your computer
or paste the data in the **text area** to the right.

Adhere to the following format:
* Every Row is a unique PSM
* Every row contains:
 + **id**: Can be any text or number
 + **score**: Score given to the PSM, higher scores are better
 + **decoy**: 0 or 1; 1 indicates that the PSM matches a decoy peptide sequence
 + **subset**:  0 or 1: 1 indicates that the PSM matches a subset peptide sequence
* The first row are the column names
* All columns should be comma separated

Example input:
<pre>
id,score,decoy,subset
1,7.67,1,0
2,10.99,1,1
3,75.10,0,0
4,73.83,0,1
</pre>
**Warning:** use this tool only with PSMs from a competitive target-decoy search!
The minimal required input are PSMs from subset targets and decoys (row 2 and 4 in example input).
Additional decoy PSMs (row 1) are used for a more stable FDR calculation.
Non subset targets (row 3) are ignored in the analysis.

**Load an example dataset:**
')),
                          actionButton("action", label = "Load example")
             ),
             mainPanel(width = 5,h3(textOutput('errormessage',container = span)),
                       tags$textarea(id = 'rawinput',label = '',
                                     rows = '20',cols = 40))
           )),



  tabPanel('Calculate FDR',
           sidebarLayout(
             sidebarPanel(width = 12,
                          withMathJax(),
                           HTML(markdownToHTML(text =
'The following columns were calculated:
* **pi_0_cons**: A conservative estimation of pi_0. (pi_0_cons = (#decoys+1) / (#targets+1})
* **FDR**: The estimated FDR at this score cutoff for subset PSMs.
Missing for non-subset or decoy PSMs.
* **FDR_stable**: The estimated stable FDR at this score cutoff for subset PSMs.
This FDR is estimated from the complete decoy set and pi_0_cons.
Missing for non-subset or decoy PSMs.
Please check the diagnostics on the next tab.
* **FDR_BH**: The estimated stable FDR at this score cutoff for subset PSMs.
This FDR is estimated from the complete decoy set and according the Benjamini-Hochberg FDR procedure.
Please check the diagnostics on the next tab.
Missing for non-subset or decoy PSMs.

Download the results to your computer:
')),
                          downloadButton('downloadData', 'Download')),

             mainPanel(width = 12,dataTableOutput('contents'))
           )),
  tabPanel('Check diagnostic plots',
           sidebarLayout(
             sidebarPanel(width = 12,
                          HTML(markdownToHTML(text =  '
These are diagnostic plots to evaluate the quality of the decoy set and the estimation of the fraction of incorrect target PSMs (pi0).
This allows an informed choice on the use of the stable all-sub FDR estimator and the large decoy set.

**Panel a** shows the posterior distribution of pi_0 given the observed number of target and decoy PSMs in the subset.
The vertical line indicates our conservative estimation of pi_0 used in the calculations.
At very high pi_0 uncertainty (broad peak), you can opt to use the BH procedure to minimize sample to sample variability.
However, this will come at the expense of too conservative PSM lists.

Our improved TDA for subsets relies on the assumption that incorrect subset PSMs and the complete set of decoys are following the same distribution.
This distributional assumption can be verified through a PP-plot where the empirical Cumulative Distribution Function (eCDF) of the decoys is plotted against the eCDF of the subset target PSMs.
The PP-plots in **panel b - d** display the target subset PSMs plotted against all decoy PSMs from the complete search, the decoy subset PSMs plotted against all decoy PSMs from the complete search, and the target subset PSMs plotted against the decoy PSMs from the complete search, respectively.
The full line in panel **b** and **d** indicates a line with a slope of pi_0.
The full line in panel **c** indicates the identity line.
The first part of the plot in **a** and **b** is linear with a slope that equals pi_0.
This indicates that the decoy distribution and the mixture component for incorrect PSMs of the target mixture distribution coincide.
The second part of the plot deviates from the line towards higher percentiles and will ultimately become vertical (decoy percentile = 1).
If we see this profile in panel **a**, we have a good indication that the set of decoys from the complete search is representative for the mixture component for incorrect PSMs of the target mixture distribution.
Deviations from this pattern might be subtle, therefore we provide the PP plots in **c** and **d** to support the conclusion drawn from panel **b**.
The PP-plot in panel **c**  shows the subset decoy PSMs plotted against all decoy PSMs.
The whole plot should follow the identity line, indicating that the complete set of decoys is a good representation of the subset decoys.
To verify that the subset decoys (and thus also the complete set of decoys) are representative for the mixture component for incorrect PSMs of the target mixture distribution, we look at the PP-plot of the subset decoys against the subset targets in panel **d**.
The profile should look as described for panel **b**.
If the profile matches in panel **d** but does not for panel **b**, then we suggest to not use the extra decoy set and use only the subset decoys for FDR estimation.

When you are not sure how the diagnostic plots should look like, you can simulate your own data under various (erratic) settings in the simulation tab.
'))
             ),
              mainPanel(width = 12,plotOutput('plot1'))
           ))
,
tabPanel('simulation',
         sidebarLayout(
           sidebarPanel(width = 12,
                        HTML(markdownToHTML(text =  '
In this tab, you can simulate your own data and look at examples of diagnostic plots.
Random datasets are generated based on a observed number of target and decoy subset PSMs.
Optionally, you can also generate a random number of subset decoy PSMs based on the observed number of subset target PSMs and a theoretic pi0 that you can choose.
In the default setting is the decoy distribution equal to the incorrect subset target distribution.
This means that the diagnostic plots from the simulated datasets are exemplary for experimental settings where the assumptions are grounded and you can safely use the decoys for FDR estimation.
Optionally, you can change the mean and standard deviation of the subset decoys and/or the large set of extra decoys, violating the assumption that they are representative for the incorrect subset targets.
Plots generated by these simulated datasets are examples of diagnostic plots when the assumptions are violated an you should not use these decoys for FDR estimation.
')))
,
mainPanel(width = 12,
          column(6,
                 fluidRow(column(6, '# subset target  PSMs'),
                          column(6, tags$input(id = 'subset_n',type = "number",
                                               value = 200,min = 1,step = 10))
                          ),
                 fluidRow(column(6, '# subset decoy PSMs'),
                          column(6, tags$input(id = 'decoy_n',type = "number",
                                               value = 20,min = 1,step = 10))),
                 fluidRow(column(6, '# extra decoys'),
                          column(6, tags$input(id = 'decoy_extra_n',type = "number",
                                               value = 2000,min = 1,step = 10))))
         ,
          column(1,actionButton("action_simulate", label = "GO")),
          column(5,checkboxInput("check_also_decoys", label = "Simulate also decoys"),
                        conditionalPanel(
                          condition = "input.check_also_decoys == true",
                          column(4, 'with pi0:'),
                          column(7, tags$input(id = 'pi0',type = "number",
                                               value = .2,min = 0,max = 1,step = .1))
                        )
                 )
         ,
          column(12,br(),
                 column(12,
                        HTML(markdownToHTML(text ='__Subset target mixture distribution__')))
                ,
                 column(6, 'incorrect targets',
                        wellPanel(
                          fluidRow(column(6, numericInput('H0_mean', 'mean', 2.75,step = .1)),
                                   column(6, numericInput('H0_sd', 'sd', 0.13,step = .1)))))
                ,
                 column(6, 'correct targets',
                        wellPanel(
                          fluidRow(column(6, numericInput('H1_mean', 'mean', 3.31,step = .1)),
                                   column(6, numericInput('H1_sd', 'sd', 0.28,step = .1))))))
         ,
          column(6,
                 HTML(markdownToHTML(text ='__Subset decoys distribution__')),
                 checkboxInput('check_subset_decoys_same','Same as incorrect target PSMs?',value = TRUE)
                ,
                        conditionalPanel(
                          condition = "input.check_subset_decoys_same == false",
                          wellPanel(
                            fluidRow(column(6, numericInput('decoy_mean', 'mean', 0,step = .1)),
                                     column(6, numericInput('decoy_sd', 'sd', 0,step = .1))))
                        )
                ,
                 HTML(markdownToHTML(text ='__Extra decoys distribution__')),
                 checkboxInput('check_extra_decoys_same','Same as subset decoy PSMs?',value = TRUE)
                ,
                        conditionalPanel(condition = "input.check_extra_decoys_same == false",
                                         wellPanel(
                                           fluidRow(column(6, numericInput('decoy_extra_mean', 'mean', 0,step = .1)),
                                                    column(6, numericInput('decoy_extra_sd', 'sd', 0,step = .1))))
                                         )
                 )
                 ,
                 column(6,plotOutput('plot_theo'))
          ,
                             fluidRow(column(12,plotOutput('plot_sim_diag')))
          )
))

))


#' Launches the GUI version of saas.
#'
#' To easily launch the GUI outside an R session (eg. on a server),
#' you can run R -e "library(saas);saas_gui()" from the terminal (on linux/mac).
#'
#' see ?shiny::runApp for help all parameters.
#'
#' @param options See help of shiny::shinyApp for details
#' @return
#' @export
#' @import shiny
#' @import markdown
#' @examples
saas_gui = function(options = list(port = 3320, host  = "0.0.0.0"))
  shinyApp(ui = ui(), server = server,options = options)
