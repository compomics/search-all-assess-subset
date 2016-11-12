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

library(shiny)
library(dplyr)
library(markdown)
library(ggplot2)
source('./functions.r')

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

server = function(input, output, session) {
  ## logic data input tab
  #######################
  
  ## read data from given path
  rawfiledata = reactive({readLines(req(input$file1$datapath))})
  # ## or read example data
  rawexamplefiledata = reactive({
    req(input$action)
    readLines('./cytoplasm.csv')})
  
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
              vars(score,pi_0,pi_0_cons, FDR, FDR_BH, FDR_stable),
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
    if(input$subsetdecoy){
      updateTextInput(session,inputId = 'decoy_mean',value = input$H0_mean)      
      updateTextInput(session,inputId = 'decoy_sd',value = input$H0_sd)
    }
    if(input$largedecoy){
      updateTextInput(session,inputId = 'decoy_large_mean',value = input$decoy_mean)      
      updateTextInput(session,inputId = 'decoy_large_sd',value = input$decoy_sd)
    }
  })

## check if there are enough PSMs for diagnostics
#enough_targets_decoys = reactive({input$H1_n + input$H0_n) > 0) & (input$decoy_n > 0)})
#output$not_enough_decoys_targets = reactive({ifelse(enough_targets_decoys,'',
#                                                    'You need at least 1 decoy or target to visualize the diagnostic plots!')})

output$plot_theo = renderPlot({
  plot_theo_dist(H1_n = input$H1_n,
                 decoy_n = input$decoy_n,
                 decoy_large_n = input$decoy_large_n,
                 H0_n = input$H0_n,
                 H0_mean = input$H0_mean,
                 H1_mean = input$H1_mean,
                 H0_sd = input$H0_sd,
                 H1_sd = input$H1_sd,
                 decoy_mean = input$decoy_mean,
                 decoy_sd = input$decoy_sd,
                 decoy_large_mean = input$decoy_large_mean,
                 decoy_large_sd = input$decoy_large_sd)$plot 
}, width = 400, height = 400)  


output$plot_sim_diag = renderPlot({
  input$simulate_action
  d = isolate(sample_dataset(H1_n = input$H1_n,
                             decoy_n = input$decoy_n,
                             decoy_large_n = input$decoy_large_n,
                             H0_n = input$H0_n,
                             H0_mean = input$H0_mean,
                             H1_mean = input$H1_mean,
                             H0_sd = input$H0_sd,
                             H1_sd = input$H1_sd,
                             decoy_mean = input$decoy_mean,
                             decoy_sd = input$decoy_sd,
                             decoy_large_mean = input$decoy_large_mean,
                             decoy_large_sd = input$decoy_large_sd))
  plot_diag(d)$all
}, width = 940,height = 410)
}

ui = fluidPage(navbarPage(
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
                          HTML(markdownToHTML(text = 
'The following columns were calculated:
* **pi_0**: The estimated fraction of false positives in the subset.
* **pi_0_cons**: A conservative estimation of pi_0.
* **FDR**: The estimated FDR at this score cutoff for subset PSMs.
Missing for non-subset or decoy PSMs.
* **FDR_stable**: The estimated stable FDR at this score cutoff for subset PSMs.
This FDR is estimated from the complete decoy set and a conservative pi_0.
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
  tabPanel('Check diagnostics plots',
           sidebarLayout(
             sidebarPanel(width = 12,
                          HTML(markdownToHTML(text =  '
**TODO: explain diagnostics according paper**
The reweighted subset FDR relies on the assumption that the properties of 
the PSMs in the subset and the complete set are similar.
You can use use the histograms of the scores below to visually check this 
assumption.
Grey is the distribution of the target scores.
Blue is the distribution of the decoy scores. This is the estimated distribution
of incorrect PSMs.
Red is the estimated distribution of the correct PSMs.

* The location of the decoy distribution should be similar in the subset and 
complete set
* The location of the estimated correct PSM distribution should be similar 
in the subset and complete set
* The shape of the decoy distribution should be similar in the subset and 
complete set
* The shape of the estimated correct PSM distribution should be similar in 
the subset and complete set
'))
             ),
              mainPanel(width = 12,plotOutput('plot1'))
           ))
,
tabPanel('simulation',
         column(12,
                HTML(markdownToHTML(text =  'explaination adn more bla bla bla bla bla bla'))
             , 'Generate',
             tags$input(id = 'subset_n',type = "number", value = 200,min = 1,step = 10),
             'subset PSMs with pi0 =',
                tags$input(id = 'subset_pi0',type = "number", value = .1,step = .1,min = 0,max = 1),
             actionButton("make_subset_action", label = "Go"),
                         fluidRow(column(6,
                                # fluidRow(column(2,'Generate'),
                                #          column(3,numericInput('subset_n','', value = 100)),
                                #          column(4,'subset PSMs with pi0 ='),
                                #          column(3,numericInput('test', '',value = .1))),
                                "Correct target distribution",
                                fluidRow(column(4, numericInput('H1_n', 'number', 160)),
                                         column(4, numericInput('H1_mean', 'mean', 3.31,step = .1)),
                                         column(4, numericInput('H1_sd', 'sd', .28,step = .1)))
                                ,"Incorrect target distribution",
                                fluidRow(column(4, numericInput('H0_n', 'number', 40)),
                                         column(4, numericInput('H0_mean', 'mean', 2.75,step = .1)),
                                         column(4, numericInput('H0_sd', 'sd', .13,step = .1)))
                                ,"Subset decoy distribution",
                                fluidRow(column(4, numericInput('decoy_n', 'number', 40)),
                                         column(4, numericInput('decoy_mean', 'mean', 2.75,step = .1)),
                                         column(4, numericInput('decoy_sd', 'sd', .13,step = .1)))
                                ,checkboxInput('subsetdecoy','mean and sd of incorrect targets',value = TRUE)
                                ,"Large decoy distribution",
                                fluidRow(column(4, numericInput('decoy_large_n', 'number', 2000)),
                                         column(4, numericInput('decoy_large_mean', 'mean', 2.75,step = .1)),
                                         column(4, numericInput('decoy_large_sd', 'sd', .13,step = .1)))
                                ,checkboxInput('largedecoy','mean and sd of subset decoys',value = TRUE))
                         ,
                         column(6,plotOutput('plot_theo'))
                )
                ,
                actionButton("simulate_action", label = "Simulate"),
                
                fluidRow(column(12,plotOutput('plot_sim_diag')))
         )
)))


shinyApp(ui = ui, server = server,options = list(port = 3320, host  = "0.0.0.0"))
