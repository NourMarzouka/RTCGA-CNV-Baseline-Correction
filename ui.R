library(shiny)
library(markdown)
library(shinythemes)

shinyUI(
  navbarPage("Baseline Correction for Copy Number Data from Cancer Samples", 
             id="baseCN",
             inverse = F, 
             theme = "bootstrap.min.css",
             tabPanel("Description",
                      fluidPage(titlePanel("Description"),
                                sidebarLayout(
                                  sidebarPanel(
                                    
                                    wellPanel(tags$p("Cancer samples face some difficulties regarding the determination of the copy number status due to the false sample centering and baseline shifting. Without solving this issue the copy number calling will be inaccurate. TCGA.CNV.BaselineCorrection package was designed to correct the baseline in cancer samples using the Maximum Density Peak Estimation (MDPE) method.
"),
                                              tags$p("The main advantages for TCGA.CNV.BaselineCorrection package are: Fast (few seconds per sample), high accuracy rate, in-sample correction, no input parameters needed, low computer sources required, and adaptable for different technologies.")),
                                    #tags$hr(),
                                    wellPanel(fluidRow(uiOutput('LoadDataButtons'))),
                                    wellPanel(tags$strong("Author: Nour-al-dain Marzouka"),
                                              tags$p(),
                                              tags$strong("Maintainer: Nour-al-dain Marzouka <nour-al-dain.marzouka@med.lu.se>"),
                                              tags$p(),
                                              tags$p("License: GPL (>= 2)"),
                                              tags$p(),
                                              tags$hr(),
                                              tags$p("Check the FAQs tab before using the tool"),
                                              tags$p("Note: For large datasets (>200), it is much faster to run this tool on your own computer, you can download it from the like below:"),
                                              downloadButton('downloadtool3', 'Download TCGA.CNV.BaselineCorrection.zip'))),
                                  
                                  mainPanel(fluidRow(img(src='photo_plan.png', align = "left")),
                                            tags$hr()
                                  )))),      
             
             tabPanel("Upload Your Own Data",
                      fluidPage(titlePanel("Upload your own segmentation data"),
                                sidebarLayout(
                                  sidebarPanel(
                                    fileInput('file1', 'Choose Regions CSV File', 
                                              accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
                                    tags$hr(),
                                    tags$div(class="header", checked=NA,
                                             tags$p("CSV manipulation")),
                                    
                                    checkboxInput('header', 'Header', TRUE),
                                    radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                                    radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"'),
                                    tags$div(tags$p("Select which column is:")),
                                    selectInput("RegionSample", "Sample:",NULL),
                                    selectInput("RegionChromosome", "Chromosome:",NULL),
                                    selectInput("Regionbpstart", "bp.Start:",NULL),
                                    selectInput("Regionbpend", "bp.End:",NULL),
                                    selectInput("RegionNumMark", "Num.of.Markers:",NULL),
                                    selectInput("RegionMean", "Mean:",NULL)),
                                  
                                  mainPanel(tableOutput('csvtableRegions'),
                                            uiOutput('regionsbuttonsGo2Sample'),
                                            uiOutput('regionsbuttonsGo2PlotRaw'))))),
             
             tabPanel("Upload sample list",
                      fluidPage(titlePanel("Optional Upload: Sample List for the comments on the plots"),
                                sidebarLayout(
                                  sidebarPanel(
                                    fileInput('file2', '*Optional* Choose Sample List File',
                                              accept=c('text/csv', 'text/comma-separated-values,text/plain','.csv')),
                                    tags$hr(),
                                    tags$div(class="headersamp", checked=NA,
                                             tags$p("CSV manipulation")),
                                    
                                    checkboxInput('headersamp', 'Header', TRUE),
                                    radioButtons('sepsamp', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
                                    radioButtons('quotesamp', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"'),
                                    tags$div(tags$p("Select which column is:")),
                                    selectInput("SampleNumber", "Number:",NULL),
                                    selectInput("SampleSample", "Sample:",NULL),
                                    selectInput("SampleComment", "Comment:",NULL)),
                                  
                                  mainPanel(tableOutput('csvtableSample'),
                                            uiOutput('sampleButtonG2Raw'))))),
             
             tabPanel("Load TCGA data",
                      fluidPage(titlePanel("TCGA copy number variants data"),
                                sidebarLayout(
                                  sidebarPanel(
                                    wellPanel(uiOutput('tcga'), 
                                              uiOutput('tcgaSamplenumber'),
                                              actionButton("load.this.TCGA","Load these samples")),
                                    wellPanel(tags$h5("TCGA data is loaded from the Bioconductor package RTCGA.CNV"),
                                              tags$h4(helpText( a("RTCGA.CNV website", target="_blank",href="http://www.bioconductor.org/packages/release/data/experiment/html/RTCGA.CNV.html")))),
                                    wellPanel(tags$h3("IMPORTANT"),
                                              tags$h5("It is much faster for large datasets like TCGA to run this tool on your own computer, you can download it from the like below:"),
                                              downloadButton('downloadtool', 'Download TCGA.CNV.BaselineCorrection.zip'))),
                                  
                                  mainPanel(tags$strong("First 10 lines from the selected data:"),
                                            tableOutput('tableTCGA'),
                                            actionButton("SampleActionButton2", label = "Plot samples?"))))),
             
             tabPanel("Plot raw data",
                      fluidPage(titlePanel("Raw plot"),
                                sidebarLayout(
                                  sidebarPanel(sliderInput("NumberSampleSlider",
                                                           tagList(p("Number of Samples:"),p("(Note: you can select all samples without plotting them)")),
                                                           min=0,
                                                           max=0,
                                                           value=1,
                                                           ticks=FALSE ),
                                               sliderInput("NumberMarkerSlider",
                                                           "Minimum number of markers per segment:",
                                                           min = 0,
                                                           max = 1000,
                                                           value = 20,
                                                           step=10),
                                               sliderInput("NumberCutoffSlider",
                                                           "Cutoff:",
                                                           min = 0,
                                                           max = 1,
                                                           value = 0.1,
                                                           step = 0.05),
                                               checkboxInput('ShowComments', 'Show Comments?', FALSE),
                                               actionButton("PlotActionButtonGo2Autocorrect", label = "Autocorrect selected samples")),
                                  
                                  mainPanel(uiOutput("plotraw"))))),
             
             tabPanel("Baseline Correction & QC",
                      fluidPage(titlePanel("Baseline Correction & Download"),       
                                wellPanel(wellPanel(
                                  p("Download panel (Baseline autocorrection were performed only for the selected samples)"),                                                
                                  uiOutput("downloadButtons")),
                                  wellPanel(tags$h4("Please check the autocorrected plots below. Use the slider to choose another baseline"),
                                            tags$hr(),
                                            p("Each page contains 25 sample maximum...click next to see the next plots"),
                                            pageruiInput('pager', page_current = 1, pages_total = 1))),
                                wellPanel(fluidRow(column(12,uiOutput("autocorrection"))),
                                          p("If you do not see the plots ... please wait few moments to complete ...")
                                ))),
             
             tabPanel("FAQs",
                      fluidPage(titlePanel("Help"),
                                sidebarLayout(
                                  sidebarPanel(
                                    wellPanel(tags$strong("For any questions please contact: nour-al-dain.marzouka@med.lu.se"),
                                              #tags$hr(),
                                              #tags$br(),
                                              #tags$br(),
                                              tags$br()#,
                                              #HTML('<script type="text/javascript" id="clustrmaps" src="//cdn.clustrmaps.com/map_v2.js?u=eOcy&d=WlzyZjY8j4Uojtqc_XjZCbgfu8A4vpcUfJYO2_Itpnk"></script>')
                                  )),
                                  mainPanel(wellPanel(tags$h3("IMPORTANT & TIPS"),
                                                      tags$h4("1) Restart the application if you want to use different set of samples."),
                                                      tags$h4("2) To avoid any problem with the sliders, let the tool finish the calculations before changing the slider again."),
                                                      tags$h4("3) If the slider started to jump between two values repeatedly then just go to another tab then go back to the baseline correction tab.")),
                                            wellPanel(tags$h3("How baseline correction works?"),
                                                      tags$h4(helpText( a("Click Here to check the paper", target="_blank",href="http://bioinformatics.oxfordjournals.org/content/32/7/1080.full")))),
                                            wellPanel(tags$h3("TCGA data"),
                                                      tags$h5("TCGA data is loaded from the Bioconductor package RTCGA.CNV"),
                                                      tags$h4(helpText( a("RTCGA.CNV website", target="_blank",href="http://www.bioconductor.org/packages/release/data/experiment/html/RTCGA.CNV.html"))),
                                                      tags$h4(helpText( a("RTCGA.CNV on github", target="_blank",href="https://github.com/RTCGA")))),
                                                              
                                                              wellPanel(tags$h3("User input data"),
                                                                        tags$h4("Here you can find example data for 4 samples (segments file & sample list with comments)"),
                                                                        downloadButton('downloadexample', 'Download example input data')
                                                              ),
                                                              
                                                              wellPanel(tags$h3("Quality Control (QC)"),
                                                                        tags$h4("The QC measurements includes the number of the segments, Standard Deviation (SD), InterQuartile Range (IQR), Median Absolute Pairwise Difference (MAPD), the highest segments density peak sharpness, the area under the curve between the cut-offs, and finally the heights (i.e. values) of the density function at the upper and lower cut-offs.  Only the last three measurements are sensitive to the baseline position (i.e. sample center) and the cut-offs. Lower values of the density function at the cut-offs indicate fewer segments located near the cut-off and lower possibility to have segments with log values higher than the cut-off by chance or due to the noise in the data.")
                                                              ),
                                                              
                                                              wellPanel(tags$h3("Output (Downloadable files)"),
                                                                        tags$h4("The corrected CN data can be downloaded as tables and plots. To facilitate the tracking of the modifications in the data, TCBC provides a file that contains the shift value for each sample that was performed by TCBC and the user."),
                                                                        tags$h4("TCBC generates 4 output files:"),
                                                                        tags$h4("1) Segmentation file, which is similar to the input segmentation file but with shifted log intensity values based on the new baselines in the samples."),
                                                                        tags$h4("2) Corrected segmentation plots."),
                                                                        tags$h4("3) QC measurements."),
                                                                        tags$h4("4) Record for all the shifts that were applied on the samples by the user or by the tool.")
                                                                        
                                                              ),
                                                              
                                                              
                                                              wellPanel(tags$h3("Download the tool"),
                                                                        tags$h4(helpText( a("The tool online: https://copynumber.shinyapps.io/RTCGA-CNV-BaselineCorrection/", target="_blank",href="https://copynumber.shinyapps.io/RTCGA-CNV-BaselineCorrection/"))),
                                                                        tags$h4("Note: For large datasets (>200), it is much faster to run this tool on your own computer, you can download it from the like below:"),
                                                                        downloadButton('downloadtool2', 'Download TCGA.CNV.BaselineCorrection.zip'))))))
                                
                      )) # END of UI 
             