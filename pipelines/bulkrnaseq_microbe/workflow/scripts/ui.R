library(shiny)
library(dplyr)


topTables<-readRDS("data/topTables.rds")
allComparisons<-names(topTables)[1:(length(topTables)-1)]
 for (i in 1:length(allComparisons)){
   names(allComparisons)[i]<-paste(strsplit(allComparisons[i], "[.]")[[1]][1], "versus", strsplit(allComparisons[i], "[.]")[[1]][2], sep=" ")
 }



shinyUI(navbarPage("Select tab:", 
  tabPanel("Differential gene expression", 
  titlePanel("Visualisation of results from differential gene expression analysis"),
  fluidRow(
    column(3,
      wellPanel(
            radioButtons(inputId="plotType", label = "Select the type of plot", 
                 choices=c("Volcano plot", "Heatmap with top DE genes", "Heatmap with user-defined genes"), selected = "Volcano plot")
            ),
      conditionalPanel(
        condition="input.plotType=='Volcano plot' | input.plotType=='Heatmap with top DE genes'",
       br(),
        selectInput(inputId="comparison", label="Select a comparison",
                  choices=names(allComparisons), selected=allComparisons[1]),
        br()
       ),
      
      ## Panels to show for volcano plots:
      
      conditionalPanel(
        condition="input.plotType=='Volcano plot'",
        textInput(inputId="threshold", 
                 label="Significance threshold", 
                 value = "0.05"),

        p("Genes with an adjusted P-value < threshold will be shown in purple",
          style="font-size:9pt"),
        br(),
        textInput(inputId="gene", 
                label="Highlight this gene in the plot (Ensembl ID):", 
                value=NULL),
        p("Ensembl IDs look like this: ENSMUSG00000050936. You can find these IDs in the excel tables or on the ", 
          a("Ensembl website.", href="http://www.ensembl.org/index.html"),
          "Make sure to select the correct species!",
          style="font-size:9pt"),
        hr(style="border-color: black"),
        downloadButton("downloadPDF", "Download plot")
      ),
      
      # panels to show for heatmap with top DE genes
      conditionalPanel(
        condition="input.plotType=='Heatmap with top DE genes'",
        numericInput(inputId="GeneNumber", label = "Display the top N genes:", value = 30, 
                     min = 1, max = 200, step = 1),
        p("Shows N genes with lowest adjusted P-values", style="font-size:9pt"),
        br(),
        hr(style="border-color: black"),
        downloadButton("downloadHeatmap", "Download heatmap")
       ),
  
    # panels to show for heatmap with user-defined genes
    conditionalPanel(
      condition="input.plotType=='Heatmap with user-defined genes'",
      br(),
      p(strong("Select data to be displayed")),
      uiOutput("resetable_input"),
      tags$hr(),
      actionButton("reset_input", "Reset"),
      tags$script('
                  Shiny.addCustomMessageHandler("myCallbackHandler",
                      function(color){
                        document.getElementById("submitButton").style.backgroundColor = color;
                      });
                  '),
      #div(style="text-align:center", 
      #    actionButton("submitButton", "Submit choices!")),
      actionButton("submitButton", "Create Heatmap!"),
      br(),
      hr(style="border-color: black"),
      downloadButton("downloadCustomHeatmap", "Download heatmap")
     )
    ),
    
  column(9,
    conditionalPanel(
      condition="input.plotType=='Volcano plot'",
      hr(style="border-color: black"),
      htmlOutput("figureLegend"),
      hr(style="border-color: black"),
      plotOutput("p", click = "plot_click"),
      hr(style="border-color: black"),
      tableOutput("click_info")),
    
      conditionalPanel(
      condition="input.plotType=='Heatmap with top DE genes'",
      mainPanel(plotOutput("heatmap"))),
      
      conditionalPanel(
        condition="input.plotType=='Heatmap with user-defined genes'",
        mainPanel(plotOutput("customHeatmap"),
                 strong(textOutput("title")),
                 textOutput("unknownIds")))
    )
  )),
  
  tabPanel("GO term enrichment", 
    titlePanel("Visualise genes within Gene Ontology (GO) terms"), 
    fluidRow(
      column(3,
             wellPanel(
               radioButtons(inputId="plotTypeGO", label = "Select the type of plot", 
                            choices=c("MA Plot", "Heatmap"), 
                            selected = "MA Plot")
               ),
           
               selectInput(inputId="comparisonGO", label="Select a comparison",
                           choices=names(allComparisons), selected=allComparisons[1]),
               br(),
               radioButtons(inputId="ontology", label = "Select a GO ontology", 
                            choices=c("BP", "MF", "CC"), selected = "BP"),
               br(),
               uiOutput("goterms"),
               p("The top 30 GO terms are listed for each ontology ranked by increasing weight01.Fisher P-values (no P-value cutoff is applied!)", style="font-size:9pt"),
               br(),
               
               
               conditionalPanel(
                 condition="input.plotTypeGO=='MA Plot'",
                 textInput(inputId="thresholdMAPlot", 
                              label="Significance threshold", 
                              value = "0.05"),
                 hr(style="border-color: black"),
                 downloadButton("downloadPDF_MAPlot", "Download plot"),
                 downloadButton("downloadGeneset", "Download data table")
               ),
       
              conditionalPanel(
                 condition="input.plotTypeGO=='Heatmap'",
                 hr(style="border-color: black"),
                 downloadButton("downloadPDF_GOHeatmap", "Download plot"),
                 downloadButton("downloadGeneset2", "Download data table")
              )
            ),
        column(9,
               conditionalPanel(
                 condition="input.plotTypeGO=='MA Plot'",
                 hr(style="border-color: black"),
                 htmlOutput("figureLegendMAPlot"),
                 hr(style="border-color: black"),
                 plotOutput("MAplot"),
                 hr(style="border-color: black")
               ),
               conditionalPanel(
                 condition="input.plotTypeGO=='Heatmap'",
                 ## Suppress error message that shows for <1sec when updating comparison or ontology (but not GO terms)
                 ## Error in matrix: 'data' must be of a vector type, was 'NULL'
                 tags$style(type="text/css",
                            ".shiny-output-error { visibility: hidden; }",
                            ".shiny-output-error:before { visibility: hidden; }"
                 ),
                 plotOutput("heatmapGO"),
                 hr(style="border-color: black")

               ),
               htmlOutput("tableHeaderTopGo"),
               tableOutput("topGoTable")

       )
    )
  )

))