#### README ####

### Requires

## 1) topTables:
# this is a list where elements 1 to N-1 are tables with the results from DESeq2 (produced by running Deseq2_Contrast_Rlog())
# and element N is a table with rlog values for all samples in the dataset


## DESeq2 results must contain the following columns with exactly these colnames:
# log2FoldChange
# pvalue
# padj
# gene_id (with ensembl id)
# gene_symbol
# heatmap.id (with gene id that will be displayed as rownames in heatmap. Typically, this is gene symbol when
#             available, ensembl id otherwise)
# <sample>.rlog (with regularised log of counts to display in heatmap)
# baseMean1 with mean of normalized counts for the two experimental groups + 1
# Arbitrary number of columns with 0/1 entries that indicate if a particular gene is in a GO term or not
# Column headers must be of format ontology_GO:XXXXX, e.g. MF_GO:0043236


# The individual tables must have names of format condition1.condition2.anything

## Correct table format can be produced by running Deseq2_Contrast_Rlog() to create DE tables, then
## CreateTopTables(). 

# 2) table with experimental groups
# expGroups.txt contains the sample id exactly as it is in counts table and dds
# and the experimental group of each sample (column headers: id and group)
# Careful: If the column names in counts table start with a number, R will introduce an
# X at the beginning. This must then also be added to expGroups$id


# 3) A list containing 1 table with topGO results per comparison. The individual tables 
# must have names of format condition1.condition2.anything.
# The GO terms will be shown in this order (if the table is produced by DESeq_TopGo_Workflow.R, it will be 
# by increasing weight01.Fisher P-value

#########################################################################################################


topTables<-readRDS("data/topTables.rds")
topGoResults<-readRDS("data/topGoResults.rds")

# add -log10(pvalue) in new column
DESeqTables<-lapply(topTables[1:(length(topTables)-1)], function(x){mutate(x, pval_log=(-log10(pvalue)))})
allComparisons<-names(DESeqTables)
for (i in 1:length(allComparisons)){
  names(allComparisons)[i]<-paste(strsplit(allComparisons[i], "[.]")[[1]][1], "versus", strsplit(allComparisons[i], "[.]")[[1]][2], sep=" ")
}

rlog<-topTables[[length(topTables)]]


temp_groups<-read.table("data/expGroups.txt", header=TRUE)
expGroups<-temp_groups$group
names(expGroups)<-temp_groups$id


allGroups<-rep(NA, 2*length(allComparisons))
j<-1
for (i in 1:length(allComparisons)){
  allGroups[j]<-strsplit(allComparisons[i], "[.]")[[1]][1]
  j<-j+1
  allGroups[j]<-strsplit(allComparisons[i], "[.]")[[1]][2]
  j<-j+1
}
allGroups<-unique(allGroups)
allGroups<-allGroups[order(allGroups)]




#########################################################################################################


library(dplyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(gplots)



shinyServer(function(input, output, session) {

#### DE analysis ######################################################################################  
    
  output$figureLegend<-renderUI({
    HTML("<b>Volcano plot:</b>", 
         "<br/>", 
         "<ul>
         <li>Click on points to display additional information below the plot (will work only if there are fewer than 10 nearby points)</li>
         <li>Positive log2(","<abbr title='Fold-change'> FC</abbr>",") = Higher expression in first condition</li>
         </ul>")
  })

  thresholdVolcano<-reactive({as.numeric(input$threshold)})
  
  ## Volcano plot
  output$p <- renderPlot({volcanoPlot(dataframe=DESeqTables[[allComparisons[input$comparison]]], 
                              threshold=thresholdVolcano(), 
                              gene=input$gene, 
                              title=input$comparison)})

  output$click_info<-renderTable({
    req(input$plot_click)
    # for some reason, the following will not work properly in browser (table is immediately updated to contain nothing)
   # temp<-nearPoints(DESeqTables[[allComparisons[input$comparison]]], input$plot_click, maxpoints = 10)
    temp<-nearPoints(DESeqTables[[allComparisons[input$comparison]]], input$plot_click)
    temp<-dplyr::select(temp, gene_id, gene_symbol, baseMean1, log2FoldChange, padj)
    colnames(temp)<-c("Ensembl ID", "Gene symbol", "Mean counts", "log2(FC)",
                      "Adjusted P")
    if(nrow(temp)==0 | nrow(temp)>10){
      NULL
    } else {
      temp
     }
    }, 
  include.rownames=FALSE)
  
  output$downloadPDF<-downloadHandler(
    filename=function() {
	paste(input$comparison,"pdf",  sep='.')
	},
    content=function(file){
      ggsave(file, 
             volcanoPlot(DESeqTables[[allComparisons[input$comparison]]], 
                         gene=input$gene, 
                         title=input$comparison,
                         threshold=thresholdVolcano()), 
             width=10 , height=7)
    }
  )
  

   
  # Heatmap with top N genes ################################################
  
   output$heatmap<-renderPlot(makeHeatmap2(DESeqTables[[allComparisons[input$comparison]]], input$GeneNumber, expGroups),
                              width="auto", height="auto")
 
   output$downloadHeatmap<-downloadHandler(
     filename=function() {
	 paste(input$comparison, input$GeneNumber, "pdf",  sep='.')
	 },
     content=function(file){
       pdf(file)
       makeHeatmap2(DESeqTables[[allComparisons[input$comparison]]], input$GeneNumber, expGroups)
       dev.off()
     }
   )
   
#### heatmap with userdefined genes ################################################


## Create resetable user input fields 
   output$resetable_input <- renderUI({
     times<-input$reset_input
     div(id=letters[(times %% length(letters)) + 1],
         fluidRow(column(5,
                         checkboxGroupInput("SelectedComp", label="Samples:", 
                                            choices=allGroups, selected=NULL)),
                  column(7,
                         p(strong("Genes:")),
                         p("Paste at least 2 Ensembl IDs (one per line), for example:", style="font-size:8pt"),
                         HTML('<textarea id="idFile" rows="5" cols="30", style="font-size:6pt">ENSMUSG00000050936\nENSMUSG00000096108</textarea>'),
                         
                         # Make example text appear/disappear
                         HTML('<script>
                              function addEvents(id) {
                              var field = document.getElementById(id);
                              field.onfocus = function () {
                              if (this.value == "ENSMUSG00000050936\\nENSMUSG00000096108") {
                              this.value = "";
                              } 
                              };
                              field.onblur = function () {
                              if (this.value == "") {
                              this.value = "ENSMUSG00000050936\\nENSMUSG00000096108";
                              this.style.color="black"; 
                              } 
                              };
                              }
                              addEvents("idFile");
                              </script>'),
                         
                         p("You can look up these IDs in the excel tables or on the ", a("Ensembl website.", href="http://www.ensembl.org/index.html"), style="font-size:8pt")
                         ))
     )
   })
   
   
# Figure out what colour the submit button should have: Orange when it needs to be pressed, green otherwise   
   button<-reactiveValues(colour="salmon")
   
   observeEvent(input$SelectedComp, {button$colour<-"salmon"})
   observeEvent(input$idFile, {button$colour<-"salmon"})
   observeEvent(input$reset_input, {button$colour<-"salmon"})
   observeEvent(input$submitButton, {button$colour<-"limegreen"})
   
   observe({
     button$colour
     session$sendCustomMessage(type="myCallbackHandler", button$colour)
   })
   
   
   samples_genes<-eventReactive(input$submitButton, {
     
     req(input$idFile)
     temp<-unlist(strsplit(input$idFile, "\n"))
     
     return(list(input$SelectedComp, temp))
   })

   
   output$customHeatmap<-renderPlot(makeHeatmap2_selectedGenes(rlog, samples_genes()[[2]], expGroups, samples_genes()[[1]]))


    
    # output ensembl ids that could not be matched to DESeq2 tables
   title<-reactive({
     title<-"Unknown Ensembl IDs (do not match results table):"
     
     inFile<-input$idFile
     if(is.null(inFile) | length(unknownIds(samples_genes()[[2]], rlog))==0) {
       return(NULL)} else {
      return(title)
       }
     })
   
   output$title<-renderText(title())
   output$unknownIds<-renderText(unknownIds(samples_genes()[[2]], rlog))

   output$downloadCustomHeatmap<-downloadHandler(
     filename=function() {
	 paste(input$comparison, "Heatmap", "pdf",  sep='.')
	 },
     content=function(file){
       pdf(file)
       makeHeatmap2_selectedGenes(rlog, samples_genes()[[2]], expGroups, samples_genes()[[1]])
       dev.off()
     }
   )

# GO enrichment  ##################################################################################
 
gotermList<-reactive({
  req(input$ontology, input$comparisonGO)
  temp<-colnames(topTables[[allComparisons[input$comparisonGO]]])
  temp<-temp[grep(paste(input$ontology, "GO", sep="_"), temp)]
  temp<-sub(paste(input$ontology, "_", sep=""), "", temp)
#  temp<-temp[order(temp)] 
  return(temp)
})

   
output$goterms<-renderUI({
    selectInput("selectedSet", "Select a GO term for plotting", choices=gotermList())
})


# create dataframe that will be used to make plot
df<-reactive({
    req(input$selectedSet, input$ontology, input$comparisonGO) 
    allGenes<-topTables[[allComparisons[input$comparisonGO]]]
    goterm<-paste(input$ontology, input$selectedSet, sep="_")
    inSet<-allGenes[allGenes[,grep(goterm, colnames(topTables[[allComparisons[input$comparisonGO]]]))]==1,]
    list(allGenes, inSet)
})


# Create a clean version of table with all genes in GO term for download
cleanTable<-reactive({
  temp<-df()[[2]]
  temp<-dplyr::select(temp, -baseMean1)
  temp<-dplyr::select(temp, gene_id:padj)
  temp<-temp[, grep(".rlog", colnames(temp), invert=TRUE)]
  return(temp)
})

## table with topGO results for the selected GO term
output$topGoTable<-renderTable({
  temp<-topGoResults[[allComparisons[input$comparisonGO]]][topGoResults[[allComparisons[input$comparisonGO]]]$GO.ID==input$selectedSet,]
  temp<-temp[,grep("Rank", colnames(temp), invert=TRUE)]
  temp<-dplyr::select(temp, -Expected, -ontology)
  },
  include.rownames=FALSE
  )
  
# Header for the topGo table
output$tableHeaderTopGo<-renderUI({
  HTML("<b>TopGo Results:</b>", 
       "Annotated=total number of genes in term; Significant=Total number of ","<abbr title='differentially expressed'> DE</abbr>",
       "genes in term (adjusted P<0.05); Remaining columns: P-values for three different gene set enrichment tests performed by TopGo",
       "<br/>"
  )
})




### MA Plot
output$figureLegendMAPlot<-renderUI({
  HTML("<b>MA plot:</b>", 
       "Plot of the mean expression level vs log2 fold-change",
       "<br/>", 
       "<ul>
       <li>Grey=all expressed genes; Red/black=genes in selected GO term. The colour indicates if the adjusted P-value from the test for differential expression is below (red) or above (black) the selected significance threshold</li>
       <li>Positive log2(","<abbr title='Fold-change'> FC</abbr>",") = Higher expression in first condition</li>
       </ul>"
       )
})


thresholdMAPlot<-reactive({as.numeric(input$thresholdMAPlot)})
                      
# draw scatterplot
output$MAplot<-renderPlot(
     tryCatch(
       {
         draw_ggplot(df()[[1]], df()[[2]], threshold=thresholdMAPlot(), title=input$selectedSet)
       },
      error=function(e) {
        return(NULL)
      }
     )
    )

# download MA plot
output$downloadPDF_MAPlot<-downloadHandler(
  filename=function() {
  paste(input$comparisonGO,"MAPlot.pdf",  sep='.')
  },
  content=function(file){
    ggsave(file, 
           draw_ggplot(df()[[1]], df()[[2]], threshold=thresholdMAPlot(), title=input$selectedSet), 
           width=10 , height=7)
  }
)

# download data for genes in GO term
output$downloadGeneset<-downloadHandler(
  filename=function() { paste(input$comparisonGO,"-", input$selectedSet,".txt",  sep='') },
  content = function(file) {write.table(cleanTable(), file, row.names = FALSE, quote=FALSE, sep='\t')}
)



### Heatmap

# Draw heatmap for genes within GO term
output$heatmapGO<-
  tryCatch(
    {
      renderPlot(makeHeatmap2(df()[[2]], nrow(df()[[2]]), expGroups, input$selectedSet), width="auto", height="auto")
    },
    error=function(e) {
      return(NULL)
    }
  )


# download Heatmap
output$downloadPDF_GOHeatmap<-downloadHandler(
  filename=function() {
  paste(input$comparisonGO, "GOHeatmap.pdf",  sep='.')
  },
  content=function(file){
    pdf(file)
    makeHeatmap2(df()[[2]], nrow(df()[[2]]), expGroups, input$selectedSet)
    dev.off()
  }
)

# download data for genes in GO term (same as above but somehow conditional panels are messed up 
# if the same output is used 2x)
output$downloadGeneset2<-downloadHandler(
  filename=function() { paste(input$comparisonGO,"-", input$selectedSet,".txt",  sep='') },
  content = function(file) {write.table(cleanTable(), file, row.names = FALSE, quote=FALSE, sep='\t')}
)


})
##################################################################################
# Functions
##################################################################################

## function to make ggplot volcano plot
volcanoPlot<-function(dataframe, gene="", title=NULL, threshold){

g<-ggplot(data=dataframe, aes(x=log2FoldChange, 
                              y=pval_log, 
                              colour=as.factor(padj<threshold))) + geom_point(alpha=0.5, size=2)
g<-g + scale_color_manual(paste("Differentially expressed\n(adjusted P<", threshold, ")", sep=''), 
                          labels=c("yes", "no"), breaks = c("TRUE", "FALSE"), values=c("TRUE" = "maroon4", "FALSE" = "black"))
g<-g + ylab("-log10(p-value)") + xlab("log2(fold-change)")

# highlight gene (if applicable):
if (gene!=""){
  if(is.element(gene, dataframe$gene_id)){  # check if the ensembl id exists in table
    geneInfo<-filter(dataframe, gene_id==gene)
    g<-g + geom_point(aes(x=log2FoldChange, 
                      y=pval_log),
                      data=geneInfo,
                      size=5,
                      colour="green")
    g<-g + ggtitle(paste(title, "  (highlighted gene: ", 
                         geneInfo$heatmap.id, ")", sep=""))} else {
    g<-g + geom_text(aes(x=0,
                     y=min(pval_log)+(max(pval_log)-min(pval_log))/2,
                     label="ENSEMBL ID NOT FOUND"),
                     show.legend=FALSE,
                     hjust=0.5, vjust=0, size=8, colour="red")
    g<-g + ggtitle(title)} 
} else {
 g<-g + ggtitle(title) 
}
return(g)
}


### Function to create heatmap with top N genes (i.e. with lowest padj)
makeHeatmap2<-function(df, NGenes, expGroups, Goterm=NULL){
  select <- order(df$padj)[1:NGenes] 
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  y<-as.matrix(df[select, grep("rlog",colnames(df))])
  rownames(y)<-df$heatmap.id[select]
  sampleNames<-sub(".rlog", "", colnames(y))
  colnames(y)<-sampleNames
  
  expGroups2<-as.character(expGroups[colnames(y)])
  groups<-unique(expGroups2)
  if (is.null(Goterm)){
    title<-paste(groups[1], " versus ", groups[2], " (top ", NGenes, " genes)", sep='')
    } else {
    title<-paste(groups[1], " versus ", groups[2], " (genes in ", Goterm, ")", sep='')
    }
  
  expGroups2<-replace(expGroups2, which(expGroups2==unique(expGroups2)[1]), "snow3")
  expGroups2<-replace(expGroups2, which(expGroups2==unique(expGroups2)[2]), "slategrey")
  
  par(cex.main=0.8)
  heatmap.2(y, col = hmcol, trace="none", margin=c(10,6),
            labCol=sampleNames, key=TRUE, density.info="none", ColSideColors = expGroups2,
            cexRow = 1/log10(NGenes), offsetRow = 0, main=title,
			key.xlab="rlog normalized counts")
}



### Function to create heatmap with user-defined genes
makeHeatmap2_selectedGenes<-function(df, GeneIds=NULL, expGroups, samples){
  req(samples)
  if(!is.null(GeneIds)){
    select <- which(is.element(rownames(df), GeneIds))
    
	
	sampleColumns<-grep(paste(paste("^", samples, "$", sep=""), collapse="|"), expGroups, value=FALSE)
    
    if(length(select)!=0){
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    y<-as.matrix(df[select, sampleColumns])
    rownames(y)<-df$heatmap.id[select]   

   expGroups2<-as.character(expGroups[colnames(y)]) #this works because expGroups is a named vector where names are sample IDs
    # create a vector with one colour per experimental group:
   # take four greys (more are difficult to distinguish)
   # and if there are more groups, different yellows
   # Maximum number of experimental groups is currently 22
   colours<-gray.colors(4, start=0.1, end=0.9)   # 0=black, 1=white
   colours<-c(colours, brewer.pal(9, "YlOrBr"))
   colours<-c(colours, brewer.pal(9, "PRGn"))

   
   comparisons<-unique(expGroups2)
   for(i in 1:length(unique(expGroups2))){
      expGroups2<-replace(expGroups2, which(expGroups2==unique(expGroups2)[i]), colours[i])
   }
    
   
   title<-"Heatmap with user-specified genes for groups: "
   for (j in 1:length(comparisons)){
     title<-paste(title, comparisons[j], sep=" ")
   }

   
   par(cex.main=0.8)
    heatmap.2(y, col = hmcol, trace="none", margin=c(10, 6),
              labCol=colnames(y), key=TRUE, density.info="none",
              ColSideColors = expGroups2,
              cexRow = 1/log10(length(select)), offsetRow = 0, main=title,
			  key.xlab="rlog normalized counts")
    } else {
      plot.new()
    }
    
  } else {
    plot.new()
  }
    
}


# return ids which cannot be matched to DE results table

unknownIds<-function(GeneIds, df){
  unknown<-GeneIds[!is.element(GeneIds, rownames(df))]
  return(unknown)
}



## Create MA plot with ggplot2
## Function can plot data from two dataframes
draw_ggplot<-function(fullData, geneset, threshold, title) { 
  
  # add colour variable to genes within gene set
  geneset$sign<-'no'
  geneset$sign[geneset$padj<=threshold]<-'yes'
  geneset$sign<-as.factor(geneset$sign)
  colours<-c(no='black', yes='red')
  # find out the range of the y-axis and the ticks to be used:
  maxFC<-ceiling(max(abs(fullData$log2FoldChange)))
  
  p<-ggplot(fullData, aes(x=baseMean1, y=log2FoldChange))
  p<-p + geom_point(colour="darkgrey", alpha=0.3, size=3)
  p<-p + ylab("log2 fold-change")
  p<-p + xlab("mean number of reads (log10 scale)")
  p<-p + ylim(-maxFC, maxFC)
  p<-p + geom_hline(aes(yintercept=0), linetype='dashed')
  p<-p + geom_point(aes(x=baseMean,
                        y=log2FoldChange,
                        colour=sign),
                    alpha=0.5,
                    size=3,
                    data=geneset)
  p<-p + scale_color_manual(values=colours,
                            name=paste("P-adj<", threshold, sep=""))
  p<-p + ggtitle(title)
  
  p<-p + scale_x_continuous(trans=log10_trans(),
                             breaks=trans_breaks("log10", function(x) 10^x),
                             labels=trans_format("log10", math_format(10^.x)))
  return(p)
}  



  
  