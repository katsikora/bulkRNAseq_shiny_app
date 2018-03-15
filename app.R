## app.R ##
Rlib="/data/manke/repository/scripts/DNA_methylation/Rlibs.3.3.1"
require(shiny,lib.loc=Rlib)
require(shinydashboard,lib.loc=Rlib)
require(rhandsontable,lib.loc=Rlib)
ui <- function(request) {dashboardPage(
    dashboardHeader(title = "Dataset selection"),
    ## Sidebar content
    dashboardSidebar(

        textInput(inputId="group", label="Group", value = "", width = NULL, placeholder = NULL),
        textInput(inputId="projectid", label="ProjectID", value = "", width = NULL, placeholder = NULL),
        fileInput(inputId="countfile",label="CountFile",multiple=FALSE,accept=NULL,buttonLabel = "Browse...", placeholder = "No file selected"),
        actionButton("submitinput", "Retrieve dataset"),
        fileInput(inputId="file1", label="UserFile", multiple = FALSE, accept = NULL, width = NULL,buttonLabel = "Browse...", placeholder = "No file selected"),
        textOutput("fileDescription"),
        bookmarkButton()
        ),
        
    dashboardBody(
        h2("Bulk RNAseq analysis"),
        uiOutput("resultPanels")
               
    )

 )}

server <- function(input, output, session) {
    
    output$walkThrough<-renderUI(HTML("<ul><li>1.Complete the Group and ProjectID fields (currently only for Boehm department) OR provide path to feature counts count table to CountFile. Click on retrieve dataset. Your data will appear in the InputData tab.</li><li>2.Provide a UserFile with ensembl gene IDs to analyze.</li><li>3.If necessary, relabel sample IDs to desired plotting IDs in the corresponding interactive sampleInfo table column, in the InputData tab. Provide group information in the group column.Click on Run analysis. Your results will appear in the corresponding tabs.</li><li>The order of providing the information matters!</li></ul>"))
    output$FAQ<-renderText("Currently supported organisms are Homo sapiens, Mus musculus, Danio rerio, Drosophila melanogaster.\n Merging data from multiple datasets or batch effect removal are currenlty not supported.\n For questions, bug reports or feature requests, contact sikora@ie-freiburg.mpg.de.\n For reporting issues or pull requests on GitHub, go to https://github.com/katsikora/bulkRNAseq_shiny_app .")
    
    output$fileDescription<-renderText("UserFile: Please provide a one-column headerless txt file with ensemble gene IDs you would like to analyze.\n CountFile: Please provide an integer count table with ensembl gene IDs as rownames and sample IDs as colnames.")

 
##########################read/load processed data from a database##################
    #Rlib="/data/manke/repository/scripts/DNA_methylation/Rlibs.3.3.1"
    require("data.table",lib.loc=Rlib)
    require("limma",lib.loc=Rlib)
    require("edgeR",lib.loc=Rlib)
    require("car",lib.loc=Rlib)
    require(gplots,lib.loc=Rlib)
    require(RColorBrewer,lib.loc=Rlib)
    require(ggplot2,lib.loc=Rlib)
    require(reshape2,lib.loc=Rlib)
    require(biomaRt,lib.loc=Rlib)

    output$sessionInfo <- renderPrint({capture.output(sessionInfo())})

    ##barplots mean +/- stdev
    summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
        require(plyr,lib.loc=Rlib)

    # New version of length which can handle NA's: if na.rm==T, don't count them
        length2 <- function (x, na.rm=FALSE) {
            if (na.rm) sum(!is.na(x))
            else       length(x)
        }

        # This does the summary. For each group's data frame, return a vector with
        # N, mean, and sd
        datac <- ddply(data, groupvars, .drop=.drop,
          .fun = function(xx, col) {
            c(N    = length2(xx[[col]], na.rm=na.rm),
              mean = mean   (xx[[col]], na.rm=na.rm),
              sd   = sd     (xx[[col]], na.rm=na.rm)
            )
          },
          measurevar
        )

        # Rename the "mean" column    
        datac <- rename(datac, c("mean" = measurevar))

        datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

        # Confidence interval multiplier for standard error
        # Calculate t-statistic for confidence interval: 
        # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
        ciMult <- qt(conf.interval/2 + .5, datac$N-1)
        datac$ci <- datac$se * ciMult
        datac$ci.sd <- datac$sd * ciMult

        return(datac)
    }




    ##query database using group/projectid fields

    dbFile<-read.table("/data/manke/group/shiny/sikora/aux_files/RNAseq.DB.csv",header=TRUE,sep="\t",quote="",as.is=TRUE)

    observeEvent(input$submitinput, {
        
        if((input$group!="")&(input$projectid!="")){
            inGroup<-isolate(input$group)
            inProjectID<-isolate(input$projectid)
  
            datPath<-dbFile$PathToData[(dbFile$Group %in% inGroup) & (dbFile$ProjectID %in% inProjectID)]
            print(datPath)
            ##get genome and annotation version from a featureCounts log file
            logPath<-dbFile$PathToLog[(dbFile$Group %in% inGroup) & (dbFile$ProjectID %in% inProjectID)]
            org<-sub("_ensembl/.+.bam","",sub("/package.+/organisms/","",system(paste0('grep -m 1 \"/data/repository/organisms\" ',logPath),intern=TRUE)))
            output$organism<-renderText(paste("Detected genome is ",org))
            annot<-sub("/genes.gtf.+.bam","",sub("/package.+/data/repository/organisms/.+/ensembl/","",system(paste0('grep -m 1 \"/data/repository/organisms\" ',logPath),intern=TRUE)))}
        else if (!is.null(input$countfile)){
             countFile<-isolate(input$countfile)
             datPath<-countFile$datapath  
             print(datPath)
             annot<-"not available"}
        inDat<-fread(input=datPath,header=TRUE,sep="\t")
        ##render the head
        output$datHead<-renderTable(head(inDat),caption="Original unnormalized input data",caption.placement = getOption("xtable.caption.placement", "top"))

        output$annotation<-renderText(paste("Detected annotation is ",annot))
        ##or grep for organism from ensembl gene ids:
        emv<-c("ENSDARG"="drerio","ENSMUS"="mmusculus","ENSG"="hsapiens","FBgn"="dmelanogaster")
        ems<-emv[grep(gsub("[0-9].+","",inDat$V1[5]),names(emv))]
        ensembl.xx<-useMart(biomart="ensembl",dataset=paste0(ems,"_gene_ensembl"))#,host="aug2017.archive.ensembl.org"

###############initiate reactive table to collect sample information ###############

        values <- reactiveValues()
        DF<-data.frame(SampleID=colnames(inDat)[2:ncol(inDat)],PlottingID=colnames(inDat)[2:ncol(inDat)],Group=(rep("NA",(ncol(inDat)-1))),stringsAsFactors = F)
        observe({
          if (!is.null(input$hot)) {
            DF = hot_to_r(input$hot)
          } else {
            if (is.null(values[["DF"]]))
              DF <- DF
            else
              DF <- values[["DF"]]
          }
          values[["DF"]] <- DF
        })

        output$hot <- renderRHandsontable({
          DF <- values[["DF"]]
          if (!is.null(DF))
            rhandsontable(DF, stretchH = "all")
       })

     
##########################read user input file with gene ids#######################################
        observeEvent(input$file1,{
        inFile<-isolate(input$file1)
        inTab<-read.table(inFile$datapath, header = FALSE,sep="\t",quote="",as.is=TRUE)
        output$inTabHead<-renderTable(head(inTab),caption="Input gene IDs",caption.placement = getOption("xtable.caption.placement", "top"))
    
            observeEvent(input$savetable, {
                sampleInfo<-isolate(values[["DF"]])
                sampleInfo<-sampleInfo[!sampleInfo$Group %in% "NA",]
                rownames(sampleInfo)<-sampleInfo$PlottingID

                countdata<-as.data.frame(inDat[,match(sampleInfo$SampleID,colnames(inDat)),with=FALSE],stringsAsFactors=FALSE)
                colnames(countdata)<-sampleInfo$PlottingID
                rownames(countdata)<-inDat$V1
                output$countDatHead<-renderTable(head(countdata),include.rownames=TRUE,caption="Relabeled unnormalized input data",caption.placement = getOption("xtable.caption.placement", "top"),width=500)

                design<-model.matrix(~1+Group,data=sampleInfo)
                rownames(design)<-sampleInfo$Plotting_ID
                dge <- DGEList(counts=countdata)
                dge <- calcNormFactors(dge)
                v <- voom(dge,design,plot=FALSE)
                normCounts<-v$E
                selNormCounts<-normCounts[match(inTab$V1,rownames(normCounts)),]
                output$NormCounts<-renderTable(selNormCounts,include.rownames=TRUE,caption="Normalized Log2CPM",caption.placement = getOption("xtable.caption.placement", "top"))
                output$downloadNormCounts <- downloadHandler(filename = "NormalizedLog2CPM.xls",content=function(file) {write.table(selNormCounts,file,row.names = TRUE,sep="\t",quote=FALSE,dec = ",")})
                ##annotate with external gene symbol
                bmk<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),filters="ensembl_gene_id",values=rownames(selNormCounts),mart=ensembl.xx)
                plotdata<-as.data.frame(selNormCounts,stringsAsFactors=FALSE)
                plotdata$GeneID<-rownames(selNormCounts)
                plotdata$GeneSymbol<-bmk$external_gene_name[match(plotdata$GeneID,bmk$ensembl_gene_id)]
                output$heatmap<-renderPlot({heatmap.2(selNormCounts, scale="row", trace="none", dendrogram="column",col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),main=input$projectid, keysize=1,
                    margins = c(10,18),labRow=plotdata[,colnames(plotdata) %in% input$XlabelChoiceHeatmap]) })
                x_choice <- reactive({switch(input$XlabelChoiceBarplot,"GeneID"="GeneID","GeneSymbol"="GeneSymbol")})
                output$barplot<-renderPlot({
                    plotdataL<-melt(plotdata,id.vars=c("GeneID","GeneSymbol"),value.name="Log2CPM",variable.name="SampleID")
                    plotdataL$Log2CPM<-as.numeric(plotdataL$Log2CPM)
                    plotdataL$Group<-sampleInfo$Group[match(plotdataL$SampleID,sampleInfo$PlottingID)]
                    plotdata.SE<-summarySE(data=plotdataL,measurevar="Log2CPM",groupvars=c("Group","GeneID"))
                    plotdata.SE$GeneSymbol<-plotdata$GeneSymbol[match(plotdata.SE$GeneID,plotdata$GeneID)]
                
                    ggplot(data=plotdata.SE,aes(x=reorder(eval(as.name(x_choice())),Log2CPM),y=Log2CPM,fill=Group))+geom_bar(position=position_dodge(.9),colour="black",stat="identity")+geom_errorbar(position=position_dodge(.9),width=.25,aes(ymin=Log2CPM-sd, ymax=Log2CPM+sd))+scale_fill_manual(values=c("#FFFFFF","#CCCCCC"))+theme(text = element_text(size=16),axis.text = element_text(size=14),axis.text.x=element_text(angle=90,vjust=0),axis.title = element_text(size=14)) +xlab(x_choice())
                                          })

            
                        })#end of observe input$file1

                    })#end of observe input$savetable
        

        },ignoreInit=TRUE)#end of observe input$submitinput
   
    
###################################################################################   
    

       
    output$resultPanels<-renderUI({myTabs<-list(tabPanel(title="WalkThrough",
                                                      fluidPage(
                                                          box(title="Walkthrough",uiOutput("walkThrough")),
                                                          box(title="Miscellaneous information",textOutput("FAQ"))                                                          
                                                               )
                                                          ),
                                                  tabPanel(title="InputData",
                                                      fluidPage(
                                                          fluidRow(
                                                              tableOutput("datHead"),
                                                              tableOutput("countDatHead")
                                                                   ),
                                                          #fluidRow(
                                                              rHandsontableOutput("hot"),
                                                              actionButton(inputId="savetable",label="Run analysis"),
                                                            #      ),
                                                          fluidRow(
                                                              box(renderText("Please complete missing sample information. Plotting ID will be used to replace Sample ID on plots, Group will be used to construct a design table. Any samples with missing annotation will be removed from analysis.")),
                                                              tableOutput("inTabHead"),
                                                              box(textOutput("organism"),
                                                                  textOutput("annotation"))
                                                                  ) 
                                                              )
                                                          ),
                                                   tabPanel(title="NormCounts.Table",
                                                      fluidPage(
                                                          tableOutput("NormCounts"),
                                                          box(title="Method Description",renderText("Depth-normalized gene expression counts are listed for selected genes as log2 counts per million.")),
                                                          downloadButton("downloadNormCounts", label="Download Normalized Counts")
                                                               )
                                                          ),
                                                   tabPanel(title="NormCounts.Heatmap",
                                                      fluidPage(
                                                          box(plotOutput("heatmap")),
                                                          box(title = "Plot controls",selectInput("XlabelChoiceHeatmap", "Gene label",choices=c("GeneID","GeneSymbol"),selected="GeneID")),
                                                          box(title="Method Description",renderText("Genewise Log2 counts per million were mean-centered and scaled across the samples. Average linkage clustering on euclidean distances was performed."))
                                                               )
                                                          ),
                                                  tabPanel(title="NormCounts.Barplot",
                                                      fluidPage(
                                                          box(plotOutput("barplot")),
                                                          box(title = "Plot controls",selectInput("XlabelChoiceBarplot", "Gene label",choices=c("GeneID","GeneSymbol"),selected="GeneID")),
                                                          box(title="Method Description",renderText("Mean Log2 counts per million are plotted with standard deviation as error bars."))
                                                          
                                                               )
                                                          ),
                                                  tabPanel(title="sessionInfo",
                                                      fluidPage(
                                                          verbatimTextOutput("sessionInfo")                                                          
                                                               )
                                                          )

                                                 )

            do.call(tabsetPanel, myTabs)})


}

shinyApp(ui, server,enableBookmarking="url")
