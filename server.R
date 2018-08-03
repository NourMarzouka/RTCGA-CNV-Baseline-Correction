library(shiny)

# options(warn=-1)

#for updating RTCGA.CNV package
options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")
library("RTCGA.CNV.20160128")

#the limit of the file upload (50MB)
options(shiny.maxRequestSize=50*1024^2)


shinyServer(function(input, output, session) {
  
  #create reactive list
  reactive           <- reactiveValues()
  reactive$data_type <- ""
  reactive$do.plot   <- FALSE
  
  #create new environment 
  senv               <- new.env()
  senv$object        <- NULL
  senv$is_autocorrect <- FALSE
  
  #New CopyNumber Functions
  #region_file fileinput path from Shiny
  #regions_colnames colnames vector from shiny input
  #sample_list file input path from Shiny input
  #sample_colnames colnames vector from shiny inputs
  ReadData <- function(session,regions_file, regions_colnames, sample_list,sample_colnames) {
    
    regions           <- read.csv(regions_file, stringsAsFactors = FALSE,
                                  header=input$header, sep=input$sep)
    #regions          <- replaceChr(read.csv(regions_file, stringsAsFactors = FALSE))
    regions           <- regions[,regions_colnames]
    colnames(regions) <- c("Sample", 
                           "Chromosome", 
                           "bp.Start", 
                           "bp.End", 
                           "Num.of.Markers", 
                           "Mean")
    regions           <- replaceChr(regions)
    
    #save the number of samples (plots) 
    #this number will be used a lot in this code
    senv$max_plots    <- length(unique(regions$Sample))
    
    #update the slider in the plotraw tab based on the number of samples
    updateSliderInput(session, "NumberSampleSlider", 
                      max = senv$max_plots, min = 1,step = 1 )
    
    #If the user uploaded sample list
    if(!missing(sample_list)){
      SL           <- read.csv(sample_list, stringsAsFactors = FALSE,
                               header=input$headersamp, sep=input$sepsamp)
      SL           <- SL[,sample_colnames]
      colnames(SL) <- c("Number", "Sample", "Comment")
    }
    
    object <- list(regions      = regions,
                   regions_save = regions,
                   regions_auto = regions)
    class(object) <- "CopyNumber450kCancer_data"
    
    #If the user did not upload sample list
    if (missing(sample_list)) {
      Number  <- c(1:length(unique(regions$Sample)))
      Sample  <- unique(regions$Sample)
      Comment <- c(rep(" ",length(unique(regions$Sample))))
      SL      <- data.frame(Number, Sample, Comment, stringsAsFactors = F)
    }
    
    #store the sample list in the object
    object$SL <- SL
    
    # to store the autocorrection modification in it
    object$mod_auto <- data.frame("Sample" = object$SL$Sample,
                                  "Density_Maximum_Peak"= 0,
                                  "Shifting"=0,
                                  "Auto_Corrected"="No",
                                  stringsAsFactors =FALSE)
    
    # to store the manual modification in it
    object$mod <- data.frame("Sample" = object$SL$Sample, 
                             "Shifting"=0,
                             "Using_slider"="No", 
                             stringsAsFactors =FALSE)
    
    object
  }
  
  #function to calculate the area under the curve to be used in the QC function 
  auc <- function (x, y, thresh = NULL, dens = 100, sort.x = TRUE) {
    x <- x[!is.na(x)]
    y <- y[!is.na(x)]
    x <- x[!is.na(y)]
    y <- y[!is.na(y)]
    if (sort.x) {
      ord <- order(x)
      x <- x[ord]
      y <- y[ord]
    }
    idx = 2:length(x)
    x <- as.vector(apply(cbind(x[idx - 1], x[idx]),
                         1, 
                         function(x) seq(x[1], x[2], length.out = dens)))
    y <- as.vector(apply(cbind(y[idx - 1], y[idx]),
                         1, 
                         function(x) seq(x[1], x[2], length.out = dens)))
    if (!is.null(thresh)) {
      y.0    <- y <= thresh
      y[y.0] <- thresh
    }
    idx = 2:length(x)
    integral <- as.double((x[idx] - x[idx - 1]) %*% 
                            (y[idx] + y[idx - 1]))/2
    integral
  }  
  
  # QC function
  #the usage:
  #the user selects cutoff and markers
  QC.function <- function(object, cutoff=0.1, markers=0, ...) {
    
    #copy the sample names and the comments to prepare the QC file 
    QC <- object$SL[, c("Sample", "Comment")]
    QC[, 3:10] <- 0
    colnames(QC) <- c("Sample", "Comment", "peak.sharpness",
                      "number.of.segments", "IQR", "SD", "MAPD",
                      "Area.between.cutoffs",
                      "Density.height.at.lower.cutoff",
                      "Density.height.at.upper.cutoff")
    
    
    for (i in 1:length(object$SL[, "Sample"])) {
      
      sam <- object$regions[which(object$regions$Sample %in% as.character(object$SL[i, "Sample"])), ]
      
      forDen <- sam[which(sam$Chromosome != "chrX" & sam$Chromosome != "chrY"& sam$Chromosome!="chr23"& sam$Chromosome!="chr24"),
                    c("Num.of.Markers", "Mean") ]
      d      <- density(forDen$Mean,
                        weights = forDen$Num.of.Markers / sum(forDen$Num.of.Markers),
                        na.rm = TRUE, kernel = "gaussian",
                        adjust = 0.15, n = 1024, from = -1, to = 1)
      
      max.peak.value    <- d$x[which.max(d$y)]
      point             <- which(d$x == (max.peak.value))
      QC.peak.sharpness <- ( (d$y[point + 20] + d$y[point - 20]) / 2 ) /
        d$y[which(d$x == max.peak.value)]
      
      QC[i, "peak.sharpness"]                 <- as.numeric(QC.peak.sharpness)
      QC[i, "number.of.segments"]              <- length(sam[, 1])
      QC[i, "IQR"]                            <- IQR(forDen[, "Mean"], na.rm = TRUE, type = 7)
      QC[i, "SD"]                             <- sd(forDen[, "Mean"], na.rm = TRUE)
      QC[i, "MAPD"]                           <- median(abs(diff(forDen[, "Mean"], na.rm = TRUE)), na.rm = TRUE)
      QC[i, "Area.between.cutoffs"]           <- (1-(  auc(d$x[d$x>(cutoff)],  d$y[d$x>(cutoff)])
                                                       + auc(d$x[d$x<(-cutoff)], d$y[d$x<(-cutoff)])
      ))
      QC[i, "Density.height.at.lower.cutoff"] <- d$y[which.min(abs(d$x-(-cutoff)))]
      QC[i, "Density.height.at.upper.cutoff"] <- d$y[which.min(abs(d$x-cutoff))]
    }
    
    object$QC <- QC
    invisible(object)
  }
  
  #Function to change the chromosomes names by adding "chr"
  replaceChr <- function(object){
    object$Chromosome <-
      ifelse( grepl('^chr', object$Chromosome),    # If the name start with chr
              object$Chromosome,                   # It's the right name
              paste0('chr', object$Chromosome))    # Else prepend 'chr'
    
    invisible(object)
  }
  
  #Function to store the TCGA dataset in the regions of the object
  #Also translates the chromosome names
  tcgaToObject <- function(session,tcganame) {
    reactive$data_type <- "tcga"
    regions            <- replaceChr(get(tcganame))
    
    colnames(regions) <- c("Sample","Chromosome","bp.Start","bp.End","Num.of.Markers","Mean")	
    
    #find the number of samples
    senv$max_plots    <- length(unique(regions$Sample))
    
    #update the slider in plot raw data tab
    updateSliderInput(session, "NumberSampleSlider", max = senv$max_plots, min = 1, step = 1)
    
    
    object        <- list(regions      = regions,
                          regions_save = regions,
                          regions_auto = regions)
    
    class(object) <- "CopyNumber450kCancer_data"
    
    #prepare the sample list (SL)
    Number        <- c(1:length(unique(regions$Sample)))
    Sample        <- unique(regions$Sample)
    Comment       <- c(rep(" ", length(unique(regions$Sample))))
    SL            <- data.frame(Number, Sample, Comment, stringsAsFactors = F)
    
    #store the sample list in the object
    object$SL <- SL
    
    # copy to store the autocorrection modification in it
    object$mod_auto <- data.frame("Sample"            = object$SL$Sample,
                                  "Density_Maximum_Peak" = 0,
                                  "Shifting"          = 0,
                                  "Auto_Corrected"    = "No",
                                  stringsAsFactors    = FALSE)
    
    # copy to store the manual modification in it
    object$mod <- data.frame("Sample"         = object$SL$Sample,
                             "Shifting"       = 0,
                             "Using_slider"   = "No",
                             stringsAsFactors = FALSE)
    
    object
  }
  
  #---this function to plot the ticks on the right side
  tick.tick <- function (nx = 2, ny = 2, tick.ratio = 0.5) {
    ax <- function(w, n, tick.ratio) {
      range <- par("usr")[if (w == "x") 1:2 else 3:4]
      tick.pos <- if (w == "x") {
        par("xaxp")
      } else {
        par("yaxp")
      }
      distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
      possible.minors        <- tick.pos[1] - (0:100) * distance.between.minor
      low.minor              <- min(possible.minors[possible.minors >= range[1]])
      if (is.na(low.minor))
        low.minor <- tick.pos[1]
      
      possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
      hi.minor        <- max(possible.minors[possible.minors <= range[2]])
      if (is.na(hi.minor))
        hi.minor <- tick.pos[2]
      axis(if (w == "x") 1 else 4,
           seq(low.minor, hi.minor, by = distance.between.minor),
           labels = FALSE, tcl = par("tcl") * tick.ratio)
    }
    if (nx > 1)
      ax("x", nx, tick.ratio = tick.ratio)
    if (ny > 1)
      ax("y", ny, tick.ratio = tick.ratio)
    invisible()
  }
  
  #To plot the segmentation data
  #Chromosome should be in this format: "chr1"
  plotRegions <- function(object, cutoff=0.1,markers=20, ...) {
    sample_segments <- object
    
    sample_segments <- sample_segments[which(sample_segments$Num.of.Markers>markers),]
    
    segment_values <- as.numeric(sample_segments[,"Mean"])
    segment_colors <- rep("black", nrow(sample_segments))
    
    segment_colors[as.numeric(segment_values) >= cutoff] <- "green"
    segment_colors[as.numeric(segment_values) <= -cutoff] <- "red"
    
    # Plotting the whole genome
    chromosomes <- unique(sample_segments[, "Chromosome"])
    site_per_chr <- cumsum(c(0, sapply(chromosomes, function(chr) max(as.numeric(sample_segments[sample_segments[,"Chromosome"] == chr, "bp.End"])))))
    # 1 instead of "chr1" #as.numeric(gsub("\\D", "", x))
    offset <- site_per_chr - min(as.numeric(sample_segments[sample_segments[, "Chromosome"] == "chr1", "bp.Start"])) 
    start <- 0
    end <- as.numeric(max(site_per_chr))
    x_axis_type <- "n"
    
    #put the Y axis limits
    yMin <- (-1) #min(c(-1, as.numeric(sample_segments[significant_segments, "Mean"])))
    yMax <- 1    #max(c(1, as.numeric(sample_segments[significant_segments, "Mean"])))
    
    myPlot <- plot(range(start, end), range(yMin, yMax), type = "n",axes=FALSE, xaxt = x_axis_type, xlab="",ylab="", ...)
    
    #to get the poistion of the labels (chromosomes names) on X axis
    xlabs <- sapply(2:length(site_per_chr), function(j) {
      ((site_per_chr[j] - site_per_chr[j - 1])/2) + site_per_chr[j - 1]
    })
    
    axis(1, at = xlabs, labels = chromosomes, lty = 0, las = 2, ...)
    axis(4)
    tick.tick(nx=0,ny=2, tick.ratio=1.6)
    tick.tick(nx=0,ny=10, tick.ratio=0.6)
    mtext("L-value", side = 4, line = 2, cex = par("cex.lab"))
    box()
    
    #lines to seperate the choromosomes
    abline(v = site_per_chr, lty = 3)
    #lines for the cutoffs
    abline(h = c(0,-cutoff,cutoff), lty = 3)
    
    #Draw the segments
    lapply(1:length(chromosomes), function(i) {
      used_segments <- sample_segments[, "Chromosome"] == chromosomes[i]
      colors        <- segment_colors[used_segments]
      starts        <- as.numeric(sample_segments[used_segments, "bp.Start"]) + offset[i]
      ends          <- as.numeric(sample_segments[used_segments, "bp.End"]) + offset[i]
      y             <- as.numeric(sample_segments[used_segments, "Mean"])
      graphics::segments(starts, y, ends, y, col = colors, lwd = 2, lty = 1)
    })
    
    myPlot
  }
  
  #NEW PLOT MODIFICATION (for correction tab)
  Plot.Manual <- function(object, select=1, cutoff=0.1,markers=20, comments =FALSE, slider_value=0,...){

    # get the sample name
    name         <- object$SL[select,"Sample"]
    
    # get the sample segments
    #sam          <- object$regions[which(object$regions$Sample %in% as.character(name)),]
    sam <- object$regions_auto[which(object$regions_auto$Sample %in% as.character(name)),]
    original.sam <- sam  	
    
    sam <- sam[which(sam$Num.of.Markers>markers),]
    
    #modify by the value from the slider
    sam$Mean <- sam$Mean + slider_value
    
    #to prepare the spaces for the plots
    par(mfrow=c(1,2), mar=c(0,0,2,0), oma=c(0,0,0,4))
    layout(matrix(c(1,2), 1, 2, byrow = TRUE),
           widths=c(3,21),
           heights=c(10),
           TRUE)
    
    #calculate the density
    forDen <- sam[which(sam$Chromosome!="chrX" & sam$Chromosome!="chrY" & sam$Chromosome!="chr23"& sam$Chromosome!="chr24"), c("Num.of.Markers", "Mean")]
    d      <- density(forDen$Mean,
                      weights = forDen$Num.of.Markers/sum(forDen$Num.of.Markers),
                      na.rm   = TRUE,
                      kernel  = "gaussian",
                      adjust  = 0.15,
                      n       = 512)
    #plot the density
    plot(d$y, d$x,
         type = 'l',
         ylim = c(-1,1),
         xlim = rev(range(d$y)),
         ylab = "",
         xlab = "",
         axes = FALSE)
    abline(h = c(0,-cutoff,cutoff), lty = 3)
    box()
    legend("bottomleft", legend="Density", cex=1)
    
    #plot the segments
    title=c(paste("Sample ",
                  select,
                  ":",
                  object$SL[which(object$SL$Sample %in% as.character(name)),"Sample"]), ...)
    
    plotRegions(sam,
                cutoff=cutoff,
                markers=markers,
                main=title)
    
    if (comments){
      legend("topleft", legend=paste("Comment:",object$SL[which(object$SL$Sample %in% as.character(name)),"Comment"]),cex=0.75)
    }
    
    #save the modifications
    object$mod[which(object$mod$Sample %in% as.character(name)), 2:3] <- c(0, "No")
    if (slider_value!=0) {
      object$mod[which(object$mod$Sample %in% as.character(name)), 2:3] <- c(slider_value, "Yes")
    }
    
    #save the new Means in object$regions
    object$regions[which(object$regions$Sample %in% as.character(name)), "Mean"] <- original.sam$Mean + slider_value
    #object$regions[which(object$regions$Sample %in% as.character(name)), "Mean"] <- sam$Mean
    
    invisible(object)
  }
  
  #plot the raw plots (for plot raw data tab)
  PlotRawData <- function(object, original.data= F, select=1, plots=TRUE,cutoff=0.1,markers=20, comments =FALSE,...){

    # get the sample name
    name <- object$SL[select,"Sample"]
    
    #get the sample segments
    if(original.data){
      sam  <- object$regions_save[which(object$regions_save$Sample %in% as.character(name)),]   
    } else {
      sam  <- object$regions[which(object$regions$Sample %in% as.character(name)),]   
    }
    
    #keep only the regions with more than minimum number of markers
    sam <- sam[which(sam$Num.of.Markers>markers),]
    
    #to prepare the spaces for the plots
    par(mfrow=c(1,2), mar=c(0,0,2,0), oma=c(0,0,0,4))
    layout(matrix(c(1,2),1,2,byrow=TRUE), widths=c(3,21), heights=c(10), TRUE)
    
    #calculate the density
    forDen <- sam[which(sam$Chromosome!="chrX" & 
                          sam$Chromosome!="chrY"& 
                          sam$Chromosome!="chr23"& 
                          sam$Chromosome!="chr24"), c("Num.of.Markers","Mean")]
  
    if (sum(forDen$Num.of.Markers) == 0) {
      return(NULL)
    }
    
    d      <- density(forDen$Mean ,
                      weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),
                      na.rm=TRUE, 
                      kernel = "gaussian", 
                      adjust=0.15, 
                      n=512)
    
    #plot the density
    plot(d$y,d$x,ylim=c(-1,1),type='l',ylab="",xlab="",axes=FALSE,xlim=rev(range(d$y)))
    abline(h = c(0,-cutoff,cutoff), lty = 3)
    box()
    legend("bottomleft", legend="Density",cex=1)
    
    #plot the regions
    plotRegions(sam,cutoff=cutoff,markers=markers,main=c(paste("Sample ",select,":",object$SL[which(object$SL$Sample %in% as.character(name)),"Sample"]), ...))
    if (comments){
      legend("topleft", legend=paste("Comment:",object$SL[which(object$SL$Sample %in% as.character(name)),"Comment"]),cex=0.75)
    }
  }
  
  #The Auto correction based on the highest density peak -------------------
  #this function correct the mean of the regions based on the highest peak and plot them
  AutoCorrectPeak <- function(object){
    
    #correction
    for (i in 1:length(object$SL[,"Sample"])){ 
      sam            <- object$regions_auto[which(object$regions_auto$Sample %in% as.character(object$SL[i,"Sample"])),]
      forDen         <- sam[which(sam$Chromosome!="chrX" & sam$Chromosome!="chrY"& sam$Chromosome!="chr23"& sam$Chromosome!="chr24"),c("Num.of.Markers","Mean")]
      sam.original   <- sam
      d              <- density(forDen$Mean,
                                weights=forDen$Num.of.Markers/sum(forDen$Num.of.Markers),
                                na.rm=TRUE,kernel = "gaussian",adjust=0.15,n=512)
      max.peak.value <- d$x[which.max(d$y)]
      sam$Mean       <- sam$Mean-max.peak.value
      
      #store the modifications
      object$regions_auto[which(object$regions_auto$Sample %in% as.character(object$SL[i,"Sample"])),"Mean"] <- sam$Mean 
      object$mod_auto[which(object$mod_auto$Sample %in% as.character(object$SL[i,"Sample"])),2:4] <- c(round(max.peak.value, 3),round(-max.peak.value, 3),"Auto")
    }
    
    #added for the new version
    object$regions <- object$regions_auto
    
    object
  }
  
  #############################################
  ######### Observers for pushbuttons #########
  ##### Watch buttons and update views etc ####
  #############################################
  observeEvent(input$RegionsActionButtonGo2Sample, {
    updateNavbarPage(session, "baseCN", selected = "Upload sample list")
  })
  
  observeEvent(input$PlotActionButtonGo2Autocorrect, {
    
    senv$object         <- AutoCorrectPeak(senv$object) 
    senv$is_autocorrect <- TRUE
    
    updateNavbarPage(session, "baseCN", selected = "Baseline Correction & QC")
    
    
  })
  
  observeEvent(input$SampleActionButton1, {
    updateNavbarPage(session, "baseCN", selected = "Plot raw data")
  })
  
  observeEvent(input$LoadFromCancerAtlasData, {
    updateNavbarPage(session, "baseCN", selected = "Load TCGA data")
  })
  
  observeEvent(input$UpLoadData, {
    updateNavbarPage(session, "baseCN", selected = "Upload Your Own Data")
  })
  
  observeEvent(input$RegionsActionButtonGo2PlotRaw, {
    updateNavbarPage(session, "baseCN", selected = "Plot raw data")
  })
  
  #check boxes: select and unselect samples 
  observeEvent(input$SelectAllSamples, {
    if(input$SelectAllSamples){
      
      for (i in 1:senv$max_plots) {
        updateCheckboxInput(session, paste("PlotRawSamplecheckbox", i, sep=""), value =TRUE)
        updateCheckboxInput(session, "SelectPageSamples", value =FALSE)
      }} else {
        for (i in 1:senv$max_plots) { 
          updateCheckboxInput(session, paste("PlotRawSamplecheckbox", i, sep=""), value=FALSE)
        }
      }
  })
  
  observeEvent(input$SelectPageSamples, {
    if (input$SelectPageSamples){
      for (i in 1:input$NumberSampleSlider) {
        updateCheckboxInput(session, paste("PlotRawSamplecheckbox", i, sep=""), value =TRUE)
        updateCheckboxInput(session, "SelectAllSamples", value=FALSE)
      }} else {
        for (i in 1:senv$max_plots) { 
          updateCheckboxInput(session, paste("PlotRawSamplecheckbox", i, sep=""), value=FALSE)
        }
      }})
  
  #button for load example samples
  observeEvent(input$LoadSampleData, {
    reactive$regions <- "DATA/regions.csv"
    reactive$sample_list <- "DATA/sample_list.csv"
    reactive$data_type <- "sampledata"
    updateNavbarPage(session, "baseCN", selected = "Plot raw data")
  })
  
  #buttons for TCGA tab
  observeEvent(input$SampleActionButton2, {
    senv$object <- tcgaToObject(session,input$cancerdata)
    updateNavbarPage(session, "baseCN", selected = "Plot raw data")
  })
  
  observeEvent(input$load.this.TCGA, {
    senv$object <- tcgaToObject(session,input$cancerdata)
    updateNavbarPage(session, "baseCN", selected = "Plot raw data")
  })
  
  ##############################################
  ######### Rendering UI for the tabs ##########
  ##############################################
  ##User input: segmentation file
  output$csvtableRegions <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    regions <- input$file1
    if (is.null(regions))
      return(NULL)
    reactive$data_type <- "RegionsUploaded"
    RegionInput=read.csv(regions$datapath, header=input$header, sep=input$sep,
                         quote=input$quote,nrows=10)
    RegionVariables=names(RegionInput)
    updateSelectInput(session, "RegionSample",     choices = RegionVariables, selected=grep("sample|name|sample([:blank:]|[:punct:])name|code|id",RegionVariables , value=TRUE,ignore.case =TRUE))
    updateSelectInput(session, "RegionChromosome", choices = RegionVariables, selected=grep("chromosome|chr|chromo",RegionVariables , value=TRUE,ignore.case =TRUE))
    updateSelectInput(session, "Regionbpstart",    choices = RegionVariables, selected=grep("bp([:blank:]|[:punct:])start|start|chromStart|from,",RegionVariables , value=TRUE,ignore.case =TRUE))
    updateSelectInput(session, "Regionbpend",      choices = RegionVariables, selected=grep("bp([:blank:]|[:punct:])End|ends|end|chromoEnd|to",RegionVariables , value=TRUE,ignore.case =TRUE))
    updateSelectInput(session, "RegionNumMark",    choices = RegionVariables, selected=grep("Num|num([:blank:]|[:punct:])of([:blank:]|[:punct:])Markers|markers|probes|number([:blank:]|[:punct:])of([:blank:]|[:punct:])probes|number([:blank:]|[:punct:])of([:blank:]|[:punct:])markers",RegionVariables , value=TRUE,ignore.case =TRUE))
    updateSelectInput(session, "RegionMean",       choices = RegionVariables, selected=grep("Mean|log|value|meanlog|L-value",RegionVariables , value=TRUE,ignore.case =TRUE))
    
    RegionInput
  })
  
  ##User input: sample list file
  output$csvtableSample <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    sample <- input$file2
    if (is.null(sample))
      return(NULL)
    reactive$data_type <- "RegionsAndSampleUploaded"
    SampleInput=read.csv(sample$datapath, header=input$headersamp, sep=input$sepsamp,
                         quote=input$quotesamp,nrows=10)
    SampleVariables=names(SampleInput)
    updateSelectInput(session, "SampleNumber",  choices = SampleVariables, selected=grep("number|sample([:blank:]|[:punct:])number", SampleVariables, value=TRUE,ignore.case =TRUE))
    updateSelectInput(session, "SampleSample",  choices = SampleVariables, selected=grep("sample|name|sample([:blank:]|[:punct:])name|code|id",SampleVariables , value=TRUE,ignore.case =TRUE))
    updateSelectInput(session, "SampleComment", choices = SampleVariables, selected=grep("comment", SampleVariables, value=TRUE,ignore.case =TRUE))
    SampleInput
  })
  
  #Select output from RTCGA.CNV package
  #Filling the select input from the results Item collumn.
  output$tcga <- renderUI({
    data.cancer <- data(package = "RTCGA.CNV.20160128")$results[,"Item"]
    selectInput("cancerdata", "Select the cancer type:", c(data.cancer))
  })
  
  #Output how many samples there is in the dataset.
  output$tcgaSamplenumber <- renderUI({
    if (is.null(input$cancerdata))
      return(NULL)
    if (input$cancerdata!=""){
      tags$p(paste("Number of Sample for this cancer type: ",nrow(unique(get(input$cancerdata)[1]))))
    }
  })
  
  #Render head table depending on the select input.
  output$tableTCGA <- renderTable({
    if (is.null(input$cancerdata))
      return(NULL)
    if (input$cancerdata!=""){
      head(get(input$cancerdata),10)}
  })
  
  # Render UI with raw unmodified plots of the samples.
  # Lets the user select input parameters and number of plots to autocorrect.
  output$plotraw <- renderUI({
    
    switch(reactive$data_type,
           tcga        = {reactive$do.plot <- TRUE
           senv$is_autocorrect <- FALSE
           },
           
           sampledata  = {senv$object      <- ReadData(session,
                                                       reactive$regions,
                                                       c("Sample","Chromosome",
                                                         "bp.Start","bp.End","Num.of.Markers","Mean"),
                                                       reactive$sample_list,
                                                       c("Number","Sample","Comment"))
           reactive$do.plot    <- TRUE
           senv$is_autocorrect <- FALSE
           },
           
           RegionsUploaded = {regions <- input$file1
           region_colnames  <- c(input$RegionSample,
                                 input$RegionChromosome,
                                 input$Regionbpstart,
                                 input$Regionbpend,
                                 input$RegionNumMark,
                                 input$RegionMean)
           senv$object      <- ReadData(session,regions$datapath,region_colnames)
           reactive$do.plot <- TRUE
           senv$is_autocorrect <- FALSE
           },
           
           RegionsAndSampleUploaded = {regions    <- input$file1
                             region_colnames      <- c(input$RegionSample,
                                                       input$RegionChromosome,
                                                       input$Regionbpstart,
                                                       input$Regionbpend,
                                                       input$RegionNumMark,
                                                       input$RegionMean)
                             sample_list          <- input$file2
                             sample_list_colnames <- c(input$SampleNumber, 
                                                       input$SampleSample, 
                                                       input$SampleComment)
                             
                             senv$object          <- ReadData(session,
                                                              regions$datapath,
                                                              region_colnames, 
                                                              sample_list$datapath, 
                                                              sample_list_colnames)
                             reactive$do.plot     <- TRUE
                             senv$is_autocorrect  <- FALSE
           },
           
           stop(return(tags$b("Please load some test data, upload your own data or select a dataset from TCGA tab")))
    )
    
    if(reactive$do.plot){
      for (i in 1:senv$max_plots) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i         <- i
          plotcheckbox <- paste("plotcheckbox", my_i, sep="")
          
          output[[plotcheckbox]] <- renderUI({checkboxInput(paste("PlotRawSamplecheckbox", my_i, sep=""),
                                                            paste("Select Sample", my_i, sep="") , 
                                                            FALSE)
          })
        })
      }
      
      # ##################
      # for (i in 1:input$NumberSampleSlider) {
      #   # Need local so that each item gets its own number. Without it, the value
      #   # of i in the renderPlot() will be the same across all instances, because
      #   # of when the expression is evaluated.
      #   local({
      #     my_i         <- i
      #     plotname     <- paste("Sample", my_i, sep="")
      #     
      #     output[[plotname]] <- renderPlot({ PlotRawData(senv$object, 
      #                                                    select=my_i, 
      #                                                    plots=TRUE,
      #                                                    original.data = T,
      #                                                    cutoff=input$NumberCutoffSlider,
      #                                                    markers=input$NumberMarkerSlider, 
      #                                                    comments=input$ShowComments)
      #     })
      #     
      #   })
      # }
      # 
      plot_output_list1 <- tagList(checkboxInput("SelectAllSamples", 
                                                 label = tags$strong("Select the complete dataset (all samples)")),
                                   checkboxInput("SelectPageSamples", 
                                                 label = paste("Select the first",input$NumberSampleSlider,"samples"),value = F))
      
      plot_output_list  <- lapply(1:input$NumberSampleSlider, function(i) {
        plotname <- paste("Sample", i, sep="")
        
        output[[plotname]] <- renderPlot({ PlotRawData(senv$object, 
                                                       select=i, 
                                                       plots=TRUE,
                                                       original.data = T,
                                                       cutoff=input$NumberCutoffSlider,
                                                       markers=input$NumberMarkerSlider, 
                                                       comments=input$ShowComments)
        })
        
        plotcheckbox <- paste("plotcheckbox", i, sep="")
        tags$div(class = "group-output",
                 uiOutput(plotcheckbox),
                 plotOutput(plotname))
      })
      
      # Convert the list to a tagList - this is necessary for the list of items
      # to display properly.
      append(plot_output_list1,do.call(tagList, plot_output_list))
    }
    
  })
  
  #Correction tab: plot the autocorrected plots
  output$autocorrection <- renderUI({
    
    if (is.null(input$cancerdata)) {
      tags$b("tcga is null")
    }
    
    if(is.null(senv$max_plots)){
      tags$h2("Please select samples for baseline autocorrection")} else {
        
        # find the selected samples
        selsample <- c(rep(FALSE,times=senv$max_plots))
        
        if(input$SelectAllSamples){
          selsample <- c(rep(TRUE,times=senv$max_plots))
        } else {
          for (g in 1:input$NumberSampleSlider) {
            if (input[[paste("PlotRawSamplecheckbox", g, sep="")]]) {
              selsample[g] <- TRUE
            }
          }
        }
        
        senv$selsample <- selsample
        
        NumbCorrectedPlots <- 0
        #make plots and store them 
        for (i in 1:senv$max_plots) {
          NumbCorrectedPlots <- local({
            my_i <- i
            
            if(input$SelectAllSamples || senv$selsample[my_i]){
              
              NumbCorrectedPlots   <- NumbCorrectedPlots+1
              my_corr              <- NumbCorrectedPlots
              
              plotname             <- paste("SampleCorrect", my_corr, sep="")
              plotslider           <- paste("correctplotSlider", my_corr, sep="")
              
              
              output[[plotslider]] <- renderUI({if(is.null(input[[plotslider]])){sliderInput(plotslider,
                                                                                             paste("Correct baseline: ",my_i, sep=""),
                                                                                             max=2, min=-2,
                                                                                             width='800px',
                                                                                             value=0,step=0.01)} else {
                                                                                               sliderInput(plotslider,
                                                                                                           paste("Correct baseline: ",my_i, sep=""),
                                                                                                           max=2, min=-2,
                                                                                                           width='800px',
                                                                                                           value=input[[plotslider]],step=0.01)
                                                                                             }
              })
              
              output[[plotname]]   <- renderPlot({if(!is.null(input[[plotslider]])){
                senv$object <- Plot.Manual(senv$object, 
                                           select=my_i,
                                           cutoff=input$NumberCutoffSlider,
                                           markers=input$NumberMarkerSlider, 
                                           comments=input$ShowComments,
                                           slider_value=input[[plotslider]])} else {
                                             senv$object <- Plot.Manual(senv$object, 
                                                                        select=my_i,
                                                                        cutoff=input$NumberCutoffSlider,
                                                                        markers=input$NumberMarkerSlider, 
                                                                        comments=input$ShowComments,
                                                                        slider_value=0)
                                           }
              }) 
              
            }
            return(NumbCorrectedPlots) })
        }
        
        if (NumbCorrectedPlots==0){
          tags$div(tags$p("Please select at Least one sample"))
        } else {
          #process bar
          withProgress(message = 'Plotting...', value = 0, {
            
            #update the pages
            updatePageruiInput(session,"pager",pages_total = ceiling((sum(senv$selsample)/25)))
            
            #make the list of the plots to show on the page 
            plot_output_list_corr <- lapply(
              ((25*(input$pager$page_current-1))+1) : min(sum(senv$selsample),(25*input$pager$page_current)), 
              function(s) {
                
                #update the process bar
                incProgress(1/sum(senv$selsample), detail = paste("Processing...."))
                
                #take the plot and the slider
                plotname <- paste("SampleCorrect", s, sep="")
                plotslider <- paste("correctplotSlider", s, sep="")
                tags$div(class = "group-output",
                         wellPanel(output$code <- renderUI({}),
                                   plotOutput(plotname),
                                   uiOutput(plotslider,align = "center"),
                                   output$code <- renderUI({}),
                                   output$code <- renderUI({})))
              })
            
            do.call(tagList, plot_output_list_corr)
            
          })
        }
      }
  })
  
  #######################
  ### Go to buttons #####
  #######################
  output$sampleButtonG2Raw <- renderUI({
    if (is.null(input$file2))
      return(NULL)
    actionButton("SampleActionButton1", label = "Plot samples?")
  })
  
  output$regionsbuttonsGo2Sample <- renderUI({
    if (is.null(input$file1))
      return(NULL)
    actionButton("RegionsActionButtonGo2Sample", label = "Load sample sheet? (Optional)")
  })
  
  output$regionsbuttonsGo2PlotRaw <- renderUI({
    if (is.null(input$file1))
      return(NULL)
    actionButton("RegionsActionButtonGo2PlotRaw", label = "Plot samples directly?")
  })
  
  output$LoadDataButtons <- renderUI({
    
    tagList(tags$h3("Select one of these options:"),
            #img(src='bar.gif', align = "left"),
            p(),
            actionButton("UpLoadData", label = "Upload Your Data") ,
            p(),
            actionButton("LoadSampleData", label = "Load Example Data (4 samples)"),
            
            p(),
            actionButton("LoadFromCancerAtlasData", label = tags$strong("Load Data from The Cancer Genome Atlas (TCGA)"))
            
    )
  })
  
  ####################################
  ###   Download panel buttond #######
  ####################################
  output$downloadButtons <- renderUI({tagList(
    downloadButton("downloadQC", "Download QC File"),
    downloadButton("downloadPlot", "Download Corrected Plot(s)"),
    downloadButton("downloadPlot.original", "Download Original Plot(s)"),
    downloadButton("downloadRegions", "Download Corrected Segments Table"),
    downloadButton("downloadShifts", "Download Shifts (modification by tool & user)")
  )
  })
  
  #######################################################
  #############  Output for download ####################
  ##### Creates and downloadHandler and plots and QC ####
  #######################################################
  output$downloadPlot <- downloadHandler(
    filename = 'Corrected_plots.zip',
    content  = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())
      SelPlots <- 0
      
      #########
      for (i in 1:senv$max_plots) {
        my_i <- i
        
        if(input$SelectAllSamples || senv$selsample[my_i]){
          
          SelPlots <- SelPlots+1
          plotname <- paste0("Sample_",my_i,"_",senv$object$SL[my_i,"Sample"],".png",sep="")
          
          if (SelPlots==1) {
            fs <- c(plotname)
          } else {
            fs <- append(fs,plotname)
          }
          
          png(filename=plotname,width=1920,height=1200)
          PlotRawData(senv$object,
                      select=my_i,
                      plots=TRUE,
                      cutoff=input$NumberCutoffSlider,
                      markers=input$NumberMarkerSlider,
                      comments=input$ShowComments)
          dev.off()
        }
      }
      
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip")
  
  output$downloadPlot.original <- downloadHandler(
    filename = 'Original_plots.zip',
    content  = function(fname) {
      tmpdir <- tempdir()
      setwd(tempdir())
      SelPlots <- 0
      
      #########
      for (i in 1:senv$max_plots) {
        my_i <- i
        
        if(input$SelectAllSamples || senv$selsample[my_i]){
          
          SelPlots <- SelPlots+1
          #plotname <- paste("Sample",my_i,".png",sep="")
          plotname <- paste0("Sample_",my_i,"_",senv$object$SL[my_i,"Sample"],".png",sep="")
          
          if (SelPlots==1) {
            fs <- c(plotname)
          } else {
            fs <- append(fs,plotname)
          }
          
          png(filename=plotname,width=1920,height=1200)
          PlotRawData(senv$object,
                      original.data = T,
                      select=my_i,
                      plots=TRUE,
                      cutoff=input$NumberCutoffSlider,
                      markers=input$NumberMarkerSlider,
                      comments=input$ShowComments)
          dev.off()
        }
      }
      
      zip(zipfile=fname, files=fs)
    },
    contentType = "application/zip")
  
  output$downloadRegions <- downloadHandler(
    filename = "corrected_regions.csv", 
    content = function(file) {
      names         <- senv$object$SL[senv$selsample,"Sample"]
      # get the segments for the selected samples
      selregions          <- senv$object$regions[which(senv$object$regions$Sample %in% as.character(names)),]
      write.csv(selregions,file)}) 
  
  output$downloadQC <- downloadHandler(
    filename = "QC.csv",
    content = function(file) {
      senv$object <- QC.function(senv$object)
      
      names         <- senv$object$SL[senv$selsample,"Sample"]
      # get the QC for the selected samples only
      selQC          <- senv$object$QC[which(senv$object$QC$Sample %in% as.character(names)),]
      
      write.csv(selQC,file,row.names = F)
    })
  
  output$downloadShifts <- downloadHandler(
    filename = "Shifts.csv",
    content = function(file) {
      
      #Prepare the output file
      tab <- senv$object$mod_auto
      colnames(tab) <- c("Sample", "Original_Maximum_Peak", "Shifting", "Corrected_by")
      
      #calculate total shifts
      tab_save <- senv$object$regions_save[!duplicated(senv$object$regions_save$Sample),c("Sample","Mean")]
      tab_auto <- senv$object$regions_auto[!duplicated(senv$object$regions_auto$Sample),c("Sample","Mean")]
      tab_last <- senv$object$regions[!duplicated(senv$object$regions$Sample),c("Sample","Mean")]
      
      #order the tab_save and tab_auto same order of tab
      tab_save <- tab_save[match(tab_save$Sample,tab$Sample),]
      tab_auto <- tab_auto[match(tab_auto$Sample,tab$Sample),]
      tab_last <- tab_last[match(tab_last$Sample,tab$Sample),]
      
      #calculate tool shifts
      diff_auto <- tab_auto$Mean - tab_save$Mean
      diff_last <- tab_last$Mean - tab_save$Mean
      diff_user <- diff_last - diff_auto
      
      df <- data.frame(Sample=tab_last$Sample, diff_last=diff_last, diff_user=diff_user)
      
      #find the sample with user manual modification
      #get sample name for samples with user modification
      diff_user_samples <- unique(df$Sample[which(df$diff_user != 0)])
      
      # get the correction by: User or Tool?
      tab$Corrected_by <- "Tool"
      tab$Corrected_by[which(tab$Sample %in% diff_user_samples)] <- "User"
      
      #get the Shifts to the output file
      tab$Shifting <- df$diff_last[!is.na(df$diff_last)]#[df$Sample %in% tab$Sample]
      
      #keep the selected samples only
      names  <- senv$object$SL[senv$selsample,"Sample"]
      tab    <- tab[which(tab$Sample %in% as.character(names)),]
      
      write.csv(tab,file)
    })
  
  output$downloadtool <- downloadHandler( 
    filename = "TCGA.CNV.BaselineCorrection.zip",
    content = function(file) {
      
      zip(zipfile=file, files=file.path("www/TCGA.CNV.BaselineCorrection.zip"))
    },
    contentType = "application/zip")
  
  output$downloadtool2 <- downloadHandler(
    filename = "TCGA.CNV.BaselineCorrection.zip",
    content = function(file) {
      
      zip(zipfile=file, files=file.path("www/TCGA.CNV.BaselineCorrection.zip"))
    },
    contentType = "application/zip")
  
  output$downloadtool3 <- downloadHandler(
    filename = "TCGA.CNV.BaselineCorrection.zip",
    content = function(file) {
      
      zip(zipfile=file, files=file.path("www/TCGA.CNV.BaselineCorrection.zip"))
    },
    contentType = "application/zip")
  
  output$downloadexample <- downloadHandler(
    filename = "Example input data.zip",
    content = function(file) {
      
      zip(zipfile=file, files=file.path("DATA/DATA.zip"))
    },
    contentType = "application/zip")
  
  
}) # END of server
