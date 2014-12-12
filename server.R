# Author: Matt Watts
# Date: 12 Dec 2014
# Purpose: RunMarxanLPILPTas web app server.R

require(shiny)
require(sp)
require(maptools)
require(PBSmapping)
require(foreign)
require(sqldf)
require(vegan)
require(labdsv)
require(xtable)
require(lpSolveAPI)

cat(paste0("\n",sMarxanDir,"\n"))

PrepareDisplay <- function()
{
    # prepare the map: pulayer object
    sPulayer <- paste(sMarxanDir,"/pulayer/pulayer.shp",sep="")
    pulayer <- readShapePoly(sPulayer)
    pulayer <<- SpatialPolygons2PolySet(pulayer)
    pu_table <<- read.dbf(paste(sMarxanDir,"/pulayer/pulayer.dbf",sep=""))

    # prepare the cluster analysis objects and map objects for ILP & LP
    solutions_raw<-read.table(paste0(sMarxanDir,"/output/output_solutionsmatrix.csv"),header=TRUE, row.name=1, sep=",")
    ilp_soln <- read.table(paste0(sMarxanDir,"/output/output_ilp_puidsel.csv"),header=TRUE, row.name=1, sep=",")
    lp_soln <- read.table(paste0(sMarxanDir,"/output/output_lp_puidsel.csv"),header=TRUE, row.name=1, sep=",")
    lp_ssoln <- read.table(paste0(sMarxanDir,"/output/output_lp_variables.csv"),header=TRUE, row.name=1, sep=",")
    lp_ssoln <- sqldf("SELECT x from lp_ssoln where x > 0")
    # prepare summary table data
    ILPpucount <<- nrow(ilp_soln)
    LPpucount <<- nrow(lp_soln)
    LPcost <<- read.csv(paste0(sMarxanDir,"/output/output_lp_objective.csv"))[1,2]
    ILPcost <<- read.csv(paste0(sMarxanDir,"/output/output_ilp_objective.csv"))[1,2]
    # join LP and ILP results to the solutions matrix
    ilp_row <- c()
    lp_row <- c()
    lp_ssolnrow <- c()
    for (i in 1:ncol(solutions_raw))
    {
      iILPValue <- 0
      iLPValue <- 0
      iLPssolnValue <- 0
  
      iPUID <- as.integer(substr(colnames(solutions_raw)[i],2,nchar(colnames(solutions_raw)[i])))
      iWhich <- which(ilp_soln[,1]==iPUID)
      if (length(iWhich > 0))
      { 
        iILPValue <- 1
      } 
      iWhich <- which(lp_soln[,1]==iPUID)
      if (length(iWhich > 0))
      {
        iLPValue <- 1
        iLPssolnValue <- lp_ssoln[iWhich,1]
      }

      ilp_row <- c(ilp_row,iILPValue)
      lp_row <- c(lp_row,iLPValue)
      lp_ssolnrow <- c(lp_ssolnrow,iLPssolnValue)
    }
    # move best solution to start of table
    thetable <- read.csv(paste0(sMarxanDir,"/output/output_sum.csv"))
    thetable <- round(sqldf("SELECT Score, Cost, Planning_Units, Penalty, Shortfall, Missing_Values, MPM from thetable"))
    iBest <- which.min(thetable[,1])
    Best <- solutions_raw[iBest,]
    solutions_raw <- solutions_raw[-iBest,]
    solutions_join <- rbind(ilp_row,lp_row,Best,solutions_raw)
    rownames(solutions_join) <- c("ILP","LP",paste0("S",iBest," (Best)"),row.names(solutions_raw))
    plotlabels <<- c("ILP","LP",paste0("S",iBest," (Best)"),row.names(solutions_raw))
    
    # prepare map fields
    ILPvalues <<- ilp_row+1
    LPvalues <<- lp_ssolnrow

    solutions <- unique(solutions_join)
    iUniqueSolutions <- dim(solutions)[1]
    nmdscolours <- rep("black",each = iUniqueSolutions)  
    nmdscolours[1] <- "blue"
    nmdscolours[2] <- "green"
    nmdscolours[3] <- "#40E0D0"
    nmdscolours <<- nmdscolours
    soldist<-vegdist(solutions,distance="bray")
    sol.mds<<-nmds(soldist,2)
    h<<-hclust(soldist, method="complete")
}

PrepareDisplay()

GetOutputFileext <- function(sMarxanDir,sParam)
# For the specified Marxan output file, return the file extension (.csv or .txt)
# Scan input.dat for the parameter,
# if value = 1, .dat, tab delimited, no header
# if value = 2, .txt, comma delimited (Qmarxan sets this value)
# if value = 3, .csv, comma delimited
{
  inputdat <- readLines(paste(sMarxanDir,"/input.dat",sep=""))
  iParam <- which(regexpr(sParam,inputdat)==1)
  
  iValue <- as.integer(unlist(strsplit(inputdat[iParam], split=" "))[2])
  
  if (iValue == 1)
  {
    return(".dat")
  }
  if (iValue == 2)
  {
    return(".txt")
  }
  if (iValue == 3)
  {
    return(".csv")
  }
}

GenerateSolnFilename <- function(iRunNumber,sMarxanDir)
{
  sFilename <- paste(sMarxanDir,"/output/output_r",sep="")  
  iPadding <- 5 - nchar(as.character(iRunNumber))
  if (iPadding > 0)
  {
    for (i in 1:iPadding)
    {
      sFilename <- paste(sFilename,"0",sep="")
    }
  }
  sFilename <- paste(sFilename,iRunNumber,GetOutputFileext(sMarxanDir,"SAVERUN"),sep="")  
}

ImportOutputsCsvToShpDbf <- function(sPuShapeFileDbf, sMarxanDir, iNumberOfRuns, sPUID)
# Imports the relevant contents of output files to the planning unit shape file dbf.
{
  # load and prepare pu_table
  pu_table <- read.dbf(sPuShapeFileDbf)
  pu_table <- sqldf(paste("SELECT ", sPUID, " from pu_table",sep=""))
  colnames(pu_table)[1] <- "PUID"
                    
  pu_table$PUID <- as.integer(pu_table$PUID)
  
  # load and prepare ssoln_table
  ssoln_table <- read.csv(paste(sMarxanDir,"/output/output_ssoln",GetOutputFileext(sMarxanDir,"SAVESUMSOLN"),sep=""))
  colnames(ssoln_table)[1] <- "PUID"
  colnames(ssoln_table)[2] <- "SSOLN2"
  ssoln_table$SSOLN1 <- as.integer(iNumberOfRuns - ssoln_table$SSOLN2)
  ssoln_table$SSOLN2 <- as.integer(ssoln_table$SSOLN2)
  
  # join pu_table and ssoln_table
  pu_table <- sqldf("SELECT * from pu_table LEFT JOIN ssoln_table USING(PUID)")
  
  # load and prepare best_table
  best_table <- read.csv(paste(sMarxanDir,"/output/output_best",GetOutputFileext(sMarxanDir,"SAVEBEST"),sep=""))
  best_table$BESTSOLN <- as.integer(best_table$SOLUTION + 1)
  best_table <- sqldf("SELECT PUID, BESTSOLN from best_table")
  
  # join pu_table and best_table
  pu_table <- sqldf("SELECT * from pu_table LEFT JOIN best_table USING(PUID)")
  
  for (i in 1:iNumberOfRuns)
  {
    sFieldName <- paste("SOLN",i,sep="")
    
    # load and prepare solnX_table
    solnX_table <- read.csv(GenerateSolnFilename(i,sMarxanDir))
    solnX_table[sFieldName] <- as.integer(solnX_table$SOLUTION + 1)
    solnX_table <- sqldf(paste("SELECT PUID, ",sFieldName," from solnX_table",sep=""))
  
    # join pu_table and solnX_table
    pu_table <- sqldf("SELECT * from pu_table LEFT JOIN solnX_table USING(PUID)")
    
    rm(solnX_table)
  }
  
  # save the new pu_table
  colnames(pu_table)[1] <- sPUID
  write.dbf(pu_table,sPuShapeFileDbf)  
}

labelCol <- function(x)
{
  thetable <- read.csv(paste0(sMarxanDir,"/output/output_sum.csv"))
  thetable <- round(sqldf("SELECT Score from thetable"))
  iBest <- which.min(thetable[,1])
  if (is.leaf(x))
  {
    ## fetch label
    a <- attributes(x)
    label <- attr(x, "label") 
    colour <- "black"
    if (label == "ILP") { colour <- "blue"}
    if (label == "LP") { colour <- "green"}
    if (label == paste0("S",iBest," (Best)")) { colour <- "#40E0D0" }
    cat(paste0("label ",label,"\n"))
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- c(a$nodePar, lab.col = colour)
  }
  return(x)
}

Run_ILP_LP <- function()
{
  cat("Run_ILP_LP\n")
  
  # prepare data
  load(paste0(sMarxanDir,"/input/lpsmtx_.Rda"))
  load(paste0(sMarxanDir,"/input/puid_filter_.Rda"))
  specdat <- read.csv(paste0(sMarxanDir,"/input/spec.dat"))
  targ_ <- sqldf("SELECT target from specdat")
  colnames(targ_)[1] <- "TARG"
  rownames(targ_) <- c("SPEC.10","SPEC.11","SPEC.12","SPEC.13","SPEC.14","SPEC.15","SPEC.16",
                       "SPEC.17","SPEC.18","SPEC.19","SPEC.20","SPEC.21","SPEC.22","SPEC.23",
                       "SPEC.24","SPEC.25","SPEC.26")
  targ_ <- as.matrix(targ_)

  # prepare LP
  n <- length(lpsmtx_[,1])#no. of variables
  m <- length(targ_[,1])#no. of constraints
  msp1coln <- as.vector(colnames(lpsmtx_)) #vector containing column names of msp1 (for oview-matrix)
  #puid <- as.vector(rownames(lpsmtx_))
  puid <- as.vector(puid_filter_)
  f <- make.lp(0,n)
  set.objfn(f, c(lpsmtx_[,1]))
  for(i in 1:m){
    add.constraint(f,c(lpsmtx_[,i+1]), ">=", targ_[i])
  }
  set.bounds(f, upper = rep(1,n), columns = seq(1,n))

  # run LP solve
  solve(f)
  cat("LP solve\n")

  # save results
  obj <- get.objective(f)
  var <- get.variables(f)
  sel <- which(var>0)
  puidsel <- c()
  for(i in 1:length(sel)){
    puidsel[i] <- puid[sel[i]]
  }
  mcost <- mean(c(lpsmtx_[,1]))
  sdcost <- sd(c(lpsmtx_[,1]))
  write.csv(obj, paste0(sMarxanDir,"/output/output_lp_objective.csv"))
  write.csv(var, paste0(sMarxanDir,"/output/output_lp_variables.csv"))
  write.csv(puidsel, paste0(sMarxanDir,"/output/output_lp_puidsel.csv"))
  write.csv(c(mcost, sdcost), paste0(sMarxanDir,"/output/output_lp_costana.csv"))

  # prepare ILP
  set.type(f, columns = seq(1,n), type = "binary")
  lp.control(f, timeout = 1, verbose = "normal")

  # run ILP solve
  solve(f)
  cat("ILP solve\n")
  
  # save results
  obj <- get.objective(f)
  var <- get.variables(f)
  sel <- which(var>0)
  puidsel <- c()
  for(i in 1:length(sel)){
    puidsel[i] <- puid[sel[i]]
  }
  mcost <- mean(c(lpsmtx_[,1]))
  sdcost <- sd(c(lpsmtx_[,1]))
  write.csv(obj, paste0(sMarxanDir,"/output/output_ilp_objective.csv"))
  write.csv(var, paste0(sMarxanDir,"/output/output_ilp_variables.csv"))
  write.csv(puidsel, paste0(sMarxanDir,"/output/output_ilp_puidsel.csv"))
  write.csv(c(mcost, sdcost), paste0(sMarxanDir,"/output/output_ilp_costana.csv"))
}

shinyServer(function(input, output, session) {

    system("touch /var/shiny-server/www/Marxan_LP_ILP/Tas_Activity_rev3/restart.txt")

    observe({
        input$feature
        # change the target text control for the selected feature
        specdat <- read.csv(paste(sMarxanDir,"/input/spec.dat",sep=""),stringsAsFactors=FALSE)
        for (j in 1:nrow(specdat))
        {
            if (specdat[j,1] == input$feature)
            {
                updateNumericInput(session, "target", value = specdat[j,2])
                updateNumericInput(session, "spf", value = specdat[j,3])
            }
        }
    })

    runmarxan <- reactive({
        if (input$mrun == 0)
        {
            imrun <<- 0
            cat("init mrun\n")
        }
        else
        {
            if (input$mrun > imrun)
            {
                imrun <<- input$mrun
                cat("mrun incremented\n")
                
                # run Marxan
                system("./MarOpt_v243_Linux64 -s")
                #system("./MarOpt_v243_Mac64 -s")
                # run ILP & LP
                Run_ILP_LP()
                
                # write results to pu dbf
                ImportOutputsCsvToShpDbf(paste0(sMarxanDir,"/pulayer/pulayer.dbf"),
                                         sMarxanDir,iNUMREPS,"PUID")
                
                # fetch the results
                PrepareDisplay()

                # trigger a refresh of the UI
                irefreshinput <<- irefreshinput + 1
                updateNumericInput(session, "refreshinput", value = irefreshinput)
            }
        }
    
        return(as.character(input$mrun))
    })

    outputmap <- reactive({

        input$refreshinput
        
        if (input$map == "ilpmap")
        {
            # like bestmap
            greenramp <- colorRampPalette(c("white","green"))(2)
            colours <- rep(greenramp[1],length(ILPvalues))
            for (j in 1:length(ILPvalues))
            {
                colours[j] <- greenramp[ILPvalues[j]]
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }
        if (input$map == "lpmap")
        {
            # continuous between 0 & 1
            blueramp <- colorRampPalette(c("white","blue"))(16)
            colours <- rep(blueramp[1],length(LPvalues))
            for (j in 1:length(LPvalues))
            {
                colours[j] <- blueramp[round(15 * LPvalues[j])+1]
                if (LPvalues[j] > 0)
                {
                    cat(paste0("j ",j," LP ",LPvalues[j]," colour ",round(15 * LPvalues[j])+1,"\n"))
                }
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }
        if (input$map == "ssolnNmap")
        {
            values <- sqldf(paste("SELECT SSOLN2 from pu_table",sep=""))
            blueramp <- colorRampPalette(c("white","blue"))(16)
            colours <- rep(blueramp[1],nrow(values))
            for (j in 1:nrow(values))
            {
                colours[j] <- blueramp[round(15 / iNUMREPS * values[j,])+1]
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }

        if (input$map == "bestmap")
        {
            values <- sqldf("SELECT BESTSOLN from pu_table")
            greenramp <- colorRampPalette(c("white","green"))(2)
            colours <- rep(greenramp[1],nrow(values))
            for (j in 1:nrow(values))
            {
                colours[j] <- greenramp[values[j,]]
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }

        if (input$map == "runMmap")
        {
            values <- sqldf(paste("SELECT SOLN",input$m," from pu_table",sep=""))
            greenramp <- colorRampPalette(c("white","green"))(2)
            colours <- rep(greenramp[1],nrow(values))
            for (j in 1:nrow(values))
            {
                colours[j] <- greenramp[values[j,]]
            }
            plotPolys(pulayer,col=colours,axes=FALSE,border=NA,cex.lab=0.1,cex.axis=0.1)
        }
    })

    outputtable <- reactive({

        input$refreshinput
        
        if (input$savetargetspf == 0)
        {
            isavetargetspf <<- 0
            cat("init savetargetspf\n")
        }
        else
        {
            if (input$savetargetspf > isavetargetspf)
            {
                isavetargetspf <<- input$savetargetspf
                cat("savetargetspf incremented\n")
                
                rtarget <- input$target
                if (rtarget < 0)
                {
                    rtarget <- 0
                }
                rspf <- input$spf
                if (rtarget < 0)
                {
                    rtarget <- 0
                }

                # save the target table
                # save prop to spec.dat
                specdat <- read.csv(paste(sMarxanDir,"/input/spec.dat",sep=""),stringsAsFactors=FALSE)
                # change the value only for the row with name == input$feature
                for (j in 1:nrow(specdat))
                {
                    if (specdat[j,1] == input$feature)
                    {
                        specdat[j,2] <- rtarget
                        specdat[j,3] <- rspf
                    }
                }
                write.csv(specdat,paste0(sMarxanDir,"/input/spec.dat"),quote=FALSE,row.names=FALSE)
            }
        }

        if (input$table == "spec")
        {
            thetable <- read.csv(paste0(sMarxanDir,"/input/spec.dat"))
        }
        if (input$table == "sumtable")
        {
            sFilename <- paste(sMarxanDir,"/output/output_sum.csv",sep="")
            thetable <- read.csv(sFilename)
            thetable <- round(sqldf("SELECT Score, Cost, Planning_Units, Penalty, Shortfall from thetable"))
            iBest <- which.min(thetable[,1])
            Run <- c()
            for (j in 1:nrow(thetable))
            {
                if (j == iBest)
                {
                    Run <- c(Run,"Best")
                } else {
                    Run <- c(Run,j)
                }
            }        
                        
            thetable <- cbind(Run,thetable)#,
            thetable$Run <- as.character(thetable$Run)
            thetable <- rbind(thetable,
                              c("ILP","",round(ILPcost),ILPpucount,"",""),
                              c("LP","",round(LPcost),LPpucount,"",""))
            thetable$Run <- as.character(thetable$Run)
            for (j in 1:6)
            {
                thetable[iBest,j] <- HTML(paste0("<FONT COLOR='#40E0D0'>",thetable[iBest,j],"</FONT>"))
                thetable[11,j] <- HTML(paste0("<FONT COLOR='blue'>",thetable[11,j],"</FONT>"))
                thetable[12,j] <- HTML(paste0("<FONT COLOR='green'>",thetable[12,j],"</FONT>"))
            }
        }

        return(thetable)
    })

    output2ds <- reactive({
            input$refreshinput

            #return(plot(sol.mds$points, xlab='', ylab='', main='NMDS of solutions', col=nmdscolours))
            plot(sol.mds$points, xlab='', ylab='', main='NMDS of solutions', col=nmdscolours)
            text(sol.mds$points,labels=plotlabels,pos=4, col=nmdscolours)
    })

    outputdendogram <- reactive({
            input$refreshinput

            d <- dendrapply(as.dendrogram(h), labelCol)
            return(plot(d, xlab="Solutions", ylab="Disimilarity", main="Bray-Curtis dissimilarity of solutions"))
    })

    output$marxanmap <- renderPlot({
            print(outputmap())
    }, height=450,width=450)#height="auto", width="auto") #

    output$marxantable <- renderTable({
            dat <- data.frame(outputtable())
            dat
    }, sanitize.text.function = function(x) x)

    output$textfeedback = renderText({
        runmarxan()
        sprintf("Finished")
    })

    output$plot2ds <- renderPlot({ 
        print(output2ds())
    })
    
    output$plotdendogram <- renderPlot({ 
        print(outputdendogram())
    })
})
