library(shiny)
library(pracma) #provides error function
library(ggplot2)
library(reshape2)
library(gridExtra)
library(directlabels)



    ui = fluidPage(
      
         
      sidebarPanel(
        checkboxInput("header", "Check if your difficulty file has a header", FALSE),
        fileInput("file1", "CSV file of item difficulties (as a column vector)",accept =".csv"),
        checkboxInput('lookup',"Do you want to use a lookup table?",FALSE),
        
        #conditional panel for if they have a lookup table
        conditionalPanel(condition='input.lookup==true',
           checkboxInput("header2", "Check if your conversion file has a header", FALSE),
           fileInput('file2','Lookup matrix with two columns: Col 1=%Correct, Col 2=desired scale',accept='.csv')
        )#end lookup conditoinal
        
        ), #end sidebar panel
      
      
      mainPanel(
         checkboxInput('rantau','Check if you want to treat the true cut-score as a known constant',FALSE),
         
         #conditional panel for if they have a lookup table
         conditionalPanel(condition='input.rantau==true',
          numericInput('taustar1','True cut-score value (theta scale)',0,min=-4,max=4,step=.1)
         )#end conditional panel for rantau
         
       ), #end main panel
      
      
      #small conditional panel for inputs that goes away if tau star is constant
    conditionalPanel(condition='input.rantau==false',    # for some weird reason, shiny uses lower case TRUE/FALSE
       numericInput("meanang2", "Mean of angoff ratings (theta)",0,min=-4,max=4,step=.1),    
       numericInput("sdang2", "SD of angoff ratings (theta)", 1, min = .01, max = 100, step=.1)
    ), #end conditional panel

      #rest of inputs
      numericInput("meantruescore2", "Mean of examinee thetas", 0,min=-4,max=4,step=.1),  
      numericInput("sdobs2", "SD of examinee thetas",1, min = .01, max = 100, step=.1),
      numericInput("reliability2","Reliability of test",.5,min=0.01,max=1,step=.1),
      checkboxInput("penalty", "Check if you want penalty function values", FALSE),
      actionButton('go','Calculate and Plot',icon("arrow-right"),style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    downloadButton("report", "Generate report"),
    
    #outputs
    verbatimTextOutput("a"), #error display in case of error or processing time display if no error
    verbatimTextOutput("value"), #text for abs error min
    plotOutput(outputId="fnandfpoverlay"), #plot for abs error min
    
    #conditional panel for conditional which goes away if using penalty
    conditionalPanel(condition='input.penalty==false',
      verbatimTextOutput("value2"),
      plotOutput(outputId="rfnandrfpoverlay")
    ), #end conditional
    
    verbatimTextOutput("value3"),
    plotOutput(outputId="totalerror")

    )#end UI
    
    
    
    
    

        server = function(input, output, session) {
          
         observeEvent(input$go, {
          
            if (input$go == 0)   #this ingenious little bit of code stops the app from running prior to the first click of button
              return()
          

      output$a <- renderText({    #### This function does the math: creates a callable matrix that is the values of the WCE, called with outmat2
                                  #it also gives us a nice place for error messages, which, if all goes right, we let be processing time
        inFile <- isolate(input$file2)

        #handling the case when we have a lookup table indicated by user
        if (isolate(input$lookup) == TRUE) 
        {
          if (is.null(inFile)) #building my own halt function. I use the built-in version called validate below
            {return('You indicated you had a lookup table, please select it')} #return ends the 'a' function and gives error message
          
            convto3digit=function(catmin){
              inFile <- isolate(input$file2)
              
              convtable=isolate(read.csv(inFile$datapath, header = input$header2))
              end=nrow(convtable)+1
              
              k=1
              kbest=1
              while(k<end)
              {
                best=(abs(convtable[kbest,1]-catmin))
                if(best>abs(convtable[k,1]-catmin))
                {
                  kbest=k
                }
                
                k=k+1
              }#end loop over possible matches
              
              cconvto3digit=convtable[kbest,2]
              output=round(cconvto3digit,digits=3)
              output2=paste('Converted score at min is',output)
              return(output2)
            }#end of conversion function  
          
          
        }#end of if (isolate(input$lookup)==TRUE
        
       
        
        #case when there is no lookup table (just provide output saying they didnt' specify one
        if (isolate(input$lookup) == FALSE) 
        {
          convto3digit=function(catmin){
          return('no conversion indicated')
          }
          
        }#end of if (isolate(input$lookup)==FALSE
        
        
        #reading in difficulties file
        inFile <- isolate(input$file1)
        
        if (is.null(inFile)) #error out if they didn't include the file
          return('You must provide a difficulty file')
        
       dat=isolate(read.csv(inFile$datapath, header = input$header))
       difficulties=dat[,1]

       

       validate(    #this fun bit of code is shiny specific language which tests a condition, and returns an error message if it fails
         need(isolate(input$sdobs2>=0), "Please provide a Examinee SD which is positive"),
         need(isolate(input$sdang2>=0), 'Please provide a Angoff SD which is positive'),
         need(isolate(input$reliability2>=0), 'Please provide a Reliability which is positive')
       )
         

       #now we build all of our functions which will be called to do the math (there's a long list here)
        pitemcorrect=function(tau, bi)
        {(1/(1+exp(bi-tau)))}
        
        vectorofitemp=function(tau,dif)
        {
          leng=length(dif)
          newmat=vector(length=leng)
          i=1
          while(i<leng+1)
          {
            newmat[i]=pitemcorrect(tau, dif[i])
            i=i+1
            
          }
          return(newmat)
        }#end vector of item p
        
        pbarf=function(tau,dif){
          difficulties=dif
          itemps=vectorofitemp(tau,difficulties)
          pbar=sum(itemps)/length(difficulties)*100
          return(pbar)
        }#generates the pbar for a given value of tau
        
        sdtauf=function(tau,dif)
        {
          difficulties=dif
          itemps=vectorofitemp(tau,difficulties)
          out=sum((itemps)*((1-itemps)))
          n=length(difficulties)
          outval=(1/n)*sqrt(out)*100
          return(outval)
        } #generates sd of tau star
        
        phi=function(x)
        {(1/2)*(1+erf(x/sqrt(2)))}
        
        gofxfunc=function(x)
        {
          sqrt(exp(1))*exp(x)*phi(x+1)-phi(x)
        }
        
        
        if(isolate(input$rantau==FALSE))
        {
          ptaustarlesstauf=function(tau)
          {
            output=phi((tau-meanang)/sdang)
            return(output)
          }
          
          
          penaltyfn=function(tau)
          {
            sdtau=sdtauf(tau,difficulties)
            
            out=  phi((c-pbarf(tau,difficulties))/sdtau)*gofxfunc((tau-meanang)/sdang) 
            return(out) 
          }
          
          penaltyfp=function(tau)
          {
            sdtau=sdtauf(tau,difficulties)
            
            
            out=phi((pbarf(tau,difficulties)-c)/sdtau)*gofxfunc((meanang-tau)/sdang)
            return(out)
          }
          
         
          
        }
        
        if(isolate(input$rantau==TRUE))
        {
          ptaustarlesstauf=function(tau)
          {
            output=1
            return(output)
          }
          
          
          penaltyfn=function(tau)
          {
            sdtau=sdtauf(tau,difficulties)
            
            
            if(tau>taustar)
            {
            out=  phi((c-pbarf(tau,difficulties))/sdtau)*(exp(sqrt((taustar-tau)^2)/stderror)-1)
            }
            else
            {
              out=0
            }
            return(out) 
          }
          
          penaltyfp=function(tau)
          {
            sdtau=sdtauf(tau,difficulties)
            if(tau<taustar)
            {
              out=phi((pbarf(tau,difficulties)-c)/sdtau)*(exp(sqrt((taustar-tau)^2)/stderror)-1)
            }
            else
            {
              out=0
            }
          }
          
          
          
        }
          
        
        ptaustarlesstauf=function(tau)
        {
          output=phi((tau-meanang)/sdang)
          return(output)
        }
        
        pTnlesscf=function(tau)
        {
          pbar=pbarf(tau,difficulties)
          sdtau=sdtauf(tau,difficulties)
          output=phi((c-pbar)/sdtau)
          return(output)
        }
        
        
        falsenegf=function(tau)
        {
           (ptaustarlesstauf(tau)*pTnlesscf(tau))*(1/(sqrt(2*pi)*sdtruescore))*exp(-((tau-meantruescore)^2)/(2*(sdtruescore^2)))
        }
        
        
        falseposf=function(tau)
        {
          ((1-ptaustarlesstauf(tau))*(1-pTnlesscf(tau)))*(1/(sqrt(2*pi)*sdtruescore))*exp(-((tau-meantruescore)^2)/(2*(sdtruescore^2)))
        }

        pshouldpass=function(tau)
        {
          ptaustarlesstauf(tau)*(1/(sqrt(2*pi)*sdtruescore))*exp(-((tau-meantruescore)^2)/(2*(sdtruescore^2)))
        }
        

        
        penaltyfalsenegf=function(tau)
        {
          penaltyfn(tau)*exp(-((tau-meantruescore)^2)/(2*(sdtruescore^2)))
        }
        
        
        penaltyfalseposf=function(tau)
        {
          penaltyfp(tau)*exp(-((tau-meantruescore)^2)/(2*(sdtruescore^2)))
        }
        
    #end function building
        
      #begin calling values from the user inputs
        sdang=isolate(input$sdang2)
        meanang=isolate(input$meanang2)
        meantruescore=isolate(input$meantruescore2)
        rel=isolate(input$reliability2)
        sdobs=isolate(input$sdobs2)
        sdtruescore=sqrt(rel*(sdobs^2))
        taustar=isolate(input$taustar1)
        stderror=sqrt((sdobs^2)*(1-rel))
      #end calling user values
        
        
      #depending on if penalty is selected, we use a different function  
        if(isolate(input$penalty==FALSE))
        {
        falsevec=Vectorize(falsenegf,'tau')
        posvec=Vectorize(falseposf,'tau')
        }
        if(isolate(input$penalty==TRUE))
        {
          falsevec=Vectorize(penaltyfalsenegf,'tau')
          posvec=Vectorize(penaltyfalseposf,'tau')
        }
        ppassvec=Vectorize(pshouldpass)
        
        
      #creating vectors to save output, since we have variable length (due to optimization) we make them huge, set all to -1, then at end delete all below zero  
        fnvec=vector(length=1000) ; fnvec[1:1000]=-1
        fpvec=vector(length=1000); fpvec[1:1000]=-1
        rfnvec=vector(length=1000); rfnvec[1:1000]=-1
        rfpvec=vector(length=1000); rfpvec[1:1000]=-1
        cvec=vector(length=1000); cvec[1:1000]=-1
        
        #finding smart start value by testing 3 points, then we brach + from highest min and - 25 from the lowest
        startsvec=matrix(nrow=3,ncol=4)
          
        starts=c(25,50,75)
        k=1
          
          
        withProgress(message = 'Preparing to run calculations...', min=0,max=100, value = 0, {
        while(k<4)
        {
         c=starts[k] 
          
         
          if(isolate(input$rantau==FALSE)) 
          { 
            fn2=integrate(falsevec,-4,4)#\\,rel.tol=1e-13)#subdivisions=1000000000)#,rel.tol=.0000000000001)
            fp2=integrate(posvec,-4,4)
          }
          if(isolate(input$rantau==TRUE))
          {
            fn2=integrate(falsevec,taustar,4)#\\,rel.tol=1e-13)#subdivisions=1000000000)#,rel.tol=.0000000000001)
            fp2=integrate(posvec,-4,taustar)
          }
            
            
          ppass=integrate(ppassvec,-4,4)
          rfn2=fn2$value/ppass$value  
          rfp2=fp2$value/(1-ppass$value)
          
          
          maxve=max(rfn2,rfp2)

          
          maxve2=max(fp2$value,fn2$value)

          
          
          totalp=fp2$value+fn2$value # use this for lower
          
          startsvec[k,1]=starts[k]
          startsvec[k,2]=totalp
          startsvec[k,3]=maxve
          startsvec[k,4]=maxve2
          
          k=k+1
        }#end of startsloop
          
          
  
          })#end withprogress
          
          
          
          
          
          
          #here we identify which metric has the smallest min c (and largest) and use those points -25 (and +25) as starting and end points
          totalminloc=which(startsvec[,2]==min(startsvec[,2]))
          totalmin=starts[totalminloc]
          
          absminloc=which(startsvec[,4]==min(startsvec[,4]))
          absmin=starts[absminloc]
          
          condminloc=which(startsvec[,3]==min(startsvec[,3]))
          condmin=starts[condminloc]
          
          cminvec=c(totalmin,absmin,condmin)
          cactmin=min(cminvec)-25
          cactmax=max(cminvec)+25
          
          cmax=cactmax
          cmin=cactmin

          #just in case we've found values below zero or above 100
          if(cmin<=0)
          {
            cactmin=0
          }
          if(cmax>=100)
          {
            cactmax=100
          }
          
          dist=cactmax-cactmin+57 #This is a precise statement telling the withProgress statment when to end
                                  #57 is from the 3 .1 incriment runs with abs,cond,and total
          
        #progress bar
          time1=Sys.time()
          withProgress(message = 'Calculating... please wait', min=0,max=dist, value = 0, {
          
          
          c=cactmin
          end=cactmax
          i=1
          
          
          
          #first we do a loop over the possible ranges, with increments of 1
          while(c<end)
          {
            
            
            
            
            fn2=integrate(falsevec,-4,4)#\\,rel.tol=1e-13)#subdivisions=1000000000)#,rel.tol=.0000000000001)
            fp2=integrate(posvec,-4,4)
            ppass=integrate(ppassvec,-4,4)
            rfn2=fn2$value/ppass$value  # us this for upper
            rfp2=fp2$value/(1-ppass$value)
            
            fnvec[i]=round(fn2$value,digits=5)
            fpvec[i]=round(fp2$value,digits=5)
            
            totalp=fnvec+fpvec # use this for lower
            
            rfnvec[i]=round(rfn2,digits=5)
            rfpvec[i]=round(rfp2,digits=5)
            cvec[i]=c
            
            
          
          i=i+1

          c=c+1
          
          
          
          incProgress(1)
          }#end loop with big increments
          


          outmat=cbind(cvec,fnvec,fpvec,rfnvec,rfpvec)
          outmat2=subset(outmat,outmat[,2]>=0)
          
          
          #now we go more precise for each error,  (<-1->) with steps of .1
          
            #absolute error loop
            outmat=outmat2
            fnvec2=outmat[,2]
            fpvec2=outmat[,3]
            cvec2=outmat[,1]
            maxvec=pmax(fnvec2,fpvec2)
            bound=cbind(cvec2,maxvec)
            minofmax=min(maxvec)
            loc=which(bound[,2]==min(bound[,2]))
            catmin=bound[loc,1]
            threedig=convto3digit(catmin)
            
            
  
            c=catmin-.9
            end=catmin+1
            while(c<end)
              {
              
                fn2=integrate(falsevec,-4,4)#\\,rel.tol=1e-13)#subdivisions=1000000000)#,rel.tol=.0000000000001)
                fp2=integrate(posvec,-4,4)
                ppass=integrate(ppassvec,-4,4)
                rfn2=fn2$value/ppass$value  # us this for upper
                rfp2=fp2$value/(1-ppass$value)
                
                fnvec[i]=round(fn2$value,digits=5)
                fpvec[i]=round(fp2$value,digits=5)
                
                totalp=fnvec+fpvec # use this for lower
                
                rfnvec[i]=round(rfn2,digits=5)
                rfpvec[i]=round(rfp2,digits=5)
                cvec[i]=c
                
                
                
                i=i+1
                
                c=c+.1  
                incProgress(1)
              }
          
            #conditional error loop
            outmat=outmat2
            fnvec2=outmat[,4]
            fpvec2=outmat[,5]
            cvec2=outmat[,1]
            maxvec=pmax(fnvec2,fpvec2)
            bound=cbind(cvec2,maxvec)
            minofmax=min(maxvec)
            loc=which(bound[,2]==min(bound[,2]))
            catmin=bound[loc,1]
            threedig=convto3digit(catmin)
        
            c=catmin-.9
            end=catmin+1
            while(c<end)
            {
              fn2=integrate(falsevec,-4,4)#\\,rel.tol=1e-13)#subdivisions=1000000000)#,rel.tol=.0000000000001)
              fp2=integrate(posvec,-4,4)
              ppass=integrate(ppassvec,-4,4)
              rfn2=fn2$value/ppass$value  # us this for upper
              rfp2=fp2$value/(1-ppass$value)
              
              fnvec[i]=round(fn2$value,digits=5)
              fpvec[i]=round(fp2$value,digits=5)
              
              totalp=fnvec+fpvec # use this for lower
              
              rfnvec[i]=round(rfn2,digits=5)
              rfpvec[i]=round(rfp2,digits=5)
              cvec[i]=c
              
              
              
              i=i+1
              
              c=c+.1  
              incProgress(1)
            }
          
            #total error loop
            outmat=outmat2
            fnvec2=outmat[,3]
            fpvec2=outmat[,2]
            cvec2=outmat[,1]
            
            totalp=fnvec2+fpvec2
            bound=cbind(cvec2,totalp)
            minofmax=min(totalp)
            #minnie2=minnie[1]
            loc=which(bound[,2]==minofmax)
            loc2=min(loc) #neccessary becauase loc returns a vector in this case because there are two adjacent points with same error
            catmin=bound[loc2,1]
            
            
            c=catmin-.9
            end=catmin+1
            while(c<end)
            {
              fn2=integrate(falsevec,-4,4)#\\,rel.tol=1e-13)#subdivisions=1000000000)#,rel.tol=.0000000000001)
              fp2=integrate(posvec,-4,4)
              ppass=integrate(ppassvec,-4,4)
              rfn2=fn2$value/ppass$value  # us this for upper
              rfp2=fp2$value/(1-ppass$value)
              
              fnvec[i]=round(fn2$value,digits=5)
              fpvec[i]=round(fp2$value,digits=5)
              
              totalp=fnvec+fpvec # use this for lower
              
              rfnvec[i]=round(rfn2,digits=5)
              rfpvec[i]=round(rfp2,digits=5)
              cvec[i]=c
              
              
              
              i=i+1
              
              c=c+.1  
              incProgress(1)
            }
          
          
        #output  
          outmat=cbind(round(cvec,digits=1),fnvec,fpvec,rfnvec,rfpvec)
          outmat1=subset(outmat,outmat[,2]>=0)
          q=duplicated(outmat1[,1])
          outmat2=outmat1[!q,]
       
        })#end withprogress
        
      

  #####*****---------------------------------------------------
          #All calculations have been done above. Everything below is 
          #part of the output text or plots, all statements are
          #calls to the data generated above, or manipulations for
          #pretty output
      

      
      
 
      
      
      
      
      
      
      
      
      output$value=renderText({   #makes the text output

      #absolute error text  
                
        outmat=outmat2
        fnvec=outmat[,2]
        fpvec=outmat[,3]
        cvec=outmat[,1]
        maxvec=pmax(fnvec,fpvec)
        bound=cbind(cvec,maxvec)
        minofmax=min(maxvec)
        loc=which(bound[,2]==min(bound[,2]))
        catmin=bound[loc,1]
        threedig=convto3digit(catmin)
        

  
        if(isolate(input$penalty==FALSE))
        {
          textout="the optimal absolute error value is"
        }
        
        if(isolate(input$penalty==TRUE))
        {
          textout="the optimal absolute penalty error value is"
        }
        paste(paste(textout, round(min(minofmax),digits=3),"{fp = fn =",round(min(minofmax),digits=3),'}'), paste("%C at min is",round(catmin,digits=3),',',threedig),sep='\n')
        
        
        
     # })#end with prog  
        
         }) #end absolute error text  
      
      
        
     #absolute error plot
      output$fnandfpoverlay=renderPlot({
        outmat=outmat2
        #sdobs=input$sdobs
        plotmat1=outmat[,1]
        plotmat2=outmat[,3]
        plotmat23=outmat[,2]
        plotmat3=cbind(plotmat1,plotmat2,plotmat23)
        
        leng=length(plotmat1)
        fname=vector(length=leng)
        fname[1:leng]=1#'FN'
        plottest=cbind(plotmat1,plotmat2,fname)
        fname[1:leng]=2#'FP'
        plottest2=cbind(plotmat1,plotmat23,fname)
        plottest3=rbind(plottest,plottest2)
        
        colnames(plottest3)=c('Cvalue','value','fname')
        plottest4=as.data.frame(plottest3)
        
        colnames(plotmat3)=c('Cvalue', 'FPcond','FNcond')
        
        melty=melt(plottest4,id=c('Cvalue','fname'))
        
        if(isolate(input$penalty==FALSE))
        {
          titles='Absolute Errors'
          lab=c('FP','FN')
        }
        
        if(isolate(input$penalty==TRUE))
        {
          titles='Absolute Penalty Errors'
          lab=c('PFP','PFN')
        }
        
        
        
        fnvec=outmat[,2];fpvec=outmat[,3];cvec=outmat[,1]
        maxvec=pmax(fnvec,fpvec)
        bound=cbind(cvec,maxvec)
        minofmax=min(maxvec)
        ymax=min(minofmax)*4
        
        ggplot(melty,aes(x=Cvalue,y=value,group=factor(fname),color=factor(fname)))+xlab('%Correct cut value')+labs(color='Error Type')+geom_line()+scale_color_manual(labels = lab, values = c("blue", "red"))+ggtitle(titles)+coord_cartesian(ylim = c(0,ymax))
        
      })# end absolute error plot  
      
    #conditional error text   
      output$value2=renderText({   #makes the text output
        
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        outmat=outmat2
        fnvec=outmat[,4]
        fpvec=outmat[,5]
        cvec=outmat[,1]
        maxvec=pmax(fnvec,fpvec)
        bound=cbind(cvec,maxvec)
        minofmax=min(maxvec)
        loc=which(bound[,2]==min(bound[,2]))
        catmin=bound[loc,1]
        threedig=convto3digit(catmin)
        

        
        if(isolate(input$penalty==TRUE))
        {
          textout='no conditional calculations done for penalty function'
          paste(textout)
        }
        else
        {
          textout="the optimal conditional error value is"
          paste(    paste(textout, round(min(minofmax),digits=3),"{rfp = rfn =",round(min(minofmax),digits=3),'}') ,paste("%C at min is",round(catmin,digits=3),',',threedig), sep='\n'     )
          
        }
        
        
        
        
      })#end conditional error text   
      
      
      
    #conditional error plot
      output$rfnandrfpoverlay=renderPlot({
        outmat=outmat2
        #sdobs=input$sdobs
        plotmat1=outmat[,1]
        plotmat2=outmat[,5]
        plotmat23=outmat[,4]
        plotmat3=cbind(plotmat1,plotmat2,plotmat23)
        
        leng=length(plotmat1)
        fname=vector(length=leng)
        fname[1:leng]=1#'FN'
        plottest=cbind(plotmat1,plotmat2,fname)
        fname[1:leng]=2#'FP'
        plottest2=cbind(plotmat1,plotmat23,fname)
        plottest3=rbind(plottest,plottest2)
        
        colnames(plottest3)=c('Cvalue','value','fname')
        plottest4=as.data.frame(plottest3)
        
        colnames(plotmat3)=c('Cvalue', 'FPcond','FNcond')
        
        melty=melt(plottest4,id=c('Cvalue','fname'))
        
        
        fnvec=outmat[,5];fpvec=outmat[,4];cvec=outmat[,1]
        maxvec=pmax(fnvec,fpvec)
        bound=cbind(cvec,maxvec)
        minofmax=min(maxvec)
        ymax=min(minofmax)*4
        
      

        if(isolate(input$penalty==TRUE))
        {
          ggplot()+geom_blank()
        }
        else
        {
          ggplot(melty,aes(x=Cvalue,y=value,group=factor(fname),color=factor(fname)))+xlab('%Correct cut value')+labs(color='Error Type')+geom_line()+scale_color_manual(labels = c("RFP", "RFN"), values = c("blue", "red"))+ggtitle("Relative Errors")+coord_cartesian(ylim = c(0,ymax))       
        }
        
        
        #explaining everything will take a moment...
        #melty makes a long format,we could have made that ourselves by using rbind
        #group=factor(fname) tells r that we have a discrete category (fname)
        #color=factor(fname) tells r to give us different line colors and a discrete legend
        #scale_color_manual lets us change the labels in the legend
        
      }) # end conditional error plot
      
      
    #total error text
      output$value3=renderText({   #makes the text output
        
        inFile <- input$file1
        
        if (is.null(inFile))
          return(NULL)
        
        outmat=outmat2
        fnvec=outmat[,3]
        fpvec=outmat[,2]
        cvec=outmat[,1]
        
        totalp=fnvec+fpvec
        bound=cbind(cvec,totalp)
        minofmax=min(totalp)
        #minnie2=minnie[1]
        loc=which(bound[,2]==minofmax)
        loc2=min(loc) #neccessary becauase loc returns a vector in this case because there are two adjacent points with same error
        catmin=bound[loc2,1]
        
        threedig=convto3digit(catmin)
        
        if(isolate(input$penalty==FALSE))
        {
          textout="the optimal total error value is"
          
        }
        
        if(isolate(input$penalty==TRUE))
        {
          textout="the optimal penalty total error value is"
        }
        
       #paste(paste(textout, round(min(minofmax),digits=3),"{fp=",round(outmat[loc,3],digits=3)," fn=",round(outmat[loc,2],digits=3),"}"), paste("%C at min is",round(catmin,digits=3),',',threedig),sep='\n')
       paste(paste(textout, round(min(minofmax),digits=3),"{fp = ",round(outmat[loc2,3],digits=3),'fn =',round(outmat[loc2,2],digits=3), '}'), paste("%C at min is",round(catmin,digits=3),',',threedig),sep='\n')
       
       
       
        
      })  # end total error text
      
    #total error plot 
      output$totalerror=renderPlot({
        outmat=outmat2
        #sdobs=input$sdobs
        plotmat1=outmat[,1]
        plotmat2=outmat[,5]
        plotmat23=outmat[,4]
        plotmat3=cbind(plotmat1,plotmat2,plotmat23)
        

        fnvec=outmat[,3]
        fpvec=outmat[,2]
        cvec=outmat[,1]
        
        totalp=fnvec+fpvec
        bound=cbind(cvec,totalp)
        
        minofmax=min(totalp)
        #minnie2=minnie[1]
        loc=which(bound[,2]==minofmax)
        loc2=min(loc) #neccessary becauase loc returns a vector in this case because there are two adjacent points with same error
        catmin=bound[loc2,1]
        
        
        colnames(bound)=c('Cvalue','value')
        plottest4=as.data.frame(bound)
        
        
        if(isolate(input$penalty==FALSE))
        {
          lab='Total Errors'
        }
        
        if(isolate(input$penalty==TRUE))
        {
          lab='Total Penalty Errors'
        }
        
        

        ggplot(plottest4,aes(x=Cvalue,y=value))+geom_line()+ggtitle(lab)+xlab('%Correct cut value')+coord_cartesian(ylim = c(0,minofmax*4)) 
      })   # end total error text
      

      time2=Sys.time()
      timedif=difftime(time2,time1)
      unittime=units(timedif)
      paste('processing completed in',round(timedif,digits=3),unittime)
      })#end of with progress (loading bar)
      

    })#end observe (the button that made it all go)
      
          output$report <- downloadHandler(
            # For PDF output, change this to "report.pdf"
            filename = "report.html",
            content = function(file) {
              # Copy the report file to a temporary directory before processing it, in
              # case we don't have write permissions to the current working dir (which
              # can happen when deployed).
              tempReport <- file.path(tempdir(), "report.Rmd")
              file.copy("report.Rmd", tempReport, overwrite = TRUE)
              
              # Set up parameters to pass to Rmd document
              params <- list(n = input$slider)
              
              # Knit the document, passing in the `params` list, and eval it in a
              # child of the global environment (this isolates the code in the document
              # from the code in this app).
              rmarkdown::render(tempReport, output_file = file,
                                params = params,
                                envir = new.env(parent = globalenv())
              )
            }
          )
          
          
  session$onSessionEnded(function() {   #this code closes r when the user exits the browser window, crucial for rinno version
    stopApp()
  })
      
      
}#end server
        
        
  

        


shinyApp(ui = ui, server = server)
