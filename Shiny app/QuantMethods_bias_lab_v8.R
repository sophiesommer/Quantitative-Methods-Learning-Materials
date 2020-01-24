#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#download/load shiny
if(!require(shiny)) {
    install.packages("shiny")
}
library(shiny)
require(MASS)
require(xtable)

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Regression Simulation"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Choose the sample size, betas, error variance ----
      numericInput(inputId="sampsize", label ="Sample Size", value = 10, min=5),
      sliderInput("corr", label = "Correlation of X1 and X2",value = .1, min = -1,max = 1, step=.05),
      numericInput(inputId="beta1", label = "X1 Coefficient", value = 1),
      numericInput(inputId="beta2", label = "X2 Coefficient", value = 1),
      numericInput(inputId="sigma2", label = "Error Variance", value = 1, min=0),
      selectInput(inputId="whichx", label="Include in Model", choices=c("Include X1 Only",
                                                                        "Include X2 Only",
                                                                        "Include X1 and X2")),
      fluidRow(column(3, verbatimTextOutput("value"))),
      actionButton("resample", "Sample"),
      actionButton("sampmany", "Sample Many"),
      helpText("Model: Y=1+beta1*X1+beta2*X2+error")
      
    ),
    
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Sample Once", 
                           h4("Regression Summary"),
                           tableOutput("regsummary"),
                           h4("Model Sum of Squares"),
                           splitLayout(cellWidths = c("40%", "60%"), 
                                       tableOutput("msstab"), plotOutput("piechart"))
                           ),
                  tabPanel("Sample Many",
                           h4("Simulation Summary"),
                           tableOutput("simsummary"),
                           h4("Distribution of Betas"),
                           plotOutput("betasdist"))
      )
    )
  )
)
  
  
  

# Define server logic to summarize and view selected dataset ----
server = function(input, output) {
    
  Simdata <- eventReactive(input$resample, {
    data = mvrnorm(n=input$sampsize, mu=c(0, 0),
                   Sigma=matrix(c(1, input$corr, input$corr, 1), nrow=2),
                   empirical=TRUE)
    X1 = data[, 1]
    X2 = data[, 2]
    Y = 1 + input$beta1 * X1 + input$beta2 * X2 + rnorm(input$sampsize, 0, input$sigma2^.5)
    return(data.frame(X1 = X1,X2 = X2,Y = Y))
  }
  )
  
  simdatafunc = function(){
    data = mvrnorm(n=input$sampsize, mu=c(0, 0),
                   Sigma=matrix(c(1, input$corr, input$corr, 1), nrow=2),
                   empirical=TRUE)
    X1 = data[, 1]
    X2 = data[, 2]
    Y = 1 + input$beta1 * X1 + input$beta2 * X2 + rnorm(input$sampsize, 0, input$sigma2^.5)
    return(data.frame(X1 = X1,X2 = X2,Y = Y))
  }
  
  calcSS = function(lmobj){
    r2 <- summary(lmobj)$r.squared
    res <- lmobj$residuals
    RSS <- sum(res^2)
    TSS <- (RSS)/(-(r2-1))
    MSS <- TSS - RSS
    return(data.frame(MSS=MSS, TSS=TSS, RSS=RSS))
  }
  
      observeEvent(input$resample, {
      simdata <- Simdata()
      output$regsummary <- renderTable({
        if(input$whichx=="Include X1 Only"){
          lm1 = lm(Y~X1, data=simdata)
          truevals = c(1, input$beta1)
        } else if(input$whichx=="Include X2 Only"){
          lm1 = lm(Y~X2, data=simdata)
          truevals = c(1, input$beta2)
        } else{
          lm1 = lm(Y~X1+X2, data=simdata)
          truevals = c(1, input$beta1, input$beta2)
        }
        table=summary(lm1)$coef
        table=cbind(truevals, table)
        colnames(table)[1] = "True Value"
        xtable(table)},
        rownames=TRUE)
      
      
      output$piechart <- renderPlot({
        if(input$whichx=="Include X1 Only"){
          lm1 = lm(Y~X1, data=simdata)
          MSS = calcSS(lm1)$MSS
          RSS = calcSS(lm1)$RSS
          slices <- c(MSS, RSS)
          lbls <- c("MSS","RSS")
          pie(slices, labels = lbls, main="TSS Breakdown", col=c("blue", "red"), radius=1)
          
        } else if(input$whichx=="Include X2 Only"){
          lm1 = lm(Y~X2, data=simdata)
          MSS = calcSS(lm1)$MSS
          RSS = calcSS(lm1)$RSS
          slices <- c(MSS, RSS)
          lbls <- c("MSS","RSS")
          pie(slices, labels = lbls, main="TSS Breakdown", col=c("blue", "red"), radius=1)
        } else{
          lmx1x2 = lm(Y~X1+X2, data=simdata)
          lmx1 = lm(Y~X1, data=simdata)
          lmx2 = lm(Y~X2, data=simdata)
          MSSx1 = calcSS(lmx1)$MSS
          MSSx2 = calcSS(lmx2)$MSS
          MSSx1x2 = calcSS(lmx1x2)$MSS
          RSSx1x2 = calcSS(lmx1x2)$RSS
          TSS = ceiling(calcSS(lmx1x2)$TSS) 
          middle = MSSx1 + MSSx2 - MSSx1x2
          if(middle > 0){
            slices <- c((MSSx1 - middle), RSSx1x2, (MSSx2 - middle), middle)
            lbls <- c("MSS X1","RSS", "MSS X2", "MSS X1 & X2")
            pie(slices, labels = lbls, main="TSS Breakdown", 
                col=c("blue", "red", "white", "light blue"), clockwise = TRUE)
          } else {
            middle = abs(middle)
            slices <- c(MSSx1, RSSx1x2, MSSx2, middle)
            lbls <- c("MSS X1","RSS", "MSS X2", "MSS X1 & X2")
            pie(slices, labels = lbls, main="TSS Breakdown", clockwise = TRUE,
                col=c("light green", "red", "light blue", "white"))
          }
        
        }
      })
      
      output$msstab <- renderTable({
        lmx1x2 = lm(Y~X1+X2, data=simdata)
        lmx1 = lm(Y~X1, data=simdata)
        lmx2 = lm(Y~X2, data=simdata)
        MSSx1 = calcSS(lmx1)$MSS
        MSSx2 = calcSS(lmx2)$MSS
        MSSx1x2 = calcSS(lmx1x2)$MSS
        table=rbind(c(MSSx1),c(MSSx2), c(MSSx1 + MSSx2), c(MSSx1x2))
        colnames(table) = c("MSS")
        rownames(table) = c("x1 only", "x2 only", 
                            "x1 only + x2 only", "model with x1 and x2")
        xtable(table)},
        rownames=TRUE)
      
    })
      
      observeEvent(input$sampmany, {
        if(input$whichx=="Include X1 Only"){
          estb1 = rep(NA, 500)
          for(i in 1:500){
            simdata2 = simdatafunc()
            lm1 = lm(Y~X1, data=simdata2)
            estb1[i] = summary(lm1)$coefficients[2,1]
          }
          minx = min(c(estb1-.1, input$beta1-.1))
          maxx = max(c(estb1+.1, input$beta1+.1))
          output$betasdist <- renderPlot({
            hist(estb1, xlab="beta1", main="", xlim=c(minx,maxx), breaks=20)
            abline(v=mean(estb1), col=2, lwd=2)
            abline(v=input$beta1, col=4, lwd=2)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1), lwd=2)
          })
          output$simsummary <- renderTable({
            table=cbind(c(input$beta1),c(mean(estb1)), c(sd(estb1)))
            colnames(table) = c("True Value", "Sim Mean", "Sim SD")
            rownames(table) = c("Beta 1")
            xtable(table)},
            rownames=TRUE)
          
        } else if(input$whichx=="Include X2 Only"){
          estb2 = rep(NA, 500)
          for(i in 1:500){
            simdata2 = simdatafunc()
            lm1 = lm(Y~X2, data=simdata2)
            estb2[i] = summary(lm1)$coefficients[2,1]
          }
          minx = min(c(estb2-.1, input$beta2-.1))
          maxx = max(c(estb2+.1, input$beta2+.1))
          
          output$betasdist <- renderPlot({
            hist(estb2, xlab="beta2", main="", xlim=c(minx,maxx), breaks=20)
            abline(v=mean(estb2), col=2, lwd=2)
            abline(v=input$beta2, col=4, lwd=2)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1), lwd=2)
          })
          
          output$simsummary <- renderTable({
            table=cbind(c(input$beta2),c(mean(estb2)), c(sd(estb2)))
            colnames(table) = c("True Value", "Sim Mean", "Sim SD")
            rownames(table) = c("Beta 2")
            xtable(table)},
            rownames=TRUE)
          
        } else{
          estb1 = rep(NA, 500)
          estb2 = rep(NA, 500)
          for(i in 1:500){
            simdata2 = simdatafunc()
            lm1 = lm(Y~X1+X2, data=simdata2)
            estb1[i] = summary(lm1)$coefficients[2,1]
            estb2[i] = summary(lm1)$coefficients[3,1]
          }
          minx1 = min(c(estb1-.1, input$beta1-.1))
          maxx1 = max(c(estb1+.1, input$beta1+.1))
          minx2 = min(c(estb2-.1, input$beta2-.1))
          maxx2 = max(c(estb2+.1, input$beta2+.1))
          
          output$betasdist <- renderPlot({
            par(mfrow=c(1,2))
            hist(estb1, xlab="beta1", main="", xlim=c(minx1,maxx1), breaks=20)
            abline(v=mean(estb1), col=2, lwd=2)
            abline(v=input$beta1, col=4, lwd=2)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1), cex=.7, lwd=2)
            hist(estb2, xlab="beta2", main="", xlim=c(minx2,maxx2), breaks=20)
            abline(v=mean(estb2), col=2, lwd=2)
            abline(v=input$beta2, col=4, lwd=2)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1), cex=.7, lwd=2)
          })
          output$simsummary <- renderTable({
            table=cbind(c(input$beta1, input$beta2),
                        c(mean(estb1), mean(estb2)), 
                        c(sd(estb1), sd(estb2)))
            colnames(table) = c("True Value", "Sim Mean", "Sim SD")
            rownames(table) = c("Beta 1", "Beta 2")
            xtable(table)},
            rownames=TRUE)
          
        }
        
      })
    
}

# Create Shiny app ----
shinyApp(ui, server)




