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
require(eulerr)

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
            sliderInput("corr", label = "Correlation of X1 and X2",value = .1, min = 0,max = 1),
            numericInput(inputId="beta1", label = "X1 Coefficient", value = 1),
            numericInput(inputId="beta2", label = "X2 Coefficient", value = 1),
            numericInput(inputId="sigma2", label = "Error Variance", value = 1, min=0),
            selectInput(inputId="whichx", label="Include in Model", choices=c("Include X1 Only",
                                                                              "Include X2 Only",
                                                                              "Include X1 and X2")),
            h1(),
            fluidRow(column(3, verbatimTextOutput("value"))),
            
            helpText("X1 and X2 are randomly generated from a normal distribution"),
            helpText("Model: Y=1+beta1*X1+beta2*X2+error"),
            h1(),
            
            hr(),
            
            actionButton("resample", "Sample"),
            actionButton("sampmany", "Sample Many"),
            helpText("Click 'Sample Many' to see beta distributions across 100 samples")
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            h4("Regression Summary"),
            tableOutput("regsummary"),
            h4("Model Sum of Squares"),
            plotOutput("venndiag"),
            h4("Distribution of Betas"),
            plotOutput("betasdist")
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
  
  calcMSS = function(lmobj){
    r2 <- summary(lmobj)$r.squared
    res <- lmobj$residuals
    RSS <- sum(res^2)
    TSS <- (RSS)/(-(r2-1))
    MSS <- TSS - RSS
    return(MSS)
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
      
      
      output$venndiag <- renderPlot({
        if(input$whichx=="Include X1 Only"){
          lm1 = lm(Y~X1, data=simdata)
          MSS = calcMSS(lm1)
          fit <- euler(c(MSS = MSS))
          plot(fit, lab=paste("MSS =",round(MSS)))
        } else if(input$whichx=="Include X2 Only"){
          lm1 = lm(Y~X2, data=simdata)
          MSS = calcMSS(lm1)
          fit <- euler(c(A = MSS))
          plot(fit, lab=paste("MSS =",round(MSS)))
        } else{
          lmx1x2 = lm(Y~X1+X2, data=simdata)
          lmx1 = lm(Y~X1, data=simdata)
          lmx2 = lm(Y~X2, data=simdata)
          MSSx1 = calcMSS(lmx1)
          MSSx2 = calcMSS(lmx2)
          MSSx1x2 = calcMSS(lmx1x2)
          middle = MSSx1 + MSSx2 - MSSx1x2
          fit <- euler(c(A = MSSx1+.05, B = MSSx2+.05, "A&B" = round(middle+.0001)))
          plot(fit, labels=c(paste("MSS X1 =",round(MSSx1)),
                          paste("MSS X2 =",round(MSSx2)),
                          paste("MSS X1 & X2 =",round(MSSx1x2))))
        }
      })
      
    })
      
      observeEvent(input$sampmany, {
        output$betasdist <- renderPlot({
          if(input$whichx=="Include X1 Only"){
            estb1 = rep(NA, 100)
            for(i in 1:100){
              simdata2 = simdatafunc()
              lm1 = lm(Y~X1, data=simdata2)
              estb1[i] = summary(lm1)$coefficients[2,1]
            }
            hist(estb1, xlab="beta1", main="")
            abline(v=mean(estb1), col=2)
            abline(v=input$beta1, col=4)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1))
          } else if(input$whichx=="Include X2 Only"){
            estb2 = rep(NA, 100)
            for(i in 1:100){
              simdata2 = simdatafunc()
              lm1 = lm(Y~X2, data=simdata2)
              estb2[i] = summary(lm1)$coefficients[2,1]
            }
            hist(estb2, xlab="beta2", main="")
            abline(v=mean(estb2), col=2)
            abline(v=input$beta2, col=4)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1))
          } else{
            estb1 = rep(NA, 100)
            estb2 = rep(NA, 100)
            for(i in 1:100){
              simdata2 = simdatafunc()
              lm1 = lm(Y~X1+X2, data=simdata2)
              estb1[i] = summary(lm1)$coefficients[2,1]
              estb2[i] = summary(lm1)$coefficients[3,1]
            }
            par(mfrow=c(1,2))
            hist(estb1, xlab="beta1", main="")
            abline(v=mean(estb1), col=2)
            abline(v=input$beta1, col=4)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1), cex=.7)
            hist(estb2, xlab="beta2", main="")
            abline(v=mean(estb2), col=2)
            abline(v=input$beta2, col=4)
            legend("topleft", c("true value", "sim mean"), col=c(4,2), lty=c(1,1), cex=.7)
          }
        })
      })
    
}

# Create Shiny app ----
shinyApp(ui, server)
