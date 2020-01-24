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
            numericInput(inputId="sampsize", label ="Sample Size", value = 10, min=2),
            numericInput(inputId="corr", label ="Correlation of X1 and X2", value = 0, min=0, max=1),
            numericInput(inputId="beta1", label = "X1 Coefficient", value = 0),
            numericInput(inputId="beta2", label = "X2 Coefficient", value = 0),
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
            
            actionButton("resample", "Sample")
            
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            h4("Regression Summary"),
            tableOutput("regsummary"),
            h4("Distribution of Betas"),
            plotOutput("betasdist"),
            h4("Model Sum of Squares"),
            plotOutput("venndiag")
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
      
      output$betasdist <- renderPlot({
        if(input$whichx=="Include X1 Only"){
          lm1 = lm(Y~X1, data=simdata)
          SE = summary(lm1)$coefficients[2,2]
          beta1 = summary(lm1)$coefficients[2,1]
          xmin = beta1 - 3*SE
          xmax = beta1 + 3*SE
          curve(dnorm(x,beta1,SE), xlim=c(xmin,xmax), ylab="", xlab="beta1")
          abline(v=input$beta1, col=5)
        } else if(input$whichx=="Include X2 Only"){
          lm1 = lm(Y~X2, data=simdata)
          SE = summary(lm1)$coefficients[2,2]
          beta2 = summary(lm1)$coefficients[2,1]
          xmin = beta2 - 3*SE
          xmax = beta2 + 3*SE
          curve(dnorm(x,beta2,SE), xlim=c(xmin,xmax), ylab="", xlab="beta2")
          abline(v=input$beta2, col=5)
        } else{
          lm1 = lm(Y~X1+X2, data=simdata)
          SE = summary(lm1)$coefficients[2,2]
          beta1 = summary(lm1)$coefficients[2,1]
          beta2 = summary(lm1)$coefficients[3,1]
          xmin1 = beta1 - 3*SE
          xmax1 = beta1 + 3*SE
          xmin2 = beta2 - 3*SE
          xmax2 = beta2 + 3*SE
          par(mfrow=c(1,2))
          curve(dnorm(x,beta1,SE), xlim=c(xmin1,xmax1), ylab="", xlab="beta1")
          abline(v=input$beta1, col=5)
          curve(dnorm(x,beta2,SE), xlim=c(xmin2,xmax2), ylab="", xlab="beta2")
          abline(v=input$beta2, col=5)
        }
        })
      
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
          fit <- euler(c(A = MSSx1, B = MSSx2, "A&B" = middle))
          plot(fit, labels=c(paste("MSS X1 =",round(MSSx1)),
                          paste("MSS X2 =",round(MSSx2)),
                          paste("MSS X1 & X2 =",round(MSSx1x2))))
        }
      })
      
    })
    
}

# Create Shiny app ----
shinyApp(ui, server)
