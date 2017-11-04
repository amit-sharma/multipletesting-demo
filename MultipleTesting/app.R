#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(MASS)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Demystifying Multiple Testing: Calculating FDR, FNR and Power"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("t",
                     "Nominal p value threshold for each test:",
                     min = 0,
                     max = 1,
                     value = 0.05),
         numericInput(inputId = "m",
                      label = "Number of hypothesis tests:",
                      value = 100),
         numericInput(inputId="small_n",
                      label="Number of samples in each test",
                      value=15),
         sliderInput("pi",
                     "Actual fraction of null hypotheses out of all hypotheses tested",
                     min=0,
                     max=1,
                      value=0.5),
         selectInput(inputId = "null_data_gen_distr",
                     label = "Choose a distribution for generating data for the null hypothesis:",
                     choices = c("Normal", "Log-normal", "fisher")),
         numericInput(inputId="mean_null",
                      label="Mean of the null gaussian",
                      value=0),
         selectInput(inputId = "alter_data_gen_distr",
                     label =  "Choose a distribution for generating data for the alternative hypothesis:",
                     choices = c("Normal", "Log-normal", "fisher")),
         numericInput(inputId="mean_alter",
                      label="Mean of the alternative gaussian",
                      value=0.1),
         selectInput("range_yd",
                     "Choose the number of distinct empirical values of Yd",
                     choices=c(1,2,3,4),
                     selected=2)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot"),
         verbatimTextOutput("fdr")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  getData <- reactive({
    var_null_distr=1
    var_alter_distr=1
    sd_null_distr=sqrt(var_null_distr)
    sd_alter_distr=sqrt(var_alter_distr)
    mean_true_null_distr=input$mean_null
    
    if(input$null_data_gen_distr=="Normal"){
      null_data  <- replicate(input$pi*input$m, rnorm(1000, mean=input$mean_null, sd=sd_null_distr))
    } else if (input$null_data_gen_distr=="Log-normal"){
      null_data_temp  <- replicate(input$pi*input$m, rnorm(1000, mean=input$mean_null, sd=sd_null_distr))
      null_data <- exp(null_data_temp)
      mean_true_null_distr=exp(var_null_distr/2)
    } else if (input$null_data_gen_distr=="fisher"){
      tt=15
      data1 = replicate(input$pi*input$m, c(sample(1:4, input$small_n, replace=T))) 
      data2= replicate(input$pi*input$m, c(sample(1:input$range_yd, input$small_n, replace=T))) 
    }
    
    if(input$alter_data_gen_distr=="Normal"){
      alter_data <- replicate(input$m - input$pi*input$m, rnorm(1000, mean=input$mean_alter, sd=sd_alter_distr))
    } else if (input$alter_data_gen_distr=="Log-normal"){
      alter_data_temp <- replicate(input$m - input$pi*input$m, rnorm(1000, mean=input$mean_alter, sd=sd_alter_distr))
      alter_data = exp(alter_data_temp)
    } else if (input$alter_data_gen_distr=="fisher"){
      alter_data1 = replicate(input$m-input$pi*input$m, sample(1:4, input$small_n, replace=T))
      noise= replicate(input$m-input$pi*input$m, round(rnorm(input$small_n, mean=0,sd=2)))
      data2_temp = alter_data1+noise
      alter_data2 = ifelse(data2_temp>4, 4, ifelse(data2_temp<1, 1, data2_temp))
    }
    if(input$null_data_gen_distr=="fisher"){ 
      j=0
      test_ret_null=apply(data1,2, function(x) {j<<-j+1; fisher.test(x, data2[,j]) })
      print(head(test_ret_null))
      j=0
      test_ret_alter=apply(alter_data1,2, function(x) {j<<-j+1;fisher.test(x, alter_data2[,j])})
    } else if (input$null_data_gen_distr=="Normal"){
      test_ret_null=apply(null_data,2, t.test, mu=mean_true_null_distr)
      test_ret_alter = apply(alter_data, 2, t.test, mu=mean_true_null_distr)
    }
    
    pvals_null = sapply(test_ret_null, function(x) x$p.value)
    pvals_alter = sapply(test_ret_alter, function(x) x$p.value) 
    pvals_all = c(pvals_null, pvals_alter)
    
    list(pvals_null=pvals_null, pvals_alter=pvals_alter, pvals_all=pvals_all)
  })
  
  output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
        
        plot_df = rbind(
          #data.frame(hyp="all", pval=getData()$pvals_all),
          data.frame(hyp="null", pval=getData()$pvals_null),
          data.frame(hyp="alter", pval=getData()$pvals_alter)
        ) %>%
          mutate(hyp=factor(hyp, levels=c("alter", "null")))
        ggplot(plot_df, aes(x=pval, color=hyp,fill=hyp)) + geom_histogram(alpha=0.2)+
          geom_vline(xintercept=input$t, linetype="dotted")
   })
   
   output$fdr <- renderPrint({
    all_pvals=getData()$pvals_all
    P = length(all_pvals[all_pvals<input$t])
    N = input$m - P
    max_fdr = min(1, input$t*input$m/P)
    max_fnr = min(1, (1-input$t)*input$m/N)
    
    min_power=P/input$m
    min_powerdash=N/input$m
    paste("MAX FDR", max_fdr, "MAX FNR", max_fnr, "MIN POWER", min_power, "MIN POWER'", min_powerdash, P, N)
    
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

