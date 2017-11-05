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
library(RColorBrewer)
library(MASS)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Demystifying Multiple Testing: Calculating FDR, FNR and Power"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("t",
                     "Nominal p-value threshold for each test",
                     min = 0,
                     max = 1,
                     value = 0.05),
         numericInput(inputId = "m",
                      label = "Number of hypothesis tests",
                      value = 1000),
         numericInput(inputId="small_n",
                      label="Number of samples in each test",
                      value=15),
         sliderInput(inputId="pi",
                     label="Actual fraction of null hypotheses out of all hypotheses tested",
                     min=0.1,
                     max=0.9,
                     value=0.5),
         selectInput(inputId = "null_data_gen_distr",
                     label = "Choose a distribution for generating data for the null hypothesis:",
                     choices = c("Normal", "Log-normal", "Poisson")),
         conditionalPanel(
           condition = "input.null_data_gen_distr=='Normal'",
           em("There is a single null distribution."),
           numericInput(inputId="mean_null",
                        label="Mean of the data distribution from null hypothesis",
                        value=0),
           numericInput(inputId="var_null",
                        label="Variance of the data distribution from null hypothesis",
                        value=1),
           #selectInput(inputId = "alter_data_gen_distr",
           #           label =  "Choose a distribution for generating data for the alternative hypothesis:",
           #             choices = c("Normal", "Log-normal", "Poisson")),
           em("There can be multiple alternative distributions."),
           numericInput(inputId="mean_alter",
                        label="Average Mean of the individual data distributions from alternative hypothesis",
                        value=1),
           numericInput(inputId="noise_means_alter",
                        label="Variance in the means from the alternative hypothesis",
                        value=1),
           numericInput(inputId="var_alter",
                        label="(Fixed) variance of each individual data distribution from alternative hypothesis",
                        value=1)
           
         ),
         checkboxInput(inputId="is_discretized",
                       label="Discretize data",
                       value=FALSE),
        conditionalPanel(
          condition="input.is_discretized == 1",
          numericInput(inputId="num_discrete_bins",
                       label="Number of bins for discretization",
                       value=4)
        ) 
         #selectInput("range_yd",
         #             "Choose the number of distinct empirical values of Yd",
         #           choices=c(1,2,3,4),
         #            selected=2)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         tabsetPanel(
           tabPanel("All p-values", 
                    plotOutput("distPlot"),
                    h3("FDR=False Positives/Positives = FP/P"),
                    plotOutput("distPlotFDR"),
                    verbatimTextOutput("compute_fdr")
                    ),
           tabPanel("Null and Alternative", 
                    plotOutput("distPlotSeparated"),
                    h3("FDR=False Positives/Positives = FP/P"),
                    plotOutput("distPlotFDR_withm0"),
                    verbatimTextOutput("compute_fdr_withm0")),
           tabPanel("FDR/FNR",verbatimTextOutput("fdr") )
         )
         
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  getData <- reactive({
    sd_null_distr=sqrt(input$var_null)
    sd_alter_distr=sqrt(input$var_alter)
    mean_true_null_distr=input$mean_null
    if(input$pi=="Unknown"){
      input_pi=runif(1)
    } else {
      input_pi = as.numeric(input$pi)
    }
    if(input$null_data_gen_distr=="Normal"){
      null_data  <- replicate(input_pi*input$m, rnorm(input$small_n, mean=input$mean_null, sd=sd_null_distr))
    } else if (input$null_data_gen_distr=="Log-normal"){
      null_data_temp  <- replicate(input_pi*input$m, rnorm(input$small_n, mean=input$mean_null, sd=sd_null_distr))
      null_data <- exp(null_data_temp)
      mean_true_null_distr=exp(var_null_distr/2)
    } else if (input$null_data_gen_distr=="Poisson"){
      null_data = replicate(input_pi*input$m, rpois(input$small_n, lambda=input$lambda)) #c(sample(1:4, input$small_n, replace=T))) 
      #data2= replicate(input_pi*input$m, rpois(1000, lambda=input$lambda))#c(sample(1:input$range_yd, input$small_n, replace=T))) 
    }
    
    if(input$null_data_gen_distr=="Normal"){
      alter_data <- replicate(input$m - input_pi*input$m, rnorm(input$small_n, mean=rnorm(1,mean=input$mean_alter,sd=input$noise_means_alter), sd=sd_alter_distr))
    } else if (input$alter_data_gen_distr=="Log-normal"){
      alter_data_temp <- replicate(input$m - input_pi*input$m, rnorm(input$small_n, mean=input$mean_alter, sd=sd_alter_distr))
      alter_data = exp(alter_data_temp)
    } else if (input$alter_data_gen_distr=="Poisson"){
      alter_data = replicate(input$m-input$pi*input$m, rpois(input$small_n, lambda=input$lambda))#sample(1:4, input$small_n, replace=T))
      #alter_data2 = ifelse(data2_temp>4, 4, ifelse(data2_temp<1, 1, data2_temp))
    }
    
    if(input$is_discretized){
      null_data = apply(null_data, 2, function(x){binwidth=(max(x)-min(x))/input$num_discrete_bins; min(x)+floor((x-min(x))/binwidth)*binwidth} )
      alter_data = apply(alter_data, 2, function(x){binwidth=(max(x)-min(x))/input$num_discrete_bins;min(x)+floor((x-min(x))/binwidth)*binwidth} )
    }
    if(input$null_data_gen_distr=="Poisson"){ 
      #j=0
      #test_ret_null=apply(data1,2, function(x) {j<<-j+1; fisher.test(x, data2[,j]) })
      #print(head(test_ret_null))
      #j=0
      test_ret_alter=apply(alter_data1,2, function(x) {j<<-j+1;fisher.test(x, alter_data2[,j])})
    } else if (input$null_data_gen_distr=="Normal"){
      NUM_SIMS=10000
      all_data = cbind(null_data, alter_data)
      bootstrap_vals=apply(all_data, 2, function(x) x - mean(x)+mean_true_null_distr)
      if(input$is_discretized){
        binwidth=(max(bootstrap_vals)-min(bootstrap_vals))/input$num_discrete_bins
        bootstrap_vals = min(bootstrap_vals)+floor((bootstrap_vals-min(bootstrap_vals))/binwidth)*binwidth
      }
      sim_diff_means = replicate(NUM_SIMS, 
                                 mean(sample(bootstrap_vals, input$small_n, replace=T))-
                                 mean_true_null_distr)
      
      sorted_sim_diff_means = sort(sim_diff_means)
      rev_sorted_sim_diff_means = rev(sorted_sim_diff_means)
      print(summary(sorted_sim_diff_means))
      obs_diff_means_null = apply(null_data, 2, function(x) mean(x)-mean_true_null_distr )
      print(summary(obs_diff_means_null))
      obs_diff_means_alter= apply(alter_data, 2, function(x) mean(x)-mean_true_null_distr)
      print(summary(obs_diff_means_alter))
      ranks_null=sapply(obs_diff_means_null, function(x) min(which(rev_sorted_sim_diff_means<x)[1],which(sorted_sim_diff_means>x)[1], na.rm=T))
      #ranks_null = ifelse(is.na(ranks_null), NUM_SIMS, ranks_null)
      
      #print(ranks_null)
      pvals_null=ranks_null/NUM_SIMS*2
      ranks_alter=sapply(obs_diff_means_alter, function(x) min(which(rev_sorted_sim_diff_means<x)[1], which(sorted_sim_diff_means>x)[1], na.rm=T))
      #ranks_alter = ifelse(is.na(ranks_alter), NUM_SIMS, ranks_alter)
      pvals_alter=ranks_alter/NUM_SIMS*2
 
      
      #test_ret_null=apply(null_data,2, t.test, mu=mean_true_null_distr)
      #test_ret_alter = apply(alter_data, 2, t.test, mu=mean_true_null_distr)
    }
    
    #pvals_null = sapply(test_ret_null, function(x) x$p.value)
    #pvals_alter = sapply(test_ret_alter, function(x) x$p.value) 
    pvals_all = c(pvals_null, pvals_alter)
    
    list(pvals_null=pvals_null, pvals_alter=pvals_alter, pvals_all=pvals_all)
  })
  
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    
    plot_df = data.frame(hyp="all", pval=getData()$pvals_all)
    ggplot(plot_df, aes(x=pval), color="grey") + geom_histogram(alpha=0.9, binwidth=0.05, center=0.025)+
      geom_vline(xintercept=input$t, linetype="dotted")
  })
  
  output$distPlotFDR <- renderPlot({
    # generate bins based on input$bins from ui.R
    pval_bindwidth=0.05
    num_bins = 1/pval_bindwidth
    max_num_null=input$m
    all_pvals=getData()$pvals_all
    positive_pvals=all_pvals[all_pvals<=input$t]
    other_pvals=all_pvals[all_pvals>input$t]
    plot_df = rbind(data.frame(hyp="Others", pval=other_pvals),
                    data.frame(hyp="Positives (P)", pval=positive_pvals))
    ggplot(plot_df, aes(x=pval, color=hyp, fill=hyp)) + 
      geom_histogram(alpha=0.4,binwidth=0.05, center=0.025)+
      geom_vline(xintercept=input$t, linetype="longdash", size=2) +
      geom_hline(yintercept=max_num_null/num_bins, linetype="longdash", size=2) + 
      annotate("text", x=0.4, y=max_num_null/num_bins*1.15, label="Upper bound on expected number of null hypotheses ")+
      annotate("text", x=0, y=max_num_null/num_bins*1.15, label="Max. FP")+
      scale_fill_brewer(palette="Dark2")+
      scale_colour_brewer(palette="Dark2")
  })
  
  output$distPlotFDR_withm0 <- renderPlot({
    # generate bins based on input$bins from ui.R
    pval_binwidth=0.05
    num_bins = 1/pval_binwidth
    num_null=input$pi*input$m
    all_pvals=getData()$pvals_all
    positive_pvals=all_pvals[all_pvals<=input$t]
    other_pvals=all_pvals[all_pvals>input$t]
    plot_df = rbind(data.frame(hyp="Others", pval=other_pvals),
                    data.frame(hyp="Positives (P)", pval=positive_pvals))
    ggplot(plot_df, aes(x=pval, color=hyp, fill=hyp)) + 
      geom_histogram(alpha=0.4,binwidth=pval_binwidth, center=0.025)+
      geom_vline(xintercept=input$t, linetype="longdash", size=2) +
      geom_hline(yintercept=num_null/num_bins, linetype="longdash", size=2) + 
      annotate("text", x=0.4, y=num_null/num_bins*1.15, label="Upper bound on expected number of null hypotheses ")+
      annotate("text", x=0, y=num_null/num_bins*1.15, label="Max. FP")+
      scale_fill_brewer(palette="Dark2")+
      scale_colour_brewer(palette="Dark2")
  })
  
  output$distPlotSeparated <- renderPlot({
      # generate bins based on input$bins from ui.R
        
        plot_df = rbind(
          #data.frame(hyp="all", pval=getData()$pvals_all),
          data.frame(hyp="null", pval=getData()$pvals_null),
          data.frame(hyp="alter", pval=getData()$pvals_alter)
        ) %>%
          mutate(hyp=factor(hyp, levels=c("alter", "null")))
        ggplot(plot_df, aes(x=pval, color=hyp,fill=hyp)) + geom_histogram(alpha=0.2, binwidth=0.05, center=0.025)+
          geom_vline(xintercept=input$t, linetype="dotted")
   })
  
   output$compute_fdr <- renderPrint({
    all_pvals=getData()$pvals_all
    P = length(all_pvals[all_pvals<=input$t])
    max_FP = input$m*input$t
    #N = input$m - P
    max_fdr = min(1, max_FP/P)
    print(paste("FDR (assuming all hypotheses are null hypotheses) <=", max_FP, "/", P, "=", max_fdr))
    
    #Now using Mosig's method
    num_bins=20
    m0=input$m
    pvals_df = data.frame(pvalue=all_pvals) %>%
      mutate(rounded_pvalue=round(pvalue*num_bins)/num_bins)
    mosig_counts=count(pvals_df, rounded_pvalue) %>% arrange(rounded_pvalue)
    
    for(i in 1:30){
      index=which(mosig_counts$n<=m0/num_bins)[1]
      P_index = sum(mosig_counts$n[1:(index-1)])
      diff=P_index-m0/num_bins*(index-1)
      m0 = input$m-diff
      #print(paste(index, P_index, diff, m0))
    }
    
    est_FP = m0*input$t
    #N = input$m - P
    max_fdr = min(1, est_FP/P)
    print(paste("FDR (using an estimate for number of null hypotheses) <=", est_FP, "/", P, "=", max_fdr))
    #return(max_fdr)
   })
   
   output$compute_fdr_withm0 <- renderPrint({
    all_pvals=getData()$pvals_all
    P = length(all_pvals[all_pvals<=input$t])
    expected_FP = input$pi*input$m*input$t
    N = input$m - P
    max_fdr = min(1, expected_FP/P)
    print(paste("FDR =", expected_FP, "/", P, "=", max_fdr))
    #return(max_fdr)
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

