library(shiny)
library(magrittr)
library(ggplot2)
library(ls)
library(gridExtra)

ui = fluidPage(
  titlePanel("Latent Stratification Calculator"),
  tags$a(href="https://arxiv.org/abs/1911.08438", "Latent Stratification in Advertising Experiments"), 
  br(),
  "By",
  tags$a(href="https://www.ron-berman.com/", "Ron Berman"), 
  "and",
  tags$a(href="https://www.lebow.drexel.edu/people/eleafeit", "Elea McDonnell Feit"),
  br(),
  "Citation (APA): Berman, R., & Feit, E. M. (2019). Principal stratification for advertising experiments. 
   The Wharton School Research Paper, Wharton Customer Analytics Initiative Research Paper.",
  br(),
  "Developed by Zhen Huang",
  
  tabsetPanel(id = "tabsetpanel",
    tabPanel("Upload Data",
             fluidRow(
               column(
                 12, align = "center", fileInput("file1", "Upload Data", accept = c(".csv",".xlsx", ".xls")), style="font-size:300%"
               )
             ),
             
             fluidRow(
               column(
                 12, align = "center", actionButton("analyze", "Analyze Data", style="font-size:300%")
               )
             ),
             
             br(),
             
             fluidRow(
               column(
                 12, align = "center", downloadButton("downloadData", "Download Sample Data", style="font-size:300%")
               )
             ),
             tags$h2("Instructions"),
             tags$h4("1. To use the calculator, upload a csv or excel file containing an outcome column named y and a
             treatment column named z. The treatment column z must only consist of 0 and 1 (1 = received treatment, 
             0 = no treatment)."),
             # br(),
             tags$h4("2. To see an example of the data that is accepted by the calculator, click the Download button."),
             # br(),
             tags$h4("3. After uploading the file, click the Analyze button to view the ATE estimates."),
             # br(),
             tags$h4("4. To see an example of the analysis, download and upload the sample data to the 
             calculator, then click the Analyze button.")
    ),
    tabPanel("Results",
             fluidRow(
               column(12, align = "center", uiOutput("UI")
                      ),
             actionButton("back", "Back")
             )
    )
  ),
  
  br(),
  br(),
  
  tags$footer("Copyright 2021 Elea McDonnell Feit, Ron Berman, and Zhen Huang"),
  tags$footer("Licensed under the Apache License, Version 2.0 (the 'License'); 
              you may not use this software except in compliance with the License."),
  tags$footer("Unless required by applicable law or agreed to in writing, 
              software distributed under the License is distributed on an 'AS IS' BASIS, 
              WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
              See the License for the specific language governing permissions and limitations under the License.")
  
)

server = function(input, output, session) {
  
  data = reactive({
    req(input$file1)
    path = (input$file1)$datapath
    ext = tools::file_ext(path)
    
    if (ext == "csv") file = read.csv(path)
    else if (ext == "xlsx" || ext == "xls") file = readxl::read_excel(path)
    else file = NULL
    
    return(file)
  })
  
  observeEvent(input$analyze, {
    updateTabsetPanel(session, "tabsetpanel",
                      selected = "Results")
  })
  
  observeEvent(input$back, {
    updateTabsetPanel(session, "tabsetpanel",
                      selected = "Upload Data")
  })
  
  ##########################################################################################
  
  set.seed(20030601)
  sampleData = ls::sim_latent_strat(n=10000, piA=0.2, piB=0.1, muA1=5, muA0=4.5, muB1=3, sigma=0.3)$data[,1:2]
  
  mlePars = reactive({
    ls::mle_ls(data = data())$pars
  })
  ttest = reactive({
    t.test(data()$y[data()$z==1], data()$y[data()$z==0])
  })
  
  testPars = reactive({
    mleATE = mlePars()[1, 2]
    mleSE = mlePars()[1, 3]
    mleCI = c(mleATE - 1.96 * mleSE, mleATE + 1.96 * mleSE)
    mlepValue = 2*pnorm(q=mleATE / mleSE, lower.tail=FALSE)
    mle.result = c("Latent Stratification", round(c(mleATE, mleSE, mleCI, mlepValue),3))
    
    
    ttestATE = unname(ttest()$estimate[1] - ttest()$estimate[2])
    ttestSE = ttest()$stderr
    ttestCI = ttest()$conf.int[1:2]
    ttestpValue = 2*pnorm(q=ttestATE / ttestSE, lower.tail=FALSE)
    ttest.result = c("Diff in Means t-test", round(c(ttestATE, ttestSE, ttestCI, ttestpValue),3))
    # add t test p value
    
    df = data.frame(t(cbind(ttest.result, mle.result))) %>%
      dplyr::rename(Method=X1, ATE=X2, SE=X3, "95% CI Lower"=X4, "95% CI Upper"=X5, "p-value"=X6)
  })
  
  CIgraph = reactive({
    ggplot(testPars(), aes(x = Method, y = as.numeric(ATE))) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin=as.numeric(testPars()[,4]), ymax=as.numeric(testPars()[,5])), width = 0.3) +
      ylab("ATE Estimate Value") + xlab("Method to Estimate ATE") + coord_fixed(ratio=8) +
      theme(axis.text=element_text(size=15), axis.title=element_text(size=20,face="bold")) +
      ylim(as.numeric(testPars()[1,2]) - as.numeric(testPars()[1,3])*3, 
           as.numeric(testPars()[1,2]) + as.numeric(testPars()[1,3])*3)
  })
  
  
  mleExtra = reactive({
    table = mlePars()[-1,]
    table["par"] = c("Exponential ATE", "Proportion of always buyer", "Proportion of influenced buyer",
    "mean of always buyer under treatment", "mean of always buyer under control", 
    "mean of influenced buyer under treatment", "sigma")
    table
  })
  ##########################################################################################
  
  output$downloadData = downloadHandler(
    filename = function() {
      paste("ls-sample-data.csv", sep="")
    },
    content = function(file) {
      write.csv(sampleData, file)
    }
  )
  
  output$ATETable = renderTable(
    testPars()
  )
  
  output$CIplot = renderPlot(
    CIgraph()
  )
  
  output$extraVar = renderTable(
    mleExtra()
  )
  
  output$UI = renderUI({
    validate(need(!is.null(data()), "Please upload a csv or excel file"))
    validate(need(c("y","z") %in% colnames(data()), "Column Names y and z don't exist"))
    validate(
      need(nrow(data()[data()$y == 0 & data()$z == 0,]) >= 30, "z = 0, y = 0 samples are fewer than 30"),
      need(nrow(data()[data()$y == 0 & data()$z == 1,]) >= 30, "z = 1, y = 0 samples are fewer than 30"),
      need(nrow(data()[data()$y > 0 & data()$z == 0,]) >= 30, "z = 0, y > 0 samples are fewer than 30"),
      need(nrow(data()[data()$y > 0 & data()$z == 1,]) >= 30, "z = 1, y > 0 samples are fewer than 30"))
    
    tagList(
      tags$h2("Estimates of Average Treatment Effect", align="center"),
      tableOutput("ATETable"),
      tags$h2("Graph of The Confidence Interval of ATE", align="center"),
      plotOutput("CIplot"),
      tags$h2("Estimates of Extra Parameters Via Latent Stratification", align="center"),
      tableOutput("extraVar")
    )
  })
  

}

shinyApp(ui = ui, server = server)

# To download package
# devtools::install_github(" eleafeit/latent_strat ", ref="main", auth_token = "tokenstring")
# tokenstring is generated via Github