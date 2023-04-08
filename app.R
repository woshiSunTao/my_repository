

library(shiny)
library(purrr)
library(tidymodels)

load("traindata.Rdata")
# table(traindata$Lung_mets)
load("XGB.Rdata")

datarecipe <- recipe(Lung_mets ~ ., traindata) %>%
  step_naomit(all_predictors(), skip = F) %>%
  step_dummy(all_nominal_predictors()) %>%
  prep()


names2 <- colnames(traindata)[-9]


newinput <- function(xname, xname2, data){
  datai <- data[, xname]
  if (is.numeric(datai)) {
    numericInput(inputId = xname, 
                 label = xname2, 
                 value = median(datai))
  } else {
    selectInput(inputId = xname, 
                label = xname2, 
                choices = levels(datai))
  }
}

# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("   "),
    hr(),
    sidebarLayout(
        sidebarPanel(
          map2(names2, 
               names2, 
               newinput, data = traindata)
        ),

        mainPanel(
          h2("XGBoost for distant metastasis"),
          h2("New Sample"),
          dataTableOutput("newdf"),
          hr(),
          h2("Prediction"),
          actionButton("pred", "Predict"),
          br(),
          br(),
          textOutput("newlab"),
          br(),
          textOutput("newprob")
        )
    )
)

# Define server logic
server <- function(input, output) {
  newsample <- reactive({
    newsamplei <- 
      data.frame(matrix(NA, 
                        nrow = 1, 
                        ncol = length(names2)))
    for (i in seq_along(names2)) {
      valuei <- input[[names2[i]]]
      newsamplei[1, i] <- valuei
    }
    names(newsamplei) <- names2
    newsamplei
  })
  
  output$newdf <- renderDataTable({
    newsample()
  }, options = list(scrollX = TRUE))
  
  observeEvent(input$pred, {
    newsamplei1 <- newsample()
    newsamplei1$Lung_mets <- "No"
    newsamplei2 <- rbind(traindata, newsamplei1)
    newsamplei3 <- bake(datarecipe, new_data = newsamplei2)
    
    predprob <- 
      round(predict(final_xgboost, 
                    newsamplei3, 
                    type = "prob")[nrow(newsamplei3), 2], 4)
    
    output$newprob <- renderText({
      
      paste("预测Lung_mets=Yes的概???:", predprob)
      
    })
    
    output$newlab <- renderText({
      
      paste("预测Lung_mets:", 
            ifelse(predprob > 0.5, "Yes", "No"))
      
    })
    
  })
  

  
}

# Run the application 
shinyApp(ui = ui, server = server)
