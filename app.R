library(shiny)
library(DT)
library(tidyverse)
library(shinythemes)
library(broom)
library(cowplot)
library(colourpicker)
library(plotly)

floor_ceil <- function(a){
    choose <- sample(c(TRUE, FALSE), 1)
    if (choose) return(floor(a))
    if (!choose) return(ceiling(a))
}

# CURRENTLY ASSUMES SINGLE OFFSPRING PER MATING PAIR
herit_choices <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, "Random")
default_herit <- 0.1
min_spots <- 2
max_spots <- 14
theme_set(theme_cowplot())

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("flatly"),

    # Application title
    titlePanel("BIOL01104 Heritability Lab"),
    p("Written by", a("Stephanie J. Spielman, Ph.D.",href="https://spielmanlab.github.io"), "and released under GPL-3 License."),
    p("Source code:", a("https://github.com/sjspielman/biol01104_heritability_lab",href="https://github.com/sjspielman/biol01104_heritability_lab")), 
    p(strong("NOTE: CURRENTLY ASSUMES A SINGLE OFFSPRING PER MATING PAIR.")),
    
    sidebarLayout(
        sidebarPanel( width=3,
            selectInput("heritability",
                        "Choose true heritability:",
                        herit_choices, 
                        selected = default_herit),
            
            actionButton("new_mating",
                         "Generate new mating."),
            
           
            
            colourpicker::colourInput("point_color", "Point color:", value = "seagreen"),
            
           br(),br(),
           actionButton("reveal_herit",
                        "Reveal true heritability!"),
           br(),br(),
            actionButton("reset",
                         "Reset app.")            

        ),

        # Show a plot of the generated distribution
        mainPanel(
            span(textOutput("print_herit"), style="color:red;font-weight:600;font-size:20px"),
            
            #div(style="display:inline-block;vertical-align:top;",
            #    radioButtons("show_smooth", "Show regression line?", c("Yes", "No"), selected = "No")
            #),
            #div(style="display:inline-block;vertical-align:top;",
            #    radioButtons("show_se", "Show confidence interval?", c("Yes", "No"), selected = "No")
            #),
            #div(style="display:inline-block;vertical-align:top;",
            #    radioButtons("show_formula", "Show line formula?", c("Yes", "No"), selected = "No")
            #),
            #plotlyOutput("parent_offspring_plot", width = "800px", height = "500px"),
            fluidRow(
                column(9, plotlyOutput("parent_offspring_plot", width = "800px", height = "600px")),
                column(3, br(),br(),
                       radioButtons("show_smooth", "Show regression line?", c("Yes", "No"), selected = "No"),
                            radioButtons("show_se", "Show confidence interval?", c("Yes", "No"), selected = "No"),
                            radioButtons("show_formula", "Show line formula?", c("Yes", "No"), selected = "No"))
            ),
            br(),
            dataTableOutput("parent_offspring_table")
        )
    )
)

server <- function(input, output, session) {

    
    target_intercept <- sample(3:7, 1)
    
    true_herit <- reactive({
        if (input$heritability == "Random"){
            return( round(runif(1), 2) )
        } else {
            return(as.numeric(input$heritability))
        }
        
        
    })
    mating_replicate <- reactiveVal(0)
    herit_table <- reactiveVal(tibble("Mating replicate" = numeric(),
                       "Parent 1"= numeric(),
                       "Parent 2" = numeric(),
                       "Mid-parent Value" = numeric(),
                       "Mid-offspring Value" = numeric()))  
    
    # Reset app ----------------------------------------------------
    observeEvent(input$reset,{
        updateSelectInput(session, "heritability", 
                          "Choose true heritability:",
                          herit_choices, 
                          selected = default_herit)
        mating_replicate(0)
        herit_table(tibble("Mating replicate" = numeric(),
                           "Parent 1"= numeric(),
                           "Parent 2" = numeric(),
                           "Mid-parent Value" = numeric(),
                           "Mid-offspring Value" = numeric()))  
        output$print_herit <- renderText({""})
    })
    

    observeEvent(input$new_mating, {
        last_mating <- mating_replicate()
        mating_replicate(last_mating + 1)
        p1_draw <- floor_ceil(sample(min_spots:max_spots, 1))
        p2_draw <- floor_ceil(sample(min_spots:max_spots, 1))
        
        mean_parent <- mean(c(p1_draw, p2_draw))
        # ASSUMES A SINGLE OFFSPRING
        offspring <- floor_ceil( mean_parent * true_herit() ) + target_intercept
        if (offspring <= 0) offspring <- 1
        
        herit_table(bind_rows(herit_table(),
                              tibble("Mating replicate" = mating_replicate(),
                                     "Parent 1" = p1_draw,
                                     "Parent 2" = p2_draw,
                                     "Mid-parent Value" = mean_parent,
                                     "Mid-offspring Value" = offspring))
        )
    })
    

    observeEvent(input$reveal_herit, {
        output$print_herit <- renderText({paste0("The true heritability is ", true_herit())})
    })
    
    
    output$parent_offspring_table <- renderDataTable({
       herit_table() 
            
    })
        
        
    output$parent_offspring_plot <- renderPlotly({
        req(nrow(herit_table()) >= 1)
        ggplot(herit_table(), aes(x = `Mid-parent Value`, y = `Mid-offspring Value`)) + 
            geom_point(size = 3, pch = 21, fill = input$point_color) +
            scale_x_continuous(limits=c(min_spots-1, max_spots+1), breaks=seq(min_spots-1, max_spots +1, 2)) + 
            scale_y_continuous(limits=c(min_spots-1, max_spots+1), breaks=seq(min_spots-1, max_spots +1, 2)) + 
            # NOTE: FOR CONSISTENT FORMATTING WE ADD THE TITLE VIA PLOTLY
            xlab("Mid-parent trait value") + 
            ylab("Mid-offspring trait value") -> herit_plot
        
        if (input$show_smooth == "Yes"){
            if (input$show_se == "Yes"){
                herit_plot <- herit_plot + geom_smooth(method = "lm")
            } else {
                herit_plot <- herit_plot + geom_smooth(method = "lm", se=FALSE)
            }
        }
        if ((input$show_formula == "Yes") & (nrow(herit_table()) >= 2)){
            
            tidy(lm(`Mid-offspring Value` ~ `Mid-parent Value`, data = herit_table())) %>%
                pull(estimate) -> coeffs
            
            intercept <- round(coeffs[1], 2)
            slope <- round(coeffs[2], 2)
            if (intercept < 0){
                form_string <- paste0("Y = ", slope, "X - ", abs(intercept))
            } else {
                form_string <- paste0("Y = ", slope, "X + ", intercept)
            }
            
            
            ggplotly(herit_plot) %>%
                layout(title = list(text = paste0("Midparent-midoffspring regression",
                                                  '<br>',
                                                  '<sup>',
                                                  form_string,
                                                  '</sup>')),
                       margin = list(t = 90))
        } else{
            ggplotly(herit_plot) %>%
                layout(title = list(text = paste0("Midparent-midoffspring regression")),
                       margin = list(t = 90))
        }
        
        
    })
    
}







# Run the application 
shinyApp(ui = ui, server = server)
