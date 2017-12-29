
############################################################
#This is the ui file for the Macroparasite model by Dobson and Hudson
#written by Andreas Handel and modified by Ania Majewska
#last updated 12/28/2017
############################################################

library(shiny)

#This is the UI part of the shiny App
ui <- fluidPage(
  
  
  
  
  div(align = "center"),
  h1('Macroparasite Model by Dobson and Hudson', align = "center", style = "background-color:#123c66; color:#fff"),
  
  #section to add buttons
  fluidRow(
    column(6,
           actionButton("submitBtn", "Run Simulation", class="submitbutton")  
    ),
    column(6,
           actionButton("exitBtn", "Exit App", class="exitbutton")
    ),
    align = "center"
  ), #end section to add buttons
  
  tags$hr(),
  
  
  #################################
  #Split screen with input on left, output on right
  fluidRow(
    #all the inputs in here
    column(6,
           #################################
           # Inputs section
           h2('Simulation Settings'),
           fluidRow(
             column(6,
                    numericInput("H0", "initial number of hosts (H0)", min = 0, max = 5000, value = 50, step = 1)
             ),
             column(6,
                    numericInput("W0", "initial number of free living parasites (W0)", min = 0, max = 10000, value = 200, step = 1)
             )
           ), #close fluidRow structure for input
           fluidRow(
             column(6,
                    numericInput("P0", "initial number of adult parasites (P0)", min = 0, max = 10000, value = 200, step = 1)
             ),
             column(6,
                    numericInput("tmax", "Time", min = 1, max = 500, value = 50, step = 1)
             )
           ), #close fluidRow structure for input 
           fluidRow(
             column(6,
                    numericInput("a", "Host birth rate (a)", min = 0.1, max = 5, value = 1.1, step = 0.1)
             
                    ),
             column(6,
                    numericInput("b", "Host natural death rate (b)", min = 0.01, max = 5, value = 0.75, step = 0.1  )
             )
           ), #close fluidRow structure for input
           fluidRow(
             column(6,
                    numericInput("alpha", HTML("&alpha;"), min = 0, max = 1, value = 0.04, step = 0.001)
             ),
             column(6,
                    numericInput("delta", HTML("&delta;"), min = 0, max = 1, value = 0.001, step = 0.0001  )
             )
           ), #close fluidRow structure for input 
           fluidRow(
             column(4,
                    numericInput("lambda", HTML("&lambda;"), min = 0, max = 50, value = 11, step = 1)
             ),
             column(4,
                    numericInput("gamma", HTML("&gamma;"), min = 0, max = 100, value = 6, step = 1)
             ),
             column(4,
                    numericInput("beta", HTML("&beta;"), min = 0, max = 2, value = 0.4, step = 0.001 )
             )
           ), #close fluidRow structure for input
           fluidRow(
             column(4,
                    numericInput("u", HTML("&mu;"), min = 0, max = 5, value = 1, step = 0.01)
             ),
             column(4,
                    numericInput("k", "Parameter of the negative binomial distribution (k)", min = 0, max = 15, value = 1, step = 0.1 )
             )
           ) #close fluidRow structure for input
    ), #end sidebar column for inputs
    
    #all the outcomes here
    column(6,
           
           #################################
           #Start with results on top
           h2('Simulation Results'),
           plotOutput(outputId = "plot"),
           # PLaceholder for results of type text
           htmlOutput(outputId = "text"),
           #Placeholder for any possible warning or error messages (this will be shown in red)
           htmlOutput(outputId = "warn"),
           
           tags$head(tags$style("#warn{color: red;
                                font-style: italic;
                                }")),
           tags$hr()
           
           ) #end main panel column with outcomes
  ), #end layout with side and main panel
  #################################
  #Ends the 2 column structure with inputs on left and outputs on right
  
  
  
  #################################
  #Instructions section at bottom as tabs
  h2('Instructions'),
  
  #use external function to generate all tabs with instruction content
  
  div(align="center", style="font-size:small") #footer
  ) #end fluidpage

