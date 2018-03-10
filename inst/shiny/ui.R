library(shiny)
# Define UI for miles per gallon application
library(shinydashboard)
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Input", tabName = "input", icon = icon("file-text")),
    menuItem("Estimation Settings", icon = icon("rocket"), tabName = "est"),
    menuItem("Estimation Summary", icon = icon("check-square-o"), tabName = "summary"),
    menuItem("Parameter Estimates", icon = icon("superscript"), tabName = "par"),
    menuItem("Q-matrix Validation Outputs", icon = icon("th"), tabName = "Qval"),
    menuItem("Model selection Outputs", icon = icon("list"), tabName = "ms"),
    menuItem("Plots", icon = icon("bar-chart"), tabName = "plot"),
    menuItem("About", icon = icon("users"), tabName = "about")
  )
)

body <- dashboardBody(
  tabItems(
    # input files
    tabItem(tabName = "input",
            h2("Response and Q-matrix files"),
            fluidRow(
              shinydashboard::box(
              title = "Input Files", width = 6, solidHeader = TRUE, status = "primary",
            fileInput('file1', 'Choose Response File',
                      accept = c('text/csv','text/comma-separated-values',
                                 'text/tab-separated-values', 'text/plain',
                                 '.csv', '.tsv')
            ),
            checkboxInput('header', 'Header', FALSE),
            radioButtons('sep', 'Separator',
                         c(Tab='\t', Comma=',',Semicolon=';',
                           Space=''),
                         '\t',inline = TRUE)
            ),
            shinydashboard::box(
              title = "Input Files", width = 6, solidHeader = TRUE, status = "primary",
              fileInput('file2', 'Choose Q-matrix',
                        accept = c('text/csv','text/comma-separated-values',
                                   'text/tab-separated-values','text/plain',
                                   '.csv','.tsv')
              ),
              checkboxInput('header2', 'Header', FALSE),
              radioButtons('sep2', 'Separator',
                                               c(Tab='\t', Comma=',',Semicolon=';',
                                                 Space=''),
                                               '\t',inline = TRUE)
              )),
            fluidRow(
              shinydashboard::box(title = "First 6 observations of the responses", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                  tableOutput('contents1')
              )
            ),
            fluidRow(
              shinydashboard::box(title = "First 6 items of the Q-matrix", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                  tableOutput('contents2')
              )
            )
    ),

    tabItem(tabName = "est",
            h2("Estimation Specifications"),
            fluidRow(
              shinydashboard::box(
                title = "Models", width = 3, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                selectInput("type", label = "Fitted CDMs",
                            choices = list("GDINA" = "GDINA", "DINA" = "DINA","DINO" = "DINO",
                                           "ACDM" = "ACDM", "LLM" = "LLM", "RRUM" = "RRUM"), selected = "GDINA")
              ),
              shinydashboard::box(
                title = "Attribute Distribution", width = 5, collapsible = TRUE, solidHeader = TRUE,status = "primary",
                selectInput("attdis", label = "Attribute distribution",
                            choices = list("Saturated" = 0, "Higher-order" = 1,"Fixed"=2), selected = 0),
                selectInput("hom", label = "Applicable only when higher-order model is selected",
                            choices = list("Rasch" = "Rasch", "1PL" = "1PL","2PL" = "2PL"), selected = "1PL")
              ),
              shinydashboard::box(
                title = "Other Settings", width = 4, solidHeader = TRUE,collapsible = TRUE, status = "primary",
                checkboxInput("seq", label = "Sequential models?", value = FALSE),
                checkboxInput("mono", label = "Monotonic Constraints?", value = FALSE),
                checkboxInput("qvalcheck", label = "Q-matrix validation?", value = FALSE),
                checkboxInput("modelsel", label = "Item-level model selection?", value = FALSE)
              ),
            shinydashboard::box(actionButton("goButton", "Click to Estimate!"),width = 3)
    )
),
tabItem(tabName = "summary",
        h2("Summary Information"),

        fluidRow(
          shinydashboard::box(
            title = "Estimation Summary", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            verbatimTextOutput('iter.info')
          )),
        fluidRow(
          shinydashboard::box(
            title = "Relative test fit", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            verbatimTextOutput('info')
          )),
        fluidRow(
          shinydashboard::box(
            title = "Absolute test fit", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            verbatimTextOutput('itfit')
          ))),
tabItem(tabName = "par",
        h2("Parameter Estimation"),


        fluidRow(
          shinydashboard::box(
            title = "Item Parameter Estimation Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            selectInput("ips", label = "Item Parameters",
                        choices = list("catprob" = "catprob", "gs" = "gs","delta" = "delta","LCprob"="LCprob"), selected = "catprob"),
            checkboxInput("ipse", label = "Estimate S.E.?", value = TRUE)
            ),shinydashboard::box(
            title = "Item parameters (Category Probability of Success)", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            verbatimTextOutput('ip')
          )),
        fluidRow(
          shinydashboard::box(
          title = "Person Parameter Estimation Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          selectInput("pp", label = "Person Parameter Estimation Method:",
                      choices = list("EAP" = "EAP", "MAP" = "MAP","MLE" = "MLE","MasteryProb"="mp"), selected = "EAP"),
          downloadButton('downloadpp', 'Download'),
          radioButtons("ppfiletype", "File type:",
                       choices = c("csv", "tsv"))
        ),
        shinydashboard::box(
            title = "Person Parameter Estimates of first 10 observations", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            verbatimTextOutput('pparm')
          ))
),
tabItem(tabName = "Qval",
          h2("Q-matrix validation"),
          fluidRow(
            shinydashboard::box(
              title = "Q-matrix Validation Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
              sliderInput("PVAFcutoff", label = h3("PVAF cutoff"), min = 0,
                          max = 1, value = 0.95)
            ),
            shinydashboard::box(
              title = "Suggested Q-matrix", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
              verbatimTextOutput('sugQ')
            )),
        fluidRow(shinydashboard::box(
          title = "Mesa Plot Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          numericInput("item.mesaplot",label = "Item #",
                       value = 1),
          radioButtons("mesatype", "Plot type:",
                       choices = c("best", "all")),
          checkboxInput("datalabel", label = "Data Labels?", value = FALSE),
          downloadButton('downloadMesaplot', 'Download Plot as PDF file')
        ),
        shinydashboard::box(
          title = "Mesa plot", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          plotOutput("mesaplot")
        ))
),
tabItem(tabName = "ms",
        h2("Item-level model selection outputs"),
        fluidRow(
          shinydashboard::box(
            title = "Wald statistics", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            verbatimTextOutput('ws')
          ),
          shinydashboard::box(
            title = "p-values", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            verbatimTextOutput('pv')
          ))
),
tabItem(tabName = "plot",
        h2("Item Response Function Plots"),
        fluidRow(box(
          title = "Item #", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          numericInput("item.plot",label = "",
                       value = 1),
          checkboxInput("IRFplotse", label = "Errorbars?", value = FALSE),
          downloadButton('downloadIRFplot', 'Download Plot as PDF file')
        ),
        shinydashboard::box(
            title = "IRF plot", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            plotOutput("plot")
          ))
),
tabItem(tabName = "about",
        h2("About the GDINA R package and this GUI"),
        fluidRow(box(
          title = "", width = 12,

          p("This package is developed by Wenchao Ma and Jimmy de la Torre. The development of the GDINA R package aims to help students and researchers conduct CDM analyses as easily as possible."),
          p('This GUI application is developed with',
            a("Shiny", href="http://www.rstudio.com/shiny/", target="_blank"), 'and',
            a("shinydashboard.", href="http://rstudio.github.io/shinydashboard/", target="_blank"),
            ''),

          p(" This program is free software, so you can redistribute it and or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation version 3 of the License."),

          p("This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details."),

          p("This package is still under development. If you find bugs, please email Wenchao Ma at wenchao.ma@ua.edu.")
        ))
)
)
)

dashboardPage(
  dashboardHeader(title = "CDM Analysis"),
  sidebar,
  body
)

