# suppressMessages(require(shiny))
# suppressMessages(require(shinydashboard))
sidebar <- shinydashboard::dashboardSidebar(
  shinydashboard::sidebarMenu(
    shinydashboard::menuItem("Input", tabName = "input", icon = icon("file-text")),
    shinydashboard::menuItem("Estimation Settings", icon = icon("rocket"), tabName = "est"),
    shinydashboard::menuItem("Estimation Summary", icon = icon("check-square-o"), tabName = "summary"),
    shinydashboard::menuItem("Parameter Estimates", icon = icon("superscript"), tabName = "par"),
    shinydashboard::menuItem("Q-matrix Validation Outputs", icon = icon("th"), tabName = "Qval"),
    shinydashboard::menuItem("Model selection Outputs", icon = icon("list"), tabName = "ms"),
    shinydashboard::menuItem("Item Response Function Plots", icon = icon("bar-chart"), tabName = "plot"),
    shinydashboard::menuItem("About", icon = icon("users"), tabName = "about")
  )
)

body <- shinydashboard::dashboardBody(
  shinydashboard::tabItems(
    # input files
    shinydashboard::tabItem(tabName = "input",
            shiny::h2("Response and Q-matrix files"),
            shiny::fluidRow(
              shinydashboard::box(
              title = "Input Files", width = 6, solidHeader = TRUE, status = "primary",
            shiny::fileInput('file1', 'Choose Response File',
                      accept = c('text/csv','text/comma-separated-values',
                                 'text/tab-separated-values', 'text/plain',
                                 '.csv', '.tsv')
            ),
            shiny::checkboxInput('header', 'Header', FALSE),
            shiny::radioButtons('sep', 'Separator',
                         c(Tab='\t', Comma=',',Semicolon=';',
                           Space=''),
                         '\t',inline = TRUE)
            ),
            shinydashboard::box(
              title = "Input Files", width = 6, solidHeader = TRUE, status = "primary",
              shiny::fileInput('file2', 'Choose Q-matrix',
                        accept = c('text/csv','text/comma-separated-values',
                                   'text/tab-separated-values','text/plain',
                                   '.csv','.tsv')
              ),
              shiny::checkboxInput('header2', 'Header', FALSE),
              shiny::radioButtons('sep2', 'Separator',
                                               c(Tab='\t', Comma=',',Semicolon=';',
                                                 Space=''),
                                               '\t',inline = TRUE)
              )),
            shiny::fluidRow(
              shinydashboard::box(title = "First 6 observations of the responses", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                                  shiny::tableOutput('contents1')
              )
            ),
            shiny::fluidRow(
              shinydashboard::box(title = "First 6 items of the Q-matrix", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                                  shiny::tableOutput('contents2')
              )
            )
    ),

    shinydashboard::tabItem(tabName = "est",
            shiny::h2("Estimation Specifications"),
            shiny::fluidRow(
              shinydashboard::box(
                title = "Models", width = 3, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                shiny::selectInput("type", label = "Fitted CDMs",
                            choices = list("GDINA" = "GDINA", "DINA" = "DINA","DINO" = "DINO",
                                           "ACDM" = "ACDM", "LLM" = "LLM", "RRUM" = "RRUM"), selected = "GDINA")
              ),
              shinydashboard::box(
                title = "Attribute Distribution", width = 5, collapsible = TRUE, solidHeader = TRUE,status = "primary",
                shiny::selectInput("attdis", label = "Attribute distribution",
                            choices = list("Saturated" = 0, "Higher-order" = 1,"Fixed"=2), selected = 0),
                shiny::selectInput("hom", label = "Applicable only when higher-order model is selected",
                            choices = list("Rasch" = "Rasch", "1PL" = "1PL","2PL" = "2PL"), selected = "1PL")
              ),
              shinydashboard::box(
                title = "Other Settings", width = 4, solidHeader = TRUE,collapsible = TRUE, status = "primary",
                shiny::checkboxInput("seq", label = "Sequential models?", value = FALSE),
                shiny::checkboxInput("mono", label = "Monotonic Constraints?", value = FALSE),
                shiny::checkboxInput("qvalcheck", label = "Q-matrix validation?", value = FALSE),
                shiny::checkboxInput("modelsel", label = "Item-level model selection?", value = FALSE)
              ),
            shinydashboard::box(shiny::actionButton("goButton", "Click to Estimate!"),width = 3)
    )
),
shinydashboard::tabItem(tabName = "summary",
        h2("Summary Information"),

        shiny::fluidRow(
          shinydashboard::box(
            title = "Estimation Summary", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('iter.info')
          )),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Relative test fit", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('info')
          )),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Absolute test fit", width = 12, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('itfit')
          )),
        shiny::fluidRow(shinydashboard::box(
          title = "Heatmap Plot Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          shiny::radioButtons("heatmap.type", "Plot type:",
                       choices = c("log odds ratio", "transformed correlation")),
          shiny::checkboxInput("heatmap.adjust", label = "Bonferroni adjusted?", value = TRUE)
        ),
          shinydashboard::box(
            title = "Heatmap plots", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("heatplot1")
          ))),
shinydashboard::tabItem(tabName = "par",
        shiny::h2("Parameter Estimation"),


        shiny::fluidRow(
          shinydashboard::box(
            title = "Item Parameter Estimation Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::selectInput("ips", label = "Item Parameters",
                        choices = list("catprob" = "catprob", "gs" = "gs","delta" = "delta","LCprob"="LCprob"), selected = "catprob"),
            shiny::checkboxInput("ipse", label = "Estimate S.E.?", value = TRUE)
            ),
          shinydashboard::box(
            title = "Item parameter Estimates", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('ip')
          )),
        shiny::fluidRow(
          shinydashboard::box(
          title = "Person Parameter Estimation Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          shiny::selectInput("pp", label = "Person Parameter Estimation Method:",
                      choices = list("EAP" = "EAP", "MAP" = "MAP","MLE" = "MLE","MasteryProb"="mp"), selected = "EAP"),
          shiny::downloadButton('downloadpp', 'Download'),
          shiny::radioButtons("ppfiletype", "File type:",
                       choices = c("csv", "tsv"))
        ),
        shinydashboard::box(
            title = "Person Parameter Estimates of first 10 observations", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('pparm')
          ))
),
shinydashboard::tabItem(tabName = "Qval",
        shiny::h2("Q-matrix validation"),
        shiny::fluidRow(
            shinydashboard::box(
              title = "Q-matrix Validation Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
              shiny::sliderInput("PVAFcutoff", label = h3("PVAF cutoff"), min = 0,
                          max = 1, value = 0.95)
            ),
            shinydashboard::box(
              title = "Suggested Q-matrix", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
              shiny::verbatimTextOutput('sugQ')
            )),
        shiny::fluidRow(shinydashboard::box(
          title = "Mesa Plot Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          shiny::numericInput("item.mesaplot",label = "Item #",
                       value = 1),
          shiny::radioButtons("mesatype", "Plot type:",
                       choices = c("best", "all")),
          shiny::checkboxInput("datalabel", label = "Data Labels?", value = FALSE),
          shiny::downloadButton('downloadMesaplot', 'Download Plot as PDF file')
        ),
        shinydashboard::box(
          title = "Mesa plot", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          shiny::plotOutput("mesaplot")
        ))
),
shinydashboard::tabItem(tabName = "ms",
        shiny::h2("Item-level model selection outputs"),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Wald statistics", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('ws')
          ),
          shinydashboard::box(
            title = "p-values", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('pv')
          ))
),
shinydashboard::tabItem(tabName = "plot",
        shiny::h2("Item Response Function Plots"),
        shiny::fluidRow(
          shinydashboard::box(
          title = "Item #", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          shiny::numericInput("item.plot",label = "", value = 1),
          shiny::checkboxInput("IRFplotse", label = "Errorbars?", value = FALSE),
          shiny::downloadButton('downloadIRFplot', 'Download Plot as PDF file')
        ),
        shinydashboard::box(
            title = "IRF plot", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("plot")
          ))
),
shinydashboard::tabItem(tabName = "about",
        shiny::h2("About the GDINA R package and this GUI"),
        shiny::fluidRow(
          shinydashboard::box(
          title = "", width = 12,

          shiny::p("This package is developed by Wenchao Ma and Jimmy de la Torre. The development of the GDINA R package aims to help students and researchers conduct CDM analyses as easily as possible."),
          shiny::p('This GUI application is developed with',
                   shiny::a("Shiny", href="http://www.rstudio.com/shiny/", target="_blank"), 'and',
                   shiny::a("shinydashboard.", href="http://rstudio.github.io/shinydashboard/", target="_blank"),
            ''),

          shiny::p(" This program is free software, so you can redistribute it and or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation version 3 of the License."),

          shiny::p("This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details."),

          shiny::p("This package is still under development. If you find bugs, please email Wenchao Ma at wenchao.ma@ua.edu.")
        ))
)
)
)

shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "GDINA GUI"),
  sidebar,
  body
)

