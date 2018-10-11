
  library(GDINA)
sidebar <- shinydashboard::dashboardSidebar(
  shinydashboard::sidebarMenu(
    shinydashboard::menuItem("Input", tabName = "input", icon = icon("file-text")),
    shinydashboard::menuItem("Estimation Settings", icon = icon("rocket"), tabName = "est"),
    shinydashboard::sidebarMenuOutput("summary"),
    shinydashboard::sidebarMenuOutput("fit"),
    shinydashboard::sidebarMenuOutput("par"),
    shinydashboard::sidebarMenuOutput("qv"),
    shinydashboard::sidebarMenuOutput("msec"),
    shinydashboard::sidebarMenuOutput("menuplot"),
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
            shiny::fileInput('file1', 'Response matrix',
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
              shiny::fileInput('file2', 'Q-matrix',
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
                title = "Models", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                shiny::selectInput("type", label = "Select a single CDM for all items",
                            choices = list("GDINA [Generalized deterministic inputs, noisy and gate] model" = "GDINA",
                                           "logit GDINA model [loglinear CDM]" = "logitGDINA",
                                           "log GDINA model" = "logGDINA",
                                           "DINA [Deterministic inputs, noisy and gate] model" = "DINA",
                                           "DINO [Deterministic inputs, noisy or gate] model" = "DINO",
                                           "ACDM [Additive cognitive diagnosis model]" = "ACDM",
                                           "R-RUM [Reduced reparameterized unified model]" = "RRUM",
                                           "LLM/C-RUM [Linear logistic model/Compensatory RUM]" = "LLM",
                                           "To be specified..."="UM"), selected = "GDINA"),
                shiny::textInput('mv', 'Or enter a vector of models (comma delimited without quotation marks)', 'GDINA,DINA,LLM,...')
                ),
              shinydashboard::box(
                title = "Attribute Distribution", width = 6, collapsible = TRUE, solidHeader = TRUE,status = "primary",
                shiny::selectInput("attdis", label = "Attribute distribution",
                            choices = list("Saturated" = 0, "Higher-order" = 1,"Fixed"=2), selected = 0),
                shiny::selectInput("hom", label = "Applicable only when higher-order model is selected",
                            choices = list("Rasch" = "Rasch", "1PL" = "1PL","2PL" = "2PL"), selected = "1PL")
              )

    ),
    shiny::fluidRow(
      shinydashboard::box(
        title = "Other Settings", width = 6, solidHeader = TRUE,collapsible = TRUE, status = "primary",
        shiny::checkboxInput("seq", label = "Sequential models?", value = FALSE),
        shiny::checkboxInput("mono", label = "Monotonic Constraints?", value = FALSE),
        shiny::checkboxInput("qvalcheck", label = "Q-matrix validation?", value = FALSE),
        shiny::checkboxInput("modelsel", label = "Item-level model selection?", value = FALSE)
      ),
      shinydashboard::box(shiny::actionButton("goButton", "Click to Estimate!"),width = 6))
),
shinydashboard::tabItem(tabName = "summary",
                        shiny::h2("Model Estimation Summary"),
                        hr(),
                        shiny::fluidRow(
                          shinydashboard::box(
                            title = "Estimation Summary", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                            shiny::verbatimTextOutput('iter.info')
                          ),
                          shinydashboard::box(
                            title = "Classification Summary", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
                            shiny::verbatimTextOutput('iter.info2')
                          )
                        )
),
shinydashboard::tabItem(tabName = "fit",
        shiny::fluidRow(
          shinydashboard::box(
            title = "Relative test fit", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('info')
          ),
          shinydashboard::box(
            title = "Absolute test fit", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('itfit')
          )),
        shiny::fluidRow(shinydashboard::box(
          title = "Item-fit Plot Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          shiny::radioButtons("heatmap.type", "Plot type:",
                       choices = c("log odds ratio", "transformed correlation")),
          shiny::checkboxInput("heatmap.adjust", label = "Bonferroni adjusted?", value = TRUE),
          shiny::downloadButton('downloadHeatPlot', 'Download Plot as PDF file')
        ),
          shinydashboard::box(
            title = "Heatmap plots", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("heatplot")
          ))),
shinydashboard::tabItem(tabName = "par",
        shiny::h2("Parameter Estimation"),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Item Parameter Estimation Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::selectInput("ips", label = "Item Parameters",
                        choices = list("Success probabilities of reduced latent classes" = "catprob",
                                       "Guessing and slip parameters" = "gs",
                                       "Delta parameters" = "delta",
                                       "Success probabilities of all latent classes"="LCprob"), selected = "catprob"),
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
                      choices = list("EAP" = "EAP", "MAP" = "MAP","MLE" = "MLE","Probabilities of mastering each attribute"="mp"), selected = "EAP"),
          shiny::downloadButton('downloadpp', 'Download'),
          shiny::radioButtons("ppfiletype", "File type:",
                       choices = c("csv", "tsv"))
        ),
        shinydashboard::box(
            title = "Person Parameter Estimates of first 10 observations", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('pparm')
          )),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Estimated Proportions of Latent Classes Specifications", width = 4, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::selectInput("plc", label = "Sorted by:",
                               choices = list("default" = "default", "decreasing" = "decreasing","increasing" = "increasing"), selected = "default"),
            shiny::downloadButton('downloadplc', 'Download'),
            shiny::radioButtons("plcfiletype", "File type:",
                                choices = c("csv", "tsv"))
          ),
          shinydashboard::box(
            title = "Estimated Proportions of first 10 Latent Classes", width = 8, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::verbatimTextOutput('plc.output')
          ))
),
shinydashboard::tabItem(tabName = "Qval",
        shiny::h2("Q-matrix validation"),
        shiny::fluidRow(
            shinydashboard::box(
              title = "Q-matrix Validation Specifications", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
              shiny::sliderInput("PVAFcutoff", label = h3("PVAF cutoff"), min = 0,
                          max = 1, value = 0.95)
            ),
            shinydashboard::box(
              title = "Suggested Q-matrix", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
              shiny::verbatimTextOutput('sugQ')
            )),
        shiny::fluidRow(shinydashboard::box(
          title = "Mesa Plot Specifications", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
          shiny::numericInput("item.mesaplot",label = "Item #",
                       value = 1, min = 1),
          shiny::radioButtons("mesatype", "Plot type:",
                       choices = c("best", "all")),
          shiny::checkboxInput("datalabel", label = "Data Labels?", value = FALSE),
          shiny::downloadButton('downloadMesaplot', 'Download Plot as PDF file')
        ),
        shinydashboard::box(
          title = "Mesa plot", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
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
        shiny::h2("Plots for Individual Statistics"),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Specifications for Individuals' mastery plots", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::textInput('personid', 'Enter a vector of individuals (comma delimited)', "1,2,5"),
            shiny::checkboxInput("HPlot", label = "Horizontal?", value = FALSE),
            shiny::downloadButton('downloadmpplot', 'Download Plot as PDF file')
          ),
          shinydashboard::box(
            title = "Plot of probability of mastery for individuals", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("Mplot")
          )),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Specifications for Individual posterior probability plot", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::numericInput("ippid",label = "Specify an individual:", value = 1, min = 1),
            shiny::selectInput("ippplc", label = "Sorted by:",
                               choices = list("default" = "default", "decreasing" = "decreasing","increasing" = "increasing"), selected = "default"),
            shiny::textInput('inlc', 'Enter the maximum number of latent classes:', "10"),
            shiny::checkboxInput("ippAdir", label = "Horizontal?", value = FALSE),
            shiny::downloadButton('downloadiPPplot', 'Download Plot as PDF file')
          ),
          shinydashboard::box(
            title = "Individual posterior probability plot", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("iPostProbplot")
          )),
        shiny::h2("Plots for Group Statistics"),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Specifications for Proportions of Latent Classes Plot", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::selectInput("ppplc", label = "Sorted by:",
                               choices = list("default" = "default", "decreasing" = "decreasing","increasing" = "increasing"), selected = "default"),
            shiny::textInput('nlc', 'Enter the number of latent classes:', "10"),
            shiny::checkboxInput("ppAdir", label = "Horizontal?", value = FALSE),
            shiny::downloadButton('downloadPPplot', 'Download Plot as PDF file')
          ),
          shinydashboard::box(
            title = "Plot of Proportions of Latent Classes", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("PostProbplot")
          )),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Specifications for Attribute Prevalence plot", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::selectInput("palette", label = "RColorBrewer palette",
                               choices = list("Set2" = "Set2", "Greys" = "Greys","Paired" = "Paired","Accent"="Accent"), selected = "default"),
            shiny::checkboxInput("Adir", label = "Horizontal?", value = FALSE),
            shiny::downloadButton('downloadAPplot', 'Download Plot as PDF file')
          ),
          shinydashboard::box(
            title = "Plot of Attribute Prevalence", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("APplot")
          )),
        shiny::h2("Plots for Item Statistics"),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Specifications for IRF plots", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::numericInput("item.plot",label = "Specify the item for IRF plot:", value = 1,min = 1),
            shiny::checkboxInput("IRFplotse", label = "Errorbars?", value = FALSE),
            shiny::downloadButton('downloadIRFplot', 'Download Plot as PDF file')
          ),
          shinydashboard::box(
            title = "IRF plot", width = 6, solidHeader = TRUE, collapsible = TRUE, status = "primary",
            shiny::plotOutput("plot")
          ))
),
shinydashboard::tabItem(tabName = "about",
        shiny::h2("About the GDINA R package and this GUI"),
        shiny::fluidRow(
          shinydashboard::box(
          title = "", width = 12,

          shiny::p('This GUI application is developed using',
                   shiny::a("Shiny", href="http://www.rstudio.com/shiny/", target="_blank"), 'and',
                   shiny::a("shinydashboard", href="http://rstudio.github.io/shinydashboard/", target="_blank"),
                   'and distributed as part of',
            ' the',shiny::a("GDINA", href="https://wenchao-ma.github.io/GDINA/", target="_blank"),"R package by Wenchao Ma and Jimmy de la Torre."),

          shiny::p("The GUI application, as well as the GDINA R package, is free, and you can redistribute it and or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation version 3 of the License."),

          shiny::p("The program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details."),

          shiny::p("Should you have any comments or suggestions, please email Wenchao Ma at wenchao.ma@ua.edu.")
        ))
)
)
)

shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "GDINA GUI"),
  sidebar,
  body
)
