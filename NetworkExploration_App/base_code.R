library('igraph')
library('shiny')
library('tidyverse')
library('ggplot2')
library('stringi')
library('readxl')
library('shinythemes')
library('shinyWidgets')




rm(list=ls())


ui <- fluidPage(
  theme = shinytheme("flatly"),
  div(style = "padding: 1px 0px; width: '100%'",
      titlePanel(
        title = "",
        windowTitle = "IDEANET NETWORK VISUALIZER"
      )
  ),
  navbarPage(
    title = "NETWORK VISUALIZER WITH IDEANET",
    tabPanel(
      "Upload",
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "Node Data",
          sidebarPanel(
            fileInput(
              'raw_nodes', "Upload Node Data (CSV)", multiple = FALSE, 
              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"),
              buttonLabel = "Browse...", placeholder = "No file selected"
              ),
            checkboxInput("node_header", tags$b("Does the file have a header?"), TRUE),
            #pickerInput("node_group", "Group Assignments (Optional)", choices = 'output.node_raw_colnames'),
            tags$p(span("Large datasets may take a few seconds to render.", style = "color:red")),
            tags$p(HTML("<b>Upload</b> the node data. The application only accepts .csv files.")),
            tags$p(HTML("The node data <b>may</b> contain extra columns for the node characteristics.")),
            tags$p(HTML("<b>Continue</b> on to upload the edge list."))
          ),
          mainPanel(
            dataTableOutput('node_raw_upload')
          )
        ),
        tabPanel(
          "Edge List",
          sidebarPanel(
            fileInput(
              'raw_edges', "Upload Edge Data (CSV)", multiple = FALSE, 
              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"),
              buttonLabel = "Browse...", placeholder = "No file selected"
            ),
            checkboxInput("edge_header", tags$b("Does the file have a header?"), TRUE),
            tags$p(span("Large datasets may take a few seconds to render.", style = "color:red")),
            tags$p(HTML("<b>Upload</b> the edge list. The application only accepts .csv files.")),
            tags$p(HTML("The edge list <b>may</b> contain an extra column for the edge weight.")),
            tags$p(HTML("<b>Continue</b> on to process the data before visualizing it."))
          ),
          mainPanel(
            dataTableOutput('edge_raw_upload')
          )
        )
      )
    ),
    tabPanel(
      "Process",
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "Process Node Data",
          sidebarPanel(
            uiOutput("node_ids"),
            uiOutput("node_factor"),
            uiOutput("node_numeric"),
            tags$p(span("Questions with an asterisk are required.", style = "color:red")),
            tags$p(HTML("<b>Process</b> the node data by assigning the columns to their function.")),
            tags$p(HTML("The <b>node</b> <b>identifier</b> contains the name or number of each node.")),
            tags$p(HTML("The <b>factor</b> and <b>numeric</b> columns contain either qualitative (grouping) or quantitative (numerical) information about each node."))
          ),
          mainPanel(
            dataTableOutput('node_processed')
          )
        ),
        tabPanel(
          "Process Edge Data ",
          sidebarPanel(
            uiOutput("edge_in"),
            uiOutput("edge_out"),
            uiOutput("edge_weight"),
            tags$p(span("Questions with an asterisk are required.", style = "color:red")),
            tags$p(HTML("<b>Process</b> the edge data by assigning the columns to their function.")),
            tags$p(HTML("The edge <b>\"in\" and \"out\"</b> columns contain information for the start and end of each edge. If the graph is undirected, the order of columns doesn't matter.")),
            tags$p(HTML("<b>Edge weight</b> contains information on the strength of the connection between two nodes."))
          ),
          mainPanel(
            dataTableOutput('edge_processed')
          )
        )
      )
      
    ),
    tabPanel(
      "Visualize",
      sidebarPanel(
        checkboxInput("outlier_toggle", tags$b("Remove outliers?"), FALSE),
        checkboxInput("simplify_toggle", tags$b("Remove self-loops and duplicate edges?"), FALSE),
        uiOutput("layout_picker"),
        tags$p(HTML("<u>Node Features</u>")),
        #color
        #size
        tags$p(HTML("<b>Process</b> the edge data by assigning the columns to their function.")),
        tags$p(HTML("The edge <b>\"in\" and \"out\"</b> columns contain information for the start and end of each edge. If the graph is undirected, the order of columns doesn't matter.")),
        tags$p(HTML("<b>Edge weight</b> contains information on the strength of the connection between two nodes."))
      ),
      mainPanel(
        plotOutput('network')
      )
    )
  )
)


server <- function(input, output) {
  
  #### UPLOADS #####
  
  ## Node Data
  
  #Upload Node Data
  node_data <- reactive({
    req(input$raw_nodes)
    read.csv(input$raw_nodes$datapath, header = input$node_header)
  })
  
  #Display Node Data
  output$node_raw_upload <- renderDataTable({
    node_data()
  })
  
  
  ## Edge Data
  
  #Upload Edge Data
  edge_data <- reactive({
    req(input$raw_edges)
    read.csv(input$raw_edges$datapath, header = input$edge_header)
  })
  
  #Display Edge Data
  output$edge_raw_upload <- renderDataTable({
    edge_data()
  })

  
  #### PROCESS #####
  # Redisplay Datatables
  output$node_processed <- renderDataTable({
    node_data()
  })
  output$edge_processed <- renderDataTable({
    edge_data()
  })
  
  #Node Processing Options
  output$node_ids <- renderUI({
    selectInput(inputId = "node_id_col", label = "Select the column with node identifiers*", choices = append("N/A",colnames(node_data())), selected = "No Selection", multiple = FALSE)
  })
  output$node_factor <- renderUI({
    selectInput(inputId = "node_factor_col", label = "Select the columns with factors or groups", choices = append("N/A",colnames(node_data())), multiple = TRUE)
  })
  output$node_numeric <- renderUI({
    selectInput(inputId = "node_numeric_col", label = "Select the columns with a numeric node trait", choices = append("N/A",colnames(node_data())), multiple = TRUE)
  })
  
  #Edge Processing Options
  output$edge_in <- renderUI({
    selectInput(inputId = "edge_in_col", label = "Select the column with the \"in\" nodes*", choices = append("N/A",colnames(edge_data())), selected = "No Selection", multiple = FALSE)
  })
  output$edge_out <- renderUI({
    selectInput(inputId = "edge_out_col", label = "Select the column with \"out\" nodes*", choices = append("N/A",colnames(edge_data())), selected = "No Selection", multiple = FALSE)
  })
  output$edge_weight <- renderUI({
    selectInput(inputId = "edge_weight_col", label = "Select the column with the edge weight", choices = append("N/A",colnames(edge_data())), selected = "No Selection", multiple = FALSE)
  })
  
  #### Visualize####
  
  edge_list <- reactive({
    req(input$edge_in_col != "N/A")
    req(input$edge_out_col != "N/A")
    edge_data <- edge_data()
    edge_data[,input$edge_in_col] <- as.character(edge_data[,input$edge_in_col])
    edge_data[,input$edge_out_col] <- as.character(edge_data[,input$edge_out_col])
    edge_data <- edge_data[, c(input$edge_in_col, input$edge_out_col)]
    colnames(edge_data) <- c('in_node', 'out_node')
    edge_data
  })
  
  node_list <- reactive({
    req(input$node_id_col != "N/A")
    node_data <- node_data()
    node_data[,input$node_id_col] <- as.character(node_data[,input$node_id_col])
    node_data <- node_data[, input$node_id_col, drop=FALSE]
    colnames(node_data) <- c('node_ids')
    node_data
  })
  
  ## Visualizing Pipeline
  # 1. Outliers
  # 2. Self Loops / Duplicate Edges
  # 3. Layout
  
  # 1. Outliers
  
  net1 <- reactive({
    if (input$outlier_toggle == TRUE) {
      net <- graph_from_data_frame(d=edge_list(), vertices=node_list(), directed=F)
      net <- delete.vertices(net, degree(net)==0)
      net
    } else {
      net <- graph_from_data_frame(d=edge_list(), vertices=node_list(), directed=F)
      net
    }
    
  })
  
  
  
  # 2. Simplify (Self Loops and Repeating Edges)
  
  net2 <- reactive({
    if (input$simplify_toggle == TRUE) {
      net <- net1()
      net <- igraph::simplify(net)
      net
    } else {
      net <- net1()
      net
    }

  })

  # 3. Layout Picker
  
  layout_choices <- c("layout_as_star", "layout_as_tree", "layout_in_circle",
                     "layout_nicely", "layout_on_grid", "layout_on_sphere", "layout_randomly", "layout_with_dh", "layout_with_fr",
                     "layout_with_gem", "layout_with_graphopt", "layout_with_kk", "layout_with_lgl", "layout_with_mds"
                     )
  output$layout_picker <- renderUI({
    selectInput(inputId = "layout_choice", label = "Select the desired network layout", choices = layout_choices, selected = "layout_with_fr", multiple = FALSE)
  })
  
  output$network <- renderPlot({
    net <- net2()
    if(input$layout_choice == 'layout_as_star')     {layoutselected <- layout_as_star}
    if(input$layout_choice == 'layout_as_tree')     {layoutselected <- layout_as_tree}
    if(input$layout_choice == 'layout_in_circle')     {layoutselected <- layout_in_circle}
    if(input$layout_choice == 'layout_nicely')     {layoutselected <- layout_nicely}
    if(input$layout_choice == 'layout_on_grid')     {layoutselected <- layout_on_grid}
    if(input$layout_choice == 'layout_on_sphere')     {layoutselected <- layout_on_sphere}
    if(input$layout_choice == 'layout_randomly')     {layoutselected <- layout_randomly}
    if(input$layout_choice == 'layout_with_dh')     {layoutselected <- layout_with_dh}
    if(input$layout_choice == 'layout_with_fr')     {layoutselected <- layout_with_fr}
    if(input$layout_choice == 'layout_with_gem')     {layoutselected <- layout_with_gem}
    if(input$layout_choice == 'layout_with_graphopt')     {layoutselected <- layout_with_graphopt}
    if(input$layout_choice == 'layout_with_kk')     {layoutselected <- layout_with_kk}
    if(input$layout_choice == 'layout_with_lgl')     {layoutselected <- layout_with_lgl}
    if(input$layout_choice == 'layout_with_mds')     {layoutselected <- layout_with_mds}
    plot(net, edge.arrow.size=.4, edge.curved=.1, layout = layoutselected)
  })
  
}

shinyApp(ui, server)

