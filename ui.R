ui <- navbarPage("PhyloME", theme = shinytheme("yeti"),
           
           # Introduction
           tabPanel("Introduction", fluid = TRUE, icon = icon("house"),
                    #tags$style(button_color_css),
                    
                    # Page Layout
                    mainPanel(
                        h1("Welcome to PhyloME!", align = "left"),
                        imageOutput("binf"),
                        p("https://scitechdaily.com/unlimited-possibilities-new-law-of-physics-could-predict-genetic-mutations/", style = "font-size:10px")
                    )
                    
            ),
           
           # Pipeline
           tabPanel("The Pipeline", fluid = TRUE, icon = icon("tree"),
                    #titlePanel("The Pipeline"),
                    
                    h4("Customizing your Phylogenetic Tree"),
                    
                    sidebarLayout(
                      sidebar_content <- sidebarPanel(
                        
                        h5("Step 1: File Upload"),
                        fileInput("file1", "Choose Your File",
                                  multiple = TRUE),
                        
                        tags$hr(),
                    
                        h5("Step 2: Multiple Sequence Alignment"),
                        selectInput("algo_type", 
                                    label = "Select your Multiple Sequence Alignment Choice", 
                                    choices <- list("Muscle", "Clustal W - blosum", "Clustal W - pam", "Clustal W - gonnet", "Clustal W - id",
                                                    "Clustal Omega - BLOSUM30", "Clustal Omega - BLOSUM40", "Clustal Omega - BLOSUM50", 
                                                    "Clustal Omega - BLOSUM65", "Clustal Omega - BLOSUM80", "Clustal Omega - Gonnet"),
                                    selected = NULL
                                    ),
                        
                        tags$hr(),
                        
                        h5("Step 3: Constructing the Tree"),
                        # Select tree construction type 
                        radioButtons("tree_type", "Select your tree construction type",
                                     choices = list("Neighbour Joining", "Maximum Likelihood"),
                                     selected = NULL),
                        
                        # Select tree rooting type 
                        radioButtons("root_type", "Select your rooting style",
                                     choices = list("Unrooted", "Midpoint", "Outgroup"),
                                     selected = NULL),
                        
                        selectInput("outgroup_type", "Select what species will be the outgroup",
                                    choices = nj_tree$tip.label),
                        
                        
                        materialSwitch("bootstrap_switch", 
                                       label = "Bootstrap Tree",
                                       value = FALSE,
                                       status = "success"),
                        
                        tags$hr(),
                        
                        h5("Step 4: Select Your Tree"),
                        radioButtons("print_tree", "Select which tree(s) you would like to display",
                                           choices = list("Phylogenetic Tree", "Bootstrapped Phylogenetic Tree")),
                        
                        tags$hr(),
                        
                        h6("Step 5: Download"),
                        #radioButtons("download_type", "Select the file format you would like to download",
                        #               choices = list("PDF"),
                        #               selected = NULL),
                        downloadBttn("download_data", label = "Download PDF", style = "bordered", color = "primary", icon = icon("download"))
                        
                      ),
                      
                      main_content <- mainPanel(
                        #tableOutput("file_content"),
                        plotOutput("phylo_plot")
                      )
                    )
            )
           
)
