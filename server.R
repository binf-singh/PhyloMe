server <- function(input, output, session) {
  
  ## INTRODUCTION ##
  output$binf <- renderImage({
    width <- session$clientData$output_binf_width
    height <- session$clientData$output_binf_height
    list (
      src = file.path("/Users/anmol/Desktop/CAPSTONE/PhyloME Version 1/www/binf_gif.gif"),
      contentType = "image/gif",
      width = 800,
      height = 400, 
      align = "center"
    )
  }, deleteFile = F)
  
  ## PIPELINE ## 
  
  my_trees <- reactiveValues(basic_tree = NULL, bootstrap_tree = NULL)
  
  output$phylo_plot <- renderPlot({
    
    req(input$file1) 
    
    # Read in file
    location <- input$file1$datapath
    data <- readLines(location)
    
    # Function to  check sequence type 
    type_check <- function(data_file, location) {
      
      # Clean data to remove headers 
      data <- str_to_upper(data_file[-which(str_detect(data_file, ">"))])
      
      # Convert string to vector 
      data <- unlist(str_split(data, ""))
      
      # Unique letters to determine sequences 
      unique_AA <- c('R', 'N', 'D', 'Q', 'E', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'O', 'S', 'W', 'Y', 'V', 'B', 'Z', 'X', 'J')
      unique_DNA <- 'T'
      unique_RNA <- 'U'
      
      # Determine what type of sequence 
      det <- length(which(unique_AA%in%unlist(str_split(data, ""))))
      if (det > 0) {
        sequence_type <- "AA"
      } else {
        det <- length(which(unique_DNA%in%unlist(str_split(data, ""))))
        if (det > 0) {
          sequence_type <- "DNA"
        } else {
          sequence_type <- "RNA"
        }
      }
      
      return(sequence_type)
    } 
    
    # Get sequence type 
    seq_type <- type_check(data, location)
    
    
    # Prepare sequence for alignment #
    
    # Clean names 
    if (seq_type == "AA") {
      sequence <- readAAStringSet(filepath = location)
      names(sequence) <- str_remove_all(str_split_fixed(names(sequence), pattern = "\\[", n=2)[,2], "\\]")
    } else if (seq_type == "RNA") {
      RNA_string <- readRNAStringSet(filepath = location)
      sequence <- RNA2DNA(RNA_string)
      names(sequence) <- apply(str_split_fixed(names(sequence), pattern = " ",4)[,2:3], 1, base::paste, collapse = " ")
    } else {
      sequence <- readDNAStringSet(filepath = location)
      names(sequence) <- apply(str_split_fixed(names(sequence), pattern = " ",4)[,2:3], 1, base::paste, collapse = " ")
    }
    
    # Multiple Sequence Alignment #
    msa_type <- input$algo_type
    
    # Preform MSA depending on algorithm type
    if (msa_type == "Muscle") {
      msa_complete <- msaMuscle(sequence)
    } else {
      msa_algo <- strsplit(msa_type, "[-]")[[1]][1]
      sub_matrix <- strsplit(msa_type, "[- ]")[[1]][5]
      if (msa_algo == "Clustal W ") {
        msa_complete <- msaClustalW(sequence, substitutionMatrix = sub_matrix)
      } else {
        msa_complete <- msaClustalOmega(sequence, substitutionMatrix = sub_matrix)
      }
    }
    
    # Prepare for tree construction #
    
    # Convert to phyDat object
    phy <- as.phyDat(msaConvert(msa_complete, type = "seqinr::alignment"), type = seq_type)
    
    # Distance alignment based on sequence type
    if (seq_type == "AA") {
      distance <- dist.ml(phy, model = "JTT")
    } else {
      distance <- dist.ml(phy, model = "JC69")
    }
    
    # NJ distance method
    og_nj_tree <- phangorn::NJ(distance)
    tree_nodes <- og_nj_tree
    
    # Retrieve Rooting Style
    root_style <- input$root_type
    
    # Retrieve Bootstrap Value 
    bootstrap_value <- input$bootstrap_switch
    
    # Prepare Maximum Likelihood Tree
    if (input$tree_type == "Maximum Likelihood") {
      
      # Choose best model 
      mt <- modelTest(phy, model = "all", multicore = TRUE, mc.cores = 2)
      lowest <- mt$Model[mt$AIC==min(mt$AIC)]
      mod <- strsplit(lowest, "[+]")[[1]][1]
      mod <- unlist(mod)
      
      # Construct Tree
      ml_tree <- pml(og_nj_tree, phy, model = mod, k = 4, inv = .2)
      
      # Optimize Tree
      ml_tree <- optim.pml(ml_tree, optNni = TRUE, optBf = TRUE, optQ = TRUE, optInv = TRUE, optGamma = TRUE, optEdge = TRUE)
      
      tree_nodes <- ml_tree$tree
    } 
    
    # Rooting 
    root_style <- input$root_type
    if (root_style == "Outgroup") {
      species <- input$outgroup_type
      tree <- root(tree_nodes, species)
    } else if (root_style == "Midpoint") {
      tree <- midpoint(tree_nodes)
    } else {
      tree <- tree_nodes
    }
    
    # Bootstrapping
    if (bootstrap_value == TRUE) {
      if (input$tree_type == "Maximum Likelihood") {
        bs_tree <- bootstrap.pml(ml_tree, bs = 10, trees = TRUE, optNni = TRUE)
      } else {
        if (seq_type == "AA") {
          bs_tree <- bootstrap.phyDat(phy, FUN = function(x)NJ(dist.ml(x, model = "JTT")), bs = 100)
        } else {
          bs_tree <- bootstrap.phyDat(phy, FUN = function(x)NJ(dist.ml(x, model = "JC69")), bs = 500)
        }
      }
    } 
    
    # Plot Tree 
    check_print <- input$print_tree
    if (check_print == "Phylogenetic Tree") {
      plot_tree <- plot(tree_nodes)
      my_trees$basic_tree <- tree_nodes
    } else {
      plot_tree <- plotBS(tree_nodes, bs_tree, type = "phylogram")
      my_trees$bootstrap_tree <- bs_tree
    }
    
    
    return(plot_tree)
    
  })
  
  # Download file 
  output$download_data = downloadHandler (
    filename = function() {"plot.pdf"},
    content = function(file) {
      pdf(file, onefile = TRUE)
      if (input$bootstrap_switch == TRUE && input$print_tree == "Bootstrapped Phylogenetic Tree") {
        plotBS(my_trees$basic_tree, my_trees$bootstap_tree, type = "phylogram")
      } else {
        plot(my_trees$basic_tree) 
      }
      dev.off()
    }
  )
  
  
}