library(shiny)
library(tidyverse)
library(magrittr)
library(DT)
library(tibble)
library(shinythemes)
library(heatmaply)
library(plotly)

### load aBSREL data

load(file = "data/results_pval_count_all.RData")
#load(file = "Shiny_Acomys/data/results_pval_count_all.RData")

### load MEME data

load(file = "data/MEME_simple_nsites_subs_EBF.RData")
#load(file = "Shiny_Acomys/data/MEME_simple_nsites_subs_EBF.RData")

name_conversion <- read_csv("data/name_conversion.csv")
#name_conversion <- read_csv("Shiny_Acomys/data/name_conversion.csv")

### rename data

aBSREL_data <- results_tibble_count_cor
rm(results_tibble_count_cor)

MEME_data <- MEME_simple_nsites_subs_EBF
rm(MEME_simple_nsites_subs_EBF)

### get genelist and transcript IDs for aBSREL results

genenames <- aBSREL_data$genename
transcript_ids <- aBSREL_data$transcript_id
# get genelist of acomys sign genes only
aBSREL_sign <- aBSREL_data %>%
  filter(`P_value` <= 0.05)
genenames_sign <- aBSREL_sign %>%
  select(genename) %>%
  use_series(genename)
transcript_ids_sign <- aBSREL_sign %>%
  select(transcript_id) %>%
  use_series(transcript_id)
rm(aBSREL_sign)

### get genenames and transcript IDs for MEME results
MEME_genenames <- MEME_data$genename
MEME_transcript_ids <- MEME_data$transcript_id

# Define UI ----
ui <- navbarPage(
  title = HTML(paste(em("Acomys cahirinus"), "study")),
  theme = shinytheme("flatly"),
  tabPanel("aBSREL results", 
           fluidRow(
             column(3,
                    wellPanel(
                      h3("Gene selection"),
                      selectizeInput(inputId = "aBSRELGeneInput", "Select the gene for which you want to see the results", 
                                     choices = NULL),
                      checkboxInput("aBSRELSignSwitch", HTML(paste("Only allow selecting genes with significant signs of episodic",
                                                              "diversifying selection in", em("Acomys cahirinus"))), value = FALSE),
                      checkboxInput("aBSRELTranscriptIDSelect", "Choose genes based on their transcript IDs", value = FALSE)),
                    wellPanel(
                      h3("Table options"),
                      checkboxInput("aBSRELTableSignBranches", "Only show significant Branches", value = FALSE)
                    )),
             column(9,
                    uiOutput("aBSRELTableDescription"),
                    br(),
                    DT::dataTableOutput("aBSREL_table")
                    )
             )
           ),
  tabPanel("MEME results",
                    fluidRow(
                      column(3,
                             wellPanel(
                               h3("Gene selection"),
                               selectizeInput(inputId = "MEMEGeneInput", "Select the gene for which you want to see the results", 
                                              choices = NULL),
                               checkboxInput("MEMETranscriptIDSelect", "Choose genes based on their transcript IDs", value = FALSE)),
                             wellPanel(
                               h3("Table options"),
                               checkboxInput("MEMETableSignSites", "Only show significant sites", value = FALSE)
                             )),
                      column(9,
                             tabsetPanel(
                               tabPanel("Test results",
                                        br(),
                                        uiOutput("MEMETableDescription"),
                                        br(),
                                        p(paste("This HyPhy MEME analysis was a follow-up analysis after",
                                                "the HyPhy aBSREL analysis. During the aBSREL analysis,",
                                                "genes were found that exhibited signs of episodic diversifying",
                                                "selection in Acomys cahirinus. However, HyPhy aBSREL",
                                                "detects episodic diversifying selection in specific branches of",
                                                "a phylogeny without informing us which sites are responsible.",
                                                "To get an indication of which codon sites could be under episodic",
                                                "diversifying selection in the genes found by aBSREL, HyPhy MEME",
                                                "was used. HyPhy MEME is a mixed effects model of evolution (MEME)",
                                                "that can detect episodic diversifying selection at an individual",
                                                "site (among all branches). Thus, HyPhy MEME might give us additional",
                                                "insights into which codon sites might be responsible for the episodic",
                                                "diversifying selection signal in the genelist found by aBSREL.",
                                                "A word of caution is advised however, as HyPhy MEME can only detect",
                                                "if a codon site shows signs of diversifying selection among all branches",
                                                "(i.e. MEME detects selection at an individual site, not an individual",
                                                "branch-site). So, HyPhy MEME predicting that a codon site shows",
                                                "significant diversifying selection does not guarantee that any individual",
                                                "branch (like the Acomys cahirinus branch) exhibits diversifying selection",
                                                "at that codon site. nevertheless, HyPhy MEME can give us an indication", 
                                                "which sites might be involved in the episodic diversifying selection signal",
                                                "in the Acomys cahirinus genes, even when this limitation is kept in mind.")),
                                        p(paste("Every Gene was individually analyzed using HyPhy MEME. The test results are",
                                                "described in the table below, where the column names represent the following:")),
                                        tags$ul(
                                          tags$li(HTML(paste(strong(HTML("&alpha;")), ": Synonymous substitution rate at a site."))),
                                          tags$li(HTML(paste(strong(HTML("&beta;<sup>-</sup>")), ": Non-synonymous substitution rate at a site for the negative/neutral evolution component."))),
                                          tags$li(HTML(paste(strong(HTML("p<sup>-</sup>")), ": Mixture distribution weight allocated to &beta;<sup>-</sup>; loosely -- the proportion of the tree evolving neutrally or under negative selection."))),
                                          tags$li(HTML(paste(strong(HTML("&beta;<sup>+</sup>")), ": Non-synonymous substitution rate at a site for the positive/neutral evolution component."))),
                                          tags$li(HTML(paste(strong(HTML("p<sup>+</sup>")), ": Mixture distribution weight allocated to &beta;<sup>+</sup>; loosely -- the proportion of the tree evolving neutrally or under positive selection."))),
                                          tags$li(HTML(paste(strong("LRT"), ": Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;."))),
                                          tags$li(HTML(paste(strong("p-value"), ": Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> &gt; 0 <emph>and<emph> &beta;<sup>+</sup> &gt; &alpha;."))),
                                          tags$li(HTML(paste(strong("number of branches under selection"), ": The (very approximate and rough) estimate of how many branches may have been under selection at this site, i.e., had an empirical Bayes factor of 100 or more for the &beta;<sup>+</sup> rate."))),
                                          tags$li(HTML(paste(strong("Total branch length"), ": The total length of branches contributing to inference at this site, and used to scale dN-dS."))),
                                          tags$li(HTML(paste(strong("MEME LogL"), ": Site Log-likelihood under the MEME model."))),
                                          tags$li(HTML(paste(strong("FEL LogL"), ": Site Log-likelihood under the FEL model."))),
                                          tags$li(HTML(paste(strong("Variation p"), ": Asymptotic p-value for whether or not there is evidence of dN/dS variation across branches."))),
                                        ),
                                        br(),
                                        DT::dataTableOutput("MEME_table"),
                                        hr(),
                                        uiOutput("MEMEResultDescription"),
                                        hr()
                                        ),
                               tabPanel("Multiple sequence alignment",
                                        h3("Multiple sequence alignment"),
                                        p(paste("This tab displays the protein multiple sequence alignment (MSA),",
                                                "which were created by translating the TOGA codon MSAs.")),
                                        p("Description of the MSAs:"),
                                        tags$ul(
                                          tags$li("Blue columns represent columns where over 50% of the species in the MSA had an identical amino acid."),
                                          tags$li("'X' represents a residue that was a 'NNN' codon or a codon gap. Keep in mind that this MSA represents the aBSREL/MEME input sequences, where large gap regions had already been removed. Thus, the sequences represented in these MSAs are shorter than the actual sequence."),
                                          tags$li("The red arrows point at columns that have evidence of episodic diversifying selection, according to HyPhy MEME.")),
                                        uiOutput("MultipleAlignment"),
                                        hr()
                                        ),
                               tabPanel("Substitutions",
                                        h3("Predicted substitutions from the last common mammalian ancestor to Acomys cahirinus at diversifying codon sites"),
                                        p(HTML(paste("As", em("Acomys cahirinus"), "was the subject of the study, all predicted substitutions are listed in the table below."))),
                                        hr(),
                                        DT::dataTableOutput("MEME_acomys_sub"),
                                        hr(),
                                        h3("All substitutions at diversifying codon sites"),
                                        p("This table displays the predicted substitutions (by HyPhy MEME) at diversifying sites of all tested species. Keep in mind that predicted deletions and insertions in this table are not neccesarily true, as 'NNN' codons are used to label a substitution as a deletion or an insertion. It is assumed that 'NNN' codons represent gaps, but a minority of 'NNN' codons are actually masked inactivating mutations. Thus, this could lead to the misinterpretation of an 'NNN' codon."),
                                        hr(),
                                        DT::dataTableOutput("MEME_branch_subs"),
                                        hr()
                                        ),
                               tabPanel("EBF",
                                        plotlyOutput("EBF_heatmap"),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        br(),
                                        fluidRow(
                                          column(3,
                                                 checkboxInput("MEME_EBF_inf_remove", "Remove infinite EBF values", value = TRUE)),
                                          column(9,
                                                 checkboxInput("MEME_EBF_acomys_select", "Only select Acomys cahirinus EBF values", value = FALSE))),
                                        DT::dataTableOutput("EBF_table")
                                        )
                             )
                      )
                    )
           )
  )

# Define server logic ----
server <- function(input, output, session) {
  
  ### Backend for aBSREL results page
  
  observe({

    # freezes all reactive expressions until everything has been updated.
    # Without this, there is one second where it tries to treat gene symbols as
    # transcript IDs or vice versa, leading to errors.
    freezeReactiveValue(input, "aBSRELGeneInput")
    
    if (input$aBSRELTranscriptIDSelect) {
      if (input$aBSRELSignSwitch) {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = transcript_ids_sign, server = TRUE,
                             options = list(placeholder = "Transcript id"))
      } else {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = transcript_ids, server = TRUE,
                             options = list(placeholder = "Transcript id"))
      }
    } else {
      if (input$aBSRELSignSwitch) {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = genenames_sign, server = TRUE,
                             options = list(placeholder = "Gene name"))
      } else {
        updateSelectizeInput(session, "aBSRELGeneInput", choices = genenames, server = TRUE,
                             options = list(placeholder = "Gene name"))
      }
    }

  })
  
  
  # retrieve which gene has been selected
  get_genename <- reactive({
    
    aBSREL_input <- input$aBSRELGeneInput
    
    # return nothing when no gene has been chosen
    if (aBSREL_input == "" ) {
      return("")
    }
    
    aBSREL_data %>%
      {if (input$aBSRELTranscriptIDSelect) filter(., transcript_id == aBSREL_input) else filter(., genename == aBSREL_input)} %>%
     use_series(genename)
    
  })

  output$aBSRELTableDescription <- renderUI({
    
    genename <- get_genename()
    
    if (genename == "") {
      # If no gene has been selected
      strong("Select a gene to display its results")
    } else {
      # If a gene has been selected
      strong(paste("aBSREL results of", genename))
    }
  })
  
  get_aBSREL_results <- reactive({
    
    genename_chosen <- get_genename()
    
    # Makes sure nothing gets generated when no gene has been selected
    if (genename_chosen == "") {
      return(NULL)
    }
    
    # generation of datatable
    aBSREL_data %>%
      filter(genename == genename_chosen) %>%
      use_series(test_results) %>%
      extract2(1) %>%
      rename("Branch name" = "branch_name", "Likelyhood ratio test statistic" = "LRT",
             "Uncorrected P-value" = "p_value", "Corrected P-value" = "p_value_cor",
             "Uncorrected P-value significant (P-value < 0.05)" = "significant",
             "Corrected P-value significant (P-value < 0.05)" = "significant_cor") %>%
      # remove non-significant branches of input$aBSRELTableSignBranches is TRUE
      {if (input$aBSRELTableSignBranches) 
        filter(., `Corrected P-value significant (P-value < 0.05)` == TRUE) else .} %>%
      # change assembly names to species names
      mutate(`Branch name` = map_chr(`Branch name`, ~ {
        if (.x == "REFERENCE") {
          return("Homo_sapiens")
        } else if (!(.x %in% name_conversion$Assembly.name)) {
          return(.x)
        } else if (.x %in% name_conversion$Assembly.name) {
          index <- which(.x == name_conversion$Assembly.name)
          species_tree <- name_conversion$Species_Tree[index]
          return(species_tree)
        } else {
          cat("ERROR! non-recognized branch name!")
          return(NULL)
        }
      }))
  })
  
  get_aBSREL_results_length <- reactive({
    
    nrow(get_aBSREL_results())
  })
  
  output$aBSREL_table <- DT::renderDataTable({
    
    # Makes sure no table gets generated when no gene has been chosen
    if (get_genename() == "") {
      return(tibble())
    }
    
    get_aBSREL_results()
    
  }, escape = FALSE, rownames = FALSE, 
  options = list(scrollY = "500px", scrollX = TRUE, 
                 lengthMenu = c(50, 100, 250, get_aBSREL_results_length()))
  )
  
  ### Backend for MEME results page
  
  observe({
    
    freezeReactiveValue(input, "MEMEGeneInput")
    
    if (input$MEMETranscriptIDSelect) {
      
        updateSelectizeInput(session, "MEMEGeneInput", choices = MEME_transcript_ids, server = TRUE, 
                             options = list(placeholder = "Transcript id"))
    } else {
        updateSelectizeInput(session, "MEMEGeneInput", choices = MEME_genenames, server = TRUE, 
                             options = list(placeholder = "Gene name"))
    }
  })
  
  # retrieve which gene has been selected
  get_genename_MEME <- reactive({
    
    # return nothing when no gene has been chosen
    if (input$MEMEGeneInput == "") {
      return("")
    }

    MEME_data %>%
      {if (input$MEMETranscriptIDSelect) filter(., transcript_id == input$MEMEGeneInput) else filter(., genename == input$MEMEGeneInput)} %>%
      use_series(genename)
  })
  
  ### MEME backend for 'test results' tabpanel
  
  output$MEMETableDescription <- renderUI({
    
    if (get_genename_MEME() == "") {
      # If no gene has been selected
      paragraph <- strong("Select a gene to display its results")
      paragraph
    } else {
      # If a gene has been selected
      paragraph <- paste0(strong(paste("MEME results of", get_genename_MEME())))
      HTML(paragraph)
    }
  })
  
  # extract the MEME_data row from the selected gene
  get_MEME_data_row <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    MEME_data %>%
      filter(genename == get_genename_MEME())
  })
  
  get_MEME_results <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }

    MEME_row <- get_MEME_data_row()
    
    MEME_row %>%
      use_series(csv) %>%
      extract2(1) %>%
      # add an additional site column
      add_column(site = seq(nrow(.))) %>%
      relocate(site) %>%
      # remove non-significant sites if input$MEMETableSignSites is TRUE
      {if (input$MEMETableSignSites) filter(., `p-value` <= 0.05) else .} 
  })
  
  get_MEME_results_length <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    nrow(get_MEME_results())
  })
  
  output$MEME_table <- DT::renderDataTable({
    
    # Makes sure no table gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    get_MEME_results()
    
  # escape = FALSE allows HTML tags to be used in the table
  # options = list(scrollY = "500px", scrollX = TRUE) shortens the table and 
  # scrollX allows scrolling in the horizontal direction
  # lengthMenu picks the selectable lengths of the datatable
    
  }, escape = FALSE, rownames = FALSE, 
  options = list(scrollY = "500px", scrollX = TRUE, 
                 lengthMenu = c(50, 100, 250, get_MEME_results_length()))
  )
  
  # render block of summarizing text 
  output$MEMEResultDescription <- renderUI({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }

    MEME_row <- get_MEME_data_row()

    # get number of significant sites
    n_of_sites <- MEME_row %>%
      use_series(n_of_sites) %>%
      as.character()
    # get significant site location
    sites <- MEME_row %>%
      use_series(sites) %>%
      extract2(1) %>%
      as.character() %>%
      paste(collapse = ", ")
    
    paragraph <- paste0(h3(HTML(paste(em(get_genename_MEME()), "result summary:"))), 
                        p(HTML(paste(get_genename_MEME(), "was found to have", strong(n_of_sites),
                                strong("codon sites"), "with significant signs of episodic diversifying",
                                "selection"))),
                        p(HTML(paste("These codon sites are", strong(sites), ".")))
                        )
    HTML(paragraph)
  }) 
  
  ### MEME backend for 'Multiple sequence alignment' tabpanel
  
  # Render the multiple sequence alignment
  output$MultipleAlignment <- renderUI({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    # create URL
    gene_URL <- paste0("MEME_protein_alignments_v2/", get_genename_MEME(), "_alignment.pdf")
    # create pdf iframe
    tags$iframe(style="height:900px; width:100%; scrolling=yes", 
                src = gene_URL)
  })
  
  ### MEME backend for 'Substitutions' tabpanel
  
  # Get acomys substitution data from MEME_data
  get_MEME_acomys_sub <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    MEME_row <- get_MEME_data_row()
    
    MEME_row %>%
      use_series(substitutions) %>%
      extract2(1) %>%
      rename("Acomys_cahirinus" = "HLacoCah1")
  })
  
  get_MEME_acomys_sub_length <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    nrow(get_MEME_acomys_sub())
  })
  
  output$MEME_acomys_sub <- DT::renderDataTable({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    get_MEME_acomys_sub()
  }, rownames = FALSE, options = list(scrollY = "300px", scrollX = TRUE, 
                                      lengthMenu = c(get_MEME_acomys_sub_length(), 20))
  )
  
  get_MEME_branch_subs <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    MEME_row <- get_MEME_data_row()
    
    MEME_row %>%
      use_series(branch_subs) %>%
      extract2(1) %>%
      select(!node) %>%
      rename("Branch name" = "label", "Substitutions" = "subs") %>%
      # change assembly names to species names
      mutate(`Branch name` = map_chr(`Branch name`, ~ {
        if (.x == "REFERENCE") {
          return("Homo_sapiens")
        } else if (!(.x %in% name_conversion$Assembly.name)) {
          return(.x)
        } else if (.x %in% name_conversion$Assembly.name) {
          index <- which(.x == name_conversion$Assembly.name)
          species_tree <- name_conversion$Species_Tree[index]
          return(species_tree)
        } else {
          cat("ERROR! non-recognized branch name!")
          return(NULL)
        }
    }))
  })
  
  get_MEME_branch_subs_length <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    nrow(get_MEME_branch_subs())
  })
  
  output$MEME_branch_subs <- DT::renderDataTable({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    get_MEME_branch_subs()
  }, rownames = FALSE, options = list(scrollY = "500px", scrollX = TRUE, 
                                      lengthMenu = c(get_MEME_branch_subs_length(), 20))
  )
  
  ### MEME backend for 'EBF' tabpanel
  
  # retrieve EBF table
  get_EBF_data <- reactive({
    
    MEME_row <- get_MEME_data_row()
    
    MEME_row %>%
      use_series(EBF_table) %>%
      extract2(1) %>%
      rename(`Branch assembly name` = branch) %>%
      # change assembly names to species names
      mutate(`Branch name` = map_chr(`Branch assembly name`, ~ {
        if (.x == "REFERENCE") {
          return("Homo_sapiens")
        } else if (!(.x %in% name_conversion$Assembly.name)) {
          return(.x)
        } else if (.x %in% name_conversion$Assembly.name) {
          index <- which(.x == name_conversion$Assembly.name)
          species_tree <- name_conversion$Species_Tree[index]
          return(species_tree)
        } else {
          cat("ERROR! non-recognized branch name!")
          return(NULL)
        }
      })) %>%
      relocate(site, `Branch assembly name`, `Branch name`, EBF)
  })
  
  MEME_EBF_inf_remover <- reactive({
    get_EBF_data() %>%
      {if (input$MEME_EBF_inf_remove) filter(., !(is.infinite(EBF))) else .}
  })
  
  MEME_EBF_acomys_selecter <- reactive({
    MEME_EBF_inf_remover() %>%
      {if (input$MEME_EBF_acomys_select) filter(., `Branch name` == "Acomys_cahirinus") else .}
    
  })
  
  # get EBF_table length
  get_EBF_data_length <- reactive({
    nrow(MEME_EBF_acomys_selecter())
  })
  
  # render heatmap
  output$EBF_heatmap <- renderPlotly({
    
    # remove infinite rows
    EBF_heatmap <- get_EBF_data() %>%
      filter(!(is.infinite(EBF)))
    
    # limit heatmap to the following species
    included_assemblies <- c("REFERENCE", "HLpsaObe1", "HLmerUng1", "HLratNor7", "rn6", 
                             "HLmusPah1", "HLmusCar1", "mm10", "mm39", "HLacoCah1",
                             "HLonyTor1", "HLperManBai2", "HLsigHis1", "HLellLut1", 
                             "HLondZib1", "HLcriGri3", "HLmesAur2", "mesAur1")
    
    # filter for included assemblies
    EBF_heatmap %<>%
      filter(`Branch assembly name` %in% included_assemblies) %>%
      rename(name = site, variable = `Branch name`, value = EBF) %>%
      select(!`Branch assembly name`)
    
    # perform a log transformation, making the data more normal, improving the
    # visualization
    EBF_heatmap %<>%
      mutate(value = log10(value))
    
    # generate heatmap
    heatmaply(long_data = EBF_heatmap,
              dendrogram = "none",
              xlab = "", ylab = "",
              main = "",
              scale = "none",
              label_names = c("Codon", "Branch", "log10(EBF)"),
              showticklabels = c(TRUE, FALSE),
              fontsize_col = 8) %>%
      layout(width = 700, height = 1000)
  })
  
  # render datatable
  output$EBF_table <- DT::renderDataTable({
    
    MEME_EBF_acomys_selecter()
  }, rownames = FALSE, options = list(scrollY = "500px", 
                                      lengthMenu = c(20, 50, 100, 500, 1000)) 
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)