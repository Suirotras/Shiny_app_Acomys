library(shiny)
library(tidyverse)
library(magrittr)
library(DT)
library(tibble)
library(shinythemes)
library(heatmaply)
library(plotly)
library(gt)
library(ggprism)
library(scales)

### load aBSREL data
load(file = "data/aBSREL_data.RData")

### load name conversion file
name_conversion <- read_csv("data/name_conversion.csv")

### Load the STRING cluster data
acomys_clusters <- read_csv("data/acomys_bias_one2one--clustered--infl4.csv")

### Load the aBSREL dN/dS dataframe
dNdS_df <- read_delim("data/aBSREL_dNdS_df.tsv",
                      delim = "\t", escape_double = FALSE, trim_ws = TRUE)

### Load the functional enrichments
clusterprofiler_allgenes_GO_KEGG <- read_delim("data/functional_enrichments/clusterprofiler_results_GO_KEGG.tsv",
                                               delim = "\t", escape_double = FALSE,
                                               trim_ws = TRUE)

### Load the cluster-specific functional enrichments
enrichGO_aco <- read_delim("data/functional_enrichments_cluster/enrichGO_aco_all.tsv",
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)
Acomys_reduced_GO_clusters <- read_delim("data/functional_enrichments_cluster/Acomys_clusters_GO_df.tsv",
                                         delim = "\t", escape_double = FALSE,
                                         trim_ws = TRUE)

### Load the uniprot-mapped substitutions
MEME_subs <- read_delim("data/Uniprot_mapping/MEME_subs_mapped_sitefiltered_noMSAs.tsv",
                        delim = "\t", escape_double = FALSE,
                        col_types = cols(n_of_sites = col_integer(),
                                         sites = col_character()), trim_ws = TRUE)

### Create whitespace function
generate_br_tags <- function(n) {
  tagList(replicate(n, br(), simplify = FALSE))
}

### load the genenames and corresponding ensemble IDs for the significant and
### non-significant genes

cbind_genelists <- function(genelist_paths, col_names) {
  # Function that reads the genelists and binds the columns in a single dataframe.
  #   args:
  #     genelist_paths: A vector of the filepaths of the genelists.
  #     col_names: A vector of the names that the columns should have in the
  #                combined datraframe. Should be in the same order as genelist_paths.
  require(readr)
  require(dplyr)
  
  dfs <- lapply(seq_along(genelist_paths), function(i) {
    read_csv(genelist_paths[i], col_names = col_names[i], show_col_types = FALSE)
  })
  return(do.call("cbind", dfs))
}

### get function arguments (filepaths and column names)
genelist_path <- file.path("data", "genelists_bias_one2one")
genelist_sign_path <- c(file.path(genelist_path, "aBSREL_human_gene_ids_sign.txt"),
                        file.path(genelist_path, "aBSREL_human_genenames_sign.txt"),
                        file.path(genelist_path, "aBSREL_human_transcript_ids_sign.txt"))
genelist_path <- file.path("data", "genelists_one2one")
genelist_all_path <- c(file.path(genelist_path, "aBSREL_human_gene_ids.txt"),
                       file.path(genelist_path, "aBSREL_human_genenames.txt"),
                       file.path(genelist_path, "aBSREL_human_transcript_ids.txt"))
col_names <- c("ensemble", "names", "ensembleTrans")

### get genelist dataframes
sign_genes <- cbind_genelists(genelist_sign_path, col_names)
all_genes <- cbind_genelists(genelist_all_path, col_names)

### get genenames and transcript IDs for aBSREL
# all tested genes
genenames <- all_genes$names
transcript_ids <- all_genes$ensembleTrans
# sign genes
genenames_sign <- sign_genes$names
transcript_ids_sign <- sign_genes$ensembleTrans

### get genenames and transcript IDs for MEME results
MEME_genenames <- sign_genes$names
MEME_transcript_ids <- sign_genes$ensembleTrans

# Define UI ----
ui <- navbarPage(
  title = HTML(paste(em("Acomys cahirinus"), "study")),
  theme = shinytheme("flatly"),
  tabPanel("aBSREL results", 
           fluidRow(
             tabsetPanel(
               tabPanel("LRT results",
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
                               DT::dataTableOutput("aBSREL_table"),
                        )),
               tabPanel("STRING interaction network",
                        column(3,
                               wellPanel(
                                 h4("Display options for STRING clusters"),
                                 checkboxInput("STRINGplotSwitch", HTML("Display all clusters with at least three genes"), value = FALSE),
                                 br(),
                                 sliderInput("STRINGClustTableSlider", "Cluster range:",
                                             min = 1, max = max(enrichGO_aco$`\`__mclCluster\``,
                                                                na.rm = TRUE),
                                             value = c(1, 8)),
                                 checkboxInput("STRINGDisplayNonClust", HTML("Display non-clustered genes in table"), value = FALSE),
                                 actionButton("STRINGClustTableAction", "Refresh")
                               )),
                        column(9,
                               uiOutput("aBSRELSTRINGplot"),
                               br(),
                               gt_output("STRINGClusterTable"))),
               tabPanel("Functional overrepresentation analysis",
                        column(3,
                               wellPanel(
                                 h4("Display options for functionel enrichment of full aBSREL genelist"),
                               ),
                               generate_br_tags(40),
                               wellPanel(
                                 h4("Display options for functionel enrichment of STRING clusters"),
                                 checkboxInput("FunEnrichReducedSwitch", HTML("Display all GO-terms per cluster instead of only the reduced GO-terms"), value = FALSE),
                                 h4("Select which clusters to display"),
                                 sliderInput("FunEnrichSelectClust", "Cluster range:",
                                             min = 1, max = max(enrichGO_aco$`\`__mclCluster\``,
                                                                na.rm = TRUE),
                                             value = c(1, 8)),
                                 actionButton("FunEnrichClustAction", "Refresh")
                               )),
                        column(9,
                               h5("This page contains the functional overrepresentation analysis results of both complete aBSREL positively selected genelist and the cluster-specific overrepresentation results"),
                               hr(),
                               strong("Functional overrepresentation analysis of complete positively selected genelist"),
                               br(),
                               plotOutput("FunEnrichaBSRELAll", width = "100%", height = "300px"),
                               br(),
                               gt_output("FunEnrichaBSRELAllTable"),
                               hr(),
                               strong("Functional overrepresentation analysis of the individual positively selected STRING clusters"),
                               br(),
                               # The inline = TRUE is required so that the 
                               plotOutput("FunEnrichaBSRELClusters", inline = TRUE),
                               br(),
                               gt_output("FunEnrichaBSRELClustersTable")
                               ))
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
                               checkboxInput("MEMETableSignSites", "Only show significant sites", value = FALSE),
                               numericInput("p_val_select", "p-value significance cut-off", value = 0.1)
                             )),
                      column(9,
                             tabsetPanel(
                               tabPanel("LRT results",
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
                                                "described in the tables below, where the column names represent the following:")),
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
                                        hr(),
                                        p(strong(HTML(paste("This table and summary data represents the MEME likelyhood ratio test results for each site of the",
                                                "selected gene. Here, the positional site information is relative to the human gene",
                                                "that was used in the TOGA multiple sequence alignment (i.e. the multiple sequence",
                                                "alignment with mammalian orthologs of the selected gene, available at the",
                                                a("TOGA ortholog database", href="https://genome.senckenberg.de//download/TOGA/human_hg38_reference/MultipleCodonAlignments/")
                                                )))),
                                        br(),
                                        DT::dataTableOutput("MEME_table"),
                                        br(),
                                        uiOutput("MEMEResultDescription"),
                                        hr(),
                                        p(strong(HTML(paste("This table and summary data represents the data that could be mapped to a canonical Uniprot protein sequence.",
                                                            "Here, the positional site information is relative to the canonical protein sequence of the human Uniprot accession"
                                        )))),
                                        br(),
                                        DT::dataTableOutput("MEME_table_unimapped"),
                                        br(),
                                        uiOutput("MEMEResultDescription_unimapped"),
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
                                        generate_br_tags(28),
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
  
  # Selection of STRINGp cluster plot
  output$aBSRELSTRINGplot <- renderUI({
    if (input$STRINGplotSwitch) {
      tags$img(src = "./images/stringdb_all_clusters_with_three_genes.svg", width = 510 * 1.6, height = 472 * 1.6)
    } else {
      tags$img(src = "./images/Acomys_stringdb_clusters.svg", width = 510 * 1.5, height = 472 * 1.5)
    }
  })
  
  # Generate table of aBSREL info for each STRING cluster
  STRINGClustTable <- eventReactive(input$STRINGClustTableAction, {
    
    # Filter for one2one significant A. cahirinus genes
    aBSREL_table <- aBSREL_data %>%
      dplyr::filter(genename %in% sign_genes$names) %>%
      # filter for the presence of orthologs from nearest phylogenetic neighbors
      # of A. cahirinus
      dplyr::filter(MRCA_present_aco == TRUE) %>%
      dplyr::select(-all_of(c("Batch", "!_present", "JSON_exists",
                              "MRCA_present_aco", "MRCA_present_mus",
                              "significant_prop", "significant_prop",
                              "branches_tested")))
    
    # Add ensemble gene IDs to dataframe
    aBSREL_table <- all_genes %>%
      as_tibble() %>%
      select(-names) %>%
      rename(gene_id = ensemble, transcript_id = ensembleTrans) %>%
      right_join(aBSREL_table)
    
    # Join the omega (i.e. dN/dS) data with the aBSREL LRT data
    aBSREL_table_omega <- dNdS_df %>%
      dplyr::filter(species == "HLacoCah2") %>%
      dplyr::select(-all_of(c("transcript_id", "P_value", "MRCA_present_aco", "MRCA_present_mus",
                              "species"))) %>%
      right_join(aBSREL_table)
    
    # Calculate additional dN/dS values
    aBSREL_table_omega <- aBSREL_table_omega %>%
      # calculate and retrieve max omega values
      mutate(max_omega = if_else(rate_class_number == 1, omega_1,
                                 if_else(rate_class_number == 2, omega_2,
                                         if_else(rate_class_number == 3, omega_3, NA)))) %>%
      # calculate and retrieve max omega proportions
      mutate(max_omega_prop = if_else(rate_class_number == 1, omega_1_prop,
                                      if_else(rate_class_number == 2, omega_2_prop,
                                              if_else(rate_class_number == 3, omega_3_prop, NA)))) %>%
      # calculate and retrieve mean omega values
      mutate(mean_dNdS = ifelse(rate_class_number == 1, omega_1,
                                if_else(rate_class_number == 2, 
                                        omega_1*omega_1_prop+omega_2*omega_2_prop,
                                        if_else(rate_class_number == 3, 
                                                omega_1*omega_1_prop+omega_2*omega_2_prop+omega_3*omega_3_prop, NA))))
    
    # Retrieve the likelyhood ratio test statistic from the 'test_results'
    # list-column
    aBSREL_table_omega_test <- aBSREL_table_omega %>%
      rowwise() %>%
      mutate(test_results = dplyr::filter(test_results, branch_assembly == "HLacoCah2")) %>%
      unnest(test_results) %>%
      ungroup() %>%
      dplyr::select(-all_of(c("branch_name", "branch_projection", "p_value_cor",
                              "significant", "significant_cor"))) %>%
      # give Columns more representative names
      dplyr::rename(p_value_uncorrected = p_value, FDR = P_value)
    
    # add species name to data
    aBSREL_table_omega_test <- name_conversion %>%
      dplyr::select(`Assembly name`, Species_Tree) %>%
      rename(branch_assembly = `Assembly name`) %>%
      right_join(aBSREL_table_omega_test) %>%
      dplyr::select(-branch_assembly)
    
    # Add STRINGdb clusters to data
    aBSREL_table_omega_test_clusters <- acomys_clusters %>%
      rename(Cluster = `__mclCluster`, gene_id = `query term`) %>%
      dplyr::select(gene_id, Cluster) %>%
      right_join(aBSREL_table_omega_test)
    
    # Filter for selected cluster range
    if (input$STRINGDisplayNonClust) {
      aBSREL_table_omega_test_clusters <- aBSREL_table_omega_test_clusters %>%
        dplyr::filter(Cluster %in% seq(input$STRINGClustTableSlider[1],
                                       input$STRINGClustTableSlider[2]) |
                        is.na(Cluster))
    } else {
      aBSREL_table_omega_test_clusters <- aBSREL_table_omega_test_clusters %>%
        dplyr::filter(Cluster %in% seq(input$STRINGClustTableSlider[1],
                                       input$STRINGClustTableSlider[2]))
    }
    
    # Generate table for likelyhood ratio test statistics and model estimations
    # per gene and cluster
    gt_table <- aBSREL_table_omega_test_clusters %>%
      # Make custom cluster labels
      dplyr::mutate(Cluster_label = if_else(is.na(Cluster), "Not clustered",
                                            paste0("Cluster ", Cluster)),
                    Cluster_label = factor(Cluster_label,
                                           levels = c(paste0("Cluster ",
                                                             seq(max(Cluster, na.rm = TRUE))),
                                                      "Not clustered"))) %>%
      dplyr::select(all_of(c("Cluster_label", "genename", "LRT",
                             "p_value_uncorrected", "FDR", "full_adapt",
                             "full_adapt_dN", "full_adapt_dS",
                             "mean_dNdS", "max_omega", "max_omega_prop"))) %>%
      # sort by FDR
      dplyr::arrange(Cluster_label, FDR) %>%
      # Start creating the table
      gt(groupname_col = "Cluster_label", row_group_as_column = TRUE) %>%
      tab_header(title = "HyPhy aBSREL output per STRING cluster") %>%
      tab_spanner(label = html("<em><strong>Likelyhood ratio test</em></strong>"),
                  columns = c(LRT, p_value_uncorrected, FDR)) %>%
      tab_spanner(label = html("<strong><em>Model estimation</em></strong>"),
                  columns = c(full_adapt, full_adapt_dN,
                              full_adapt_dS, mean_dNdS, max_omega,
                              max_omega_prop)) %>%
      fmt_number(decimals = 2, sep_mark = "") %>%
      cols_label(contains("p_value_uncorrected") ~ "P-value uncorrected",
                 contains("LRT") ~ "likelyhood ratio test statistic",
                 contains("FDR") ~ "False discovery rate",
                 contains("full_adapt") ~ "Inferred branch length",
                 contains("full_adapt_dN") ~ "Inferred dN branch length",
                 contains("full_adapt_dS") ~ "Inferred dS branch length",
                 contains("mean_dNdS") ~ "Mean dN/dS",
                 contains("max_omega") ~ "maximum dN/dS",
                 contains("max_omega_prop") ~ "proportion of sites with maximum dN/dS") %>%
      tab_stubhead(label = "STRING cluster") %>%
      cols_align(align = "left") %>%
      opt_stylize(style = 1) %>%
      data_color(
        columns = FDR,
        palette = "Greens",
        reverse = TRUE) %>%
      data_color(
        columns = mean_dNdS,
        palette = "Blues"
      ) %>%
      tab_style(
        locations = cells_body(columns = genename),
        style = cell_fill(color = "gray95"))
    
      # Reorder row groups
    if (input$STRINGDisplayNonClust) {
      gt_table %>%
        row_group_order(groups = c(paste0("Cluster ", seq(max(aBSREL_table_omega_test_clusters$Cluster, na.rm = TRUE))), "Not clustered")) %>%
        opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
    } else {
      gt_table %>%
        row_group_order(groups = paste0("Cluster ", seq(max(aBSREL_table_omega_test_clusters$Cluster, na.rm = TRUE)))) %>%
        opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
    }
  }, ignoreNULL = FALSE)
  
  output$STRINGClusterTable <- render_gt({
    STRINGClustTable()
    })
  
  output$FunEnrichaBSRELAll <- renderPlot({
    
    ggplot(clusterprofiler_allgenes_GO_KEGG,
           aes(y = reorder(Description, p.adjust, decreasing = FALSE), 
                      x = p.adjust, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      geom_linerange(aes(xmin = 0, xmax = p.adjust), linewidth = 2) +
      ylab("GO Biological process") + xlab("False discovery rate") +
      theme_prism() +
      scale_size_continuous(name = "number of\ngenes",
                            range = c(3, 10)) +
      scale_y_discrete(labels = label_wrap(30)) +
      theme(axis.title.y = element_blank(), legend.title = element_text()) +
      scale_color_gradientn(colors = MetBrewer::met.brewer("OKeeffe2", type = "continuous"),
                            name = "number of\ngenes") +
      facet_grid(scales = "free_y", rows = vars(source),
                 space = "free_y")
  })
  
  output$FunEnrichaBSRELAllTable <- render_gt({
    gt(clusterprofiler_allgenes_GO_KEGG, groupname_col = "source",
       row_group_as_column = TRUE) %>%
      cols_align(align = "left") %>%
      opt_stylize(style = 1) %>%
      data_color(
        columns = Count,
        palette = "Blues") %>%
      data_color(
        columns = qvalue,
        palette = "Greens",
        reverse = TRUE) %>%
      opt_interactive(use_resizers = TRUE)
      
  })

  get_GOterm_plot_height <- eventReactive(input$FunEnrichClustAction, {
    # Get the number of GO-terms that will be plotted, which will
    # specify the height of the plot.
    if (input$FunEnrichReducedSwitch) {
      n_of_terms <- enrichGO_aco %>%
        dplyr::filter(Cluster %in% seq(input$FunEnrichSelectClust[1],
                                       input$FunEnrichSelectClust[2])) %>%
        nrow()
      20 * n_of_terms
    } else {
      # Get the number of reduced GO-terms that will be plotted, which will
      # specify the height of the plot.
      n_of_terms <- Acomys_reduced_GO_clusters %>%
        dplyr::filter(stringdb_cluster %in% seq(input$FunEnrichSelectClust[1],
                                                input$FunEnrichSelectClust[2])) %>%
        nrow()
      15 * n_of_terms
    }
  }, ignoreNULL = FALSE)
  
  FunEnrichClustPlot <- eventReactive(input$FunEnrichClustAction, {
    
    if (input$FunEnrichReducedSwitch) {
      ggplot(subset(enrichGO_aco, Cluster %in% seq(input$FunEnrichSelectClust[1],
                                                   input$FunEnrichSelectClust[2])),
             aes(y = reorder(Description, -log10(`p.adjust`), decreasing = FALSE), 
                 x = -log10(`p.adjust`), color = Count)) +
        geom_point(aes(size = Count)) +
        geom_linerange(aes(xmin = 0, xmax = -log10(`p.adjust`))) +
        ylab("GO-terms (Biological process)") + xlab("-log10(FDR)") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 12)) +
        scale_size_continuous(name = "number of\ngenes") +
        scale_y_discrete(labels = label_wrap(100)) +
        # reverse color palette using rev()
        scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Cross", type = "continuous")),
                              name = "number of\ngenes") +
        facet_grid(scales = "free_y", rows = vars(Cluster),
                   space = "free_y")
    } else {
      ggplot(subset(Acomys_reduced_GO_clusters,
                    stringdb_cluster %in% seq(input$FunEnrichSelectClust[1],
                                              input$FunEnrichSelectClust[2])),
             aes(y = reorder(term, score, decreasing = FALSE), 
                 x = score, color = n_of_genes)) +
        geom_point(aes(size = terms_per_parent)) +
        geom_linerange(aes(xmin = 0, xmax = score)) +
        ylab("Reduced GO-terms (Biological process)") + xlab("-log10(FDR)") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 12)) +
        scale_size_continuous(name = "number of\nGO terms") +
        scale_y_discrete(labels = label_wrap(100)) +
        # reverse color palette using rev()
        scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Cross", type = "continuous")),
                              name = "number of\ngenes") +
        facet_grid(scales = "free_y", rows = vars(stringdb_cluster),
                   space = "free_y")
    }
  }, ignoreNULL = FALSE)
    
    output$FunEnrichaBSRELClusters <- renderPlot({
      FunEnrichClustPlot()
    },
    width = 1400,
    height = function() {get_GOterm_plot_height()})
    
    FunEnrichClustTable <- eventReactive(input$FunEnrichClustAction, {
      
      # Render reduced GO-terms unless FunEnrichReducedSwitch is set to TRUE
      if (input$FunEnrichReducedSwitch) {
        enrichGO_aco %>%
          dplyr::select(-`\`__mclCluster\``) %>%
          dplyr::filter(Cluster %in% seq(input$FunEnrichSelectClust[1],
                                         input$FunEnrichSelectClust[2])) %>%
          dplyr::mutate(Cluster_label = if_else(is.na(Cluster), "Not clustered",
                                                paste0("Cluster ", Cluster)),
                        Cluster_label = factor(Cluster_label,
                                               levels = c(paste0("Cluster ",
                                                                 seq(max(Cluster, na.rm = TRUE))),
                                                          "Not clustered"))) %>%
          dplyr::select(-Cluster) %>%
          gt(groupname_col = "Cluster_label", row_group_as_column = TRUE) %>%
          opt_stylize(style = 1) %>%
          cols_align(align = "left") %>%
          data_color(
            columns = qvalue,
            palette = "Greens",
            reverse = TRUE) %>%
          data_color(
            columns = Count,
            palette = "Blues") %>%
          opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
      } else {
        Acomys_reduced_GO_clusters %>%
          dplyr::select(-all_of(c("cluster", "parent", "parentTerm", "termDispensability",
                                  "termUniquenessWithinCluster"))) %>%
          dplyr::filter(stringdb_cluster %in% seq(input$FunEnrichSelectClust[1],
                                                  input$FunEnrichSelectClust[2])) %>%
          dplyr::mutate(Cluster_label = if_else(is.na(stringdb_cluster), "Not clustered",
                                                paste0("Cluster ", stringdb_cluster)),
                        Cluster_label = factor(Cluster_label,
                                               levels = c(paste0("Cluster ",
                                                                 seq(max(stringdb_cluster, na.rm = TRUE))),
                                                          "Not clustered"))) %>%
          dplyr::select(-stringdb_cluster) %>%
          dplyr::rename("Representative GO of reduced GO-terms" = go,
                        "-log10(FDR) of representative GO" = score,
                        "representative GO-term size" = size,
                        "GO description" = term,
                        "Term Uniqueness" = termUniqueness,
                        "number of genes in reduced GO-terms" = n_of_genes,
                        "Gene IDs" = geneID,
                        "Number of reduced GO-terms" = terms_per_parent) %>%
          gt(groupname_col = "Cluster_label", row_group_as_column = TRUE) %>%
          opt_stylize(style = 1) %>%
          cols_align(align = "left") %>%
          data_color(
            columns = "-log10(FDR) of representative GO",
            palette = "Greens") %>%
          data_color(
            columns = "number of genes in reduced GO-terms",
            palette = "Blues") %>%
          opt_interactive(use_compact_mode = TRUE, use_filters = TRUE)
      }
    }, ignoreNULL = FALSE)
    
    output$FunEnrichaBSRELClustersTable <- render_gt({
      FunEnrichClustTable()
    })
    
  # retrieve which gene has been selected
  get_genename <- reactive({
    
    aBSREL_input <- input$aBSRELGeneInput
    
    # return nothing when no gene has been chosen
    if (aBSREL_input == "" ) {
      return("")
    }
    
    # select using transcript ID, if that option has been set
    if (input$aBSRELTranscriptIDSelect) {
      all_genes$names[all_genes$ensembleTrans == aBSREL_input]
    } else {
      all_genes$names[all_genes$names == aBSREL_input]
    }
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
        } else if (!(.x %in% name_conversion$`Assembly name`)) {
          return(.x)
        } else if (.x %in% name_conversion$`Assembly name`) {
          index <- which(.x == name_conversion$`Assembly name`)
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
    
    if (input$MEMETranscriptIDSelect) {
      all_genes$names[all_genes$ensembleTrans == input$MEMEGeneInput]
    } else {
      all_genes$names[all_genes$names == input$MEMEGeneInput]
    }
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
    
    # load the required RData object
    load(file.path("data", "MEME_data",
                   paste0("MEME_data_", get_genename_MEME(), ".RData")))
    return(RData_row)
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
      {if (input$MEMETableSignSites) filter(., `p-value` <= input$p_val_select) else .} 
  })
  
  get_MEME_results_length <- reactive({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    nrow(get_MEME_results())
  })
  
  get_unimapped_MEME_results <- reactive({
    
    # Makes sure no table gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    MEME_subs %>%
      dplyr::filter(genename == get_genename_MEME()) %>%
      separate_longer_delim(c(sites, branch_subs, branch_subs_ancestral, `&alpha;`,
                              `&beta;<sup>-</sup>`, `p<sup>-</sup>`,
                              `&beta;<sup>+</sup>`, `p<sup>+</sup>`, `LRT`, `p-value`,
                              pval_fdr, `#_branches_under_selection`, `Total_branch_length`,
                              `MEME_LogL`, `FEL_LogL`, `Variation_p`),
                            delim = stringr::regex("[;,]")) %>%
      # rename columns
      dplyr::rename(site = sites) %>%
      ## format some columns as doubles and integers
      mutate(site = as.integer(site),
             `#_branches_under_selection` = as.integer(`#_branches_under_selection`),
             `&alpha;` = as.double(`&alpha;`),
             `&beta;<sup>-</sup>` = as.double(`&beta;<sup>-</sup>`),
             `p<sup>-</sup>` = as.double(`p<sup>-</sup>`),
             `&beta;<sup>+</sup>` = as.double(`&beta;<sup>+</sup>`),
             `p<sup>+</sup>` = as.double(`p<sup>+</sup>`),
             `LRT` = as.double(`LRT`),
             `p-value` = as.double(`p-value`),
             `pval_fdr` = as.double(`pval_fdr`),
             `Total_branch_length` = as.double(`Total_branch_length`),
             `MEME_LogL` = as.double(`MEME_LogL`),
             `FEL_LogL` = as.double(`FEL_LogL`),
             `Variation_p` = as.double(`Variation_p`)) %>%
      # remove non-significant sites if input$MEMETableSignSites is TRUE
      {if (input$MEMETableSignSites) filter(., `p-value` <= input$p_val_select) else .} 
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

    # filter with p-value <= 0.1 if input$MEMETableSignSites is FALSE
    filtered_results <- get_MEME_results() %>%
      {if (input$MEMETableSignSites == FALSE) filter(., `p-value` <= 0.1) else .}
      
    # get number of significant sites
    n_of_sites <- filtered_results %>%  
      use_series(site) %>%
      length()
    
    # get significant site location
    sites <- filtered_results %>%
      use_series(site) %>%
      as.character() %>%
      paste(collapse = ", ")
    
    paragraph <- paste0(h4(HTML(paste(em(get_genename_MEME()), "result summary:"))), 
                        p(HTML(paste(get_genename_MEME(), "was found to have", strong(n_of_sites),
                                strong("codon sites"), "with significant signs (p <= ", input$p_val_select, ") of episodic diversifying",
                                "selection"))),
                        p(HTML(paste("These codon sites are", strong(sites), ".")))
                        )
    HTML(paragraph)
  })
  
  # Rendering table for uniprot-mapped sites
  output$MEME_table_unimapped <- DT::renderDataTable({
    
    if(nrow(get_unimapped_MEME_results()) == 0) {
      return(NULL)
    }
    
    get_unimapped_MEME_results() %>%
      dplyr::select(-all_of(c("branch_subs_ancestral", "site_included",
                              "uniprot_gn_symbol", "ensembl_gene_id",
                              "pval_fdr", "n_of_sites", "genename",
                              "transcript_id"))) %>%
      rename("Amino acid difference between Human and A. cahirinus sequence" = "branch_subs")
    }, escape = FALSE, rownames = FALSE, 
    options = list(scrollY = "500px", scrollX = TRUE, 
                   lengthMenu = c(25, 50, 100)))
  
  # render block of summarizing text for uniprot-mapped sites 
  output$MEMEResultDescription_unimapped <- renderUI({
    
    # Makes sure nothing gets generated when no gene has been selected
    if (get_genename_MEME() == "") {
      return(NULL)
    }
    
    if(nrow(get_unimapped_MEME_results()) == 0) {
      return(HTML(paste0(h4(HTML(paste("No sites could be mapped to the uniprot accession of gene",
                           em(get_genename_MEME())))))))
    }
    
    # filter with p-value <= 0.1 if input$MEMETableSignSites is FALSE
    filtered_results <- get_unimapped_MEME_results() %>%
      {if (input$MEMETableSignSites == FALSE) filter(., `p-value` <= 0.1) else .}
    
    # get number of significant sites
    n_of_sites <- filtered_results %>%  
      use_series(site) %>%
      length()
    
    # get significant site location
    sites <- filtered_results %>%
      use_series(site) %>%
      as.character() %>%
      paste(collapse = ", ")
    
    uniprot_acc <- filtered_results %>%
      use_series(uniprotswissprot) %>%
      .[1]
    
    paragraph <- paste0(h4(HTML(paste(em(get_genename_MEME()), "result summary:"))), 
                        p(HTML(paste(get_genename_MEME(), "was found to have", strong(n_of_sites),
                                     strong("codon sites"), "with significant signs (p <= ", input$p_val_select, ") of episodic diversifying",
                                     "selection that we unambigiously mapped to the canonical protein sequence of",
                                     get_genename_MEME(), paste0("(", uniprot_acc, ")")
                                     ))),
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
    
    MEME_row <- MEME_row %>%
      use_series(substitutions) %>%
      extract2(1)
    
    HLacoCah2_colname <- colnames(MEME_row)[str_detect(colnames(MEME_row), "HLacoCah2")]
    
    MEME_row %>%
      rename("Acomys_cahirinus" = all_of(HLacoCah2_colname))
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
        } else if (!(.x %in% name_conversion$`Assembly name`)) {
          return(.x)
        } else if (.x %in% name_conversion$`Assembly name`) {
          index <- which(.x == name_conversion$`Assembly name`)
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
        } else if (!(.x %in% name_conversion$`Assembly name`)) {
          return(.x)
        } else if (.x %in% name_conversion$`Assembly name`) {
          index <- which(.x == name_conversion$`Assembly name`)
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
                             "HLmusPah1", "HLmusCar1", "mm10", "mm39", "HLacoCah2",
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