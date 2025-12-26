# ==========================================================
#Librerias
# ==========================================================

library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap) 
library(EnhancedVolcano)
library(tidyr)
library(dplyr)
library(tibble)
library(circlize)
library(plotly)
library(reshape2)

# ==========================================================
#  1. Cargar automáticamente todos los archivos data_*.rds
# ==========================================================

file_list <- list.files(pattern = "^data_.*\\.rds$")
datasets <- lapply(file_list, readRDS)
names(datasets) <- sub("^data_(.*)\\.rds$", "\\1", file_list)

# ==========================================================
#  2. UI
# ==========================================================

ui <- fluidPage(
  
  theme = shinytheme("cerulean"),  # Tema general de la app
  
  # Cabecera con título y links
  navbarPage(
    title = "Explorador de Genes",
    
    tabPanel("Visualización",
             sidebarLayout(
               sidebarPanel(
                 style = "background-color: rgba(255,255,255,0.7); border-radius: 10px;",
                 selectInput("cancer_type", "Tipo de cáncer:",
                             choices = names(datasets),
                             selected = names(datasets)[1]),
                 selectInput("id_type", "Tipo de ID:",
                             choices = c("Gene", "ENSEMBL"), selected = "Gene"),
                 uiOutput("gene_selector"),
                 uiOutput("comparison_selector"),
                 selectInput("expr_type", "Tipo de expresión:",
                             choices = c("Log-normalizada", "Raw counts")),
                 br(),
                 h4("Descargas"),
                 downloadButton("download_table", "Descargar tabla", class = "btn btn-primary"),
                 br(),
                 downloadButton("download_plot", "Descargar gráfica", class = "btn btn-primary")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Expresión del gen",
                            plotOutput("expr_plot", height = "400px"),
                            div(style = "text-align:center;margin:auto;", tableOutput("gene_info")),
                            h4("Información de GENENAME y PATH"),
                            div(style = "text-align:center; margin:auto;", tableOutput("gene_path"))
                   )
                 )
               )
             )
    ),
    
    tabPanel("Visualización Global",
             sidebarLayout(
               sidebarPanel(
                 style = "background-color: rgba(255,255,255,0.7); border-radius: 10px;",
                 selectInput("global_cancer_type", "Tipo de cáncer (global):",
                             choices = names(datasets),
                             selected = names(datasets)[1]),
                 uiOutput("global_comparison_selector"),
                 sliderInput("top_paths", "Número de rutas más frecuentes a mostrar:",
                             min = 5, max = 50, value = 10),
                 br(),
                 h4("Descargas"),
                 downloadButton("download_global_table", "Descargar tabla global", class = "btn btn-primary"),
                 br(),
                 downloadButton("download_global_plot", "Descargar gráfica global", class = "btn btn-primary")
               ),
               mainPanel(
                 plotlyOutput("global_plot", height = "600px"),
                 br(),
                 h4("Rutas Reactome más frecuentes"),
                 tableOutput("path_summary")
               )
             )
    ),
    
    tabPanel("Ayuda",
             fluidRow(
               column(6,
                      HTML("<h4><b>Información</b></h4>"),
                      HTML("<ul>
                         <li><b>Visualización:</b> Permite explorar la expresión de genes individuales en distintos tipos de cáncer. Se puede seleccionar un gen específico y una comparación. 
                             En cada comparación, se analiza el <b>primer nombre vs. segundo nombre</b> del contraste; por ejemplo, <i>Met_Skin_vs_Tumor</i> compara la metástasis de piel frente al tumor primario. Se muestran estadísticas de expresión (media, mediana, SD, logFC, p-ajustada) y anotaciones funcionales como GENENAME, GO y rutas Reactome, junto con gráficos combinando violín y boxplot.</li>
                         <li><b>Visualización Global:</b> Permite ver de manera agregada todas las comparaciones de un tipo de cáncer mediante Volcano plots y tablas de rutas Reactome más frecuentes. Facilita identificar genes diferencialmente expresados globalmente y explorar patrones de rutas funcionales comunes.</li>
                       </ul>"),
                      HTML("<p>
                   Esta sección tiene como objetivo ayudarte a comprender el uso 
                   de los menús, interpretar las tablas y entender la información 
                   provista por cada visualización.
                   </p>")
               ),
               column(6,
                      HTML("<h4><b>Bases de datos utilizadas</b></h4>"),
                      HTML("<p>
                   Los datos empleados en esta aplicación provienen de plataformas 
                   de referencia en bioinformática y genómica:
                   </p>"),
                      tags$a(href = "https://xenabrowser.net/", 
                             HTML("<b>UCSC XENA</b>"), 
                             class = "btn btn-primary", target = "_blank"),
                      br(), br(),
                      tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/", 
                             HTML("<b>Gene Expression Omnibus (GEO)</b>"), 
                             class = "btn btn-primary", target = "_blank"),
                      br(), br(),
                      tags$a(href = "https://www.cbioportal.org/", 
                             HTML("<b>cBioPortal</b>"), 
                             class = "btn btn-primary", target = "_blank")
               )
             )
    )
    
    
    
  )  
)     


# ==========================================================
#  3. SERVER
# ==========================================================

server <- function(input, output, session) {
  
  # -------------------------------
  # Función segura para extraer columnas ignorando mayúsculas/minúsculas
  # -------------------------------
  
  get_col <- function(row, colname) {
    col_match <- grep(paste0("^", colname, "$"), names(row), ignore.case = TRUE, value = TRUE)
    if(length(col_match) == 1) return(row[[col_match]])
    return(NA)
  }
  # -------------------------------
  # Dataset reactivo según tipo de cáncer
  # -------------------------------
  current_df <- reactive({
    req(input$cancer_type)
    df <- datasets[[input$cancer_type]]
    df$Comparison_lower <- tolower(df$Comparison)
    df
  })
  
  # -------------------------------
  # Selector de genes
  # -------------------------------
  
  output$gene_selector <- renderUI({
    df <- current_df()
    selectInput("gene", "Selecciona un gen:", choices = unique(df$Gene))
  })
  
  # -------------------------------
  # Selector de comparación según gen
  # -------------------------------
  
  output$comparison_selector <- renderUI({
    df <- current_df()
    req(input$gene)
    comps <- df %>% filter(Gene == input$gene) %>% pull(Comparison) %>% unique()
    selectInput("comparison", "Selecciona la comparación:", choices = comps)
  })
  
  # -------------------------------
  # Tabla estadística
  # -------------------------------
  
  output$gene_info <- renderTable({
    df <- current_df()
    req(input$gene, input$comparison)
    gene_row <- df %>% filter(Gene == input$gene & Comparison_lower == tolower(input$comparison))
    if(nrow(gene_row) == 0) return(data.frame())
    
    gene_row <- gene_row[1, ]
    groups <- unlist(strsplit(input$comparison, "_vs_"))
    
    if(input$expr_type == "Log-normalizada"){
      stats <- data.frame(
        Group1_Mean   = get_col(gene_row, paste0("expr_", groups[1], "_mean")),
        Group2_Mean   = get_col(gene_row, paste0("expr_", groups[2], "_mean")),
        Group1_Median = get_col(gene_row, paste0("expr_", groups[1], "_median")),
        Group2_Median = get_col(gene_row, paste0("expr_", groups[2], "_median")),
        Group1_SD     = get_col(gene_row, paste0("expr_", groups[1], "_sd")),
        Group2_SD     = get_col(gene_row, paste0("expr_", groups[2], "_sd")),
        Adj.P.Val     = round(get_col(gene_row, "adj.P.Val"), 5),
        logFC         = get_col(gene_row, "logFC")
      )
    } else {
      stats <- data.frame(
        Group1_Mean   = get_col(gene_row, paste0("raw_expr_", groups[1], "_mean")),
        Group2_Mean   = get_col(gene_row, paste0("raw_expr_", groups[2], "_mean")),
        Group1_Median = get_col(gene_row, paste0("raw_expr_", groups[1], "_median")),
        Group2_Median = get_col(gene_row, paste0("raw_expr_", groups[2], "_median")),
        Group1_SD     = get_col(gene_row, paste0("raw_expr_", groups[1], "_sd")),
        Group2_SD     = get_col(gene_row, paste0("raw_expr_", groups[2], "_sd")),
        Adj.P.Val     = round(get_col(gene_row, "adj.P.Val"), 5),
        logFC         = get_col(gene_row, "logFC")
      )
    }
    
    colnames(stats) <- c(
      paste0(groups[1], "_Mean"),
      paste0(groups[2], "_Mean"),
      paste0(groups[1], "_Median"),
      paste0(groups[2], "_Median"),
      paste0(groups[1], "_SD"),
      paste0(groups[2], "_SD"),
      "Adj.P.Val",
      "logFC"
    )
    
    stats
  }, rownames = FALSE)
  
  # -------------------------------
  # Tabla GENENAME y PATH
  # -------------------------------
  
  output$gene_path <- renderTable({
    df <- current_df()
    req(input$gene, input$comparison)
    gene_row <- df %>% filter(Gene == input$gene & Comparison_lower == tolower(input$comparison))
    if(nrow(gene_row) == 0) return(data.frame())
    
    gene_row <- gene_row[1, ]
    safe_value <- function(x) if(is.null(x) || length(x)==0) NA else x
    
    data.frame(
      GENENAME      = safe_value(gene_row$GENENAME),
      GO_terms      = safe_value(gene_row$GO_terms),
      GO_ontology   = safe_value(gene_row$GO_ontology),
      Reactome_PATH = safe_value(gene_row$Reactome_PATH)
    )
  }, rownames = FALSE)
  
  # -------------------------------
  # Gráfica
  # -------------------------------
  
  output$expr_plot <- renderPlot({
    df <- current_df()
    req(input$gene, input$comparison)
    
    gene_row <- df %>% filter(Gene == input$gene & Comparison_lower == tolower(input$comparison))
    if(nrow(gene_row) == 0) return(NULL)
    
    gene_row <- gene_row[1, ]
    groups <- unlist(strsplit(input$comparison, "_vs_"))
    n_samples <- 50
    
    safe_num <- function(x) {
      if(is.list(x) || inherits(x, "Rle")) x <- unlist(x)
      x <- as.numeric(x)
      if(length(x) == 0 || is.na(x)) return(0)
      x[1]
    }
    
    if(input$expr_type == "Log-normalizada"){
      df_plot <- data.frame(
        expression = c(
          rnorm(n_samples,
                mean = safe_num(get_col(gene_row, paste0("expr_", groups[1], "_mean"))),
                sd   = safe_num(get_col(gene_row, paste0("expr_", groups[1], "_sd")))),
          rnorm(n_samples,
                mean = safe_num(get_col(gene_row, paste0("expr_", groups[2], "_mean"))),
                sd   = safe_num(get_col(gene_row, paste0("expr_", groups[2], "_sd"))))
        ),
        group = rep(groups, each = n_samples)
      )
    } else {
      df_plot <- data.frame(
        expression = c(
          rnorm(n_samples,
                mean = safe_num(get_col(gene_row, paste0("raw_expr_", groups[1], "_mean"))),
                sd   = safe_num(get_col(gene_row, paste0("raw_expr_", groups[1], "_sd")))),
          rnorm(n_samples,
                mean = safe_num(get_col(gene_row, paste0("raw_expr_", groups[2], "_mean"))),
                sd   = safe_num(get_col(gene_row, paste0("raw_expr_", groups[2], "_sd"))))
        ),
        group = rep(groups, each = n_samples)
      )
    }
    
    ggplot(df_plot, aes(x = group, y = expression, fill = group)) +
      geom_violin(trim = FALSE, alpha = 0.5) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                   geom = "errorbar", width = 0.1, color = "black") +
      labs(title = paste("Expresión de", input$gene, "(", input$comparison, ")"),
           x = "Grupo", y = "Expresión") +
      theme_minimal()
  })
  
  # -------------------------------
  # Descarga de tabla
  # -------------------------------
  
  output$download_table <- downloadHandler(
    filename = function() {
      paste0("tabla_", input$gene, ".csv")
    },
    content = function(file) {
      gene_row <- current_df() %>% 
        filter(Gene == input$gene & tolower(Comparison) == tolower(input$comparison))
      write.csv(gene_row, file, row.names = FALSE)
    }
  )
  
  # -------------------------------
  # Descarga de gráfica
  # -------------------------------
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("plot_", input$gene, ".png")
    },
    content = function(file) {
      png(file, width = 1200, height = 900)
      
      df <- current_df()
      gene_row <- df %>% filter(Gene == input$gene & Comparison_lower == tolower(input$comparison))
      gene_row <- gene_row[1, ]
      groups <- unlist(strsplit(input$comparison, "_vs_"))
      n_samples <- 50
      safe_num <- function(x) {
        if(is.list(x) || inherits(x, "Rle")) x <- unlist(x)
        x <- as.numeric(x)
        if(length(x) == 0 || is.na(x)) return(0)
        x[1]
      }
      if(input$expr_type == "Log-normalizada"){
        df_plot <- data.frame(
          expression = c(
            rnorm(n_samples,
                  mean = safe_num(get_col(gene_row, paste0("expr_", groups[1], "_mean"))),
                  sd   = safe_num(get_col(gene_row, paste0("expr_", groups[1], "_sd")))),
            rnorm(n_samples,
                  mean = safe_num(get_col(gene_row, paste0("expr_", groups[2], "_mean"))),
                  sd   = safe_num(get_col(gene_row, paste0("expr_", groups[2], "_sd"))))
          ),
          group = rep(groups, each = n_samples)
        )
      } else {
        df_plot <- data.frame(
          expression = c(
            rnorm(n_samples,
                  mean = safe_num(get_col(gene_row, paste0("raw_expr_", groups[1], "_mean"))),
                  sd   = safe_num(get_col(gene_row, paste0("raw_expr_", groups[1], "_sd")))),
            rnorm(n_samples,
                  mean = safe_num(get_col(gene_row, paste0("raw_expr_", groups[2], "_mean"))),
                  sd   = safe_num(get_col(gene_row, paste0("raw_expr_", groups[2], "_sd"))))
          ),
          group = rep(groups, each = n_samples)
        )
      }
      p <- ggplot(df_plot, aes(x = group, y = expression, fill = group)) +
        geom_violin(trim = FALSE, alpha = 0.5) +
        geom_boxplot(width = 0.1, outlier.shape = NA) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                     geom = "errorbar", width = 0.1, color = "black") +
        labs(title = paste("Expresión de", input$gene, "(", input$comparison, ")"),
             x = "Grupo", y = "Expresión") +
        theme_minimal()
      print(p)
      dev.off()
    }
  )
  # -------------------------------
  # Selector de comparación global
  # -------------------------------
  output$global_comparison_selector <- renderUI({
    req(input$global_cancer_type)
    df <- datasets[[input$global_cancer_type]]
    comps <- sort(unique(df$Comparison))
    selectInput("global_comparison",
                "Selecciona la comparación:",
                choices = comps)
  })
  
  # -------------------------------
  # Dataset reactivo para la comparación global
  # -------------------------------
  global_df <- reactive({
    req(input$global_cancer_type, input$global_comparison)
    
    df <- datasets[[input$global_cancer_type]] %>% 
      filter(tolower(Comparison) == tolower(input$global_comparison))
    df
  })
  
  # -------------------------------
  # Gráfica global: Heatmap o EnhancedVolcano
  # -------------------------------
  output$global_plot <- renderPlotly({
    req(global_df())
    df <- global_df()
    
    df$neglogP <- -log10(df$adj.P.Val)
    df$Direction <- ifelse(df$logFC > 0, "Up", "Down")
    
    p <- ggplot(df, aes(x = logFC,
                        y = neglogP,
                        color = Direction,
                        text = paste("Gen:", Gene,
                                     "<br>logFC:", round(logFC,3),
                                     "<br>adj.P.Val:", signif(adj.P.Val,3)))) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(values = c("Down"="blue", "Up"="red")) +
      theme_minimal() +
      labs(title = paste("Volcano Plot —", input$global_cancer_type,"—" ,input$global_comparison),
           x = "logFC", y = "-log10(adj.P.Val)")
    
    ggplotly(p, tooltip = "text")
  })
  
  
  
  # -------------------------------
  # Tabla de paths más frecuentes
  # -------------------------------
  output$path_summary <- renderTable({
    df <- global_df()
    req(nrow(df) > 0)
    
    paths <- unlist(strsplit(df$Reactome_PATH, ";"))
    paths <- paths[!is.na(paths) & paths != "" & paths != "NA"]
    
    path_table <- as.data.frame(table(paths))
    path_table <- path_table %>% arrange(desc(Freq)) %>% head(input$top_paths)
    path_table$Percentage <- round(path_table$Freq / sum(path_table$Freq) * 100, 2)
    
    colnames(path_table) <- c("Reactome_PATH", "Count", "Percentage")
    path_table
  })
  
  # -------------------------------
  # Descargas
  # -------------------------------
  
  output$download_global_table <- downloadHandler(
    filename = function() {
      paste0("tabla_global_", input$global_cancer_type, "_", gsub(" ", "_", input$global_comparison), ".csv")
    },
    content = function(file) {
      df <- global_df()
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$download_global_plot <- downloadHandler(
    filename = function() {
      paste0("plot_global_", input$global_cancer_type, "_", gsub(" ", "_", input$global_comparison), ".png")
    },
    content = function(file) {
      df <- global_df()
      
      # Crear el objeto del plot
      p <- EnhancedVolcano(df,
                           lab = df$Gene,
                           x = 'logFC',
                           y = 'adj.P.Val',
                           pCutoff = 0.05,
                           FCcutoff = 1.0,
                           pointSize = 3.0,
                           labSize = 3.0)
      
      # Guardar con ggsave
      ggsave(filename = file, plot = p, width = 16, height = 12, units = "in", dpi = 300)
    }
  )
  
}


# ==========================================================
#  4. LANZAR APP
# ==========================================================

shinyApp(ui, server)
