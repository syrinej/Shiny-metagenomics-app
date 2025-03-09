# Load necessary libraries
library(shiny)
library(phyloseq)
library(ggplot2)

# User Interface
ui <- fluidPage(
  titlePanel("Metagenomic Data Visualization with Phyloseq"),
  sidebarLayout(
    sidebarPanel(
      fileInput("otuFile", "Upload OTU table", accept = c(".csv")),
      fileInput("sampleFile", "Upload Sample Data table", accept = c(".csv")),
      fileInput("taxFile", "Upload Taxonomy table", accept = c(".csv")),
      selectInput("analysisType", "Choose analysis type:", 
                  choices = c("Bar Plot", "Heatmap", "Alpha Diversity"))
    ),
    mainPanel(
      plotOutput("mainPlot")
    )
  )
)

# Server logic
server <- function(input, output) {
  
  observeEvent(list(input$otuFile, input$sampleFile, input$taxFile), {
    
    # Ensure all files are uploaded
    req(input$otuFile, input$sampleFile, input$taxFile)
    
    otu <- read.csv(input$otuFile$datapath, row.names = 1)
    sample_data <- read.csv(input$sampleFile$datapath, row.names = 1)
    tax <- read.csv(input$taxFile$datapath, row.names = 1)
    
    otu_mat <- as.matrix(otu)
    tax_mat <- as.matrix(tax)
    
    otu_phylo <- otu_table(otu_mat, taxa_are_rows = TRUE)
    tax_phylo <- tax_table(tax_mat)
    sample_phylo <- sample_data(sample_data)
    
    physeq <- phyloseq(otu_phylo, tax_phylo, sample_phylo)
    
    output$mainPlot <- renderPlot({
      if (input$analysisType == "Bar Plot") {
        p <- plot_bar(physeq, x="SampleType") + geom_bar(aes(fill = Rank1), position="stack")
        print(p)
      } else if (input$analysisType == "Heatmap") {
        p <- plot_heatmap(physeq, sample.label="SampleType", taxa.label="Rank1")
        print(p)
      } else if (input$analysisType == "Alpha Diversity") {
        p <- plot_richness(physeq, x="SampleType", measures=c("Shannon", "Simpson"))
        print(p)
      }
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)





