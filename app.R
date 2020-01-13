library(shiny)
library(dplyr)
library(readr)
library(DT)

clinical_ann_metadata <- read_delim("clinical_ann_metadata.tsv", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)

chemical <- unique(clinical_ann_metadata$`Related Chemicals`)
chemical <- sort(chemical[-grep(",", chemical)])

vcf_M <- readRDS("vcf_1")
vcf_L <- readRDS("vcf_2")
vcf_J <- readRDS("vcf_3")

url <- "https://www.pharmgkb.org/clinicalAnnotation/"

lucky <- sample(1:length(chemical), length(chemical)/10)
lucky_drugs <- chemical[lucky]


ui <- fluidPage(
  
  imageOutput("image", height = 50),
  sidebarLayout(
    sidebarPanel(
      fileInput("vcf_file", "Upload VCF File"),
      fileInput("bam_file", "Upload BAM File"),
      actionButton("button1", "Patient Mantas"),
      actionButton("button2", "Patient Lukas"),
      actionButton("button3", "Patient Juozas"),
      selectInput(inputId = "chemical", "Select a Chemical", choices = chemical),
      htmlOutput("text"),
      uiOutput("vcf_example"),
      uiOutput("bam_example"),
      uiOutput("evidence")
      ),
    mainPanel(
      htmlOutput("text2"),
      DT::dataTableOutput("phenotype"),
      uiOutput("drug"),
      uiOutput("tag")
      
    )
  )
)

server <- function(input, output, session) {
  
  
  clinical_ann_metadata <- read_delim("clinical_ann_metadata.tsv", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
  clinical_ann <- read_delim("clinical_ann.tsv", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
  
  options(shiny.maxRequestSize = 200000*1024^2)
  
  
  output$text <- renderUI({
   HTML('<b>More Information:</b>')
  })
  link <-  a("About VCF file", href="http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_VCF.htm")
  output$vcf_example <- renderUI({
    tagList(link)
  })
  
  link2 <- a("About BAM file", href = "https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm")
  
  output$bam_example <- renderUI({
    tagList(link2)
  })
  
  link3 <- a("Levels of Evidence", href = "https://www.pharmgkb.org/page/clinAnnLevels")
  output$evidence <- renderUI({
    tagList(link3)
  })
  
  
  output$text2 <- renderUI({
    HTML(paste('This is the prototype of the decision support system, designed to help with for better drug prescription. All information provided here is for testing purposes and can’t be used in real decision making process or diagnostic procedures. Service is not released officially and it is still in development process controlled by omica.lt',

         ' ', '<b>How to use the software:</b>',
         
         '<b>1.</b>     Upload genomic data stored in vcf format or choose a sample patient',
         '<b>2.</b>    Simulate drug prescription process by choosing a drug from the list', sep = "<br/>"))
  })
  output$image <- renderImage({
    list(
      src = "linkgene_logo.png",
      contentType = "image/png",
      height = 50)},
    deleteFile = FALSE)
  
  
  output$phenotype <- DT::renderDataTable({
    req(input$vcf_file)
    vcf_dt <- data.table::fread(input$vcf_file$datapath, skip = "CHROM")
    id <- vcf_dt$ID
    allele <- vcf_dt$ALT
    clin_ann <- filter(clinical_ann_metadata, clinical_ann_metadata$`Related Chemicals` == input$chemical)
    ann_out <- NULL
    linkID <- NULL
    for ( i in 1:length(id)){
      var <- filter(clin_ann, clin_ann$Location %in% id[i])
      varID <- unlist(gsubfn::strapply(var$`Genotype-Phenotype IDs`, '([0-9]+)'))
      var <- select(var, `Clinical Annotation ID`, `Level of Evidence`,`Clinical Annotation Types`)
      var_ann <- filter(clinical_ann, clinical_ann$`Genotype-Phenotype ID` %in% varID)
      var_ann <- filter(var_ann, var_ann$Genotype == allele[i]) %>% select(c(2,3))
      var_ann <- cbind(var_ann$Genotype, var_ann$`Clinical Phenotype`, var$`Level of Evidence`, var$`Clinical Annotation Types`)
      linkID <- c(linkID, var$`Clinical Annotation ID`)
      ann_out <- rbind(ann_out, var_ann)
    }
    pa <- gsubfn::strapply(input$chemical, '(PA[0-9]+)')
    refsd <- a(paste0("More information about ", input$chemical), href = paste("https://www.pharmgkb.org/chemical/", pa, "/overview"))
    output$drug <- renderUI({
      tagList(refsd)
    })
    datatable(ann_out, colnames = c("Genotype", "Clinical Phenotype", "Level of Evidence", "Clinical Annotation Types"))
  })

  
  observeEvent(input$button1, {
    output$phenotype <-DT::renderDataTable({
      vcf_dt <- vcf_M
      id <- vcf_dt$ID
      allele <- vcf_dt$ALT
      clin_ann <- filter(clinical_ann_metadata, clinical_ann_metadata$`Related Chemicals` == input$chemical)
      ann_out <- NULL
      for ( i in 1:length(id)){
        var <- filter(clin_ann, clin_ann$Location %in% id[i])
        varID <- unlist(gsubfn::strapply(var$`Genotype-Phenotype IDs`, '([0-9]+)'))
        var <- select(var, `Clinical Annotation ID`, `Level of Evidence`,`Clinical Annotation Types`)
        var_ann <- filter(clinical_ann, clinical_ann$`Genotype-Phenotype ID` %in% varID)
        var_ann <- filter(var_ann, var_ann$Genotype == allele[i]) %>% select(c(2,3))
        var_ann <- cbind(var_ann$Genotype, var_ann$`Clinical Phenotype`, var$`Level of Evidence`, var$`Clinical Annotation Types`)
        ann_out <- rbind(ann_out, var_ann)
      }

      pa <- gsubfn::strapply(input$chemical, '(PA[0-9]+)')
      refsd <- a(paste0("More information about ", input$chemical), href = paste("https://www.pharmgkb.org/chemical/", pa, "/overview"))
      output$drug <- renderUI({
        tagList(refsd)
      })
      datatable(ann_out, colnames = c("Genotype", "Clinical Phenotype", "Level of Evidence", "Clinical Annotation Types"))
    })
  })
  
  observeEvent(input$button2, {
    output$phenotype <-DT::renderDataTable({
      vcf_dt <- vcf_J
      id <- vcf_dt$ID
      allele <- vcf_dt$ALT
      clin_ann <- filter(clinical_ann_metadata, clinical_ann_metadata$`Related Chemicals` == input$chemical)
      ann_out <- NULL
      for ( i in 1:length(id)){
        var <- filter(clin_ann, clin_ann$Location %in% id[i])
        varID <- unlist(gsubfn::strapply(var$`Genotype-Phenotype IDs`, '([0-9]+)'))
        var <- select(var, `Clinical Annotation ID`, `Level of Evidence`,`Clinical Annotation Types`)
        var_ann <- filter(clinical_ann, clinical_ann$`Genotype-Phenotype ID` %in% varID)
        var_ann <- filter(var_ann, var_ann$Genotype == allele[i]) %>% select(c(2,3))
        var_ann <- cbind(var_ann$Genotype, var_ann$`Clinical Phenotype`, var$`Level of Evidence`, var$`Clinical Annotation Types`)
        ann_out <- rbind(ann_out, var_ann)
      }

      pa <- gsubfn::strapply(input$chemical, '(PA[0-9]+)')
      refsd <- a(paste0("More information about ", input$chemical), href = paste("https://www.pharmgkb.org/chemical/", pa, "/overview"))
      output$drug <- renderUI({
        tagList(refsd)
      })
      datatable(ann_out, colnames = c("Genotype", "Clinical Phenotype", "Level of Evidence", "Clinical Annotation Types"))
    })
  })
  
  observeEvent(input$button3, {
    output$phenotype <-DT::renderDataTable({
      vcf_dt <- vcf_L
      id <- vcf_dt$ID
      allele <- vcf_dt$ALT
      clin_ann <- filter(clinical_ann_metadata, clinical_ann_metadata$`Related Chemicals` == input$chemical)
      ann_out <- NULL
      for ( i in 1:length(id)){
        var <- filter(clin_ann, clin_ann$Location %in% id[i])
        varID <- unlist(gsubfn::strapply(var$`Genotype-Phenotype IDs`, '([0-9]+)'))
        var <- select(var, `Clinical Annotation ID`, `Level of Evidence`,`Clinical Annotation Types`, `Location`)
        var_ann <- filter(clinical_ann, clinical_ann$`Genotype-Phenotype ID` %in% varID)
        var_ann <- filter(var_ann, var_ann$Genotype == allele[i]) %>% select(c(2,3))
        var_ann <- cbind(var_ann$Genotype, var_ann$`Clinical Phenotype`, var$`Level of Evidence`, var$`Clinical Annotation Types`)
        ann_out <- rbind(ann_out, var_ann)
      }
     
      pa <- gsubfn::strapply(input$chemical, '(PA[0-9]+)')
      refsd <- a(paste0("More information about ", input$chemical), href = paste("https://www.pharmgkb.org/chemical/", pa, "/overview"))
      output$drug <- renderUI({
        tagList(refsd)
      })
      datatable(ann_out, colnames = c("Genotype", "Clinical Phenotype", "Level of Evidence", "Clinical Annotation Types"))
      })
    })
}


shinyApp(ui, server)







