library(tidyverse)
library(ggplot2)
library(cowplot)
library(Rtsne)
library(RColorBrewer)
library(vegan)
load("meta2.Rdata")
box_col=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#b3e2cd")
theme_set(theme_cowplot())

ui <- fluidPage(
  titlePanel("Dimensionality reduction in Microbiome Studies"),
  sidebarLayout(
    sidebarPanel(
      h3("Parameters"),
      selectInput("tax", "Choose taxonomic rank (for Taxonomic composition plot only):",
                  list("Phylum","Family", "Genus")),
      selectInput("type", "Choose data type (for Ordination plot only)",
                  list("OTU_abundances","OTU_presence_absence")),
      selectInput("method", "Choose a dimensionality reduction method (for Ordination plot only):",
                  list("PCoA","nMDS","t-SNE")),
      hr(),
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Introduction",
                           p("Two central goals of microbiome research are to a) characterize the composition and functional capabilities of the microbiome at one or more body sites and to b) determine whether specific factors are associated with variation in the microbiome among individuals. In order to study these questions, microbiome studies usually have to use an array of statistical techniques to make sense of the data. A sequencing run will generate a table of samples and counts, whether that be counts of OTUs (for 16S rRNA gene sequencing) or counts of genes (for shotgun metagenomic sequencing). But how do you go from a table of counts, to actual biologically meaningful patterns that can be interpreted? "),
                           p("Usually, some type of dimensionality reduction technique is applied. These techniques will take our table of counts of OTUs (potentially thousands of OTUs, which means potentially thousands of dimensions), and transform it into data with fewer dimensions while still capturing as much of the variance present in the original dataset. In this R Shiny App you will explore 3 dimensionality reduction methods used in microbiome and genomic studies: Principal Coordinates Analyses (PCoA), non-metric multidimensional scaling (nMDS), and t-distributed stochastic neighbor embedding (t-SNE). Read more about each method here: https://tinyurl.com/wbnrzb4"),
                           p("For this R Shiny App, you will be visualizing the microbiome of wild spotted hyenas at multiple bodysites. Visualize the composition of the bacterial communities at the phylum, family or genera level (note that only the top 10 most abundant taxa will show, otherwise the plot gets too crowded). In this plot, each bar represents a unique hyena and each color represents a bacterial taxa. Then, visualize the clustering of communities after applying dimensionality reduction. You will also get to decide if your input data will be OTU proportions or OTU presence/absence (meaning the data matrix will only have 1s (OTU present in sample) and 0’s (OTU absent in that sample))."),
                           p("Here are a few questions to guide you along: "),
                           p("For taxonomic composition plot: What bacterial taxa are abundant or unique to each body site? Which body sites are similar? What is the variability like within and between bodysites?"),
                           p("For ordination plot: Do bacterial communities cluster by bodysite? Which bodysites overlap? How do the different reduction methods vary? How do they agree? Does this make sense mathematically given what you know about these methods?"),
                           p("Why does using OTU presence absence data lead to less-defined clusters?"),
                           p("Does the ordination plot make sense given what you observed in the taxonomic composition plots?")),
                  tabPanel("Plots",
                           h3("Taxonomic Composition"),
                           plotOutput("plot_tax", height = "300px"),
                           h3("Ordination Plot"),
                           plotOutput("plot_ord", height = "500px"))
      )
    )
  )
)

##YANN'S CODE FROM LAST APP TO HELP ME BUILD CODE FOR NEW APP
p = ggplot()

server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  
  output$plot_tax <- renderPlot({
    if (input$tax=="Phylum") {
      load("Phylum.Rdata")
      tax_meta = tax_meta
      new_col=brewer.pal(n=6, name = "Paired")}
    if (input$tax=="Family") {
      load("Family.Rdata")
      tax_meta = tax_meta
      new_col=c("#a6cee3", "blue", "#b2df8a","#33a02c","grey42", "#e31a1c", "#fdbf6f", 
                 "#ff7f00", "#cab2d6", "magenta", "#ffff99","#b15928","#e7298a")}
    if (input$tax=="Genus") {
      load("Genus.Rdata")
      tax_meta = tax_meta
      new_col=c(brewer.pal(n=12, name = "Paired"), "grey")}
    p=ggplot(tax_meta, aes(x = Group, y =abun , fill = taxon)) + 
      geom_bar(stat = "identity") + facet_grid(~bodysite, scales="free_x")+  
      scale_fill_manual(values=new_col)+
      labs(y="Relative Abundance (proportion)",x="",fill="")+
      theme(text = element_text(size=11),legend.position="right")+
      theme_bw()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
      theme(plot.title = element_text(size = 12, face = "bold"))+theme(legend.text=element_text(size=13))
    print(p)
  })
  
  output$plot_ord <- renderPlot({
    if (input$type=="OTU_abundances") {
      load("weighted.Rdata")
      dist.mat = dist.mat}
    if(input$type=="OTU_presence_absence"){
      load("unweighted.Rdata")
      dist.mat = dist.mat}
    if(input$method=="PCoA"){
      pcoa_dec=cmdscale(dist.mat, eig=TRUE); pcoa=as.data.frame(pcoa_dec$points)
      colnames(pcoa)=c("Axis1","Axis2"); pcoa=tibble::rownames_to_column(as.data.frame(pcoa), "Group")
      pcoa_met=merge(pcoa,meta_data,by="Group") 
      pcoa_per=(pcoa_dec$eig/sum(pcoa_dec$eig))*100; sum(pcoa_per) 
      
      p=ggplot(pcoa_met, aes(Axis1, Axis2))+geom_point(aes(fill=bodysite), size = 5, shape=21)+
        scale_fill_manual(values=box_col) +
        labs(colour="")+ylab(paste("PC2 (",round(pcoa_per[2],digits=2),"%)"))+xlab(paste("PC1 (",round(pcoa_per[1],digits=2),"%)"))+
        theme_bw()+ theme(legend.title=element_blank())+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(colour = "black", size=2))+
        theme(legend.position="right",legend.text=element_text(size=15))}
    if(input$method=="nMDS"){
      bod.nmds <- metaMDS(dist.mat, k=4)
      meta_ds=meta_data; rownames(meta_ds)=meta_ds$Group; meta_ds$Group=NULL
      ax1 = bod.nmds$points[,1]
      ax2 = bod.nmds$points[,2]
      nmds_met=data.frame(ax1 = ax1, ax2 = ax2, bodysite=meta_ds$bodysite)
      p=ggplot(nmds_met, aes(ax1, ax2, fill=bodysite))+geom_point(size =5, shape=21)+stat_ellipse()+
        scale_fill_manual(values=box_col) +
        labs(colour="")+ylab("NMDS2")+xlab("NMDS1")+labs(title=paste("NMDS plot, stress=",round(bod.nmds$stress, digits=2)))+
        theme_bw()+theme(legend.text=element_text(size=15),legend.title=element_blank())}
    if(input$method=="t-SNE"){
      set.seed(20)
      tsne1 = Rtsne(dist.mat, check_duplicates=FALSE, pca=TRUE, perplexity=20, theta=0.5, dims=3)
      tsne_plot <- data.frame(x = tsne1$Y[,1], y = tsne1$Y[,2], col = meta_data$bodysite)
      p=ggplot(tsne_plot, aes(x, y, fill=meta_data$bodysite))+geom_point(size =5, shape=21)+
        scale_fill_manual(values=box_col) +
        labs(colour="")+ylab("t-SNE2")+xlab("t-SNE1")+
        theme_bw()+theme(legend.text=element_text(size=15),legend.title=element_blank())}
    print(p)
  })
}

shinyApp(ui = ui, server = server)
