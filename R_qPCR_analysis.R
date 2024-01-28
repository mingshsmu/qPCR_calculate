#' @title R_qPCR_analysis
#' @author GM.W
#' @date 2023-05-16
require(dplyr,quietly = T,warn.conflicts = F)
require(optparse,quietly = T,warn.conflicts = F)
require(ggplot2,quietly = T,warn.conflicts = F)

wd <- getwd()
option_list <- list(
  make_option(opt_str = c("-i","--input"),type = "character",default = NULL,
              action = "store", help = "The name of input qPCR file."),
  make_option(opt_str = c("-o","--output"),type = "character",default = NULL,
              action = "store", help = "The name of output file."),
  make_option(opt_str = c("-c","--control"),type = "character",default = NULL,
              action = "store", help = "The name of control in experiment."),
  make_option(opt_str = c("-r","--reference"),type = "character",default = NULL,
              action = "store", help = "The reference gene name."),
  make_option(opt_str = c("-f","--fold"),type = "double",default = NULL,
              action = "store", help = "The foldchange of sd to deal with outliar values."),
  make_option(opt_str = c("-p","--palette"),type = "character",default = NULL,
              action = "store", help = "Either 'Set1', 'Set2','Set3','Dark2','Accent','Pastel1','Pastel2','Paired'"),
  make_option(opt_str = c("-e","--relativeExpression"),type = "logical",default = TRUE,
              action = "store", help = "If calculate the relative expression?")            
)

args <- parse_args(OptionParser(option_list = option_list))

input <- args$input
output <- args$output
con <- args$control
ref <- args$reference
fold <- args$fold
palette <- args$palette
rq <- args$relativeExpression

output <- stringr::str_split(output,pattern = "[.]",n = 2,simplify = T)[1]

ct_data <- readxl::read_excel(path = input,
                           sheet = "Results",skip = 24)
melt_data <- readxl::read_excel(path = input,
                        sheet = "Melt Curve Raw",skip = 24)
melt_data <- dplyr::left_join(melt_data,ct_data[,c("Well Position","Sample")],by=c("Well Position"))
amp_data <- readxl::read_excel(path = input,
                               sheet = "Amplification Data",skip = 24)

data <- ct_data[,c("Sample","Target","Cq")]
colnames(data) <- c("samples","gene","ct")


leve <- unique(data$samples)
leve <- leve[!leve %in% con]
leve <- c(con, leve)
data$samples <- factor(data$samples, levels = leve)

primer_level <- unique(data$gene)
primer_level <- primer_level[!primer_level %in% ref]
primer_level <- c(ref, primer_level)
data$gene <- factor(data$gene, levels = primer_level)

color_in <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", "#9467BDFF",
              "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
              "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
              "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
# scales::show_col(color_in)
# color_in <- sample(color_in, length(color_in), replace=FALSE)
qpcr_analysis <- function(data,ref,con,fold=NULL){
  if(is.null(fold)){
    fold <- 3
  }
  data <- data %>% group_by(gene,samples) %>% 
    mutate(ct_mean = mean(ct)) %>% 
    mutate(ct_sd = sd(ct)) %>% 
    mutate(ct_iqr = IQR(ct)) %>% 
    mutate(keep = ifelse(abs(ct-ct_mean) > fold*ct_sd, FALSE,TRUE)) %>% 
    filter(keep==TRUE)
  
  
  ref_mean <- data %>% filter(gene==ref) %>% 
    group_by(samples) %>%
    summarise(ref_mean=mean(ct))
  ref_mean <- setNames(ref_mean$ref_mean, ref_mean$samples)
  
  res <- data.frame()
  for (i in names(ref_mean)) {
    res_tmp <- data %>% filter(samples==i) %>% 
      mutate(dct = ct-ref_mean[i])
    res <- rbind(res,res_tmp)
  }
  
  con_dmean <- res %>% filter(samples==con) %>% 
    group_by(gene) %>% 
    summarise(con_dmean = mean(dct))
  con_dmean <- setNames(con_dmean$con_dmean, con_dmean$gene)
  
  result <- data.frame()
  for (j in names(con_dmean)) {
    res_tmp <- res %>% filter(gene==j) %>% 
      mutate(ddct = dct - con_dmean[j])
    result <- rbind(result, res_tmp)
  }
  
  # final
  result$exp <- 2^(-result$ddct)
  result <- result %>% group_by(gene,samples) %>% 
    mutate(exp_mean = mean(exp)) %>% 
    mutate(exp_sd = sd(exp))
  return(result)
}
qpcr_summary <- function(result){
  summ <- result %>% distinct(samples, gene,.keep_all = T)
  return(summ)
}
qpcr_plot <- function(result, summ, palette=NULL){
  if(is.null(palette)){
    palette <- color_in
  }else{
    palette <- stringr::str_to_title(palette)
  }
  p1 <- ggplot()+
    geom_bar(data=summ,mapping=aes(x=samples,y=exp_mean,fill=samples,color=samples),alpha=0.8,stat = "identity")+
    geom_errorbar(data=summ,mapping=aes(x=samples,ymin=exp_mean-exp_sd,ymax=exp_mean+exp_sd,
                                        color=samples),width=0.3)+
    geom_jitter(data=result,mapping = aes(x=samples,y=exp),width = 0.3)+
    facet_wrap(.~gene,scales = "free_y",nrow = 1)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45,hjust = 1),
          strip.background = element_rect(fill = "#eaeae0"))+
    labs(x=NULL,y="Relative expression")
  if(length(palette)==1){
    p1 <- p1+scale_fill_brewer(palette = palette,direction = -1)+
      scale_color_brewer(palette = palette,direction = -1)
  }else{
    p1 <- p1+scale_fill_manual(values = palette)+
      scale_color_manual(values = palette)
  }
    
  return(p1)
}
melt_plot <- function(data){
  sample_level <- unique(data$Sample)
  target_level <- unique(data$Target)
  data$Sample <- factor(data$Sample, levels = sample_level)
  data$Target <- factor(data$Target, levels = target_level)
  p <- ggplot(data = data, aes(x=Temperature,y=Derivative,color=Target,group=Well))+
    geom_line()+
    facet_grid(Sample~Target)+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "#eaeae0"))+
    scale_color_manual(values = color_in)
  return(p)
}

amp_plot <- function(data){
  sample_level <- unique(data$Sample)
  target_level <- unique(data$Target)
  data$Sample <- factor(data$Sample, levels = sample_level)
  data$Target <- factor(data$Target, levels = target_level)
  p <- ggplot(data = data, aes(x=`Cycle Number`,y=dRn,color=Target,group=Well))+
    geom_line()+
    facet_grid(Sample~Target)+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "#eaeae0"))+
    scale_color_manual(values = color_in)
  return(p)
}


if(rq){
  result <- qpcr_analysis(data = data,ref = ref,con = con,fold = fold)
  summ <- qpcr_summary(result = result)
  summ <- summ[,c(1,2,11,12)]
  message("qPCR analysis completed!")
  write.csv(result,file = paste0(output,".csv"),row.names = FALSE)
  message("Calculating and plotting...")
  p1 <- qpcr_plot(result = result, summ = summ, palette = palette)
  
  ggsave(filename = paste0("Rplot_expression_",output,".png"),device = "png",
         width = 2.5 * length(unique(data$gene)),height = 5,dpi = 300,plot = p1)
  ggsave(filename = paste0("Rplot_expression_",output,".pdf"),device = "pdf",
         width = 2.5 * length(unique(data$gene)),height = 5,plot = p1)
  
}

# MeltCurve & AmplicationCurve
height_ratio <- ifelse(length(unique(melt_data$Sample))==1,2,1.6)

message("MeltCurve plotting...")
p2 <- melt_plot(data = melt_data)
ggsave(filename = paste0("Rplot_MeltCurve_",output,".png"),device = "png",
       width = 2 * length(unique(melt_data$Target)),height = height_ratio*length(unique(melt_data$Sample)),
       dpi = 300,plot = p2)
ggsave(filename = paste0("Rplot_MeltCurve_",output,".pdf"),device = "pdf",
       width = 2 * length(unique(melt_data$Target)),height = height_ratio*length(unique(melt_data$Sample)),
       plot = p2)

message("AmplicationCurve plotting...")
p3 <- amp_plot(data = amp_data)
ggsave(filename = paste0("Rplot_AmplicationCurve_",output,".png"),device = "png",
       width = 2 * length(unique(amp_data$Target)),height = height_ratio*length(unique(amp_data$Sample)),
       dpi = 300,plot = p3)
ggsave(filename = paste0("Rplot_AmplicationCurve_",output,".pdf"),device = "pdf",
       width = 2 * length(unique(amp_data$Target)),height = height_ratio*length(unique(amp_data$Sample)),
       plot = p3)
message("Done! :)")  


