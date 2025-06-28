

load("./output/models_all.RData")

pacman::p_load("dplyr","tibble","flextable","ggplot2","gamlss","sf","gridExtra",
               "kableExtra","ggspatial")

TEMA<-theme(
  axis.title.x = element_text(color=1, size=8),# face="bold"),
  axis.title.y = element_text(color=1, size=8),# face="bold"),
  axis.text.x = element_text(color=1, size=8),#,margin = margin(l = 500)),
  axis.text.y = element_text(color=1, size=8),
  plot.title = element_text(size=18),
  legend.key.height = unit(1, 'cm'), #change legend key height
  legend.key.width = unit(0.2, 'cm'),
  legend.title = element_text(size=18)
) + theme_bw() +
  theme(legend.position = "bottom")


datos_st <-st_read('../datos/datos_301223.shp', quiet=TRUE) 
sabana_st <-st_read('../datos/Area_sabana_2013.shp', quiet=TRUE) 

sabana_st <- st_transform(sabana_st, crs = st_crs(datos_st))


points_within <- st_within(datos_st, sabana_st)
is_within <- lengths(points_within) > 0

datos_st <- datos_st %>% mutate(sabana=is_within)

datos_st <- datos_st %>% mutate(INC_d = case_when(INC>0~1,
                                                  .default = INC),
                                INC_d = as.factor(INC_d),
                                id = 1:n(),
                                logFA= log(FA+0.1),
                                logLSF= log(LSF+0.1),
                                logTRI= log(TRI+0.1),
                                logINC = log(INC+1))

# Lo ideal es quitar el área que queda fuera de la región?
datos_st <- datos_st %>% filter(SLOPE != 0)


datos_df <- datos_st %>% as.data.frame() %>% 
  dplyr::select(id,INC,DAH,logFA,LST,logLSF,SLOPE,logTRI,TWI,WE,WExpo,sabana)

datos_df <- datos_df %>%
  mutate_at(c(3:11), ~scale(.))

# Figure
INC_hist<- 
  ggplot(datos_df, aes(x = INC)) +
  geom_histogram(fill = "red", color = "black", bins = 30) +
  labs(title = "", x = "Wildfire", y = "Frequency") +
  theme(axis.text = element_text(size = 15),
        axis.title.x = element_text(color=1, size=20),# face="bold"),
        axis.title.y = element_text(color=1, size=20)) + 
  theme_bw()
INC_hist

ggsave(INC_hist,
       filename = paste0("./output_all/INC_hist.jpg"),
       height = 7, width = 7)

# Split data into training and testing data (0.5 and 0.5).
set.seed(4646)
rand2<-sample(2, nrow(datos_df), replace=TRUE, prob=c(0.7,0.3))

datos_training <- datos_df[rand2==1,] # training data
datos_validation <- datos_df[rand2==2,]

load("./output/models_all.RData")

model.list <- list(modPO=modPO,         #step 1
                   modNBII=modNBII,
                   modZIP=modZIP,
                   modZINBI=modZINBI,
                   modPO_red=modPO_red,     #step 2
                   modNBII_red=modNBII_red,
                   modZIP_red=modZIP_red,
                   modZINBI_red=modZINBI_red,
                   modZIPfinal=modZIPfinal,   #step 3
                   modZINBIfinal=modZINBIfinal)


for(i in 1:length(model.list)){
  
  # res <- residuals(model.list[[i]], what = "z-scores", type = c("simple"))
  # yfit <- fitted(model.list[[i]])
  # datos_training1 <- datos_training %>% mutate(res=res,
  #                                              yfit=yfit,
  #                                              index=1:nrow(datos_training))
  # 
  # fig_res1 <- ggplot(datos_training1) + geom_point(
  #   aes(x=yfit,y=res)) + TEMA +labs(title="(NQ) residual v.s. fitted values",
  #                                   x ="fitted values", y = "normalized quantile residuals")
  # 
  # fig_res2 <- ggplot(datos_training1) + geom_point(
  #   aes(x=index,y=res)) + TEMA +labs(title="index plot",
  #                                    x ="index", y = "normalized quantile residuals")
  # 
  # fig_res3 <- ggplot(data = datos_training1, aes(x = res, after_stat(density)))+
  #   geom_histogram(alpha = 0.50)+
  #   geom_density(alpha = 0.50)+ TEMA +labs(title="histogram and density",
  #                                          x ="normalized quantile residuals", y = "density")
  # 
  # fig_res4 <- ggplot(datos_training1, aes(sample = res)) + stat_qq() + stat_qq_line() + TEMA +
  #   labs(title="Q-Q plot",
  #        x ="sample quantiles", y = "theoretical quantiles")
  # 
  # 
  # # Figure
  # fig1<-gridExtra::grid.arrange(fig_res1,fig_res2,fig_res3,fig_res4,ncol=2)
  # 
  # ggsave(fig1,
  #        filename = paste0("./output_all/SuppM/",names(model.list)[i],"_diagnostics.jpg"),
  #        height = 10, width = 10)
  
  
  # Figure
  jpeg(paste0("output_all/SuppM/",names(model.list)[i],"_wp.jpg"))
  wp(model.list[[i]],ylim.all = 1.5, xlim.all = 5)
  #wp(model.list[[i]])
  dev.off()
  
}


