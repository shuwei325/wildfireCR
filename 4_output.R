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
theme_bw() +
  theme(aspect.ratio = 0.5)  # Height:Width = 0.5
INC_hist

INC_freq <- ggplot(datos_df, aes(x = INC)) +
  geom_bar(fill = "red") +
  labs(title = "", x = "Wildfire", y = "Frequency") +
  theme(axis.text = element_text(size = 15),
        axis.title.x = element_text(color=1, size=20),# face="bold"),
        axis.title.y = element_text(color=1, size=20)) + 
  theme_bw()+
  theme(aspect.ratio = 0.3)  # Height:Width = 0.5
INC_freq

ggsave(INC_freq,
       filename = paste0("./output_all/INC_hist.jpg"),
       height = 0.3*14, width = 7)

# Split data into training and testing data (0.5 and 0.5).
set.seed(4646)
rand2<-sample(2, nrow(datos_df), replace=TRUE, prob=c(0.7,0.3))

datos_training <- datos_df[rand2==1,] # training data
datos_validation <- datos_df[rand2==2,]



## Statistical model

load("./output/models_final.RData")


model.list <- list(modPO,         #step 1
                   modNBII,
                   modZIP,
                   modZINBI,
                   modPO_red,     #step 2
                   modNBII_red,
                   modZIP_red,
                   modZINBI_red,
                   modZIPfinal,   #step 3
                   modZINBIfinal)

GAIC.l <- cbind(sapply(model.list,GAIC,k=2))
BIC.l <- cbind(sapply(model.list,GAIC,k=log(nrow(datos_training))))
VGD.l <- cbind(sapply(model.list,VGD))
R2adj <- cbind(sapply(model.list,Rsq,type="both"))

CoxSnell <- round(as.numeric(R2adj[1,]),4)
CraggUhler <- round(as.numeric(R2adj[2,]),4)


IC_table <- cbind(GAIC.l,BIC.l,VGD.l,CoxSnell, CraggUhler)

colnames(IC_table) <- c("AIC","BIC","VGD",'R2-CoxSnell', 'R2-Cragg Uhler')
rownames(IC_table) <- c('PO',         #step 1
                        'NBII',
                        'ZIP',
                        'ZINBI',
                        'PO_red',     #step 2
                        'NBII_red',
                        'ZIP_red',
                        'ZINBI_red',
                        'ZIP_red_final',   #step 3
                        'ZINBI_red_final')

IC_table %>% as.data.frame() %>% 
  rownames_to_column("Model") %>% flextable() %>% 
  theme_booktabs() 


kbl(IC_table, format="latex")





load("./output/models_all.RData")
res <- residuals(modZIPfinal, what = "z-scores", type = c("simple"))
yfit <- fitted(modZIPfinal)
datos_training1 <- datos_training %>% mutate(res=res,
                                             yfit=yfit,
                                             index=1:nrow(datos_training))

fig_res1 <- ggplot(datos_training1) + geom_point(
  aes(x=yfit,y=res)) + TEMA +labs(title="(a) (NQ) residual v.s. fitted values",
                                  x ="fitted values", y = "normalized quantile residuals")

fig_res2 <- ggplot(datos_training1) + geom_point(
  aes(x=index,y=res)) + TEMA +labs(title="(b) index plot",
                                   x ="index", y = "normalized quantile residuals")

fig_res3 <- ggplot(data = datos_training1, aes(x = res, after_stat(density)))+
  geom_histogram(alpha = 0.50)+
  geom_density(alpha = 0.50)+ TEMA +labs(title="(c) histogram and density",
                                         x ="normalized quantile residuals", y = "density")

fig_res4 <- ggplot(datos_training1, aes(sample = res)) + stat_qq() + stat_qq_line() + TEMA +
  labs(title="(d) Q-Q plot",
       x ="sample quantiles", y = "theoretical quantiles")


# Figure
fig1<-gridExtra::grid.arrange(fig_res1,fig_res2,fig_res3,fig_res4,ncol=2)

ggsave(fig1,
        filename = paste0("./output_all/diagnostics.jpg"),
        height = 10, width = 10)

# Figure
jpeg("output_all/wp.jpg")
wp(modZIPfinal,ylim.all = 1.5, xlim.all = 5)
dev.off()


# Table
invisible(tabla_coef<- summary(modZIPfinal))

tabla_coef1 = tabla_coef %>% round(4) %>% as.data.frame() %>%
  rownames_to_column("Covariate")
tabla_coef1$Covariate = row.names(tabla_coef)


# Table

tabla_coef1 %>% add_row(Covariate ="Mu", .before = 1) %>% 
  add_row(Covariate ="Sigma", .before = 12) %>% flextable() %>% 
  theme_booktabs() %>%
  hline(i = c(1,11,12), border = officer::fp_border(color = "black", width = 2))


tabla_coef1 %>% add_row(Covariate ="Mu", .before = 1) %>% 
  add_row(Covariate ="Sigma", .before = 12) %>%
kbl(format="latex")



hmodZIPfinal <-predictAll(modZIPfinal)
hmodZIPfinal_val <-predictAll(modZIPfinal,newdata = datos_validation)

training_final <- datos_training %>% mutate(ZIPmu = hmodZIPfinal$mu,
                                            ZIPsigma = hmodZIPfinal$sigma)

valid_final <- datos_validation %>% mutate(ZIPmu = hmodZIPfinal_val$mu,
                                           ZIPsigma = hmodZIPfinal_val$sigma)

datos_final <- rbind(training_final,valid_final)

datos_final <- datos_final %>% dplyr::select(id, ZIPmu, ZIPsigma) %>%
  mutate(ZIPsigma_dummy = as.factor(ifelse(ZIPsigma>0.5,0,1)),
         ZIPmu_1 = ifelse(ZIPsigma>0.5,NA,ZIPmu))

datos_st_final <- datos_st %>% left_join(datos_final, by= "id")


BZone_mapa <-st_read('../datos/buffer zone/amortiguamiento_wgs.shp', quiet=TRUE) 
pila_mapa <-st_read('../datos/shp_PILA/pila_wgs.shp', quiet=TRUE) 

datos_map <- st_transform(datos_st_final, 'EPSG:4326')
BZone_mapa <- st_transform(BZone_mapa, 'EPSG:4326')
pila_mapa <- st_transform(pila_mapa, 'EPSG:4326')


ZIP_mu_map <- datos_map %>% ggplot() + 
  geom_sf(data = BZone_mapa, fill = "white") +  ###
  geom_sf(data = pila_mapa, fill = "white") +  ###
  geom_sf(mapping = aes(col=ZIPmu),size = 1)+
  coord_sf(crs=st_crs(4326))  + 
  coord_sf(xlim = c(-83.4, -83.1), ylim = c(9.1, 9.4)) + 
  scale_colour_distiller(palette = "Reds", direction = 1) +
  labs(title="(a)", color = "Fire count mean") +
  TEMA +
  annotation_north_arrow(location = "tr",       # top right
                         which_north = "true",  # true north
                         pad_x = unit(0.2, "in"),
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) 

# ZIP_sigma_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZIPsigma))+
#   coord_sf(crs=st_crs(3857))  + scale_colour_distiller(palette = "Greens", direction = 1) +
#   labs(color = "Probability of no fire") + TEMA +
#   annotation_north_arrow(location = "tr",       # top right
#                          which_north = "true",  # true north
#                          pad_x = unit(0.2, "in"),
#                          pad_y = unit(0.2, "in"),
#                          style = north_arrow_fancy_orienteering)


  
  
ZIP_sigma_map <- datos_map %>% ggplot() + 
  geom_sf(data = BZone_mapa, fill = "white") +  ###
  geom_sf(data = pila_mapa, fill = "white") +  ###
  geom_sf(mapping = aes(col=ZIPsigma),size = 1)+
  coord_sf(crs=st_crs(4326))  + 
  coord_sf(xlim = c(-83.4, -83.1), ylim = c(9.1, 9.4)) + 
  scale_colour_gradientn(
    colors = c("red", "yellow", "darkgreen"),
    values = scales::rescale(c(0, 0.5, 1))  # where 0.5 is the midpoint
  ) +
  #scale_colour_distiller(palette = "Greens", direction = 1) +
  labs(title="(b)",color = "Probability of no fire") + TEMA +
  annotation_north_arrow(location = "tr",       # top right
                         which_north = "true",  # true north
                         pad_x = unit(0.2, "in"),
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)




INC_map <- datos_map %>% ggplot() + 
  geom_sf(data = BZone_mapa, fill = "white") +  ###
  geom_sf(data = pila_mapa, fill = "white") +  ###
  geom_sf(mapping = aes(col=INC),size = 1)+
  coord_sf(crs=st_crs(4326))  + 
  coord_sf(xlim = c(-83.4, -83.1), ylim = c(9.1, 9.4)) + 
  scale_colour_distiller(palette = "Reds", direction = 1) +
  labs(title="(a)",color = "INC") +
  TEMA +
  annotation_north_arrow(location = "tr",       # top right
                         which_north = "true",  # true north
                         pad_x = unit(0.2, "in"),
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)

ZIP_mu_1_map <- datos_map %>% ggplot() + 
  geom_sf(data = BZone_mapa, fill = "white") +  ###
  geom_sf(data = pila_mapa, fill = "white") +  ###
  geom_sf(mapping = aes(col=ZIPmu_1),size = 1)+
  coord_sf(crs=st_crs(4326))  + 
  coord_sf(xlim = c(-83.4, -83.1), ylim = c(9.1, 9.4)) + 
  scale_colour_distiller(palette = "Reds", direction = 1, 
                         na.value = "gray80") +
  labs(title="(b)", color = "Fire count mean") +
  TEMA +
  annotation_north_arrow(location = "tr",       # top right
                         which_north = "true",  # true north
                         pad_x = unit(0.2, "in"),
                         pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)

# ggsave(ZIP_mu_map,
#        filename = paste0("./output_all/ZIP_mean.jpg"),
#        height = 7, width = 7)
# 
# ggsave(ZIP_sigma_map,
#        filename = paste0("./output_all/ZIP_sigma.jpg"),
#        height = 7, width = 7)
# ggsave(INC_map,
#        filename = paste0("./output_all/wildfire_count.jpg"),
#        height = 7, width = 7)
# ggsave(ZIP_mu_1_map,
#        filename = paste0("./output_all/wildfire_predict.jpg"),
#        height = 7, width = 7)


fig5<-gridExtra::grid.arrange(ZIP_mu_map,ZIP_sigma_map,ncol=2)
fig6<-gridExtra::grid.arrange(INC_map,ZIP_mu_1_map,ncol=2)

ggsave(fig5,
       filename = paste0("./output_all/fig5.jpg"),
       height = 7, width = 14)
    
ggsave(fig6,
       filename = paste0("./output_all/fig6.jpg"),
       height = 7, width = 14)

      

