library(tidyverse)
require(GGally)
library(sf)
require(ggplot2)
require(pscl)
require(MASS)
require(boot)
library(jtools)
library(gamlss)
library(gridExtra)

datos_st <-st_read('../datos/datos_301223.shp') 

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

with(datos_st, hist(FA))
with(datos_st, hist(LSF))
with(datos_st, hist(TRI))

with(datos_st, hist(logFA))
with(datos_st, hist(logLSF))
with(datos_st, hist(logTRI))


datos_df <- datos_st %>% as.data.frame() %>% 
  dplyr::select(id,INC,DAH,logFA,LST,logLSF,SLOPE,logTRI,TWI,WE,WExpo)

datos_df <- datos_df %>%
  mutate_at(c(3:11), ~scale(.))

library(corrplot)
M = cor(datos_df[,-c(1,2)])
corrplot(M, method = 'number')

#datos_df <- datos_df %>% dplyr::select(-logTRI, -logLSF)

# M1 = cor(datos_df[,-c(1,2)])
# corrplot(M1, method = 'number')


# Split data into training and testing data (0.5 and 0.5).
set.seed(4646)
rand2<-sample(2, nrow(datos_df), replace=TRUE, prob=c(0.7,0.3))

datos_training <- datos_df[rand2==1,] # training data
datos_validation <- datos_df[rand2==2,]

modPO <- gamlss(INC ~ DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo, 
                data = datos_training, 
                #newdata=datos_validation, 
                family="ZIP")

modZIP <- gamlss(INC ~ DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo, 
                 #sigma.formula = ~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                 data = datos_training, 
                 #newdata=datos_validation, 
                 family="ZIP")

modNBII <- gamlss(INC ~ DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo, 
                  #sigma.formula = ~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                  data = datos_training, 
                  #newdata=datos_validation, 
                  family="NBII")

modZINBI <- gamlss(INC ~ DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo, 
                   #sigma.formula = ~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                   #nu.formula = ~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                   data = datos_training, 
                   #newdata=datos_validation, 
                   family="ZINBI")

modPO_red <- step(modPO)
modZIP_red <- step(modZIP)
modNBII_red <- step(modNBII)
modZINBI_red <- step(modZINBI)
#modZINBI_red <- refit(modZINBI_red)

modPO_red$converged
modZIP_red$converged
modNBII_red$converged
modZINBI_red$converged

nn <- nrow(datos_training)

modZIP_red1 <- stepGAIC(modZIP_red, parameter="sigma",
                        scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                        k=log(nn))
modZIP_red1$converged

modNBII_red1 <- stepGAIC(modNBII_red, parameter="sigma",
                         scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                         k=log(nn))
modNBII_red1<- refit(modNBII_red1)
#converged


modZINBI_red1 <- stepGAIC(modZINBI_red, parameter="sigma",
                          scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                          k=log(nn))
#converged
modZINBI_red1<- refit(modZINBI_red1)
modZINBI_red1<- refit(modZINBI_red1)
modZINBI_red1$converged

modZINBI_red1a <- stepGAIC(modZINBI_red1, parameter="nu",
                           scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                           k=log(nn))
modZINBI_red1a<- refit(modZINBI_red1a)
modZINBI_red1a$converged

modZINBI_red2 <- stepGAIC(modZINBI_red, parameter="nu",
                          scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                          k=log(nn))
modZINBI_red2<- refit(modZINBI_red2)
modZINBI_red2$converged
modZINBI_red2a<- stepGAIC(modZINBI_red2, parameter="sigma",
                          scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                          k=log(nn))
modZINBI_red2a$converged
modZINBI_red2a<- refit(modZINBI_red2a)
#converged


# empezar con nu
modZINBI3 <- gamlss(INC ~ 1, 
                   nu.formula = ~1,
                   data = datos_training, 
                   #newdata=datos_validation, 
                   family="ZINBI")
modZINBI3 <- refit(modZINBI3)
modZINBI3$converged

modZINBI3_red1 <- stepGAIC(modZINBI3, parameter="nu",
                           scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                           k=log(nn))
modZINBI3_red1 <- refit(modZINBI3_red1)
modZINBI3_red1$converged

modZINBI3_red2 <- stepGAIC(modZINBI3_red1, parameter="mu",
                          scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                          k=log(nn))
modZINBI3_red2$converged
modZINBI3_red2<- refit(modZINBI3_red2)

modZINBI3_red2a <- stepGAIC(modZINBI3_red2, parameter="sigma",
                           scope=~DAH+logFA+LST+logLSF+SLOPE+logTRI+TWI+WE+WExpo,
                           k=log(nn))
modZINBI3_red2a<- refit(modZINBI3_red2a)
modZINBI3_red2a$converged

GAIC(modPO,
     modZIP,
     modNBII,
     modZINBI,
     modPO_red,
     modZIP_red,
     modNBII_red,
     modZINBI_red,
     modZIP_red1,
     modNBII_red1,
     modZINBI_red1a,
     modZINBI_red2,
     modZINBI_red2a,
     modZINBI3_red2a)

GAIC(modPO,
     modZIP,
     modNBII,
     modZINBI,
     modPO_red,
     modZIP_red,
     modNBII_red,
     modZINBI_red,
     modZIP_red1,
     modNBII_red1,
     modZINBI_red1a,
     modZINBI_red2,
     modZINBI_red2a,
     modZINBI3_red2a,k=log(nrow(datos_training)))


model.list <- list(modZIP_red1,
                   modNBII_red1,
                   modZINBI_red1a,
                   modZINBI3_red2a)

cbind(sapply(model.list,Rsq,type="both"))

GAIC.l <- cbind(sapply(model.list,GAIC,k=2))
BIC.l <- cbind(sapply(model.list,GAIC,k=log(nrow(datos_training))))
IC_table <- cbind(GAIC.l,BIC.l)

colnames(IC_table) <- c("AIC","BIC")
rownames(IC_table) <- c('modZIP_red1',
                        'modNBII_red1',
                        'modZINBI_red1a',
                        'modZINBI3_red2a')

IC_table

plot(modZIP_red1)
plot(modNBII_red1)
plot(modZINBI_red2a)
plot(modZINBI3_red2a)

par(mfrow=c(2,2))
wp(modZIP_red1,ylim.all = 1.5, xlim.all = 5)
wp(modNBII_red1,ylim.all = 1.5, xlim.all = 5)
wp(modZINBI_red2a,ylim.all = 1.5, xlim.all = 5)
wp(modZINBI3_red2a,ylim.all = 1.5, xlim.all = 10)


#fit all

modZIPfinal <- gamlss(modZIP_red1$mu.formula, 
                      sigma.formula = modZIP_red1$sigma.formula,
                      data = datos_training, 
                      newdata=datos_validation, 
                      family="ZIP")

modNBIIfinal <- gamlss(modNBII_red1$mu.formula, 
                       sigma.formula = modNBII_red1$sigma.formula,
                       data = datos_training, 
                       newdata=datos_validation, 
                       family="NBII")
modNBIIfinal$converged
modNBIIfinal <- refit(modNBIIfinal)

modZINBIfinal <- gamlss(modZINBI3_red2a$mu.formula, 
                        sigma.formula = modZINBI3_red2a$sigma.formula,
                        nu.formula = modZINBI3_red2a$nu.formula,
                        data = datos_training, 
                        newdata=datos_validation, 
                        family="ZINBI")
#modZINBIfinal <- refit(modZINBIfinal) #x3
modZINBIfinal <- refit(modZINBIfinal)
modZINBIfinal$converged


save(modPO,
     modZIP,
     modNBII,
     modZINBI,
     modPO_red,
     modZIP_red,
     modNBII_red,
     modZINBI_red,
     modZIP_red1,
     modNBII_red1,
     modZINBI_red1a,
     modZINBI_red2,
     modZINBI_red2a,
     modZINBI3_red2a, 
     # finales
     modZIPfinal,
     modNBIIfinal,
     modZINBIfinal,file = "./output/models_A4X3.RData")

load("./output/models_A4X3.RData")

model.list <- list(modZIPfinal,modNBIIfinal,modZINBIfinal)

pseudoR2 <- cbind(sapply(model.list,Rsq,type="both"))
colnames(pseudoR2) <- c('modZIPfinal','modNBIIfinal','modZINBIfinal')

pseudoR2

summary(modZIPfinal)
plot(modZIPfinal)

summary(modNBIIfinal)
plot(modNBIIfinal)

summary(modZINBIfinal)
plot(modZINBIfinal)


par(mfrow=c(2,2))
wp(modZIPfinal,ylim.all = 1.5, xlim.all = 5)
wp(modNBIIfinal,ylim.all = 1.5, xlim.all = 5)
wp(modZINBIfinal,ylim.all = 1.5, xlim.all = 5)

par(mfrow=c(1,1))

hmodZIPfinal <-predictAll(modZIPfinal)
hmodZIPfinal_val <-predictAll(modZIPfinal,newdata = datos_validation)

training_final <- datos_training %>% mutate(ZIPmu = hmodZIPfinal$mu,
                                            ZIPsigma = hmodZIPfinal$sigma)

valid_final <- datos_validation %>% mutate(ZIPmu = hmodZIPfinal_val$mu,
                                           ZIPsigma = hmodZIPfinal_val$sigma)



hmodNBIIfinal <-predictAll(modNBIIfinal)
hmodNBIIfinal_val <-predictAll(modNBIIfinal,newdata = datos_validation)

training_final <- training_final %>% mutate(NBIImu = hmodNBIIfinal$mu,
                                            NBIIsigma = hmodNBIIfinal$sigma)

valid_final <- valid_final %>% mutate(NBIImu = hmodNBIIfinal_val$mu,
                                      NBIIsigma = hmodNBIIfinal_val$sigma)


hmodZINBIfinal <-predictAll(modZINBIfinal)
hmodZINBIfinal_val <-predictAll(modZINBIfinal,newdata = datos_validation)

training_final <- training_final %>% mutate(ZINBImu = hmodZINBIfinal$mu,
                                            ZINBIsigma = hmodZINBIfinal$sigma,
                                            ZINBInu = hmodZINBIfinal$nu)

valid_final <- valid_final %>% mutate(ZINBImu = hmodZINBIfinal_val$mu,
                                      ZINBIsigma = hmodZINBIfinal_val$sigma,
                                      ZINBInu = hmodZINBIfinal_val$nu)

datos_final <- rbind(training_final,valid_final)
names(datos_final)

datos_final <- datos_final %>% select(id, ZIPmu, ZIPsigma,
                                      NBIImu,
                                      NBIIsigma,
                                      ZINBImu,
                                      ZINBIsigma,
                                      ZINBInu) %>%
  mutate(ZIPsigma_dummy = as.factor(ifelse(ZIPsigma>0.5,0,1)),
         ZIPmu_1 = ifelse(ZIPsigma>0.5,NA,ZIPmu),
         ZINBInu_dummy = as.factor(ifelse(ZINBInu>0.5,0,1)),
         ZINBImu_1 = ifelse(ZINBInu_dummy==0,NA,ZINBImu))

datos_st_final <- datos_st %>% left_join(datos_final, by= "id")

datos_map <- st_transform(datos_st_final, 'EPSG:3857')

INC_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=INC))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

ZIP_mu_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZIPmu))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

ZIP_mu_1_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZIPmu_1))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

# ZIP_mu_2_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZIPmu_10))+
#   coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)
# 
# ZIP_mu_3_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZIPmu_3))+
#   coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)


ZIP_sigma_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZIPsigma))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

ZIP_sigma1_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZIPsigma_dummy))+
  coord_sf(crs=st_crs(3857)) +scale_colour_manual(values = c("green", "red"),
                                                  labels= c("no incendio", "incendio")) 



pdf("model_ouput_ZIP.pdf")
grid.arrange(INC_map, ZIP_mu_1_map,
             #ZIP_mu_2_map,ZIP_mu_3_map,
             ZIP_mu_map, ZIP_sigma_map,ZIP_sigma1_map,
             ncol=2)
dev.off()



NBII_mu_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=NBIImu))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

NBII_sigma_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=NBIIsigma))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

pdf("model_ouput_NBII.pdf")
grid.arrange(NBII_mu_map,
             NBII_sigma_map,
             INC_map,  ncol=2)
dev.off()



ZINBI_mu_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZINBImu))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

ZINBI_sigma_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZINBIsigma))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

ZINBI_nu_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZINBInu))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

ZINBI_mu_1_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZINBImu_1))+
  coord_sf(crs=st_crs(3857))  + scale_colour_distiller(direction = 1)

ZINBI_nu1_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=ZINBInu_dummy))+
  coord_sf(crs=st_crs(3857)) +scale_colour_manual(values = c("green", "red"),
                                                  labels= c("no incendio", "incendio")) 




pdf("model_ouput_ZINBI.pdf")
grid.arrange(INC_map, ZINBI_mu_1_map, 
             ZINBI_mu_map, ZINBI_sigma_map,
             ZINBI_nu_map, ZINBI_nu1_map,
             ncol=2)
dev.off()


# covariables


DAH_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=DAH))+
  coord_sf(crs=st_crs(3857)) 
logFA_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=logFA))+
  coord_sf(crs=st_crs(3857)) 
LST_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=LST))+
  coord_sf(crs=st_crs(3857)) 
SLOPE_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=SLOPE))+
  coord_sf(crs=st_crs(3857)) 
TWI_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=TWI))+
  coord_sf(crs=st_crs(3857)) 
WE_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=WE))+
  coord_sf(crs=st_crs(3857)) 
WExpo_map <- datos_map %>% ggplot() + geom_sf(mapping = aes(col=WExpo))+
  coord_sf(crs=st_crs(3857)) 

pdf("covariates.pdf")
grid.arrange(DAH_map, LST_map, logLSF_map, 
             SLOPE_map, logTRI_map, TWI_map,
             WE_map,  WExpo_map
             , ncol=3)
dev.off()


######

names(datos_st)

plot(datos_st, max.plot=36)







