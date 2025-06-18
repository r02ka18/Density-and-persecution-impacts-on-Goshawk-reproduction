#Goshawk reproduction GAMMs

# Import data
goshawks_final <- read.csv("Data/goshawks_final_dataset_anonymised_redacted_v2.csv", header = TRUE)
goshawks_final$Year_f <- as.factor(goshawks_final$Year_f)
goshawks_final$Site.code.f <- as.factor(goshawks_final$Site.code.f)

## Rain temporal lag matrix
rain_final <- read.csv("Data/Rain_temporal_lag.csv", header = TRUE)
rain_final <- rain_final[, 2:47]
rain_final <- apply(as.matrix.noquote(rain_final), 2, as.numeric)

## Persecution spatial lag matrix
pers_final <- read.csv("Data/Persecution_spatial_lag.csv", header = TRUE)
pers_final <- pers_final[,2:50]
pers_final <- as.matrix(pers_final)

## Gamebird release pens spatial lag matrix
pens_final <- read.csv("Data/Gamebird_release_pens_spatial_lag.csv", header = TRUE)
pens_final <- pens_final[, 2:50]
pens_final <- as.matrix(pens_final)

## Goshawk nest density spatial lag matrix
nests_final <- read.csv("Data/Goshawk_nest_density_spatial_lag.csv", header = TRUE)
nests_final <- nests_final[,2:50]
nests_final <- as.matrix(nests_final)

# Create lag indexes for models
lag_rain <- t(matrix(1:20, nrow= 20, ncol= nrow(goshawks_final)))
lag_rain_long <- t(matrix(1:46, nrow= 46, ncol= nrow(goshawks_final)))
lag_spatial <- t(matrix(1:49, nrow= 49, ncol= nrow(goshawks_final)))


# Required packages
library(mgcv)
library(ggplot2)
library(arm)
library(mgcViz)


# Density only
prod_bs_ziplss_dens <- gam(list(Number.of.Chicks.Fledged ~
                                  s(Year_f, bs = "re") +
                                  s(Site.code.f, bs ="re") +
                                  s(lag_spatial, by = nests_final, k = 49), 
                                ~ s(Year_f, bs = "re") +
                                  s(Site.code.f, bs ="re") +
                                  s(lag_spatial, by = nests_final, k = 49)), 
                           family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens)
gam.check(prod_bs_ziplss_dens)
plot(prod_bs_ziplss_dens, scale = FALSE)

prod_bs_ziplss_dens$outer.info

range(predict(prod_bs_ziplss_dens,type="response")[prod_bs_ziplss_dens$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens, type="response"), residuals(prod_bs_ziplss_dens))
plot(predict(prod_bs_ziplss_dens, type="response"), prod_bs_ziplss_dens$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens$linear.predictors[,1], prod_bs_ziplss_dens$y)
qq.gam(prod_bs_ziplss_dens,rep=20,level=1)


