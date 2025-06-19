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

#### Density lag is significant on productivity and has a marginal effect on success. 

# Density interactions
## Elevation 
prod_bs_ziplss_dens_ele <- gam(list(Number.of.Chicks.Fledged ~
                                      s(Year_f, bs = "re") +
                                      s(Site.code.f, bs ="re") +
                                      Elevation +
                                      s(lag_spatial, by = nests_final, k = 49), 
                                    ~ s(Year_f, bs = "re") +
                                      s(Site.code.f, bs ="re") +
                                      Elevation +
                                      s(lag_spatial, by = nests_final, k = 49)), 
                               family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_ele)
gam.check(prod_bs_ziplss_dens_ele)
plot(prod_bs_ziplss_dens_ele, scale = FALSE)
prod_bs_ziplss_dens_ele$outer.info

range(predict(prod_bs_ziplss_dens_ele,type="response")[prod_bs_ziplss_dens_ele$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_ele, type="response"), residuals(prod_bs_ziplss_dens_ele))
plot(predict(prod_bs_ziplss_dens_ele, type="response"), prod_bs_ziplss_dens_ele$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_ele$linear.predictors[,1], prod_bs_ziplss_dens_ele$y)
qq.gam(prod_bs_ziplss_dens_ele,rep=20,level=1)

#### Elevation as a continuous linear variable has no effect on productivity or breeding success.

#### To test the interaction, we need to scale back the lag variable to a vector and categorize the elevation.
goshawks_final$rad_10km_dens <- rowSums(nests_final[,1:9])
## Chose 10km radius based in the results of the previous models


for (i in 1:nrow(goshawks_final)) {
  if (goshawks_final$Elevation[i] >= 355) {
    goshawks_final$Elevation_cat[i] <- "High"
  } else {
    goshawks_final$Elevation_cat[i] <- "Low"
  }
  
}
# 355m based on Sarah Hoy's paper

goshawks_final$Elevation_cat_f <- as.factor(goshawks_final$Elevation_cat)
summary(goshawks_final$Elevation_cat_f)

### Density included as a smooth
prod_bs_ziplss_dens_ele_int <- gam(list(Number.of.Chicks.Fledged ~
                                          s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Elevation_cat_f +
                                          s(rad_10km_dens) +
                                          s(rad_10km_dens, by = Elevation_cat_f, m = 1), 
                                        ~ s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Elevation_cat_f +
                                          s(rad_10km_dens) +
                                          s(rad_10km_dens, by = Elevation_cat_f, m = 1)), 
                                   family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_ele_int)
gam.check(prod_bs_ziplss_dens_ele_int)
plot(prod_bs_ziplss_dens_ele_int, scale = FALSE)
prod_bs_ziplss_dens_ele_int$outer.info

range(predict(prod_bs_ziplss_dens_ele_int, type="response")[prod_bs_ziplss_dens_ele_int$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_ele_int, type="response"), residuals(prod_bs_ziplss_dens_ele_int))
plot(predict(prod_bs_ziplss_dens_ele_int, type="response"), prod_bs_ziplss_dens_ele_int$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_ele_int$linear.predictors[,1], prod_bs_ziplss_dens_ele_int$y)
qq.gam(prod_bs_ziplss_dens_ele_int,rep=20,level=1)

#### No interaction between elevation and density.

## Density included without smooth
prod_bs_ziplss_dens_ele_int_no_smooth <- gam(list(Number.of.Chicks.Fledged ~
                                                    s(Year_f, bs = "re") +
                                                    s(Site.code.f, bs ="re") +
                                                    Elevation_cat_f +
                                                    rad_10km_dens +
                                                    rad_10km_dens:Elevation_cat_f, 
                                                  ~ s(Year_f, bs = "re") +
                                                    s(Site.code.f, bs ="re") +
                                                    Elevation_cat_f +
                                                    rad_10km_dens +
                                                    rad_10km_dens:Elevation_cat_f), 
                                             family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_ele_int_no_smooth)
gam.check(prod_bs_ziplss_dens_ele_int_no_smooth)
plot(prod_bs_ziplss_dens_ele_int_no_smooth, scale = FALSE)
prod_bs_ziplss_dens_ele_int_no_smooth$outer.info

range(predict(prod_bs_ziplss_dens_ele_int_no_smooth, type="response")[prod_bs_ziplss_dens_ele_int_no_smooth$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_ele_int_no_smooth, type="response"), residuals(prod_bs_ziplss_dens_ele_int_no_smooth))
plot(predict(prod_bs_ziplss_dens_ele_int_no_smooth, type="response"), prod_bs_ziplss_dens_ele_int_no_smooth$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_ele_int_no_smooth$linear.predictors[,1], prod_bs_ziplss_dens_ele_int_no_smooth$y)
qq.gam(prod_bs_ziplss_dens_ele_int_no_smooth,rep=20,level=1)


#Looks to be no interaction of elevation and density. When fitting the model with density 
#as a vector versus as a lag matrix, we need to fit as a smooth in the 0/1 part of the model, 
#but it's linear in the productivity part of the model. P-value for density smooth is only 
#marginal for the 0/1 part of the model, but that is consistent across all models.

## Density included as a smooth in 0/1, and linear in productivity
prod_bs_ziplss_dens_ele_int2 <- gam(list(Number.of.Chicks.Fledged ~
                                           s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Elevation_cat_f +
                                           rad_10km_dens +
                                           s(rad_10km_dens, by = Elevation_cat_f, m = 1), 
                                         ~ s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Elevation_cat_f +
                                           s(rad_10km_dens) +
                                           s(rad_10km_dens, by = Elevation_cat_f, m = 1)), 
                                    family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_ele_int2)
gam.check(prod_bs_ziplss_dens_ele_int2)
plot(prod_bs_ziplss_dens_ele_int2, scale = FALSE)
prod_bs_ziplss_dens_ele_int2$outer.info

range(predict(prod_bs_ziplss_dens_ele_int2, type="response")[prod_bs_ziplss_dens_ele_int2$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_ele_int2, type="response"), residuals(prod_bs_ziplss_dens_ele_int2))
plot(predict(prod_bs_ziplss_dens_ele_int2, type="response"), prod_bs_ziplss_dens_ele_int2$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_ele_int2$linear.predictors[,1], prod_bs_ziplss_dens_ele_int2$y)
qq.gam(prod_bs_ziplss_dens_ele_int2,rep=20,level=1)

#Can conclude that there is no evidence of an effect of elevation or an elevation/density 
#interaction on breeding success or productivity so can discard from the model. 
#When fitting rad_10km_dens, this is linear in the productivity model and smooth in the 
#success part of the model.

## Raptor Worker/Surveyor
### Next variable to test is raptor worker, as a way of separating areas. 
### This is a categorical to be fitted as a fixed effect, and an interaction with density.

prod_bs_ziplss_dens_rw_int <- gam(list(Number.of.Chicks.Fledged ~
                                         s(Year_f, bs = "re") +
                                         s(Site.code.f, bs ="re") +
                                         Surveyor_Af +
                                         rad_10km_dens,
                                       #rad_10km_dens:Surveyor_f,
                                       ~ s(Year_f, bs = "re") +
                                         s(Site.code.f, bs ="re") +
                                         Surveyor_Af +
                                         s(rad_10km_dens) +
                                         s(rad_10km_dens, by = Surveyor_Af, m = 1)), 
                                  family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_rw_int)
gam.check(prod_bs_ziplss_dens_rw_int)
plot(prod_bs_ziplss_dens_rw_int, scale = FALSE)
prod_bs_ziplss_dens_rw_int$outer.info

range(predict(prod_bs_ziplss_dens_rw_int, type="response")[prod_bs_ziplss_dens_rw_int$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_rw_int, type="response"), residuals(prod_bs_ziplss_dens_rw_int))
plot(predict(prod_bs_ziplss_dens_rw_int, type="response"), prod_bs_ziplss_dens_rw_int$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_rw_int$linear.predictors[,1], prod_bs_ziplss_dens_rw_int$y)
qq.gam(prod_bs_ziplss_dens_rw_int,rep=20,level=1)

## Interaction appears to be significant in the success part of the model, for S006 and S007.

## Ownership
goshawks_final$Ownership.f <- as.factor(goshawks_final$Ownership)
prod_bs_ziplss_dens_own_int <- gam(list(Number.of.Chicks.Fledged ~
                                          s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Ownership.f +
                                          rad_10km_dens +
                                          rad_10km_dens:Ownership.f,
                                        ~ s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Ownership.f +
                                          s(rad_10km_dens) +
                                          s(rad_10km_dens, by = Ownership.f, m = 1)), 
                                   family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_own_int)
gam.check(prod_bs_ziplss_dens_own_int)
plot(prod_bs_ziplss_dens_own_int, scale = FALSE)
prod_bs_ziplss_dens_own_int$outer.info

range(predict(prod_bs_ziplss_dens_own_int, type="response")[prod_bs_ziplss_dens_own_int$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_own_int, type="response"), residuals(prod_bs_ziplss_dens_own_int))
plot(predict(prod_bs_ziplss_dens_own_int, type="response"), prod_bs_ziplss_dens_own_int$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_own_int$linear.predictors[,1], prod_bs_ziplss_dens_own_int$y)
qq.gam(prod_bs_ziplss_dens_own_int,rep=20,level=1)

## Ownership on its own in the productivity part of the model is significant. 
## Interaction of density and ownership is insignificant in both parts of the model.

## Conclude
#Density is significant in the productivity part of the model as a lag and as a vector. 
#Density is marginal in the breeding success model as lag and as a vector. 
#Ownership is significant in the productivity part of the model.

#The density-raptor worker interaction is significant in the breeding success part 
#of the model, with PG and TL/ES being different. This could not be tested in the 
#productivity part of the model.

#There was no evidence of an interaction between elevation and density, or ownership 
#and density, in either part of the model.

# Persecution Interactions
#Next we'll add in persecution and explore some interactions with it. We are only 
#testing persecution effect on success as productivity assumes success which would not 
#be the case when persecution occurs.

goshawks_final$rad_10km_pers <- rowSums(pers_final[,1:9])

prod_bs_ziplss_dens_rw_pers <- gam(list(Number.of.Chicks.Fledged ~
                                          s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Ownership.f +
                                          s(lag_spatial, by = nests_final, k = 49), 
                                        ~ s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Surveyor_Af +
                                          s(lag_spatial, by = pers_final, k =  49) +
                                          s(rad_10km_dens) +
                                          s(rad_10km_dens, by = Surveyor_Af, m = 1)), 
                                   family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_rw_pers)
gam.check(prod_bs_ziplss_dens_rw_pers)
plot(prod_bs_ziplss_dens_rw_pers, scale = FALSE)
prod_bs_ziplss_dens_rw_pers$outer.info

range(predict(prod_bs_ziplss_dens_rw_pers,type="response")[prod_bs_ziplss_dens_rw_pers$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_rw_pers, type="response"), residuals(prod_bs_ziplss_dens_rw_pers))
plot(predict(prod_bs_ziplss_dens_rw_pers, type="response"), prod_bs_ziplss_dens_rw_pers$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_rw_pers$linear.predictors[,1], prod_bs_ziplss_dens_rw_pers$y)
qq.gam(prod_bs_ziplss_dens_rw_pers,rep=20,level=1)

### No spatial effect of persecution on success. Now Ownership is not significant in productivity model.

## Raptor worker/Surveyor
prod_bs_ziplss_dens_rw_pers_int_prod <- gam(list(Number.of.Chicks.Fledged ~
                                                   s(Year_f, bs = "re") +
                                                   s(Site.code.f, bs ="re") +
                                                   Ownership.f +
                                                   s(lag_spatial, by = nests_final, k = 49), 
                                                 ~ s(Year_f, bs = "re") +
                                                   s(Site.code.f, bs ="re") +
                                                   Surveyor_Af +
                                                   s(rad_10km_pers) +
                                                   s(rad_10km_dens) +
                                                   s(rad_10km_pers, by = Surveyor_Af, m = 1)), 
                                            family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_rw_pers_int_prod)
gam.check(prod_bs_ziplss_dens_rw_pers_int_prod)
plot(prod_bs_ziplss_dens_rw_pers_int_prod, scale = FALSE)
prod_bs_ziplss_dens_rw_pers_int_prod$outer.info

range(predict(prod_bs_ziplss_dens_rw_pers_int_prod, type="response")[prod_bs_ziplss_dens_rw_pers_int_prod$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_rw_pers_int_prod, type="response"), residuals(prod_bs_ziplss_dens_rw_pers_int_prod))
plot(predict(prod_bs_ziplss_dens_rw_pers_int_prod, type="response"), prod_bs_ziplss_dens_rw_pers_int_prod$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_rw_pers_int_prod$linear.predictors[,1], prod_bs_ziplss_dens_rw_pers_int_prod$y)
qq.gam(prod_bs_ziplss_dens_rw_pers_int_prod,rep=20,level=1)

## Persecution-surveyor interaction significant (S003) in success part of the model. Ownership no longer significant in productivity model.

## Elevation
prod_bs_ziplss_dens_rw_pers_int_ele <- gam(list(Number.of.Chicks.Fledged ~
                                                  s(Year_f, bs = "re") +
                                                  s(Site.code.f, bs ="re") +
                                                  Ownership.f +
                                                  s(lag_spatial, by = nests_final, k = 49), 
                                                ~ s(Year_f, bs = "re") +
                                                  s(Site.code.f, bs ="re") +
                                                  Surveyor_Af +
                                                  Elevation_cat_f +
                                                  s(rad_10km_pers) +
                                                  s(rad_10km_dens) +
                                                  s(rad_10km_pers, by = Elevation_cat_f, m = 1) +
                                                  s(rad_10km_dens, by = Surveyor_Af, m = 1)), 
                                           family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_rw_pers_int_ele)
gam.check(prod_bs_ziplss_dens_rw_pers_int_ele)
plot(prod_bs_ziplss_dens_rw_pers_int_ele, scale = FALSE)
prod_bs_ziplss_dens_rw_pers_int_ele$outer.info

range(predict(prod_bs_ziplss_dens_rw_pers_int_ele, type="response")[prod_bs_ziplss_dens_rw_pers_int_ele$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_rw_pers_int_ele, type="response"), residuals(prod_bs_ziplss_dens_rw_pers_int_ele))
plot(predict(prod_bs_ziplss_dens_rw_pers_int_ele, type="response"), prod_bs_ziplss_dens_rw_pers_int_ele$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_rw_pers_int_ele$linear.predictors[,1], prod_bs_ziplss_dens_rw_pers_int_ele$y)
qq.gam(prod_bs_ziplss_dens_rw_pers_int_ele,rep=20,level=1)


#### No interaction of elevation and persecution, elevation can again be removed from the model.

## Ownership
goshawks_final$Ownership_f <- as.factor(goshawks_final$Ownership)
prod_bs_ziplss_dens_rw_pers_int_own <- gam(list(Number.of.Chicks.Fledged ~
                                                  s(Year_f, bs = "re") +
                                                  s(Site.code.f, bs ="re") +
                                                  Ownership_f +
                                                  s(lag_spatial, by = nests_final, k = 49), 
                                                ~ s(Year_f, bs = "re") +
                                                  s(Site.code.f, bs ="re") +
                                                  Surveyor_Af +
                                                  Ownership_f +
                                                  s(rad_10km_pers) +
                                                  s(rad_10km_dens) +
                                                  s(rad_10km_pers, by = Ownership_f, m = 1) +
                                                  s(rad_10km_dens, by = Surveyor_Af, m = 1)), 
                                           family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_rw_pers_int_own)
gam.check(prod_bs_ziplss_dens_rw_pers_int_own)
plot(prod_bs_ziplss_dens_rw_pers_int_own, scale = FALSE)
prod_bs_ziplss_dens_rw_pers_int_own$outer.info

range(predict(prod_bs_ziplss_dens_rw_pers_int_own, type="response")[prod_bs_ziplss_dens_rw_pers_int_own$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_rw_pers_int_own, type="response"), residuals(prod_bs_ziplss_dens_rw_pers_int_own))
plot(predict(prod_bs_ziplss_dens_rw_pers_int_own, type="response"), prod_bs_ziplss_dens_rw_pers_int_own$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_rw_pers_int_own$linear.predictors[,1], prod_bs_ziplss_dens_rw_pers_int_own$y)
qq.gam(prod_bs_ziplss_dens_rw_pers_int_own,rep=20,level=1)

#### A marginal p-value in the state owned sites in the success part of the model. 

## Conclude
# Persecution has no effect on it's own but there is an interaction between persecution 
# and raptor worker (MD) on breeding success.
# Ownership was no longer significant in any of these models.

# Gamebird release pens
prod_bs_ziplss_dens_pens <- gam(list(Number.of.Chicks.Fledged ~
                                       s(Year_f, bs = "re") +
                                       s(Site.code.f, bs ="re") +
                                       Ownership.f +
                                       s(lag_spatial, by = pens_final, k = 49) +
                                       s(lag_spatial, by = nests_final, k = 49), 
                                     ~ s(Year_f, bs = "re") +
                                       s(Site.code.f, bs ="re") +
                                       Surveyor_Af +
                                       s(rad_10km_dens) +
                                       #s(rad_10km_pers) +  
                                       s(lag_spatial, by = pens_final, k = 49) +
                                       #s(rad_10km_pers, by = Surveyor_f, m = 1) +
                                       s(rad_10km_dens, by = Surveyor_Af, m = 1)), 
                                family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_pens)
gam.check(prod_bs_ziplss_dens_pens)
plot(prod_bs_ziplss_dens_pens, scale = FALSE)
prod_bs_ziplss_dens_pens$outer.info

range(predict(prod_bs_ziplss_dens_pens,type="response")[prod_bs_ziplss_dens_pens$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_pens, type="response"), residuals(prod_bs_ziplss_dens_pens))
plot(predict(prod_bs_ziplss_dens_pens, type="response"), prod_bs_ziplss_dens_pens$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_pens$linear.predictors[,1], prod_bs_ziplss_dens_pens$y)
qq.gam(prod_bs_ziplss_dens_pens,rep=20,level=1)

#### No effect of pens on its own on breeding success or productivity.

## Raptor worker/Surveyor
goshawks_final$rad_10km_pens <- rowSums(pens_final[,1:9])


prod_bs_ziplss_dens_pens_rw <- gam(list(Number.of.Chicks.Fledged ~
                                          s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Ownership.f +
                                          Surveyor_Af +
                                          s(rad_10km_pens) +
                                          s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                          s(lag_spatial, by = nests_final, k = 49), 
                                        ~ s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          Surveyor_Af +
                                          s(rad_10km_dens) +
                                          s(rad_10km_pens) +
                                          s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                          s(rad_10km_dens, by = Surveyor_Af, m = 1)), 
                                   family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_pens_rw)
gam.check(prod_bs_ziplss_dens_pens_rw)
plot(prod_bs_ziplss_dens_pens_rw, scale = FALSE)
prod_bs_ziplss_dens_pens_rw$outer.info

range(predict(prod_bs_ziplss_dens_pens_rw,type="response")[prod_bs_ziplss_dens_pens_rw$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_pens_rw, type="response"), residuals(prod_bs_ziplss_dens_pens_rw))
plot(predict(prod_bs_ziplss_dens_pens_rw, type="response"), prod_bs_ziplss_dens_pens_rw$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_pens_rw$linear.predictors[,1], prod_bs_ziplss_dens_pens_rw$y)
qq.gam(prod_bs_ziplss_dens_pens_rw,rep=20,level=1)

#### Model failed to converge. Repeat but with density as lag variable

prod_bs_ziplss_dens_pens_rw2 <- gam(list(Number.of.Chicks.Fledged ~
                                           s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Surveyor_Af +
                                           s(rad_10km_pens) +
                                           s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                           s(lag_spatial, by = nests_final, k = 49), 
                                         ~ s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Surveyor_Af +
                                           s(rad_10km_pens) +
                                           s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                           s(lag_spatial, by = nests_final, k = 49)), 
                                    family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_pens_rw2)
gam.check(prod_bs_ziplss_dens_pens_rw2)
plot(prod_bs_ziplss_dens_pens_rw2, scale = FALSE)
prod_bs_ziplss_dens_pens_rw2$outer.info

range(predict(prod_bs_ziplss_dens_pens_rw2,type="response")[prod_bs_ziplss_dens_pens_rw2$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_pens_rw2, type="response"), residuals(prod_bs_ziplss_dens_pens_rw2))
plot(predict(prod_bs_ziplss_dens_pens_rw2, type="response"), prod_bs_ziplss_dens_pens_rw2$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_pens_rw2$linear.predictors[,1], prod_bs_ziplss_dens_pens_rw2$y)
qq.gam(prod_bs_ziplss_dens_pens_rw2,rep=20,level=1)

#### Evidence for an interaction between raptor worker and gamebird pens in the success 
#### part of the model, with a significant difference between S006 and the rest. No evidence 
#### of an interaction between pens and rw in the productivity part of the model.

## Elevation

prod_bs_ziplss_dens_pens_ele <- gam(list(Number.of.Chicks.Fledged ~
                                           s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Ownership.f +
                                           s(rad_10km_pens) +
                                           Elevation_cat_f +
                                           s(rad_10km_pens, by = Elevation_cat_f, m=1) +
                                           s(lag_spatial, by = nests_final, k = 49), 
                                         ~ s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Surveyor_Af +
                                           s(rad_10km_pens) +
                                           Elevation_cat_f +
                                           s(rad_10km_pens, by = Elevation_cat_f, m=1) +
                                           s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                           s(lag_spatial, by = nests_final, k = 49)), 
                                    family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_pens_ele)
gam.check(prod_bs_ziplss_dens_pens_ele)
plot(prod_bs_ziplss_dens_pens_ele, scale = FALSE)
prod_bs_ziplss_dens_pens_ele$outer.info

range(predict(prod_bs_ziplss_dens_pens_ele,type="response")[prod_bs_ziplss_dens_pens_ele$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_pens_ele, type="response"), residuals(prod_bs_ziplss_dens_pens_ele))
plot(predict(prod_bs_ziplss_dens_pens_ele, type="response"), prod_bs_ziplss_dens_pens_ele$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_pens_ele$linear.predictors[,1], prod_bs_ziplss_dens_pens_ele$y)
qq.gam(prod_bs_ziplss_dens_pens_ele,rep=20,level=1)

#### No interaction between elevation and pens so elevation can again be removed from the model.

## Ownership
prod_bs_ziplss_dens_pens_own <- gam(list(Number.of.Chicks.Fledged ~
                                           s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           s(rad_10km_pens) +
                                           Ownership_f +
                                           s(rad_10km_pens, by = Ownership_f, m=1) +
                                           s(lag_spatial, by = nests_final, k = 49), 
                                         ~ s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Surveyor_Af +
                                           s(rad_10km_pens) +
                                           Ownership_f +
                                           s(rad_10km_pens, by = Ownership_f, m=1) +
                                           s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                           s(lag_spatial, by = nests_final, k = 49)), 
                                    family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_pens_own)
gam.check(prod_bs_ziplss_dens_pens_own)
plot(prod_bs_ziplss_dens_pens_own, scale = FALSE)
prod_bs_ziplss_dens_pens_own$outer.info

range(predict(prod_bs_ziplss_dens_pens_own,type="response")[prod_bs_ziplss_dens_pens_own$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_pens_own, type="response"), residuals(prod_bs_ziplss_dens_pens_own))
plot(predict(prod_bs_ziplss_dens_pens_own, type="response"), prod_bs_ziplss_dens_pens_own$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_pens_own$linear.predictors[,1], prod_bs_ziplss_dens_pens_own$y)
qq.gam(prod_bs_ziplss_dens_pens_own,rep=20,level=1)

#### No effect of ownership interaction with pens.

## Conclude:
#### Pens interact with raptor worker in the success part of the model, but have no 
#### significant interactions with elevation or ownership.
#### There is no effect of pens or any interaction on productivity.

# Rainfall
prod_bs_ziplss_dens_rain <- gam(list(Number.of.Chicks.Fledged ~
                                       s(Year_f, bs = "re") +
                                       s(Site.code.f, bs ="re") +
                                       s(lag_rain_long, by = rain_final, k = 46) +
                                       s(lag_spatial, by = nests_final, k = 49), 
                                     ~ s(Year_f, bs = "re") +
                                       s(Site.code.f, bs ="re") +
                                       Surveyor_Af +
                                       s(lag_rain_long, by = rain_final, k = 46) +
                                       s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                       s(lag_spatial, by = nests_final, k = 49)), 
                                family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_rain)
gam.check(prod_bs_ziplss_dens_rain)
plot(prod_bs_ziplss_dens_rain, scale = FALSE)
prod_bs_ziplss_dens_rain$outer.info

range(predict(prod_bs_ziplss_dens_rain,type="response")[prod_bs_ziplss_dens_rain$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_rain, type="response"), residuals(prod_bs_ziplss_dens_rain))
plot(predict(prod_bs_ziplss_dens_rain, type="response"), prod_bs_ziplss_dens_rain$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_rain$linear.predictors[,1], prod_bs_ziplss_dens_rain$y)
qq.gam(prod_bs_ziplss_dens_rain,rep=20,level=1)

### Consecutive days rain:
prod_bs_ziplss_dens_rain_con <- gam(list(Number.of.Chicks.Fledged ~
                                           s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           s(Rain_con) +
                                           s(lag_spatial, by = nests_final, k = 49), 
                                         ~ s(Year_f, bs = "re") +
                                           s(Site.code.f, bs ="re") +
                                           Surveyor_Af +
                                           s(Rain_con) +
                                           s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                           s(lag_spatial, by = nests_final, k = 49)), 
                                    family = ziplss(), data = goshawks_final)

summary(prod_bs_ziplss_dens_rain_con)
gam.check(prod_bs_ziplss_dens_rain_con)
plot(prod_bs_ziplss_dens_rain_con, scale = FALSE)
prod_bs_ziplss_dens_rain_con$outer.info

range(predict(prod_bs_ziplss_dens_rain_con,type="response")[prod_bs_ziplss_dens_rain_con$y==0])


par(mfrow=c(2,2))
plot(predict(prod_bs_ziplss_dens_rain_con, type="response"), residuals(prod_bs_ziplss_dens_rain_con))
plot(predict(prod_bs_ziplss_dens_rain_con, type="response"), prod_bs_ziplss_dens_rain_con$y);abline(0,1,col=2)
plot(prod_bs_ziplss_dens_rain_con$linear.predictors[,1], prod_bs_ziplss_dens_rain_con$y)
qq.gam(prod_bs_ziplss_dens_rain_con,rep=20,level=1)

## Conclude

#### Rainfall nor consecutive days of rain are significant on their own. The rain 
#### variable was taken from a single location so there should not be any interaction with 
#### spatial variables such as raptor worker, elevation and ownership.

# Final model
gos_prod_bs_ziplss_final <- gam(list(Number.of.Chicks.Fledged ~
                                       s(Year_f, bs = "re") +
                                       s(Site.code.f, bs ="re") +
                                       Ownership.f +
                                       s(lag_spatial, by = nests_final, k = 49), 
                                     ~ s(Year_f, bs = "re") +
                                       s(Site.code.f, bs ="re") +
                                       Surveyor_Af +
                                       s(rad_10km_pens) +
                                       s(rad_10km_pers) +
                                       s(rad_10km_dens) +
                                       s(rad_10km_pens, by = Surveyor_Af, m = 1) +
                                       s(rad_10km_pers, by = Surveyor_Af, m = 1) +
                                       s(rad_10km_dens, by = Surveyor_Af, m = 1)),
                                family = ziplss(), data = goshawks_final)

summary(gos_prod_bs_ziplss_final)
gam.check(gos_prod_bs_ziplss_final)
plot(gos_prod_bs_ziplss_final, scale = FALSE)
gos_prod_bs_ziplss_final$outer.info

range(predict(gos_prod_bs_ziplss_final,type="response")[gos_prod_bs_ziplss_final$y==0])


par(mfrow=c(2,2))
plot(predict(gos_prod_bs_ziplss_final, type="response"), residuals(gos_prod_bs_ziplss_final))
plot(predict(gos_prod_bs_ziplss_final, type="response"), gos_prod_bs_ziplss_final$y);abline(0,1,col=2)
plot(gos_prod_bs_ziplss_final$linear.predictors[,1], gos_prod_bs_ziplss_final$y)
qq.gam(gos_prod_bs_ziplss_final,rep=20,level=1)

## This model fails to converge, I think because we have raptor worker in all three interactions.
# Check confoundedness of surveyor and site
table(goshawks_final$Site.code.f, goshawks_final$Surveyor_Af)

#### Site and surveyor are confounded, and so would county/country with surveyor so we should remove raptor worker from the model.

# Removal raptor worker - Density only
gos_prod_bs_ziplss_dens <- gam(list(Number.of.Chicks.Fledged ~
                                      s(Year_f, bs = "re") +
                                      s(Site.code.f, bs ="re") +
                                      s(lag_spatial, by = nests_final, k = 49), 
                                    ~ s(Year_f, bs = "re") +
                                      s(Site.code.f, bs ="re") +
                                      s(lag_spatial, by = nests_final, k = 49)),
                               family = ziplss(), data = goshawks_final)

gam.check(gos_prod_bs_ziplss_dens)
plot(gos_prod_bs_ziplss_dens, scale = FALSE)
gos_prod_bs_ziplss_dens$outer.info
summary(gos_prod_bs_ziplss_dens)
range(predict(gos_prod_bs_ziplss_dens,type="response")[gos_prod_bs_ziplss_dens$y==0])


par(mfrow=c(2,2))
plot(predict(gos_prod_bs_ziplss_dens, type="response"), residuals(gos_prod_bs_ziplss_dens))
plot(predict(gos_prod_bs_ziplss_dens, type="response"), gos_prod_bs_ziplss_dens$y);abline(0,1,col=2)
plot(gos_prod_bs_ziplss_dens$linear.predictors[,1], gos_prod_bs_ziplss_dens$y)
qq.gam(gos_prod_bs_ziplss_dens,rep=20,level=1)

# Plot model output
library(gratia)
library(patchwork)
gos_prod_bs_ziplss_dens_draw <- gratia::draw(gos_prod_bs_ziplss_dens)

gos_prod_bs_ziplss_dens_draw_prod <- gratia::draw(gos_prod_bs_ziplss_dens, select = "s(lag_spatial):nests_final")
gos_prod_bs_ziplss_dens_draw_bs <- gratia::draw(gos_prod_bs_ziplss_dens, select = "s.1(lag_spatial):nests_final")

p1 <- gos_prod_bs_ziplss_dens_draw_prod + 
  geom_hline(yintercept = 0, colour = "red") +
  labs(x = "Distance (km)", title = "Productivity", y = "Effect Size") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"))

p2 <- gos_prod_bs_ziplss_dens_draw_bs +
  geom_hline(yintercept = 0, colour = "red") +
  labs(x = "Distance (km)", title = "Breeding Success", y = "Effect Size") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"))

library(reshape2)
apply(nests_final, 2, mean) 
apply(nests_final, 2, max) 
apply(nests_final, 2, min) 

nests_long <- melt(nests_final, value.name = "Density")

colnames(nests_long)[2] <- "Distance_in_m"

nests_long$Distance <- nests_long$Distance_in_m/1000

distances <- unique(nests_long$Distance)

nests_dens <- data.frame(Distance = distances, Density = NA, density.se = NA)

for (i in 1:length(distances)) {
  temp_gos <- nests_long[nests_long$Distance == distances[i], ]
  nests_dens$Density[i] <- mean(temp_gos$Density)
  nests_dens$density.se[i] <- sd(temp_gos$Density)/sqrt(nrow(temp_gos))
}

nests_dens$upper_ci <- nests_dens$Density + 1.96*nests_dens$density.se
nests_dens$lower_ci <- nests_dens$Density - 1.96*nests_dens$density.se

nests_dens

nests_dens_plot <- ggplot(nests_dens, aes(Distance, Density)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci)) +
  geom_point() +
  geom_line() +
  ylab("Mean Nest Density") +
  xlab("Distance (km)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"))


nests_dens_plot


gos_dens_lag_plots <- nests_dens_plot + p2 + p1 +
  plot_layout(ncol = 3, nrow = 1, widths = c(5,5,5), heights = c(3,3,3)) +
  plot_annotation(tag_levels = "A")

## To make predictions, we need to switch to the vector density variable instead of the lag. 

gos_prod_bs_ziplss_dens_vec <- gam(list(Number.of.Chicks.Fledged ~
                                          s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          rad_10km_dens, 
                                        ~ s(Year_f, bs = "re") +
                                          s(Site.code.f, bs ="re") +
                                          rad_10km_dens),
                                   family = ziplss(), data = goshawks_final)

gam.check(gos_prod_bs_ziplss_dens_vec)
plot(gos_prod_bs_ziplss_dens_vec, scale = FALSE, all.terms = TRUE)
gos_prod_bs_ziplss_dens_vec $outer.info
summary(gos_prod_bs_ziplss_dens_vec)
range(predict(gos_prod_bs_ziplss_dens_vec ,type="response")[gos_prod_bs_ziplss_dens_vec$y==0])


par(mfrow=c(2,2))
plot(predict(gos_prod_bs_ziplss_dens_vec, type="response"), residuals(gos_prod_bs_ziplss_dens_vec))
plot(predict(gos_prod_bs_ziplss_dens_vec, type="response"), gos_prod_bs_ziplss_dens_vec$y);abline(0,1,col=2)
plot(gos_prod_bs_ziplss_dens_vec$linear.predictors[,1], gos_prod_bs_ziplss_dens_vec$y)
qq.gam(gos_prod_bs_ziplss_dens_vec,rep=20,level=1)

gam.vcomp(gos_prod_bs_ziplss_dens_vec)
variance_comp(gos_prod_bs_ziplss_dens_vec)

gos_bs_prod_prediction <- predict(gos_prod_bs_ziplss_dens_vec, type = "link", terms = NULL, se = TRUE)

# Same as linear predictors above
gos_bs_prod_prediction_df <- as.data.frame(gos_bs_prod_prediction)

## Back transform the predictions
ilink <- binomial(link = "cloglog")$linkinv

newd <- transform(gos_bs_prod_prediction, fit.1.trans = exp(gos_bs_prod_prediction$fit[,1]), fit.2.trans = ilink(gos_bs_prod_prediction$fit[,2]), fitted.trans = exp(gos_bs_prod_prediction$fit[,1]) * ilink(gos_bs_prod_prediction$fit[,2]))

newd$fit.1_lower95 <- newd$fit.1 - 1.96*newd$se.fit.1
newd$fit.1_upper95 <- newd$fit.1 + 1.96*newd$se.fit.1

newd$fit.1_lower95_trans <- exp(newd$fit.1_lower95)
newd$fit.1_upper95_trans <- exp(newd$fit.1_upper95)

newd$fit.2_lower95 <- newd$fit.2 - 1.96*newd$se.fit.2
newd$fit.2_upper95 <- newd$fit.2 + 1.96*newd$se.fit.2

newd$fit.2_lower95_trans <- ilink(newd$fit.2_lower95)
newd$fit.2_upper95_trans <- ilink(newd$fit.2_upper95)

goshawks_final_2 <- goshawks_final[!is.na(goshawks_final$Number.of.Chicks.Fledged),]
goshawks_predict <- cbind(goshawks_final_2, newd)

## Plot the predictions
nest_dens <- unique(goshawks_final_2$rad_10km_dens)

gos_prod_dens_predict_df <- data.frame(rad_10km_dens = nest_dens, Site.code.f = "G001", Year_f = "1981")

gos_prod_dens_predict <- predict(gos_prod_bs_ziplss_dens_vec, gos_prod_dens_predict_df, exclude = c("s(site_f)", "s(Year_f)", "s.1(site_f)", "s.1(Year_f)"), type = "link", terms = NULL, se = TRUE)

gos_prod_dens_predict_df2 <- as.data.frame(gos_prod_dens_predict)

## Back transform the predictions
# ilink <- binomial(link = "cloglog")$linkinv
# 
# gos_prod_dens_predict_df2  <- transform(gos_prod_dens_predict_df2, fit.1.trans = exp(gos_prod_dens_predict$fit[,1]), fit.2.trans = ilink(gos_prod_dens_predict$fit[,2]))

### TC's edits
gos_prod_dens_predict_df_overall <-  (1-exp(-exp(gos_prod_dens_predict_df2[, 2]))) * (exp(gos_prod_dens_predict_df2[, 1]) / (1 - dpois(0, exp(gos_prod_dens_predict_df2[, 1]))))

gos_prod_dens_predict_count <- exp(gos_prod_dens_predict_df2[, 1]) / (1 - dpois(0, exp(gos_prod_dens_predict_df2[, 1])))

gos_prod_dens_predict_zeros <- (1-exp(-exp(gos_prod_dens_predict_df2[, 2])))

gos_prod_dens_predict_df_new <- cbind(nest_dens, gos_prod_dens_predict_df2, gos_prod_dens_predict_count, gos_prod_dens_predict_zeros)
colnames(gos_prod_dens_predict_df_new)[1] <- "Density"

## Calculate 95% CIs and back transform

gos_prod_dens_predict_df_new$fit.1_lower95 <- gos_prod_dens_predict_df_new$fit.1 - 1.96*gos_prod_dens_predict_df_new$se.fit.1
gos_prod_dens_predict_df_new$fit.1_upper95 <- gos_prod_dens_predict_df_new$fit.1 + 1.96*gos_prod_dens_predict_df_new$se.fit.1

gos_prod_dens_predict_df_new$fit.1_lower95_trans <- exp(gos_prod_dens_predict_df_new$fit.1_lower95) / (1 - dpois(0, exp(gos_prod_dens_predict_df_new$fit.1_lower95)))
gos_prod_dens_predict_df_new$fit.1_upper95_trans <- exp(gos_prod_dens_predict_df_new$fit.1_upper95) / (1 - dpois(0, exp(gos_prod_dens_predict_df_new$fit.1_upper95)))

gos_prod_dens_predict_df_new$fit.2_lower95 <- gos_prod_dens_predict_df_new$fit.2 - 1.96*gos_prod_dens_predict_df_new$se.fit.2
gos_prod_dens_predict_df_new$fit.2_upper95 <- gos_prod_dens_predict_df_new$fit.2 + 1.96*gos_prod_dens_predict_df_new$se.fit.2

gos_prod_dens_predict_df_new$fit.2_lower95_trans <- (1-exp(-exp(gos_prod_dens_predict_df_new$fit.2_lower95)))
gos_prod_dens_predict_df_new$fit.2_upper95_trans <- (1-exp(-exp(gos_prod_dens_predict_df_new$fit.2_upper95)))

gos_prod_dens_predict_df_new


### Plot raw data
## Breeding success
nest_dens

dens_raw <- data.frame(Density = nest_dens, Breeding_Success = NA)

for (i in 1:length(nest_dens)) {
  temp_gos <- goshawks_final_2[goshawks_final_2$rad_10km_dens == nest_dens[i], ]
  dens_raw$Breeding_Success[i] <- nrow(temp_gos[temp_gos$Breeding.Success == 1,]) / nrow(temp_gos)
  
}

##Productivity
nest_dens

dens_raw_prod <- data.frame(Density = nest_dens, Productivity = NA, productivity.se = NA)

for (i in 1:length(nest_dens)) {
  temp_gos <- goshawks_final_2[goshawks_final_2$rad_10km_dens == nest_dens[i], ]
  temp_gos_success <- temp_gos[temp_gos$Breeding.Success == 1,]
  dens_raw_prod$Productivity[i] <- mean(temp_gos_success$Number.of.Chicks.Fledged)
  dens_raw_prod$productivity.se[i] <- sd(temp_gos_success$Number.of.Chicks.Fledged)/sqrt(nrow(temp_gos_success))
}

dens_raw_prod$upper_ci <- dens_raw_prod$Productivity + 1.96*dens_raw_prod$productivity.se
dens_raw_prod$lower_ci <- dens_raw_prod$Productivity - 1.96*dens_raw_prod$productivity.se


prod_dens_merge2 <- merge(gos_prod_dens_predict_df_new, dens_raw_prod, by ="Density")
gos_prod_dens_plot2 <- ggplot(data = prod_dens_merge2 , aes(Density, gos_prod_dens_predict_count)) +
  geom_ribbon(aes(ymin = fit.1_lower95_trans, ymax = fit.1_upper95_trans), fill = "grey") +
  geom_line(colour = "black") +
  geom_errorbar(data = prod_dens_merge2, mapping = aes(x = Density, ymin = lower_ci, ymax = upper_ci)) +
  geom_point(data = prod_dens_merge2, mapping = aes(Density, Productivity)) +
  xlab("Nest Count in 10km Radius") +
  ylab("Number of Chicks Fledged") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black")) 

gos_prod_dens_plot2


prod_dens_merge3 <- merge(gos_prod_dens_predict_df_new, dens_raw, by ="Density")
gos_bs_dens_plot2 <- ggplot(prod_dens_merge3, aes(Density, gos_prod_dens_predict_zeros)) +
  geom_ribbon(aes(ymin = fit.2_lower95_trans, ymax = fit.2_upper95_trans), fill = "grey") +
  geom_line(colour = "black") +
  geom_point(prod_dens_merge3, mapping = aes(Density, Breeding_Success)) +
  xlab("Nest Count in 10km Radius") +
  ylab("Probability of Success") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black")) 


gos_bs_dens_plot2

library(ggpubr)

gos_bs_prod_plots2 <- ggarrange(gos_bs_dens_plot2, gos_prod_dens_plot2, ncol = 2, labels = c("A", "B"))
gos_bs_prod_plots2

# Investigate random effects
plot(gos_prod_bs_ziplss_dens_vec)
#### Actually look worst in the count part of the model (s()) compared to the zero part of the model (s.1()).

concurvity(gos_prod_bs_ziplss_dens_vec, full = TRUE)
concurvity(gos_prod_bs_ziplss_dens_vec, full = FALSE)


## Run a model without site as a random effect in the zero part of the model:
gos_prod_bs_ziplss_dens_vec_2 <- gam(list(Number.of.Chicks.Fledged ~
                                            s(Year_f, bs = "re") +
                                            s(Site.code.f, bs ="re") +
                                            rad_10km_dens, 
                                          ~ s(Year_f, bs = "re") +
                                            rad_10km_dens),
                                     family = ziplss(), data = goshawks_final_2)

gam.check(gos_prod_bs_ziplss_dens_vec_2)
plot(gos_prod_bs_ziplss_dens_vec_2, scale = FALSE, all.terms = TRUE)
gos_prod_bs_ziplss_dens_vec_2$outer.info
summary(gos_prod_bs_ziplss_dens_vec_2)
range(predict(gos_prod_bs_ziplss_dens_vec_2 ,type="response")[gos_prod_bs_ziplss_dens_vec_2$y==0])

par(mfrow=c(2,2))
plot(predict(gos_prod_bs_ziplss_dens_vec_2, type="response"), residuals(gos_prod_bs_ziplss_dens_vec_2))
plot(predict(gos_prod_bs_ziplss_dens_vec_2, type="response"), gos_prod_bs_ziplss_dens_vec_2$y);abline(0,1,col=2)
plot(gos_prod_bs_ziplss_dens_vec_2$linear.predictors[,1], gos_prod_bs_ziplss_dens_vec_2$y)
qq.gam(gos_prod_bs_ziplss_dens_vec_2,rep=20,level=1)

concurvity(gos_prod_bs_ziplss_dens_vec_2, full = TRUE)

concurvity(gos_prod_bs_ziplss_dens_vec_2, full = FALSE)

# Plot the predictions no site in zeros
nest_dens <- unique(goshawks_final_2$rad_10km_dens)

gos_prod_dens_predict_df <- data.frame(rad_10km_dens = nest_dens, Site.code.f = "G001", Year_f = "1981")

gos_prod_dens_predict_no_site <- predict(gos_prod_bs_ziplss_dens_vec_2, gos_prod_dens_predict_df, exclude = c("s(site_f)", "s(Year_f)", "s.1(site_f)", "s.1(Year_f)"), type = "link", terms = NULL, se = TRUE)

gos_prod_dens_predict_df2_no_site <- as.data.frame(gos_prod_dens_predict_no_site)

## Back transform the predictions
# ilink <- binomial(link = "cloglog")$linkinv
# 
# gos_prod_dens_predict_df2  <- transform(gos_prod_dens_predict_df2, fit.1.trans = exp(gos_prod_dens_predict$fit[,1]), fit.2.trans = ilink(gos_prod_dens_predict$fit[,2]))

### TC's edits
gos_prod_dens_predict_df_overall_no_site <-  (1-exp(-exp(gos_prod_dens_predict_df2_no_site[, 2]))) * (exp(gos_prod_dens_predict_df2_no_site[, 1]) / (1 - dpois(0, exp(gos_prod_dens_predict_df2_no_site[, 1]))))

gos_prod_dens_predict_count_no_site <- exp(gos_prod_dens_predict_df2_no_site[, 1]) / (1 - dpois(0, exp(gos_prod_dens_predict_df2_no_site[, 1])))

gos_prod_dens_predict_zeros_no_site <- (1-exp(-exp(gos_prod_dens_predict_df2_no_site[, 2])))

gos_prod_dens_predict_df_new_no_site <- cbind(nest_dens, gos_prod_dens_predict_df2_no_site, gos_prod_dens_predict_count_no_site, gos_prod_dens_predict_zeros_no_site)
colnames(gos_prod_dens_predict_df_new_no_site)[1] <- "Density"

## Calculate 95% CIs and back transform

gos_prod_dens_predict_df_new_no_site$fit.1_lower95 <- gos_prod_dens_predict_df_new_no_site$fit.1 - 1.96*gos_prod_dens_predict_df_new_no_site$se.fit.1
gos_prod_dens_predict_df_new_no_site$fit.1_upper95 <- gos_prod_dens_predict_df_new_no_site$fit.1 + 1.96*gos_prod_dens_predict_df_new_no_site$se.fit.1

gos_prod_dens_predict_df_new_no_site$fit.1_lower95_trans <- exp(gos_prod_dens_predict_df_new_no_site$fit.1_lower95) / (1 - dpois(0, exp(gos_prod_dens_predict_df_new_no_site$fit.1_lower95)))
gos_prod_dens_predict_df_new_no_site$fit.1_upper95_trans <- exp(gos_prod_dens_predict_df_new_no_site$fit.1_upper95) / (1 - dpois(0, exp(gos_prod_dens_predict_df_new_no_site$fit.1_upper95)))

gos_prod_dens_predict_df_new_no_site$fit.2_lower95 <- gos_prod_dens_predict_df_new_no_site$fit.2 - 1.96*gos_prod_dens_predict_df_new_no_site$se.fit.2
gos_prod_dens_predict_df_new_no_site$fit.2_upper95 <- gos_prod_dens_predict_df_new_no_site$fit.2 + 1.96*gos_prod_dens_predict_df_new_no_site$se.fit.2

gos_prod_dens_predict_df_new_no_site$fit.2_lower95_trans <- (1-exp(-exp(gos_prod_dens_predict_df_new_no_site$fit.2_lower95)))
gos_prod_dens_predict_df_new_no_site$fit.2_upper95_trans <- (1-exp(-exp(gos_prod_dens_predict_df_new_no_site$fit.2_upper95)))

gos_prod_dens_predict_df_new_no_site


### Plot raw data
## Breeding success
nest_dens

dens_raw_no_site <- data.frame(Density = nest_dens, Breeding_Success = NA)

for (i in 1:length(nest_dens)) {
  temp_gos <- goshawks_final_2[goshawks_final_2$rad_10km_dens == nest_dens[i], ]
  dens_raw_no_site$Breeding_Success[i] <- nrow(temp_gos[temp_gos$Breeding.Success == 1,]) / nrow(temp_gos)
  
}

##Productivity
nest_dens

dens_raw_prod_no_site <- data.frame(Density = nest_dens, Productivity = NA, productivity.se = NA)

for (i in 1:length(nest_dens)) {
  temp_gos <- goshawks_final_2[goshawks_final_2$rad_10km_dens == nest_dens[i], ]
  temp_gos_success <- temp_gos[temp_gos$Breeding.Success == 1,]
  dens_raw_prod_no_site$Productivity[i] <- mean(temp_gos_success$Number.of.Chicks.Fledged)
  dens_raw_prod_no_site$productivity.se[i] <- sd(temp_gos_success$Number.of.Chicks.Fledged)/sqrt(nrow(temp_gos_success))
}

dens_raw_prod_no_site$upper_ci <- dens_raw_prod_no_site$Productivity + 1.96*dens_raw_prod_no_site$productivity.se
dens_raw_prod_no_site$lower_ci <- dens_raw_prod_no_site$Productivity - 1.96*dens_raw_prod_no_site$productivity.se

prod_dens_merge2_no_site <- merge(gos_prod_dens_predict_df_new_no_site, dens_raw_prod_no_site, by ="Density")
gos_prod_dens_plot2_no_site <- ggplot(data = prod_dens_merge2_no_site , aes(Density, gos_prod_dens_predict_count_no_site)) +
  geom_ribbon(aes(ymin = fit.1_lower95_trans, ymax = fit.1_upper95_trans), fill = "grey") +
  geom_line(colour = "black") +
  geom_errorbar(data = prod_dens_merge2_no_site, mapping = aes(x = Density, ymin = lower_ci, ymax = upper_ci)) +
  geom_point(data = prod_dens_merge2_no_site, mapping = aes(Density, Productivity)) +
  xlab("Nest Count in 10km Radius") +
  ylab("Number of Chicks Fledged") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black")) 

gos_prod_dens_plot2_no_site

prod_dens_merge3_no_site <- merge(gos_prod_dens_predict_df_new_no_site, dens_raw_no_site, by ="Density")
gos_bs_dens_plot2_no_site <- ggplot(prod_dens_merge3_no_site, aes(Density, gos_prod_dens_predict_zeros_no_site)) +
  geom_ribbon(aes(ymin = fit.2_lower95_trans, ymax = fit.2_upper95_trans), fill = "grey") +
  geom_line(colour = "black") +
  geom_point(prod_dens_merge3_no_site, mapping = aes(Density, Breeding_Success)) +
  xlab("Nest Count in 10km Radius") +
  ylab("Probability of Success") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black")) 


gos_bs_dens_plot2_no_site

gos_bs_prod_plots2_no_site <- ggarrange(gos_bs_dens_plot2_no_site, gos_prod_dens_plot2_no_site, ncol = 2, labels = c("A", "B"))
gos_bs_prod_plots2_no_site


# Run a model with year in the zero part and site in the count part:
gos_prod_bs_ziplss_dens_vec_3 <- gam(list(Number.of.Chicks.Fledged ~
                                            s(Site.code.f, bs ="re") +
                                            rad_10km_dens, 
                                          ~ s(Year_f, bs = "re") +
                                            rad_10km_dens),
                                     family = ziplss(), data = goshawks_final_2)

gam.check(gos_prod_bs_ziplss_dens_vec_3)
plot(gos_prod_bs_ziplss_dens_vec_3, scale = FALSE, all.terms = TRUE)
gos_prod_bs_ziplss_dens_vec_3$outer.info
summary(gos_prod_bs_ziplss_dens_vec_3)
range(predict(gos_prod_bs_ziplss_dens_vec_3 ,type="response")[gos_prod_bs_ziplss_dens_vec_3$y==0])


par(mfrow=c(2,2))
plot(predict(gos_prod_bs_ziplss_dens_vec_3, type="response"), residuals(gos_prod_bs_ziplss_dens_vec_3))
plot(predict(gos_prod_bs_ziplss_dens_vec_3, type="response"), gos_prod_bs_ziplss_dens_vec_3$y);abline(0,1,col=2)
plot(gos_prod_bs_ziplss_dens_vec_3$linear.predictors[,1], gos_prod_bs_ziplss_dens_vec_3$y)
qq.gam(gos_prod_bs_ziplss_dens_vec_3,rep=20,level=1)

concurvity(gos_prod_bs_ziplss_dens_vec_3, full = TRUE)

concurvity(gos_prod_bs_ziplss_dens_vec_3, full = FALSE)

## Plot the predictions for year in zeros and site in count
nest_dens <- unique(goshawks_final_2$rad_10km_dens)

gos_prod_dens_predict_df <- data.frame(rad_10km_dens = nest_dens, Site.code.f = "G001", Year_f = "1981")

gos_prod_dens_predict_yzsc <- predict(gos_prod_bs_ziplss_dens_vec_3, gos_prod_dens_predict_df, exclude = c("s(site_f)", "s(Year_f)", "s.1(site_f)", "s.1(Year_f)"), type = "link", terms = NULL, se = TRUE)

gos_prod_dens_predict_df2_yzsc <- as.data.frame(gos_prod_dens_predict_yzsc)


### TC's edits
gos_prod_dens_predict_df_overall_yzsc <-  (1-exp(-exp(gos_prod_dens_predict_df2_yzsc[, 2]))) * (exp(gos_prod_dens_predict_df2_yzsc[, 1]) / (1 - dpois(0, exp(gos_prod_dens_predict_df2_yzsc[, 1]))))

gos_prod_dens_predict_count_yzsc <- exp(gos_prod_dens_predict_df2_yzsc[, 1]) / (1 - dpois(0, exp(gos_prod_dens_predict_df2_yzsc[, 1])))

gos_prod_dens_predict_zeros_yzsc <- (1-exp(-exp(gos_prod_dens_predict_df2_yzsc[, 2])))

gos_prod_dens_predict_df_new_yzsc <- cbind(nest_dens, gos_prod_dens_predict_df2_yzsc, gos_prod_dens_predict_count_yzsc, gos_prod_dens_predict_zeros_yzsc)
colnames(gos_prod_dens_predict_df_new_yzsc)[1] <- "Density"

## Calculate 95% CIs and back transform

gos_prod_dens_predict_df_new_yzsc$fit.1_lower95 <- gos_prod_dens_predict_df_new_yzsc$fit.1 - 1.96*gos_prod_dens_predict_df_new_yzsc$se.fit.1
gos_prod_dens_predict_df_new_yzsc$fit.1_upper95 <- gos_prod_dens_predict_df_new_yzsc$fit.1 + 1.96*gos_prod_dens_predict_df_new_yzsc$se.fit.1

gos_prod_dens_predict_df_new_yzsc$fit.1_lower95_trans <- exp(gos_prod_dens_predict_df_new_yzsc$fit.1_lower95) / (1 - dpois(0, exp(gos_prod_dens_predict_df_new_yzsc$fit.1_lower95)))
gos_prod_dens_predict_df_new_yzsc$fit.1_upper95_trans <- exp(gos_prod_dens_predict_df_new_yzsc$fit.1_upper95) / (1 - dpois(0, exp(gos_prod_dens_predict_df_new_yzsc$fit.1_upper95)))

gos_prod_dens_predict_df_new_yzsc$fit.2_lower95 <- gos_prod_dens_predict_df_new_yzsc$fit.2 - 1.96*gos_prod_dens_predict_df_new_yzsc$se.fit.2
gos_prod_dens_predict_df_new_yzsc$fit.2_upper95 <- gos_prod_dens_predict_df_new_yzsc$fit.2 + 1.96*gos_prod_dens_predict_df_new_yzsc$se.fit.2

gos_prod_dens_predict_df_new_yzsc$fit.2_lower95_trans <- (1-exp(-exp(gos_prod_dens_predict_df_new_yzsc$fit.2_lower95)))
gos_prod_dens_predict_df_new_yzsc$fit.2_upper95_trans <- (1-exp(-exp(gos_prod_dens_predict_df_new_yzsc$fit.2_upper95)))

gos_prod_dens_predict_df_new_yzsc


### Plot raw data
## Breeding success
nest_dens

dens_raw_yzsc <- data.frame(Density = nest_dens, Breeding_Success = NA)

for (i in 1:length(nest_dens)) {
  temp_gos <- goshawks_final_2[goshawks_final_2$rad_10km_dens == nest_dens[i], ]
  dens_raw_yzsc$Breeding_Success[i] <- nrow(temp_gos[temp_gos$Breeding.Success == 1,]) / nrow(temp_gos)
  
}

##Productivity
nest_dens

dens_raw_prod_yzsc <- data.frame(Density = nest_dens, Productivity = NA, productivity.se = NA)

for (i in 1:length(nest_dens)) {
  temp_gos <- goshawks_final_2[goshawks_final_2$rad_10km_dens == nest_dens[i], ]
  temp_gos_success <- temp_gos[temp_gos$Breeding.Success == 1,]
  dens_raw_prod_yzsc$Productivity[i] <- mean(temp_gos_success$Number.of.Chicks.Fledged)
  dens_raw_prod_yzsc$productivity.se[i] <- sd(temp_gos_success$Number.of.Chicks.Fledged)/sqrt(nrow(temp_gos_success))
}

dens_raw_prod_yzsc$upper_ci <- dens_raw_prod_yzsc$Productivity + 1.96*dens_raw_prod_yzsc$productivity.se
dens_raw_prod_yzsc$lower_ci <- dens_raw_prod_yzsc$Productivity - 1.96*dens_raw_prod_yzsc$productivity.se


prod_dens_merge2_yzsc <- merge(gos_prod_dens_predict_df_new_yzsc, dens_raw_prod_yzsc, by ="Density")
gos_prod_dens_plot2_yzsc <- ggplot(data = prod_dens_merge2_yzsc , aes(Density, gos_prod_dens_predict_count_yzsc)) +
  geom_ribbon(aes(ymin = fit.1_lower95_trans, ymax = fit.1_upper95_trans), fill = "grey") +
  geom_line(colour = "black") +
  geom_errorbar(data = prod_dens_merge2_yzsc, mapping = aes(x = Density, ymin = lower_ci, ymax = upper_ci)) +
  geom_point(data = prod_dens_merge2_yzsc, mapping = aes(Density, Productivity)) +
  xlab("Nest Count in 10km Radius") +
  ylab("Number of Chicks Fledged") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black")) 

gos_prod_dens_plot2_yzsc


prod_dens_merge3_yzsc <- merge(gos_prod_dens_predict_df_new_yzsc, dens_raw_yzsc, by ="Density")
gos_bs_dens_plot2_yzsc <- ggplot(prod_dens_merge3_yzsc, aes(Density, gos_prod_dens_predict_zeros_yzsc)) +
  geom_ribbon(aes(ymin = fit.2_lower95_trans, ymax = fit.2_upper95_trans), fill = "grey") +
  geom_line(colour = "black") +
  geom_point(prod_dens_merge3_yzsc, mapping = aes(Density, Breeding_Success)) +
  xlab("Nest Count in 10km Radius") +
  ylab("Probability of Success") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black")) 


gos_bs_dens_plot2_yzsc

gos_bs_prod_plots2_yzsc <- ggarrange(gos_bs_dens_plot2_yzsc, gos_prod_dens_plot2_yzsc, ncol = 2, labels = c("A", "B"))
gos_bs_prod_plots2_yzsc

ggsave("Plots/Goshawk_BS_Prod_Density_raw_data_yzsc.jpg", gos_bs_prod_plots2_yzsc, width = 10, height = 5)
