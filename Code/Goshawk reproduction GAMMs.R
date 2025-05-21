#Goshawk reproduction GAMMs

# Import data
goshawks_final_2 <- read.csv("Data/goshawks_final_dataset_anonymised.csv", header = TRUE)
### turn site code into a factor for use as random effect in models
goshawks_final_2$Site.code.f <- as.factor(goshawks_final_2$Site.Code)

## Rain temporal lag matrix
rain_final <- read.csv("Data/Rain_temporal_lag.csv", header = TRUE)
rain_final <- as.matrix(rain_final)

## Persecution spatial lag matrix
pers_final <- read.csv("Data/Persecution_spatial_lag.csv", header = TRUE)
pers_final <- as.matrix(pers_final)

## Gamebird release pens spatial lag matrix
pens_final <- read.csv("Data/Gamebird_release_pens_spatial_lag.csv", header = TRUE)
pens_final <- as.matrix(pens_final)

## Goshawk nest density spatial lag matrix
nests_final <- read.csv("Data/Goshawk_nest_density_spatial_lag.csv", header = TRUE)
nests_final <- as.matrix(nests_final)

# Create lag indexes for models
lag_rain <- t(matrix(1:20, nrow= 20, ncol= nrow(goshawks_final_2)))
lag_rain_long <- t(matrix(1:46, nrow= 46, ncol= nrow(goshawks_final_2)))
lag_spatial <- t(matrix(1:49, nrow= 49, ncol= nrow(goshawks_final_2)))

nfkjhdkgjb


