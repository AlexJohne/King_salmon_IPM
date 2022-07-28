library(nlme)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(modelr)
library(broom)
library(fields)
library(dplyr)
library(mdthemes)

# 1. data ----
## 1.1. Full data set ----
caw_data <- readRDS("caw_data1.rds")

cawdata <- caw_data %>% 
  select(z=weight, z1=Nextweight, FL, ration, timepoint, ID, time, tank, daily_FI, FCR) %>%
  mutate(ration=case_when(ration == "A" ~ "60S",ration == "B" ~ "80S",ration == "C" ~ "100S"))

d <- cawdata %>%
  select(tank,ID,ration,time,z, z1) %>%
  filter(!is.na(z)) %>%
  mutate(ration=factor(ration,levels=c("60S","80S","100S")))

## 1.2. Ration subsets ----
d60 <- d %>% filter(ration=="60S") %>% filter(time < 276)
d80 <- d %>% filter(ration=="80S") %>% filter(time < 276)
d100 <- d %>% filter(ration=="100S") %>% filter(time < 276)


# 2. Fitting simple linear regressions to each time * ration combination of the raw data ----

by_ration_time <- d %>%
  filter(time < 276) %>%
  group_by(ration, time) %>%   #### this also works with group_nest(ration, time)
  nest()

ration_time_model <- function(df) {
  lm(z1 ~ z,  data = df)
}

models <- map(by_ration_time$data, ration_time_model)

by_ration_time <- by_ration_time %>% 
  mutate(model = map(data, ration_time_model))

by_ration_time %>% 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(glance)

coefs <- by_ration_time %>%
  mutate(tidied = map(model, broom::tidy)) %>%
  unnest(tidied)

pred <- by_ration_time %>%
  mutate(preds = map(model, broom::augment)) %>%
  unnest(preds)

ggplot(pred, aes(z, z1, colour = ration)) +
  geom_point(alpha = 0.5) +
  facet_grid(ration ~ time) +
  geom_point(mapping=aes(x=z,y=.fitted),color="black",size=0.4) +
  theme_classic() +
  guides(colour = "none")
  

# 3. Forcing the linear regression to go through the origin ----

# If we force the fitted line through the origin, the model fits are much better 
# (especially for the early time steps), 
# but may only be reasonable for logistic growth (which is totally fine).

fo_model <- function(df) {
  lm(z1 ~ 0 + z,  data = df)
}

fo_models <- map(by_ration_time$data, fo_model)

forced_origin <- by_ration_time %>% 
  mutate(model = map(data, fo_model))

fo_modelfit <- forced_origin %>%                                 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(glance)     
### THIS HAS SIGNIFICANTLY IMPROVED THE MODEL FIT (according to R-squared, but not AIC)

fo_coefs <- forced_origin %>%
  mutate(tidied = map(model, broom::tidy, conf.int = TRUE)) %>%
  unnest(tidied)

fo_pred <- forced_origin %>%
  mutate(preds = map(model, broom::augment, conf.int = TRUE)) %>%
  unnest(preds) %>%
  select(ration, time, z1, z, .fitted, .resid, .hat, .sigma, .cooksd, .std.resid)

# statistic = t-value
# sigma = SD
# variance = sigma^2

## 3.1. Plotting the predicted values over the real data ----
ggplot(fo_pred, aes(z, z1, colour = ration)) +
  geom_point(alpha = 0.5) +
  facet_grid(ration ~ time) +
  geom_point(mapping=aes(x=z,y=.fitted),color="black",size=0.4) +
  theme_classic() +
  guides(colour = "none") +
  labs(x = "Body mass (g) at t",
       y = "Body mass (g) at t+1") +
  scale_x_continuous(breaks = c(250, 750, 1250))
ggsave("Figure2.jpeg")

## 3.2. Diagnostic plots ----
### Residuals vs fitted ----
ggplot(fo_pred, aes(.fitted, .resid)) +
  geom_point(alpha = 0.5) +
  facet_grid(ration ~ time) +
  geom_smooth(method = loess, colour = "red", se = FALSE) +
  guides(colour = "none") +
  theme_classic()
## these look like clouds with no trends --> perfect

ggplot(fo_pred, aes(.fitted, .std.resid)) +
  geom_point(alpha = 0.5) +
  facet_grid(ration ~ time) +
  geom_smooth(method = loess, colour = "red", se = FALSE) +
  guides(colour = "none") +
  theme_classic()
# pretty much look like trend-less clouds as well

### QQ plot of residuals ----
ggplot(fo_pred, aes(sample = .resid)) +
  stat_qq(alpha = 0.5) + 
  stat_qq_line(colour = "red") +
  facet_grid(ration ~ time) +
  guides(colour = "none") +
  theme_classic() +
  labs(y = "Quantiles of standard normal", x = "Residuals")
# looks pretty normal to me

## not enough data to add a random effect for ID

# 4. IPM ----
## 4.1. parameter list for regressions ----
params <- full_join(fo_coefs, fo_modelfit, by = c("ration", "time")) %>%
  select(ration, time, estimate, sigma) %>%
  rename(growth.slope = estimate, growth.sd = sigma)

## 4.2. growth function ----
g.yx=function(xp,x,params) {            
  dnorm(xp, mean = 0 + params$growth.slope*x,
        sd=params$growth.sd)
}

## 4.3. Setting up the matrix ----
# make a function
growth_kernel = function(data, n, params){
min.size=.9*min(c(data$z,data$z1),na.rm=T)
max.size=1.1*max(c(data$z,data$z1),na.rm=T)
# n = number of cells in the discretized kernel
# boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 
# mesh points (midpoints of the cells)
y=0.5*(b[1:n]+b[2:(n+1)])

# width of the cells
h=y[2]-y[1]

G=h*outer(y,y,g.yx,params=params)
rownames(G) <- y
colnames(G) <- y
return(G) # G is a matrix (useful for image.plot)
}

g1 <- growth_kernel(data = d60, n = 100, params = params[1,])
g2 <- growth_kernel(data = d60, n = 100, params = params[4,])
g3 <- growth_kernel(data = d60, n = 100, params = params[7,])
g4 <- growth_kernel(data = d60, n = 100, params = params[10,])
g5 <- growth_kernel(data = d60, n = 100, params = params[13,])

# setting up the ration-specific matrix
min.size=.9*min(c(d60$z,d60$z1),na.rm=T)
max.size=1.1*max(c(d60$z,d60$z1),na.rm=T)
# number of cells in the discretized kernel
n=100 
# boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 
# mesh points (midpoints of the cells)
y60=0.5*(b[1:n]+b[2:(n+1)])

# width of the cells
h=y[2]-y[1]

## 4.4. Limiting the matrix according to size for time segment ----
### 4.4.1. ration 60S ----
max(d60[ which(d60$time == 0),6], na.rm = T)
# g1 <- g1[1:22, 1:22]
g1[,22:100] <- 0
g1[22:100,] <- 0

max(d60[ which(d60$time == 91),6], na.rm = T)
min(d60[ which(d60$time == 91),5], na.rm = T)
# g2 <- g2[6:37, 6:37]
g2[1:6,] <- 0
g2[, 1:6] <- 0
g2[37:100,] <- 0
g2[, 37:100] <- 0

max(d60[ which(d60$time == 124),6], na.rm = T)
min(d60[ which(d60$time == 124),5], na.rm = T)
# g3 <- g3[5:53, 5:53]
g3[1:5,] <- 0
g3[, 1:5] <- 0
g3[53:100,] <- 0
g3[, 53:100] <- 0

max(d60[ which(d60$time == 173),6], na.rm = T)
min(d60[ which(d60$time == 173),5], na.rm = T)
# g4 <- g4[6:68, 6:68]
g4[1:6,] <- 0
g4[, 1:6] <- 0
g4[68:100,] <- 0
g4[, 68:100] <- 0

max(d60[ which(d60$time == 221),6], na.rm = T)
min(d60[ which(d60$time == 221),5], na.rm = T)
# g5 <- g5[23:100, 23:100]
g5[1:23,] <- 0
g5[, 1:23] <- 0

### 4.4.2. ration 80S ----
g6 <- growth_kernel(data = d80, n = 100, params = params[2,])
g7 <- growth_kernel(data = d80, n = 100, params = params[5,])
g8 <- growth_kernel(data = d80, n = 100, params = params[8,])
g9 <- growth_kernel(data = d80, n = 100, params = params[11,])
g10 <- growth_kernel(data = d80, n = 100, params = params[14,])

# setting up the ration-specific matrix
min.size=.9*min(c(d80$z,d80$z1),na.rm=T)
max.size=1.1*max(c(d80$z,d80$z1),na.rm=T)
# number of cells in the discretized kernel
n=100 
# boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 
# mesh points (midpoints of the cells)
y80=0.5*(b[1:n]+b[2:(n+1)])

# width of the cells
h=y[2]-y[1]


max(d80[ which(d80$time == 0),6], na.rm = T)
g6[,16:100] <- 0
g6[16:100,] <- 0

max(d80[ which(d80$time == 91),6], na.rm = T)
min(d80[ which(d80$time == 91),5], na.rm = T)
g7[1:5,] <- 0
g7[, 1:5] <- 0
g7[29:100,] <- 0
g7[, 29:100] <- 0

max(d80[ which(d80$time == 124),6], na.rm = T)
min(d80[ which(d80$time == 124),5], na.rm = T)
g8[1:5,] <- 0
g8[, 1:5] <- 0
g8[48:100,] <- 0
g8[, 48:100] <- 0

max(d80[ which(d80$time == 173),6], na.rm = T)
min(d80[ which(d80$time == 173),5], na.rm = T)
g9[1:10,] <- 0
g9[, 1:10] <- 0
g9[76:100,] <- 0
g9[, 76:100] <- 0

min(d80[ which(d80$time == 221),5], na.rm = T)
g10[1:14,] <- 0
g10[, 1:14] <- 0

### 4.4.3. ration 100S ----
g11 <- growth_kernel(data = d100, n = 100, params = params[3,])
g12 <- growth_kernel(data = d100, n = 100, params = params[6,])
g13 <- growth_kernel(data = d100, n = 100, params = params[9,])
g14 <- growth_kernel(data = d100, n = 100, params = params[12,])
g15 <- growth_kernel(data = d100, n = 100, params = params[15,])

# setting up the ration-specific matrix
min.size=.9*min(c(d100$z,d100$z1),na.rm=T)
max.size=1.1*max(c(d100$z,d100$z1),na.rm=T)
# number of cells in the discretized kernel
n=100 
# boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 
# mesh points (midpoints of the cells)
y100=0.5*(b[1:n]+b[2:(n+1)])

# width of the cells
h=y[2]-y[1]

max(d100[ which(d100$time == 0),6], na.rm = T)
g11[,18:100] <- 0
g11[18:100,] <- 0

max(d100[ which(d100$time == 91),6], na.rm = T)
min(d100[ which(d100$time == 91),5], na.rm = T)
g12[1:4,] <- 0
g12[, 1:4] <- 0
g12[31:100,] <- 0
g12[, 31:100] <- 0

max(d100[ which(d100$time == 124),6], na.rm = T)
min(d100[ which(d100$time == 124),5], na.rm = T)
g13[1:4,] <- 0
g13[, 1:4] <- 0
g13[57:100,] <- 0
g13[, 57:100] <- 0

max(d100[ which(d100$time == 173),6], na.rm = T)
min(d100[ which(d100$time == 173),5], na.rm = T)
g14[1:12,] <- 0
g14[, 1:12] <- 0
g14[77:100,] <- 0
g14[, 77:100] <- 0

min(d100[ which(d100$time == 221),5], na.rm = T)
g15[1:12,] <- 0
g15[, 1:12] <- 0

# 5. Plots ----
## 5.1. Probability density of data overlaid with predictions ----
d.pr3 <- d.pr %>% filter(time > 0) %>%
  group_by(ration, time) %>%
  bind_cols(time2 = time2) %>%
  rename(t = time, time = time2)
  

time2 <- rep(c(0, 91, 124, 173, 221), 1625)

ggplot(fo_pred, aes(z1, y=..density.., colour =  "Data")) +
  geom_histogram(alpha = 0.2) +
  geom_density(d.pr3, kernel = "gaussian", mapping = aes(pred, colour = "Logistic growth model")) +
  geom_density(fo_pred, kernel = "gaussian", mapping = aes(.fitted, colour = "Linear Regressions")) +
  facet_grid(ration~time) +
  theme_classic() +
  labs(colour="", 
       x="Body mass (g) at t+1",
       y="Probability density") +
  scale_colour_manual(values=c("black", "red", "orange")) +
  scale_x_continuous(breaks = c(250, 750, 1250)) +
  theme(legend.position="bottom",
        axis.title.y = ggtext::element_markdown())
ggsave("data_vs_pred_dist.jpeg")

## 5.2. Matrix plots ----
tiff(file="Growth_kernels.tiff",
     width=4, height=4, units="in", res=500)
#jpeg("Growth_kernels.jpeg")
par(mfrow = c(3,1), mai = c(.4,.4,.1,.1), 
    oma=c(2.5,2.5,0.5,0.5),  mar = c(2, 2, 0.5, 0.5), cex.lab = 0.5)

# 60%
image.plot(y60,y60,t(g1 + g2 + g3 + g4 + g5)) 
abline(0,1, lty=2, lwd=3)  
legend("topleft", paste0("IPM 60S"), text.font = 2, cex= 1, bg = 'white')

# 80%
image.plot(y80,y80,t(g6 + g7 + g8 + g9 + g10)) 
abline(0,1, lty=2, lwd=3)
legend("topleft", paste0("IPM 80S"), text.font = 2, cex= 1, bg = 'white')

# 100 ration
image.plot(y100,y100,t(g11 + g12 + g13 + g14 + g15)) 
abline(0,1, lty=2, lwd=3)
legend("topleft", paste0("IPM 100S"), text.font = 2, cex= 1, bg = 'white')

mtext("Body mass at t (g)",side=1,line=1,outer=TRUE,cex=0.8)
mtext("Body mass at t+1 (g)",side=2,line=1,outer=TRUE,cex=0.8,las=0)
dev.off()

# topographical color scheme -->  col = topo.colors(50)

# 6. Logistic growth model ----
## 6.1. Initial data plot ----
ggplot(d,aes(x=time,y=z,group=ID,colour=tank))+
  geom_line(alpha=0.1)+
  facet_wrap(~ration)+
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme_classic() +
  labs(colour = "Tank", 
       x="Time (days)",
       y="Body mass (g)") +
  theme(legend.position = "bottom")
ggsave("growth_by_ration_tank.jpeg")


## 6.2. Fitting logistic growth model ----
# Getting rid of correlation
d <- d%>%
  mutate(w = z/1500)

# initial model that ignores tank effects
fit0 <- nlme(w ~ wmax/(1+(wmax-w0)/w0*exp(-K*time)),
             data=d,
             fixed=wmax+K+w0 ~ ration,
             random=wmax+K+w0~1|ID,
             start=c(0.5,0.1,0.2,
                     0.01,0,0,
                     0.02,0,0))
summary(fit0)

# refitting the model with random tank effects using the previous results as starting values
fit <- nlme(w ~ wmax/(1+(wmax-w0)/w0*exp(-K*time)),
            data=d,
            fixed=wmax+K+w0 ~ ration,
            random=wmax+K+w0~1|tank/ID,
            start=fixef(fit0),
            control = nlmeControl(maxIter=80,msMaxIter=120,tolerance=1.0E-5,rel.tol=1.0E-7))
summary(fit)

# making a summary table for model (potentially extracting CIs)
fit.ci <- intervals(fit, level = 0.95, which = "fixed")
fit.ci <- as.data.frame(do.call(cbind, fit.ci))

fit.unlist <- unlist(fit$coefficients$fixed)
fit.coef <- data.frame(as.list(fit.unlist))

# plotting predicted over real data
d.pr <- unique(d[,c("ID","tank","ration")])
d.pr <- cbind(d.pr[rep(1:nrow(d.pr),each=6),],time=c(0,91,124,173,221,276))
d.pr$W <- predict(fit,d.pr)
d.pr <- d.pr %>%
  mutate(pred = W*1500)

d1 <- d %>% group_by(ration, time) %>%
  summarise(mean_z = mean(z), sd_z = sd(z), min_z = min(z), max_z = max(z))

d.pr2 <- d.pr %>% group_by(ration, time) %>%
  summarise(mean_pred = mean(pred), min_pred = min(pred), max_pred = max(pred))

fo_pred2 <- fo_pred %>% group_by(ration, time) %>%
  summarise(mean_fit = mean(.fitted), min_fit = min(.fitted), max_fit = max(.fitted)) %>%
  mutate(time2 = c(91, 124, 173, 221, 276))

ggplot(d,aes(x=time,y=z, group=ID))+
  geom_line(alpha=0.1) +
  geom_pointrange(d.pr2, mapping = aes(time, mean_pred, ymin = min_pred, ymax = max_pred, group = ration, color = "Prediction range - logistic model")) +
  geom_pointrange(fo_pred2, mapping = aes(time2-10, mean_fit, ymin = min_fit, ymax = max_fit, group = ration, color = "Prediction range - linear submodels")) +
  facet_wrap(~ration) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_manual(values=c("red", "orange")) +
  theme_classic() +
  labs(colour = "", 
       x="Time (days)",
       y="Body mass (g)") +
  theme(legend.position = "bottom")
ggsave("data_plus_pred_growth2.jpeg")

## 6.3. Calculating mean percentage deviation (MPD %) between model and data ----
# MPD = 100/n * sum((actual_val - pred_val)/actual_val)

d2 <- d.pr %>% group_by(ration, time) %>%
  summarise(mean_pred = mean(pred), sd_pred = sd(pred), min_pred = min(pred), max_pred = max(pred))

d2 <- d2 %>% 
  left_join(d1, by = c("ration", "time")) %>%
  mutate(mean_err = abs((mean_z-mean_pred)/mean_z),
         min_err = abs((min_z-min_pred)/min_z),
         max_err = abs((max_z-max_pred)/max_z),
         cv = 100*sd_z/mean_z)

MPD_nls <- d2 %>%
  group_by(ration) %>%
  summarise(mean_av_perc_err = 100*mean(mean_err),
            min_av_perc_err = 100*mean(min_err),
            max_av_perc_err = 100*mean(max_err))

# MPE for simple linear regressions
d3 <- d1%>% filter(time>1) %>%
  rename(time2=time)

d3 <- fo_pred2 %>% 
  left_join(d3, by = c("ration", "time2")) %>%
  mutate(mean_err = abs((mean_z-mean_fit)/mean_z),
         min_err = abs((min_z-min_fit)/min_z),
         max_err = abs((max_z-max_fit)/max_z))

d3 <- d3 %>%
  mutate(perc_mean_err=100*mean_err)

MPD_linear <- d3 %>%
  group_by(ration) %>%
  summarise(mean_av_perc_err = 100*mean(mean_err),
            min_av_perc_err = 100*mean(min_err),
            max_av_perc_err = 100*mean(max_err))


# 7. Potential additional figures ----
DFI <- cawdata %>%
  select(tank,ID,ration,time,z, z1, daily_FI) %>%
  filter(!is.na(z)) %>%
  filter(time > 123) %>%
  mutate(timepoint = as.factor(time))

DFI$ration <- factor(DFI$ration,
                         levels = c("60S", "80S", "100S"))

ggplot(DFI, aes(timepoint, daily_FI, fill = ration)) +
  geom_violin() +
  facet_wrap(.~ration) +
  theme_classic() +
  guides(fill = "none")

ggplot(DFI, aes(timepoint, daily_FI, colour = z)) +
  geom_point() +
  facet_wrap(.~ration) +
  theme_classic() +
  guides(fill = "none") +
  labs(x="Sampling day",
       y="Individual daily feed intake (g)")
ggsave("DFI.jpeg")

ggplot(DFI, aes(z, daily_FI, colour = timepoint)) +
  geom_point(alpha = 0.5) +
  facet_wrap(ration~.) +
  theme_classic() +
  guides(fill = "none") +
  labs(x="Body mass (g)",
      y="Individual daily feed intake (g)") +
  theme(legend.position = "bottom")
ggsave("DFI2.jpeg")

DFI2 <- DFI %>%
  group_by(ration, time) %>%
  summarise(mean_DFI = mean(daily_FI, na.rm = TRUE),
            sd_DFI = sd(daily_FI, na.rm = TRUE))

cawdata$ration <- factor(cawdata$ration,
                     levels = c("60S", "80S", "100S"))
cawdata <- cawdata %>%
  filter(!is.na(daily_FI))

ggplot(cawdata, aes(FCR, daily_FI, colour = z)) +
  geom_point() +
  scale_colour_gradient(low = "yellow", high = "blue") +
  facet_wrap(.~ration) +
  theme_classic() +
  labs(x = "Feed Conversion Ratio (FCR)",
       y = "Individual daily feed intake (g)",
       colour = "Body mass (g)") +
  theme(legend.position = "bottom")

# 8. Cohort proportion above mean weight ----
## 8.1. Last time point ----
prop <- fo_pred %>% filter(time == 221) %>%
  select(ration, time, z1, .fitted)

prop %>% filter(ration == "60S") %>%
  summarise(min = min(z1), mean = mean(z1), max = max(z1))

prop %>% filter(ration == "80S") %>%
  summarise(min = min(z1), mean = mean(z1), max = max(z1))

prop %>% filter(ration == "100S") %>%
  summarise(min = min(z1), mean = mean(z1), max = max(z1))

z1_60S <-  prop %>% filter(ration =="60S") %>% pull(z1)
z1_80S <- prop %>% filter(ration =="80S") %>% pull(z1)
z1_100S <- prop %>% filter(ration =="100S") %>% pull(z1)

breaks = seq(200, 1500, by = 50)

z1_60S_cut = cut(z1_60S, breaks, right = FALSE)
freq_60S = table(z1_60S_cut)

z1_80S_cut = cut(z1_80S, breaks, right = FALSE)
freq_80S = table(z1_80S_cut)

z1_100S_cut = cut(z1_100S, breaks, right = FALSE)
freq_100S = table(z1_100S_cut)


## 8.2. Fish who remain under mid-way mean weight by last time point ----
prop2 <- fo_pred %>% filter(time == 124) %>%
  select(ration, time, z1, .fitted) %>%
  summarise(mean = mean(z1), min = min(z1), max = max(z1))

z1_60S <-  prop %>% filter(ration =="60S") %>% filter(time == 221)  #309 individuals
z1_80S <- prop %>% filter(ration =="80S") %>% filter(time == 221)  #303 individuals
z1_100S <- prop %>% filter(ration =="100S") %>% filter(time == 221)  #256 individuals

z1_60S_dwarfs <- z1_60S %>% count(z1 < 375.6144)
z1_80S_dwarfs <- z1_80S %>% count(z1 < 463.7237)
z1_100S_dwarfs <- z1_100S %>% count(z1 < 519.0861)


