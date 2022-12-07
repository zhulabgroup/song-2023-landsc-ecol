data_df<-data_df %>% 
  rowwise() %>% 
  mutate(random=sample(c("in", "out"), 1)) %>% 
  ungroup()

library(parallel)
library(doSNOW)
cl <- makeCluster(20, outfile = "")
registerDoSNOW(cl)
data_df_new_list<-
  foreach (sp = sp_list,
           .packages = c("tidyverse", "spBayes", "gstat", "raster")) %dopar% {
             
             data_sp<-data_df %>% 
               filter(species==sp) %>% 
               drop_na() %>% 
               dplyr::select(x=lon, y=lat, temp, fl, period,  species, year, random)
             
             data_train <-data_sp %>% 
               filter(random=="in")
             
             # Preparing to run Gaussian spatial model; get model phi (spatial decay parameter) estimate
             z <- resid(lm(fl ~ temp, data_train)) # residuals of non-spatial model
             data_train$z <- z
             lm.vgm <- variogram(z ~ 1, data = data_train %>% `coordinates<-`(data_train[, c("x", "y")]))
             # plot(lm.vgm)
             lm.fit <- fit.variogram(lm.vgm, model = vgm("Exp"), fit.method = 6, fit.sills = F)
             phi <- lm.fit[2, ]$range # extract phi value from fitted variogram
             
             # Set priors: loose priors on beta and residual error variance (tausq), and on spatial variance parameter (sigma.sq), but very tight on Phi (spatial decay parameter).
             priors <- list("beta.Norm" = list(rep(0, 2), diag(100, 2)), "phi.Unif" = c(-log(0.05) / (phi * 100), -log(0.05) / (phi / 100)), "sigma.sq.IG" = c(2, 2), "tau.sq.IG" = c(2, 0.1)) # shape and scale for IG
             # Set starting and tuning values
             starting <- list("phi" = -log(0.05) / phi, "sigma.sq" = 50, "tau.sq" = 1)
             tuning <- list("phi" = (log(0.05) / (phi * 100) - log(0.05) / (phi / 100)) / 10, "sigma.sq" = 0.01, "tau.sq" = 0.01)
             
             # Knots for Gaussian models
             # knots = kmeans(coords, 20,iter.max=100)$centers
             
             # Run Gaussian spatial model
             splm <- spLM(fl ~ temp, data = data_train, coords = data_train[, c("x", "y")] %>% as.matrix(), cov.model = "exponential", priors = priors, tuning = tuning, starting = starting, n.samples = 10000, n.report = 1000
                          # , knots=knots
             )
             
             ### Prediction
             splm.pred <- spPredict(splm, pred.covars=cbind(1, data_sp[, "temp"]) %>% as.matrix(), pred.coords=data_sp[, c("x", "y")] %>% as.matrix(),
                                    start=5001, thin = 10)
             quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}
             y.hat <- apply(splm.pred$p.y.predictive.samples, 1, quant) %>% 
               t() %>% 
               as.data.frame() %>% 
               `colnames<-`(c("lower", "predict", "upper"))
             
             data_sp %>% 
               bind_cols(y.hat)
           }
data_mis<-bind_rows(data_df_new_list) %>% 
  mutate(resid=fl-predict) %>% 
  mutate(temp_bin = cut(temp, breaks=5)) %>% 
  mutate(lat_bin = cut(y, breaks=5))




data_mis %>% 
  filter(random=="out") %>% 
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975))

test<-wilcox.test(data_mis %>% 
              filter(random=="out") %>% 
              pull(resid),
            alternative = "two.sided"
)
t.test(data_mis %>% 
         filter(random=="out") %>% 
         pull(resid),
       alternative = "two.sided")

t_df<-data_mis %>% 
  filter(random=="out") %>% 
  group_by(species) %>%
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975),
            p=wilcox.test(x=resid)$p.value,
            estimate=t.test(x=resid)$estimate,
            med_fl=median(fl)) 
t_df %>% 
  filter(p<0.05) %>% 
  group_by(median<0) %>% 
  summarise(n=n(),
            lower=min(median),
            upper=max(median))

t_df %>% 
  ggplot()+
  geom_point(aes(x=med_fl, y=estimate))+
  geom_smooth(aes(x=med_fl, y=estimate), method="lm")+
  ggpubr::stat_cor(
    aes(
      x=med_fl, y=estimate,
      label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
    ),
    p.accuracy = 0.05,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = F,
    col = "blue"
  ) +
  theme_classic()

# sp_list_order<-t_df %>% 
#   arrange(med_fl) %>% 
#   pull(species)
# data_mis<-data_mis %>% 
#   mutate(species=fct_relevel(species, levels=sp_list_order))

p_mis<-ggplot()+
  # geom_violin(data=data_mis %>% filter(period=="late"),
  #             aes(x=as.factor("all species"),y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(random=="out"),
              aes(x=reorder(species, desc(species)),y=resid, col=species), draw_quantiles = 0.5)+
  geom_violin(data=data_mis %>% filter(random=="out") %>% left_join(t_df, by="species") %>% filter(p<0.05),
              aes(x=reorder(species, desc(species)),y=resid, fill=species), draw_quantiles = 0.5)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  guides(col="none", fill="none")+
  coord_flip()+
  ylab("Deviation of observed flowering time \n from predicted flowering time (day)")+
  xlab ("Species")+
  theme(axis.text.y =element_text(face="italic")) 
p_mis

data_mis %>% filter(random=="out") %>% 
  group_by(species) %>% 
  summarise(median=median(resid)) %>% 
  pull(median) %>% 
  t.test()

library(nlme)
# lme.fit <- lme(resid ~ y, random = ~ 1 + y | species, data = data_mis %>% filter(period=="late"))
# summary(lme.fit)
lme.fit <- lme((resid) ~ 1, random = ~ 1 | species,
               # corr = corSpatial(form = ~x + y, type ="exponential", nugget = T),
               data = data_mis %>% filter(random=="out") %>%
                 # left_join(t_df, by="species") %>%
                 # filter(estimate<0 & p<0.05) %>% 
                 mutate(x=jitter(x, amount=0.000001), y=jitter(y, amount=0.000001)) )
summary(lme.fit)
