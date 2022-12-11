set.seed(1)
data_df_new<-data_df %>% 
  group_by(species) %>% 
  mutate(random=rbernoulli( n())) %>% 
  mutate(random=case_when(random==1~"in-sample",
                          TRUE~"out-of-sample")) %>% 
  ungroup()

p_func_new<-ggplot(data_df_new %>% filter(species==sp_vis))+
  geom_point(aes(x=temp, y=fl, col=random), alpha=0.3)+
  geom_smooth(aes(x=temp, y=fl, group=random, col=random), method="lm")+
  theme_classic()+
  # facet_wrap(.~species)+
  theme(legend.position="bottom")+
  labs(color='Group') +
  xlab("Spring mean temperature (°C)")+
  ylab("Flowering time (day of year)")
p_func_new

library(parallel)
library(doSNOW)
cl <- makeCluster(20, outfile = "")
registerDoSNOW(cl)
data_df_new_list<-
  foreach (sp = sp_list,
           .packages = c("tidyverse", "spBayes", "gstat", "raster")) %dopar% {
             
             set.seed(42)
             data_sp<-data_df_new %>% 
               filter(species==sp) %>% 
               drop_na() %>% 
               dplyr::select(x=lon, y=lat, temp, fl, period,  species, year, random)
             
             data_train <-data_sp %>% 
               filter(random=="in-sample")
             
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
data_mis_new<-bind_rows(data_df_new_list) %>% 
  mutate(resid=fl-predict) %>% 
  mutate(temp_bin = cut(temp, breaks=5)) %>% 
  mutate(lat_bin = cut(y, breaks=5))
write_rds(data_mis_new, "./empirical/herbarium/data/mismatch_new.rds")

data_mis_new<-read_rds("./empirical/herbarium/data/mismatch_new.rds")

data_mis_new %>% 
  group_by(species, random) %>% 
  summarise (R2=cor(predict, fl)^2,
             rmse=Metrics::rmse(predict, fl)) %>% 
  gather(key="stat", value="value", -species, -random) %>% 
  spread(key="random", value="value") %>% 
  mutate(diff=`out-of-sample`-`in-sample`) %>% 
  gather(key="report", value="value", -species, -stat) %>% 
  group_by(stat, report ) %>% 
  summarise(median=quantile(value, 0.5),
            lower=quantile(value, 0.025),
            upper=quantile(value, 0.975))

p_pred_new<-ggplot(data_mis_new  %>% filter(species==sp_vis))+
  geom_point(aes(x=predict, y=fl), alpha=0.5)+
  geom_abline(intercept = 0, slope=1, col="red")+
  geom_text(data=data_mis_new  %>% 
              filter(species==sp_vis) %>% 
              group_by(random) %>% 
              summarise(median=quantile(resid, 0.5),
                        lower=quantile(resid, 0.025),
                        upper=quantile(resid, 0.975)),
            aes(x=150, y=250, label=paste0("Ymis=",round(median,2),"(", round(lower,2),", ", round(upper,2),")")), col="blue")+
theme_classic()+
  facet_wrap(.~random, nrow=1)+
  guides(col="none")+
  xlab("Predicted flowering time (day of year)")+
  ylab("Observed flowering time (day of year)")
p_pred_new


data_mis_new %>% 
  filter(random=="out-of-sample") %>% 
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975))

t.test(data_mis_new %>% 
         filter(random=="out-of-sample") %>% 
         pull(resid),
       alternative = "two.sided")

t_df_new<-data_mis_new %>% 
  filter(random=="out-of-sample") %>% 
  group_by(species) %>%
  summarise(median=quantile(resid, 0.5),
            lower=quantile(resid, 0.025),
            upper=quantile(resid, 0.975),
            p=t.test(x=resid)$p.value,
            estimate=t.test(x=resid)$estimate,
            med_fl=median(fl)) 
t_df_new %>% 
  filter(p<0.05) %>% 
  group_by(median<0) %>% 
  summarise(n=n(),
            lower=min(median),
            upper=max(median))

t_df_new %>% 
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

p_mis_new<-ggplot()+
  geom_violin(data=data_mis_new %>% filter(random=="out-of-sample"),
              aes(x=reorder(species, desc(species)),y=resid, col=species), draw_quantiles = 0.5)+
  geom_violin(data=data_mis_new %>% filter(random=="out-of-sample") %>% left_join(t_df_new, by="species") %>% filter(p<0.05),
              aes(x=reorder(species, desc(species)),y=resid, fill=species), draw_quantiles = 0.5)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  guides(col="none", fill="none")+
  coord_flip()+
  ylab("Deviation of observed flowering time \n from predicted flowering time (day)")+
  xlab ("Species")+
  theme(axis.text.y =element_text(face="italic")) 
p_mis_new

data_mis_new %>% 
  filter(random=="out-of-sample") %>% 
  group_by(species) %>% 
  summarise(median=median(resid)) %>% 
  pull(median) %>% 
  t.test()

cairo_pdf("./empirical/herbarium/figure_new.pdf", width = 10, height = 10)
grid.arrange(annotate_figure(p_map, fig.lab = "a"),
             annotate_figure(p_func_new, fig.lab = "b"),
             annotate_figure(p_pred_new, fig.lab = "c"),
             annotate_figure(p_mis_new, fig.lab = "d"),
             layout_matrix = rbind(
               c(1, 3),
               c(2, 4)
             ),
             widths=2:3,
             heights=2:3
)
dev.off()
