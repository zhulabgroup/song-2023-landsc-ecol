set.seed(1)
data_agg_df_new<-data_agg_df %>% 
  group_by(species) %>% 
  mutate(random=rbernoulli( n())) %>% 
  mutate(random=case_when(random==1~"in-sample",
                          TRUE~"out-of-sample")) %>% 
  ungroup()

p_func_new<-ggplot(data_agg_df_new%>% filter(species==sp_vis))+
  geom_point(aes(x=mat, y=bh, col=random), alpha=0.1)+
  geom_smooth(aes(x=mat, y=bh,group=random, col=random), method="lm")+
  theme_classic()+
  xlab("Mean annual temperature (°C)")+
  ylab("Nestling ringing time (day of year)")+
  labs(color='Group')
p_func_new


data_df_new_list<-
  foreach (sp = sp_list,
           .packages = c("tidyverse", "spBayes", "gstat", "raster")) %dopar% {
             set.seed(42)
             
             # spatial regression
             data_sp<-data_agg_df_new %>% 
               filter(species  == sp) %>% 
               drop_na() %>% 
               dplyr::select(x=hex.x, y=hex.y, mat, bh, period, species, mig, year, area, random)
             ### Fit nonspatial and spatial model
             
             # Prepare data frame
             data_train <-data_sp %>% 
               filter(random=="in-sample")
             
             # Preparing to run Gaussian spatial model; get model phi (spatial decay parameter) estimate
             z <- resid(lm(bh ~ mat, data_train)) # residuals of non-spatial model
             data_train$z <- z
             lm.vgm <- variogram(z ~ 1, data = data_train %>% `coordinates<-`(data_train[, c("x", "y")]))
             # plot(lm.vgm)
             lm.fit <- fit.variogram(lm.vgm, model = vgm("Exp"), fit.method = 6, fit.sills = F)
             phi <- lm.fit[2, ]$range # extract phi value from fitted variogram
             
             # Set priors: loose priors on beta and residual error variance (tausq), and on spatial variance parameter (sigma.sq), but very tight on Phi (spatial decay parameter).
             priors <- list("beta.Norm" = list(rep(0, 2), diag(100, 2)), "phi.Unif" = c(-log(0.05) / (phi * 100), -log(0.05) / (phi / 100)), "sigma.sq.IG" = c(2, 2), "tau.sq.IG" = c(2, 0.1)) # shape and scale for IG
             # Set starting and tuning values
             starting <- list("phi" = -log(0.05) / phi, "sigma.sq" = 50, "tau.sq" = 1)
             tuning <- list("phi" = (log(0.05) / (phi * 100) - log(0.05) / (phi / 100)) / 10, "sigma.sq" = 0.1, "tau.sq" = 0.1)
             
             # Knots for Gaussian models
             # knots = kmeans(data_train[, c("x", "y")] %>% as.matrix(), 10,iter.max=100)$centers
             knots = c(3,3)
             
             # Run Gaussian spatial model
             t1<-Sys.time()
             splm <- spLM(bh ~ mat, data = data_train, coords = data_train[, c("x", "y")] %>% as.matrix(), cov.model = "exponential", priors = priors, tuning = tuning, starting = starting, n.samples = 1000, n.report = 100
                          , knots=knots
             )
             t2<-Sys.time()
             t2-t1
             
             # splm <- spRecover(splm, get.beta = TRUE, get.w = TRUE, start = 1, n.report = 10)
             # library(ggmcmc)
             # splm_mcmc1 <- ggs(splm$p.beta.recover.samples)
             # ggs_traceplot(splm_mcmc1) + theme_classic() + theme(strip.background = element_rect(color = "white", size = 0))
             
             ### Prediction
             t1<-Sys.time()
             splm.pred <- spPredict(splm, pred.covars=cbind(1, data_sp[, "mat"]) %>% as.matrix(), pred.coords=data_sp[, c("x", "y")] %>% as.matrix(),
                                    start=501)
             t2<-Sys.time()
             t2-t1
             
             quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}
             y.hat <- apply(splm.pred$p.y.predictive.samples, 1, quant) %>% 
               t() %>% 
               as.data.frame() %>% 
               `colnames<-`(c("lower", "predict", "upper"))
             
             data_sp %>% 
               bind_cols(y.hat)
             
           }
data_mis_new<-bind_rows(data_df_new_list) %>% 
  mutate(resid=bh-predict) 

data_mis_new %>% 
  group_by(species, random) %>% 
  summarise (R2=cor(predict, bh)^2,
             rmse=Metrics::rmse(predict, bh)) %>% 
  gather(key="stat", value="value", -species, -random) %>% 
  spread(key="random", value="value") %>% 
  mutate(diff=`out-of-sample`-`in-sample`) %>% 
  gather(key="report", value="value", -species, -stat) %>% 
  group_by(stat, report ) %>% 
  summarise(median=quantile(value, 0.5),
            lower=quantile(value, 0.025),
            upper=quantile(value, 0.975))
# plot mismatch

p_pred<-ggplot(data_mis_new  %>% filter(species==sp_vis))+
  geom_point(aes(x=predict, y=bh), alpha=0.2)+
  # geom_errorbarh(aes(y=bh, xmin=lower, xmax=upper), alpha=0.5)+
  # geom_smooth(aes(x=predict, y=bh), method="lm")+
  geom_abline(intercept = 0, slope=1, col="red")+
  geom_text(data=data_mis_new  %>% 
              filter(species==sp_vis) %>% 
              group_by(random) %>% 
              summarise(median=quantile(resid, 0.5),
                        lower=quantile(resid, 0.025),
                        upper=quantile(resid, 0.975)),
            aes(x=182, y=220, label=paste0("Ymis=",round(median,2),"(", round(lower,2),", ", round(upper,2),")")), col="blue")+
  # ggpubr::stat_cor(
  #   aes(
  #     x=predict, y=bh,
  #     label = paste( ..rr.label..,..p.label.., sep = "*`,`~")
  #   ),
  #   p.accuracy = 0.05,
  #   label.x.npc = "left",
  #   label.y.npc = "top",
  #   show.legend = F,
  #   col = "blue"
  # ) +
theme_classic()+
  facet_wrap(.~random, nrow=1)+
  guides(col="none")+
  xlab("Predicted nestling ringing time (day of year)")+
  ylab("Observed nestling ringing time (day of year)")
p_pred

# test if mismatch is different from 0
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
            estimate=t.test(x=resid)$estimate) 

t_df_new %>% 
  filter(p<0.05) %>% 
  group_by(median<0) %>% 
  summarise(n=n(),
            lower=min(median),
            upper=max(median))

p_mis_new<-ggplot()+
  # geom_violin(data=data_mis %>% filter(period=="late"),
  #             aes(x=as.factor("all species"),y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_violin(data=data_mis_new %>% filter(random=="out-of-sample"),
              aes(x=reorder(species, desc(species)),y=resid, col=species), draw_quantiles = 0.5)+
  geom_violin(data=data_mis_new %>% filter(random=="out-of-sample") %>% left_join(t_df_new, by="species") %>% filter(p<0.05),
              aes(x=reorder(species, desc(species)),y=resid, fill=species), draw_quantiles = 0.5)+
  # geom_violin(data=data_mis %>% filter(period=="late"),
  #             aes(x=mig,y=resid), fill="grey", draw_quantiles = 0.5)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  guides(col="none", fill="none")+
  coord_flip()+
  ylab("Deviation of observed nestling ringing time \n from predicted nestling ringing time (day)")+
  xlab ("Species")+
  theme(axis.text.y =element_text(face="italic")) 
p_mis_new

library(nlme)
# lme.fit <- lme(resid ~ y, random = ~ 1 + y | species, data = data_mis %>% filter(period=="late"))
# summary(lme.fit)
lme.fit <- lme((resid) ~ 1, random = ~ 1 | species,
               # corr = corSpatial(form = ~x + y, type ="exponential", nugget = T),
               data = data_mis_new %>% filter(random=="out-of-sample") %>%
                 # left_join(t_df, by="species") %>%
                 # filter(estimate<0 & p<0.05) %>% 
                 mutate(x=jitter(x, amount=0.000001), y=jitter(y, amount=0.000001)) )
summary(lme.fit)

cairo_pdf("./empirical/bird/figure_new.pdf", width = 10, height = 10)
grid.arrange(annotate_figure(p_map, fig.lab = "a"),
             annotate_figure(p_func_new, fig.lab = "b"),
             annotate_figure(p_pred_new, fig.lab = "c"),
             annotate_figure(p_mis_new, fig.lab = "d"),
             # annotate_figure(p_area, fig.lab = "(e)"),
             layout_matrix = rbind(
               c(1, 4),
               c(2, 4),
               c(3, 4)
             ),
             widths=c(1,1),
             heights=c(3,2,2)
)
dev.off()
