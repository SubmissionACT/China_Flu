################################## Load packages & data ###################################
library(tidyverse); library(lubridate); library(dlnm); library(splines);
library(tsModel); library(gnm);library(ggpubr);library(metafor);
library(mgcv); library(imputeTS)

df <- read.csv("china_flu_season.csv")
dfA <- df %>% mutate(fluP=ifelse(fluP==0, 1e-5, fluP), year = year(date), month = month(date))

################################## State-level GLM ###################################
glm_function <- function(dat, Y='fluP', x1='value', xlag=0, cityname, variable) {
  
  df <- dat %>%
    filter(state == cityname) %>%
    mutate(loglag1 = log(Lag(fluP, 1)))%>%
    mutate(loglag2 = log(Lag(fluP, 2))) %>% 
    na.omit()
  
  formula_str <- switch(variable,
                        o3   = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(', x1, ',', xlag+1, ')+ Lag(', x1, ',', xlag+2, ') + 
                                     Lag(ah,', xlag, ')+ Lag(temp,', xlag, ')+ Lag(pm2.5,', xlag, ') + 
                                     as.factor(year) + as.factor(month) + 
                                     loglag1+loglag2', sep=''),
                        pm2.5= paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(', x1, ',', xlag+1, ')+ Lag(', x1, ',', xlag+2, ') + 
                                     Lag(ah,', xlag, ')+ Lag(temp,', xlag, ')+ Lag(o3,', xlag, ') + 
                                     as.factor(year) + as.factor(month) + 
                                     loglag1+loglag2', sep=''),
                        ah   = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(', x1, ',', xlag+1, ')+ Lag(', x1, ',', xlag+2, ') + 
                                     Lag(temp,', xlag, ')+ Lag(pm2.5,', xlag, ')+ Lag(o3,', xlag, ') + 
                                     as.factor(year) + as.factor(month) + 
                                     loglag1+loglag2', sep=''),
                        temp = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(', x1, ',', xlag+1, ')+ Lag(', x1, ',', xlag+2, ') + 
                                     Lag(ah,', xlag, ')+ Lag(pm2.5,', xlag, ')+ Lag(o3,', xlag, ') + 
                                     as.factor(year) + as.factor(month) + 
                                     loglag1+loglag2', sep='')
  )
  
  fit <- glm(as.formula(formula_str), data=df, family=quasibinomial, na.action=na.omit)
  summ <- summary(fit)
  
  outA <- data.frame(
    beta = summ$coefficients[2,1],
    se   = summ$coefficients[2,2],
    t    = summ$coefficients[2,3],
    p    = summ$coefficients[2,4],
    state = cityname,
    lag   = xlag
  )
  
  outB <- data.frame(
    beta = summ$coefficients[3,1],
    se   = summ$coefficients[3,2],
    t    = summ$coefficients[3,3],
    p    = summ$coefficients[3,4],
    state = cityname,
    lag   = xlag+1
  )

  outC <- data.frame(
    beta = summ$coefficients[4,1],
    se   = summ$coefficients[4,2],
    t    = summ$coefficients[4,3],
    p    = summ$coefficients[4,4],
    state = cityname,
    lag   = xlag+2
  )

  out <- rbind(outA, outB, outC)
  
  return(out)
}

variables <- c("o3", "ah", "temp", "pm2.5")
results <- list()

for (var in variables) {
  plist <- list(
    dat      = list(dfA),
    Y        = 'fluP',
    x1       = var,
    xlag     = 0,
    cityname = unique(dfA$state),
    variable = var
  ) %>% cross_df()
  results[[var]] <- plist %>% 
    pmap_df(glm_function)
}

# Store the results for easy access
out_o3   <- results$o3
out_pm25 <- results$pm2.5
out_ah   <- results$ah
out_temp <- results$temp

########################### Meta analysis of GLM beta results ############################
fmeta=function(dat){
  out=dat
  
  GAMresCI  =  out %>%
    mutate(betalow=beta-1.96*se, 
           betahigh=beta+1.96*se) %>% 
    mutate(state=as.factor(state),
           lag=factor(lag, levels = c(0,1,2),
                      labels = c("Lag0","Lag1","Lag2"))) %>%
    select(state,lag,beta,betalow,betahigh,p)
  
  meta.ci.result  =  as.data.frame(matrix(rep(NA,18),nrow=3))
  names(meta.ci.result)  =  names(GAMresCI)
  
  for (i in 1:3) {
    out.lag  =  out %>%
      filter(lag==i-1)
    meta  =  rma(yi=beta, sei=se, slab=state, method="REML", data=out.lag, level=99.9)
    meta.re  =  with(meta, c(b, ci.lb, ci.ub,pval))
    meta.ci.result[i,1] = "All"
    meta.ci.result[i,2] = paste('Lag',i-1)
    meta.ci.result[i,3] = meta.re[1]
    meta.ci.result[i,4] = meta.re[2]
    meta.ci.result[i,5] = meta.re[3]
    meta.ci.result[i,6] = meta.re[4]
  }
  
  meta.ci.result  =  meta.ci.result %>% as.data.frame() 
  abc_levels=str_sort(toupper(unique(GAMresCI$state)),decreasing=T)
  GAMresCI_all  =  rbind(GAMresCI, meta.ci.result) %>%
    mutate(type=factor(ifelse(state=="All",2,1))) %>%
    mutate(state=factor(toupper(state),
                        levels=c(abc_levels,"ALL"))) %>%
    as.tibble()
}

gam_o3=fmeta(out_o3)
gam_pm25=fmeta(out_pm25)
gam_ah=fmeta(out_ah)
gam_temp=fmeta(out_temp)

# Calculate overall and state-specific SD of each environmental predictor
df_SD_all=dfA %>%
  dplyr::summarize(o3=sd(o3),pm2.5=sd(pm2.5), ah=sd(ah), temp=sd(temp)) %>%
  gather(plt, SD) %>% 
  mutate(state='ALL')

df_SD_ST <- dfA %>% group_by(state) %>%
  dplyr::summarize(o3=sd(o3), pm2.5=sd(pm2.5), ah=sd(ah), temp=sd(temp)) %>%
  gather(plt, SD, -state) %>%
  mutate(state=toupper(state))

df_SD <- rbind(df_SD_ST, df_SD_all%>%select(state,plt,SD))

# Transform beta to effect size of each SD change
df1=gam_o3  %>% as.data.frame() %>% mutate(plt='o3')
df2=gam_ah  %>% as.data.frame() %>% mutate(plt='ah')
df3=gam_temp %>% as.data.frame() %>% mutate(plt='temp')
df4=gam_pm25 %>% as.data.frame() %>% mutate(plt='pm2.5')

df_gam_final=rbind(df1,df2,df3,df4) %>%
  left_join(df_SD,by=c("state",'plt')) %>%
  mutate(Size=beta*SD,SizeL=betalow*SD,SizeH=betahigh*SD) %>%
  select(state, plt, lag, Size, SizeL, SizeH, gam.p=p) %>%
  mutate(sig=ifelse(gam.p<0.05 & state!="ALL", 'sig', ifelse(gam.p<0.001 & state=="ALL", 'sig','non_sig' )),
         sig=factor(sig, levels=c("non_sig","sig"))) 

########################################## Plot ##########################################
df_gam_final$state <- ifelse(tolower(df_gam_final$state) == "all", "Overall",
                             paste0(toupper(substring(df_gam_final$state, 1, 1)),
                                    tolower(substring(df_gam_final$state, 2))))
df_gam_final[df_gam_final$state=="Inner mongolia", "state"] <- "Inner Mongolia"
df_gam_final$plt <- factor(df_gam_final$plt, levels=c( "o3","pm2.5","ah","temp"))

# Plotting tate-specific GLM results
abc_levels=str_sort(unique(df_gam_final$state)[1:30],decreasing=F)

dat_glm_ST <- df_gam_final %>% 
  filter(state!='Overall' & lag=='Lag0') %>% 
  mutate(state=factor(state, levels= abc_levels)) 

dat_glm_ST$plt <- factor(dat_glm_ST$plt, levels=c("o3","pm2.5","ah","temp"),
                         labels = c(expression(O[3]), expression(PM[2.5]), expression("AH"), expression("T")))

mytheme <- theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = NA, colour ="grey90" ),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 15, color = "black", face='bold', family='serif'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))+
  theme(
    legend.position = "None",
    legend.title=element_text(size=18,family='serif'),
    legend.key.width= unit(1.1, 'cm'),
    legend.text = element_text(size=14, color = "black",family='serif'),
    legend.spacing.y = unit(0.1, 'cm'),
    legend.background = element_rect(color = NA),
    legend.box.margin = margin(0.1,0.1,0.1,0.1,"cm")) +
  theme(axis.title.x  = element_text(size=15,family='serif'),
        axis.title.y  = element_text(size=20,family='serif'),
        axis.text.x = element_text(color="black", size=14,family='serif'),
        axis.text.y = element_text(color="black", size=14,family='serif'),
        plot.title = element_text(size = 20, hjust=0.5, family='serif'))


p_glm_ST_lag= ggplot(dat_glm_ST) + 
  geom_errorbar(aes(ymin=SizeL, ymax=SizeH, x=fct_rev(state)), color='gray',
                position = position_dodge(0.8), width=0, size=0.8) +
  geom_point(aes(y=Size, x=fct_rev(state), shape=sig),
             color='#B4020A', size=3, stroke = 0.5,
             position = position_dodge(0.8))+
  facet_wrap(~plt, ncol = 4, labeller = "label_parsed") +
  scale_shape_manual(name="Statistical significance test:", values = c(1,16),
                     labels=c("Non-significant","Significant")) +
  labs(x='',y='') +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    scale_y_continuous(limits = c(-1.0, 1.0), 
                       breaks = seq(-0.5, 0.5, by = 0.5)) + 
  coord_flip() +
  mytheme

p_glm_ST_lag

# Plotting meta-analyzed GLM results
dat_glm_ALL <- df_gam_final %>% 
  filter(state=='Overall' )

mytheme <- theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = NA, colour ="grey90" ),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 14, color = "gray10",
                                    face='bold',family='serif'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(1.1, 'cm'),
    legend.text = element_text(size = 12, color = "black", family = 'serif'),
    legend.spacing.y = unit(0.1, 'cm'),
    legend.background = element_rect(color = NA),
    legend.box.margin = margin(0, 0, 0, 0, "cm"),  
    legend.margin = margin(0, 0, 0, 0, "cm"),   
    plot.margin = margin(0.5, 0.5, 0.1, 0.5, "cm")
  )+
  theme(axis.title.x  = element_text(size=16,family='serif'),
        axis.title.y  = element_text(size=16,family='serif'),
        axis.text.x = element_text(color="black", size=16,family='serif'),
        axis.text.y = element_text(color="black", size=16,family='serif'),
        plot.title = element_text(size = 20, hjust=0.5, family='serif'))

p_glm_vs_1 <- ggplot(dat_glm_ALL, aes(x = fct_rev(plt), group = interaction(plt, lag))) +
  geom_errorbar(aes(ymin = SizeL, ymax = SizeH), 
                color = "black", 
                position = position_dodge(0.6), width = 0, linewidth = 0.6, alpha = 0.5) +
  geom_point(aes(y = Size, shape = sig, color = as.factor(lag)), 
             size = 4, stroke = 0.8, position = position_dodge(0.6)) +
  scale_shape_manual(guide = "none", values = c(1, 16)) +
  scale_color_manual(name = "Lag", values = c("gray30", "gray50", "gray70")) + 
  labs(x = '', y = expression(paste( 'Effect estimates (GLM)'))) + 
  scale_x_discrete(limits = rev(c( "temp", "ah","pm2.5", "o3")), 
                   labels = rev(c( expression("T"), expression("AH"),expression(PM[2.5]), expression(O[3])))) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  mytheme 

p_glm_vs_1

