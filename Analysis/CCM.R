################################## Load packages, data, and functions ###################################
# packages
packages=c('tidyverse','knitr','lubridate','dplyr','purrr','rEDM','metap',
           'doParallel','foreach','imputeTS', "kableExtra",
           'ggplot2', 'ggpubr', 'ggthemes', 'cowplot',
           'customLayout', 'patchwork', 'grid', 'gridExtra',
           'usmap', 'maps', 'metap', 'scales', 'ggridges',
           'ggforce', 'ggbeeswarm')
lapply(packages, require, character.only=T)

# parallel computing parameters
cores_all=detectCores()
cores=ifelse(cores_all<9,4,60)
core_type='PSOCK'

num_sample=100
num_surr=1000

## load data and functions
df <- read.csv("china_flu_season.csv")
df$date <- as.Date(df$date )

## logit transformation of fluP
dfA=df 
logitTransform <- function(p, epsilon = 1e-5) { 
  p <- ifelse(p == 0, epsilon, ifelse(p == 1, 1-epsilon, p))
  log(p / (1-p))
}

dfA$fluP <- logitTransform(dfA$fluP)

### data standardization --- dtrend, dseason, and normalizaiton
nomz <- function(x, normalization=T, dtrend=T,dseason=T, season_sd=T, sea=35,  dTtype="linear"){
  x <- as.numeric(x)
  xt <- x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt <- diff(xt)} else if (dtrend==T & dTtype=="linear"){
    lm.t <- lm(xt~c(1:length(xt)))
    xt <- xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt <- xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt <- xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt <- (xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
}

df_smapc=dfA %>% select(state,date,temp,ah,o3,pm2.5,fluP) %>% group_by(state) %>% 
  mutate_at(3:ncol(.),nomz) %>% ungroup() 

################### Calculate noise factors for surrogate data of environmental variables (the whole-year data) #######################

# whole-year environmental data
df <- read.csv("china_flu_full.csv")
df$date <- as.Date(df$date)

### data standardization --- dtrend, dseason, and normalizaiton
nomz <- function(x, normalization=T, dtrend=T,dseason=T, season_sd=T, sea=52,  dTtype="linear"){
  x <- as.numeric(x)
  xt <- x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt <- diff(xt)} else if (dtrend==T & dTtype=="linear"){
    lm.t <- lm(xt~c(1:length(xt)))
    xt <- xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt <- xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt <- xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt <- (xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
}

# data standardization
dfB=df  %>% select(state,date,temp,ah,o3,pm2.5)
dfB_smapc=dfB %>% group_by(state) %>%
  mutate_at(3:ncol(.),nomz) %>% ungroup()

#select variables
df_flu <- dfB_smapc %>%
  mutate_at(vars(temp,ah,o3,pm2.5), ~{attributes(.) <- NULL; .}) %>%
  group_by(state) %>%
  gather(key, value, -c(date, state))

# calculate spar values for independent variables
fn_state_spar=function(states,keys){
  
  splineres <- function(spar){
    res <- rep(0, length(x))
    for (i in 1:length(x)){
      mod <- smooth.spline(x[-i], y[-i], spar = spar)
      res[i] <- predict(mod, x[i])$y - y[i]
    }
    return(sum(res^2))
  }
  
  x=df_flu %>% ungroup() %>%
    filter(state==states & key==keys) %>%
    select(date) %>% pull() %>% as.numeric()
  
  y=df_flu %>% ungroup() %>%
    filter(state==states & key==keys) %>%
    select(value)  %>% pull() %>% as.numeric()
  
  spars <- seq(0, 1.5, by = 0.1)
  ss <- rep(0, length(spars))
  
  ss=foreach(i = 1:length(spars),
             .combine=rbind,
             .inorder=FALSE)  %dopar% {
               targetCol = paste("T", i, sep = "")
               ss[i] <- splineres(spars[i])
             }
  
  spar=spars[which.min(ss)]
  data.frame(state=states,plt=keys,spar)
}

plist=list(states=unique(df_flu$state),
           keys=c('o3','ah','temp',"pm2.5")) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
flu_spar=plist %>% pmap_df(fn_state_spar)
stopCluster(cl)

# calculate sd values (alpha) for independent variables (additive noise factor to produce surrogate data)
yearday_anom <- function(t,x,spars){
  doy <- as.numeric(strftime(t, format = "%j"))
  I_use <- which(!is.na(x))
  doy_sm <- rep(doy[I_use],3) + rep(c(-366,0,366),each=length(I_use))
  x_sm <- rep(x[I_use],3)
  xsp <- smooth.spline(doy_sm, y = x_sm, w = NULL,
                       spar = spars, cv = NA,
                       all.knots = TRUE,keep.data = TRUE, df.offset = 0)
  xbar <- data.frame(t=t,doy=doy) %>%
    left_join(data.frame(doy=xsp$x,xbar=xsp$y),by='doy') %>%
    select(xbar)
  out = data.frame(t=t,mean=xbar,anomaly=(x - xbar))
  names(out) <- c('date','mean','anomaly')
  return(out)
}

fn_anomaly_PNAS=function(states,plts){
  vec_t=df_flu %>%
    ungroup() %>%
    filter(state==states & key==plts) %>%
    select(date) %>% pull()
  
  vec_x=df_flu %>%
    ungroup() %>%
    filter(state==states & key==plts) %>%
    select(value) %>% pull()
  
  spars=flu_spar %>% filter(state==states & plt==plts) %>%
    select(spar) %>% pull()
  
  df_9=yearday_anom(vec_t,vec_x,spars)
  sd=sd(df_9$anomaly,na.rm=TRUE)
  data.frame(states,plt=plts,spar=spars,sd_PNAS=sd)
}

plist=list(states=unique(df_flu$state),
           plts=c('o3','ah','temp',"pm2.5")) %>% cross_df()

sd_data=plist %>% pmap_df(fn_anomaly_PNAS)

fn_surr_data=function(data,ST,plts, tp_value){
  df=data %>%
    filter(state==ST) %>%
    mutate(plt=lag(get(plts),-tp_value)) %>%
    filter(!(is.na(plt))) %>%
    select(date,"plt")
  
  alpha=sd_data %>% filter(states==ST & plt==plts) %>%
    select(sd_PNAS) %>% pull()
  
  set.seed(2019)
  surr_data=
    SurrogateData(unlist(df[,"plt"]), method = "random_shuffle",
                  num_surr = num_surr,
                  alpha=alpha) %>%
    as.data.frame()
  df=df %>% select(-"plt")
  surrA=bind_cols(df,surr_data)
}

################################### Determine optimal E for the system by each state ###########################################

jishu <- function(x){
  ifelse(x%%2 ==0,F,T)
}

make_pred_nozeroL <- function(dat) {
  dat <- dat %>% 
    mutate(year = year(date), month = month(date), day = day(date), n = 1:nrow(dat))
  
  # Identify rows for May and October
  dat1 <- dat %>% filter(month == 5) %>% group_by(year) %>% filter(day == max(day))
  dat2 <- dat %>% filter(month == 10) %>% group_by(year) %>% filter(day == min(day))
  
  # Combine indices and sort them
  I_zero_strings <- c(dat1$n, dat2$n)
  I_zero_strings <- I_zero_strings[order(I_zero_strings)]
  
  # Add first row if necessary
  if (I_zero_strings[1] != 1 | jishu(length(I_zero_strings)) == T) {
    I_zero_strings <- c(1, I_zero_strings)
  }
  
  # Remove unintended '1' duplication
  if (length(I_zero_strings) > 1 && I_zero_strings[2] == 1) {
    I_zero_strings <- I_zero_strings[-2]
  }
  
  # Add last row if necessary
  if (I_zero_strings[length(I_zero_strings)] != nrow(dat)) {
    I_zero_strings <- c(I_zero_strings, nrow(dat))
  }
  
  # Construct the output matrix
  lib_out <- matrix(I_zero_strings, ncol = 2, byrow = TRUE)
  return(lib_out)
}

fn_E_smapc=function(data,ST,y){
  
  df=data %>%
    filter(state==ST) %>%
    select(date,y) %>%
    na.omit()
  
  M <- nrow(df)
  
  lib <- make_pred_nozeroL(df)
  pred <- make_pred_nozeroL(df)

  E=EmbedDimension(dataFrame=df,
                   lib=lib,
                   pred=pred,
                   columns=y, target=y,
                   maxE=20,
                   showPlot=F)
  temp=data.frame(dis=y,E,ST)
  temp}

plist=list(data=list(df_smapc),y=c('fluP'),ST=unique(dfA$state)) %>%
  cross_df()

E_smapc_out=plist %>% pmap_df(fn_E_smapc)

E_smapc =E_smapc_out %>% filter(E %in% 2:6) %>%
  group_by(dis,ST) %>%
  filter(rho==max(rho))  %>%
  as.data.frame() 

################################### Determine optimal theta for S-map by each state ###########################################

fn_theta_justY=function(data, ST, dis, theta){
  
  E=E_smapc[E_smapc[,'ST']==ST & E_smapc[,'dis']==dis,'E']
  
  df=data %>%
    filter(state==ST) %>%
    select(date,dis)
  
  M <- nrow(df)
  
  lib <- make_pred_nozeroL(df)
  pred <- make_pred_nozeroL(df)
  
  rho_theta = PredictNonlinear(dataFrame = df,
                               embedded = FALSE,
                               columns = dis,
                               target = dis,
                               Tp=1,
                               theta=theta,
                               lib=lib,
                               pred=pred,
                               showPlot = FALSE,
                               E = E)
  
  best_theta_df=rho_theta %>%
    mutate(dis=dis, state=ST, theta=theta)
}

plist=list(data=list(df_smapc),ST=unique(dfA$state),
           dis=c('fluP'),
           theta=c(0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4,
                   5, 6, 7, 8, 9)) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
theta_out=foreach(i = 1:nrow(plist),
                  .packages = c("rEDM","tidyverse",'lubridate'),
                  .combine=rbind,
                  .inorder=FALSE)  %dopar% {
                    theta_out=plist[i,] %>% pmap_df(fn_theta_justY)
                  }
stopCluster(cl)

best_theta=theta_out %>%
  group_by(state,dis) %>%
  filter(rho==max(rho)) %>%
  as.data.frame()

############################################# surrogate CCM #####################################################
fn_season_ccm=function(data,ST,x,y,tp_value){
  
  df=data %>%
    filter(state==ST) %>%
    select(date,y,x)
  
  E=E_smapc[E_smapc[,'ST']==ST & E_smapc[,'dis']==y,'E']
  
  alpha=sd_data %>% filter(states==ST & plt==x) %>%
    select(sd_PNAS) %>% pull()
  
  surr_data <- fn_surr_data(dfB_smapc,ST,x,0)
  
  all_data <- df %>% left_join(surr_data,by="date")
  
  names(all_data) = c("date", y, 'T1',paste0("T", 2:(num_surr+1)))
  
  m=nrow(all_data) %>% as.data.frame()
  
  libSize =c(E+2,m-E-2)
  
  rho_surr <- NULL
  
  for (i in 1:(num_surr+1)) {
    targetCol = paste("T", i, sep = "")
    ccm_out = CCM(dataFrame = all_data, E = E, Tp = tp_value,
                  columns = y,
                  target = targetCol,
                  libSizes = libSize,
                  random=T,
                  sample = num_sample,
                  seed=2019)
    col = paste(y, ":", targetCol, sep = "")
    dat=ccm_out %>% select(LibSize,col)
    names(dat)=c("lib","rho")
    test1=mutate(dat,i=i,dis=y,plt=x,E=E,
                 tp_value=tp_value,state=ST)
    rho_surr <- rbind(rho_surr,test1)
    
  }
  rho_surr
}

plist=list(data=list(df_smapc),
           ST=unique(df_smapc$state),
           y=c('fluP'),
           x=c('o3',"temp","ah","pm2.5"),
           tp_value=-2:0) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
ccm_out=foreach(j = 1:nrow(plist),
                .packages = c("rEDM","tidyverse"),
                .combine=rbind,
                .inorder=FALSE)  %dopar% {
                  ccm_out=plist[j,] %>% pmap_df(fn_season_ccm)
                }
stopCluster(cl)

# Calculate the difference in cross-mapping skills obtained by the maximum and the minimum library to test convergence property
dat_min <- ccm_out %>% filter(lib<50)
dat_min <- dat_min[order(dat_min$state,dat_min$tp_value,dat_min$plt,dat_min$dis,dat_min$i),]
dat_max <- ccm_out %>% filter(lib>50)
dat_max <- dat_max[order(dat_max$state,dat_max$tp_value,dat_max$plt,dat_max$dis,dat_max$i),]

dat <- cbind(dat_max,dat_min[,"rho"])
names(dat) <- c("lib","rho_max","i","dis","plt","E","tp_value","ST","rho_min")

dat1 <- dat %>% mutate(rho=rho_max-rho_min)
ccm_out <- dat1

ccm_out_raw = ccm_out  %>%
  filter(i==1) %>%
  select(dis,plt, ST,tp_value, rho)

# Calculate the P value of significance test by comparing the original CCM skill against the null distribution of surrogate ones.
ccm_p=ccm_out %>%
  group_by(dis,plt,E,ST,tp_value) %>%
  summarise(p=1-ecdf(rho[i != 1])(rho[i == 1])) %>%
  left_join(ccm_out_raw, by=c("dis","plt", "ST", "tp_value")) %>%
  rename(rho_raw=rho) %>%
  arrange(dis,plt)

# Adjust extreme P values before meta-significance test
# If P is extremely small approximating 0, then P is deemed as 0.005 allowing for Fisher's meta-significance test
# If original CCM skill is <0, then P is deemed as 1, that is accepting the null hypothesis exactly.

ccm_p=ccm_p %>%
  mutate(p=ifelse(p==0, 0.005, p)) %>%
  mutate(p=ifelse(rho_raw<=0, 1, p))

# Calculate meta-significance test using Fisher's method.
fn_metap=function(var,tp_values,plts,diss){
  df=filter(var,tp_value==tp_values & plt==plts & dis==diss)
  out=allmetap(df$p, method = c("sumlog")) %>% as.data.frame()
  mutate(out,tp_value=tp_values,plt=plts,dis=diss)
}

plist=list(var=list(ccm_p),
           tp_values=-2:0,
           plts=c("temp","ah",'o3',"pm2.5"),
           diss=c('fluP')) %>% cross_df()

meta_ccm_p_out=plist %>% pmap_df(fn_metap)

meta_ccm_p_out

################################################ Effect strength #################################################
fn_smapc=function(data,ST,plt,dis,tp_value){
  
  E=E_smapc[E_smapc$ST==ST & E_smapc$dis==dis,"E"]
  
  df=data %>% filter(state==ST) %>%
    mutate(plt_tp=lag(.data[[plt]],-tp_value)) %>%
    filter(!(is.na(plt_tp))) %>%
    select(date,all_of(dis),plt_tp) %>%
    na.omit()
  
  M <- nrow(df)
  
  embed_1=Embed(dataFrame = df, E = E, tau = -1, columns = dis )
  
  dataFrame = cbind(df[E:M, 'date'],df[E:M, dis],
                    embed_1[E:M, 1:(E-1)], df[E:M, 'plt_tp']  ) %>%
    as.data.frame()
  
  names(dataFrame)=c('date',dis,letters[1:(E-1)],plt)
  
  columns = paste(paste(letters[1:(E-1)],collapse =' '), plt ,
                  sep=' ')
  
  m <- nrow(dataFrame)
  
  lib <- make_pred_nozeroL(dataFrame)
  pred <- make_pred_nozeroL(dataFrame)
  
  theta = best_theta[best_theta$state==ST & best_theta$dis==dis, "theta"]
  
  smap = SMap(dataFrame = dataFrame,
              embedded = TRUE,
              columns = columns,
              target = dis,
              lib = lib,
              pred=pred,
              theta=theta,
              Tp = 1,
              E = E)
  smapc_df=smap$coefficients[c(1,2+E)]
  
  names(smapc_df)=c('date','effect')
  smapc_df=smapc_df %>%
    mutate(date=lubridate::as_date(date,
                                   origin = lubridate::origin)) %>%
    mutate(dis=dis,ST=ST, plt=plt,E=E,tp_value=tp_value)
}

plist=list(data=list(df_smapc), 
           ST=unique(df_smapc$state),
           dis=c('fluP'),
           plt=c('o3','ah','temp',"pm2.5"),
           tp_value=-2:0) %>% cross_df()


cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
C_out=foreach(i = 1:nrow(plist),
              .packages = c("rEDM","tidyverse","lubridate"),
              .combine=rbind,
              .export='best_theta',
              .inorder=FALSE)  %dopar% {
                C_out=plist[i,] %>% pmap_df(fn_smapc)
              }

stopCluster(cl)

# ############################### Plotting S-map Effect Size ###############################
# Reshape the C_out output: effect strength estimates
SEeffect<- C_out %>%
  select(date,ST,tp_value,dis,plt,effect) %>%
  mutate(tp_value=factor(tp_value,
                         levels=c(0, -1, -2)))

SEeffect$tp_value <- factor(SEeffect$tp_value, levels=c(0, -1, -2), labels=c("Lag0","Lag1","Lag2"))

# filter out the extreme values
SEeffect_ex <- SEeffect %>% group_by(ST, tp_value,plt,dis) %>%
  filter(effect < quantile(effect, probs=.95, na.rm = T),
         effect > quantile(effect, probs=.05, na.rm = T)) 


mean_effect <- SEeffect_ex %>% group_by(ST, tp_value,plt,dis) %>%
  summarise(medianE=median(effect)) %>% ungroup


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

ggplot(mean_effect, aes(x = plt, y = medianE, fill = as.factor(tp_value))) +
  geom_boxplot(position = position_dodge(0.5),alpha=1, outlier.shape = NA,width = 0.5) +  
  scale_fill_manual(name = "Lag", 
                    values = c("gray30", "gray50", "gray70"), 
                    labels = c("Lag 0", "Lag 1", "Lag 2")) +     
  scale_x_discrete(limits = rev(c( "temp", "ah","pm2.5", "o3")), 
                   labels = rev(c( expression("T"), expression("AH"), expression(PM[2.5]),expression(O[3])))) +
  labs(x = "", y = 'Effect estimates (CCM)') +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  mytheme +
  theme(legend.position = "bottom")  


