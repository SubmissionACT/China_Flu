### packages and data
packages=c('tidyverse','knitr','lubridate')
lapply(packages, require, character.only=T) 

dfA <- read.csv("china_flu_season.csv")
dfA$date <- as.Date(dfA$date)

# logit transformation
logitTransform <- function(p, epsilon = 1e-5) { 
  p <- ifelse(p == 0, epsilon, ifelse(p == 1, 1-epsilon, p))
  log(p / (1-p))
}

dfA$fluP <- logitTransform(dfA$fluP)
dfA <- dfA %>% select(date,state,fluP,temp,ah,o3,pm2.5)

# centralization
siteminus=function(x){x=x-mean(x,na.rm=T)}   

statemean=function(x){x=mean(x,na.rm=T)}

dfA=dfA %>%
  pivot_longer(cols=(temp:pm2.5)) %>%
  group_by(state,name) %>%
  mutate(siteminus=siteminus(value)) %>%
  group_by(name) %>%
  mutate(statemean=statemean(value)) %>%
  mutate(value_centr=siteminus+statemean) %>%
  select(-c(value,siteminus,statemean)) %>%
  pivot_wider(names_from='name',values_from='value_centr')

### add NA lines to separate discontinuous data 
new.row1 <- data.frame(
  date=rep(as.Date(c('2020-10-01')),30),
  state=rep(unique(dfA$state),each=1),fluP=999.,
  temp=999.,ah=999.,o3=999.,pm2.5=999.)

new.row2 <- data.frame(
  date=rep(as.Date(c('2015-07-01')),30),
  state=rep(unique(dfA$state),each=1),fluP=999.,
  temp=999.,ah=999.,o3=999.,pm2.5=999.)

new.row3 <- data.frame(
  date=rep(as.Date(c('2016-07-01')),30),
  state=rep(unique(dfA$state),each=1),fluP=999.,
  temp=999.,ah=999.,o3=999.,pm2.5=999.)

new.row4 <- data.frame(
  date=rep(as.Date(c('2017-07-01')),30),
  state=rep(unique(dfA$state),each=1),fluP=999.,
  temp=999.,ah=999.,o3=999.,pm2.5=999.)

new.row5 <- data.frame(
  date=rep(as.Date(c('2018-07-01')),30),
  state=rep(unique(dfA$state),each=1),fluP=999.,
  temp=999.,ah=999.,o3=999.,pm2.5=999.)

new.row6 <- data.frame(
  date=rep(as.Date(c('2019-07-01')),30),
  state=rep(unique(dfA$state),each=1),fluP=999.,
  temp=999.,ah=999.,o3=999.,pm2.5=999.)

### overall analysis data (combined all the states and normalized)
fn_norm=function(x){
  x=as.numeric(x)
  xt=x
  xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)
  return(xt)}

chinaflu_norm=dfA %>%
  group_by(state) %>%
  mutate_at(vars(-date,-state),fn_norm) %>%
  bind_rows(new.row1) %>%  
  bind_rows(new.row2) %>%
  bind_rows(new.row3) %>%
  bind_rows(new.row4) %>%
  bind_rows(new.row5) %>%
  bind_rows(new.row6) %>%
  arrange(state,date) %>% 
  rownames_to_column("row")

write.csv(chinaflu_norm, "chinaflu_norm.csv")
