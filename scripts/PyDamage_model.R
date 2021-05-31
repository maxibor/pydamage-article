# Required packages and custom functions
require(rstanarm)
require(sjPlot)
require(sjlabelled)
require(sjmisc)
require(ggplot2)
require(tidyverse)
require(cowplot)
require(betareg)
require(lmerTest)
require(rms)
require(data.table)

facetFix <- function(x.gg){
  gp <- ggplotGrob(x.gg)
  facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
  x.var <- sapply(ggplot_build(x.gg)$layout$panel_scales_x,
                  function(l) length(l$range$range))
  gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
  return(grid::grid.draw(gp))
}

rescale <- function(x){
  return((x-mean(x)/sd(x)))
}

# Read in data and plot

dat <- read.table('../data/PYD_simulation_results.tsv.gz',header=T,sep='\t') %>%
  janitor::clean_names() %>%
  as_tibble()

dat %>%
  dplyr::mutate(damage=factor(damage)) %>%
  dplyr::group_by(genome,damage,readlength) %>%
  dplyr::summarise(ct0=mean(cto_t_0,na.rm=T),ct1=mean(cto_t_1,na.rm=T),ct2=mean(cto_t_2,na.rm=T),ct3=mean(cto_t_3,na.rm=T),ct4=mean(cto_t_4,na.rm=T),
                   ga0=mean(gto_a_0,na.rm=T),ga1=mean(gto_a_1,na.rm=T),ga2=mean(gto_a_2,na.rm=T),ga3=mean(gto_a_3,na.rm=T),ga4=mean(gto_a_4,na.rm=T)) %>%
  dplyr::filter(!is.nan(ct0)) %>%
  gather(mutation,prop,ct0:ga4) %>%
  dplyr::mutate(position=as.numeric(substr(mutation,start=3,stop=3)),mutation=toupper(substr(mutation,start=1,stop=2))) %>%
  ggplot(aes(x=position,y=prop,colour=damage))+
  theme_bw()+
  geom_point(aes(group=readlength,shape=readlength))+
  # geom_smooth(aes(linetype=damage),method="glm",formula=y~x,se=F,
  #            method.args = list(family = gaussian(link = 'log')))+
  geom_smooth(se=F)+
  facet_grid(genome~.)+
  theme(legend.position='bottom')

# Model fitting
vs <- c('damage','contiglength','g_ccontent','median_rl','actual_cov')
dat.c.raw <- dat %>%
  dplyr::select(genome,damage,contiglength,
                g_ccontent,median_rl,actual_cov,qvalue) %>%
  dplyr::filter(!is.na(qvalue)) %>%
  dplyr::mutate(sig=qvalue<0.05,actual_cov>=0) %>%
  dplyr::mutate_at(.vars=vs[-1],rescale)

prop.table(table(ifelse(dat.c.raw$damage>0,'damaged','non-damaged'),
                 ifelse(dat.c.raw$sig,'called damaged','not called damaged')),margin=1)*100

caret::confusionMatrix(table(ifelse(dat.c.raw$damage>0,'damaged','non-damaged'),
                 ifelse(dat.c.raw$sig,'damaged','non-damaged')))
dat.fit <- sort(sample(nrow(dat.c.raw),floor(0.5*nrow(dat.c.raw))))
dat.cv <- setdiff(1:nrow(dat.c.raw),dat.fit)

dat.c <- dat.c.raw %>%
  dplyr::slice(dat.fit)
  
dat.d <- dat.c.raw %>%
  dplyr::slice(dat.cv)

n <- length(vs)
reps <- 10
repV <- matrix(0,ncol=reps,nrow=6)
combs <- as.matrix(expand.grid(as.data.frame(matrix(rep(c(0,1),times=n),nrow=2))))
trueDamage <- ifelse(dat.c$damage>0,'1','0')

for(i in 2:nrow(combs)){
  
  if(i==2){
    print('Start')
    model.tib <- tibble(vars=rep('x',nrow(combs)),acc=0,sensitivity=0,specificity=0,R2=0,F1=0,AUC=0,
                        acc.sd=0,sensitivity.sd=0,specificity.sd=0,R2.sd=0,F1.sd=0,AUC.sd=0)
    modelList <- list()
    model.tib$vars[1] <- 'null'
    
    ### Subsample data -
    for(j in 1:reps){
      #### sub-sample dat.c ----
      dat.c_j <- dat.c %>%
        dplyr::group_by(sig) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(n=min(n)) %>%
        dplyr::group_by(sig) %>%
        dplyr::sample_n(n,replace=T) %>%
        dplyr::ungroup()
      #### sub-sample dat.d ----
      dat.d_j <- dat.d %>%
        dplyr::group_by(sig) %>%
        dplyr::mutate(n=n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(n=min(n)) %>%
        dplyr::group_by(sig) %>%
        dplyr::sample_n(n,replace=T) %>%
        dplyr::ungroup()
      
      #### Calculate stats ----
      glm_j <- glm(sig ~ 1,family='binomial',data=dat.c_j,maxit=100,epsilon=1e-14)
      curr.fit_j <- table(factor(predict(glm_j,type='response')>0.5,levels=c(TRUE,FALSE)),
                          factor(dat.c_j$damage>0,levels=c(TRUE,FALSE)))
      predictedDamage_j <- factor(as.vector(ifelse(predict(glm_j,type='response')>=0.5,'1','0')),levels=c('0','1'))
      currPrediction_j <- ROCR::prediction(predict(glm_j,type='response'),ifelse(dat.c_j$damage>0,'1','0'))
      curr.sens <- caret::sensitivity(curr.fit_j)
      curr.spec <- caret::specificity(curr.fit_j)
      
      repV[2,j] <- caret::sensitivity(curr.fit_j)
      repV[3,j] <- caret::specificity(curr.fit_j)
      repV[1,j] <- mean(c(repV[2,j],repV[3,j]))
      repV[4,j] <- rcompanion::nagelkerke(glm_j)$Pseudo.R.squared.for.model.vs.null[3]
      repV[5,j] <- MLmetrics::F1_Score(y_pred=factor(ifelse(dat.c_j$damage>0,'1','0'),levels=c('0','1')),
                           y_true=factor(predictedDamage_j,levels=c('0','1')),
                           positive=1)
      repV[6,j] <- as.numeric((ROCR::performance(currPrediction_j,'auc'))@y.values)
      print(sprintf('%i of %i sub-samples.',j,reps))
    }
    print(sprintf('1 out of %i complete.',nrow(combs)))
    
    model.tib$acc[1] <- mean(repV[1,])
    model.tib$acc.sd[1] <- sd(repV[1,])
    model.tib$sensitivity[1] <- mean(repV[2,])
    model.tib$sensitivity.sd[1] <- sd(repV[2,])
    model.tib$specificity[1] <- mean(repV[3,])
    model.tib$specificity.sd[1] <- sd(repV[3,])
    model.tib$R2[1] <- mean(repV[4,])
    model.tib$R2.sd[1] <- sd(repV[4,])
    model.tib$F1[1] <- mean(repV[5,])
    model.tib$F1.sd[1] <- sd(repV[5,])
    model.tib$AUC[1] <- mean(repV[6,])
    model.tib$AUC.sd[1] <- sd(repV[6,])
    
    modelList[[1]] <- glm(sig ~ 1,family='binomial',data=dat.c,maxit=100,epsilon=1e-14)
  }
  currF <- formula(paste0('sig~',paste0(vs[which(combs[i,]==1)],collapse='+')))
  model.tib$vars[i] <- paste0(vs[which(combs[i,]==1)],collapse='\n')
  
  ### Subsample data -
  for(j in 1:reps){
    #### sub-sample dat.c ----
    dat.c_j <- dat.c %>%
      dplyr::group_by(sig) %>%
      dplyr::mutate(n=n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n=min(n)) %>%
      dplyr::group_by(sig) %>%
      dplyr::sample_n(n,replace=T) %>%
      dplyr::ungroup()
    #### sub-sample dat.d ----
    dat.d_j <- dat.d %>%
      dplyr::group_by(sig) %>%
      dplyr::mutate(n=n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n=min(n)) %>%
      dplyr::group_by(sig) %>%
      dplyr::sample_n(n,replace=T) %>%
      dplyr::ungroup()
    
    #### Calculate stats ----
    glm_j <- glm(currF,family='binomial',data=dat.c_j,maxit=100,epsilon=1e-14)
    curr.fit_j <- table(factor(predict(glm_j,type='response')>0.5,levels=c(TRUE,FALSE)),
                        factor(dat.c_j$damage>0,levels=c(TRUE,FALSE)))
    predictedDamage_j <- factor(as.vector(ifelse(predict(glm_j,type='response')>=0.5,'1','0')),levels=c('0','1'))
    currPrediction_j <- ROCR::prediction(predict(glm_j,type='response'),ifelse(dat.c_j$damage>0,'1','0'))
    curr.sens <- caret::sensitivity(curr.fit_j)
    curr.spec <- caret::specificity(curr.fit_j)
    
    repV[2,j] <- caret::sensitivity(curr.fit_j)
    repV[3,j] <- caret::specificity(curr.fit_j)
    repV[1,j] <- mean(c(repV[2,j],repV[3,j]))
    repV[4,j] <- rcompanion::nagelkerke(glm_j)$Pseudo.R.squared.for.model.vs.null[3]
    repV[5,j] <- MLmetrics::F1_Score(y_pred=factor(ifelse(dat.c_j$damage>0,'1','0'),levels=c('0','1')),
                                     y_true=factor(predictedDamage_j,levels=c('0','1')),
                                     positive=1)
    repV[6,j] <- as.numeric((ROCR::performance(currPrediction_j,'auc'))@y.values)
    print(sprintf('%i of %i sub-samples.',j,reps))
  }
  
  model.tib$acc[i] <- mean(repV[1,])
  model.tib$acc.sd[i] <- sd(repV[1,])
  model.tib$sensitivity[i] <- mean(repV[2,])
  model.tib$sensitivity.sd[i] <- sd(repV[2,])
  model.tib$specificity[i] <- mean(repV[3,])
  model.tib$specificity.sd[i] <- sd(repV[3,])
  model.tib$R2[i] <- mean(repV[4,])
  model.tib$R2.sd[i] <- sd(repV[4,])
  model.tib$F1[i] <- mean(repV[5,])
  model.tib$F1.sd[i] <- sd(repV[5,])
  model.tib$AUC[i] <- mean(repV[6,])
  model.tib$AUC.sd[i] <- sd(repV[6,])
  
  modelList[[i]] <- glm(currF,family='binomial',data=dat.c,maxit=100,epsilon=1e-14)
  print(sprintf('%i out of %i complete.',i,nrow(combs)))
}

cor.test(model.tib$AUC,model.tib$R2)
cor.test(model.tib$F1,model.tib$R2)

model.tib %>% 
  dplyr::arrange(-F1) %>%
  dplyr::mutate(F1rank=row_number()) %>%
  dplyr::arrange(-AUC) %>%
  dplyr::mutate(AUCrank=row_number()) %>%
  dplyr::arrange(-R2) %>%
  dplyr::mutate(R2rank=row_number()) %>%
  dplyr::arrange(-acc) %>%
  dplyr::mutate(accrank=row_number()) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(rank_sum=F1rank+AUCrank+accrank) %>%
  dplyr::arrange(rank_sum)

model.tib %>%
  dplyr::arrange(-AUC) %>%
  dplyr::select(-acc,-sensitivity,-specificity,-sensitivity.sd,-specificity.sd,-acc.sd,-R2.sd,-F1.sd,-AUC.sd) %>%
  dplyr::slice(1:10) %>%
  dplyr::mutate(vars=gsub('\n','/',vars)) %>%
  dplyr::mutate(vars=gsub('median_rl','readlength',vars)) %>%
  dplyr::mutate(vars=gsub('g_ccontent','GCcontent',vars)) %>%
  dplyr::mutate(vars=gsub('actual_cov','coverage',vars)) %>%
  dplyr::rename(Variables=vars) %>%
  knitr::kable("latex",digits=3,align=c('l','c','c'))
  
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,3,4,5)],applysigns = TRUE)
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,3,5)],applysigns = TRUE)
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,4,5)],applysigns = TRUE)
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,5)],applysigns = TRUE)

dat.glm <- glm(sig~actual_cov+damage+contiglength,
               data=dat.c,
               family='binomial',epsilon=1e-8)
summary(dat.glm)
fit.tab <- table(factor(predict(dat.glm,dat.c,type='response')>0.5,levels=c(TRUE,FALSE)),
                 factor(dat.c$damage>0,levels=c(TRUE,FALSE)))
caret::confusionMatrix(fit.tab)
cv.tab <- table(factor(predict(modelList[[i]],dat.d,type='response')>0.5,levels=c(TRUE,FALSE)),
                factor(dat.d$damage>0,levels=c(TRUE,FALSE)))
caret::confusionMatrix(cv.tab)
lm.beta::lm.beta(dat.glm)

dat.d.f <- dat.d %>%
  dplyr::filter(actual_cov>20)
cv.f.tab <- table(factor(predict(modelList[[i]],dat.d.f,type='response')>0.5,levels=c(TRUE,FALSE)),
                factor(dat.d.f$damage>0,levels=c(TRUE,FALSE)))
caret::confusionMatrix(cv.f.tab)







