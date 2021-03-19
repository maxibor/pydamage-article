# Rerquired packages and custom functions
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

dat <- read.table('~/Dropbox/MPI/Pydamage/continuous_data/PYD_simulation_results.csv',header=T,sep='\t') %>%
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
combs <- as.matrix(expand.grid(as.data.frame(matrix(rep(c(0,1),times=n),nrow=2))))

for(i in 2:nrow(combs)){
  if(i==2){
    model.tib <- tibble(vars=rep('x',nrow(combs)),acc=0,sensitivity=0,specificity=0,R2=0)
    modelList <- list()
    modelList[[1]] <- glm(sig ~ 1,family='binomial',data=dat.c,maxit=100,epsilon=1e-14)
    model.tib$vars[1] <- 'null'
    model.tib$R2[1] <- rcompanion::nagelkerke(modelList[[1]])$Pseudo.R.squared.for.model.vs.null[3]
    curr.fit <- table(factor(predict(modelList[[1]],type='response')>0.5,levels=c(TRUE,FALSE)),
                      factor(dat.c$damage>0,levels=c(TRUE,FALSE)))
    curr.sens <- caret::sensitivity(curr.fit)
    curr.spec <- caret::specificity(curr.fit)
    model.tib$acc[1] <- mean(c(curr.sens,curr.spec))#
  }
  currF <- formula(paste0('sig~',paste0(vs[which(combs[i,]==1)],collapse='+')))
  model.tib$vars[i] <- paste0(vs[which(combs[i,]==1)],collapse='\n')
  modelList[[i]] <- glm(currF,data=dat.c,family='binomial',maxit=100,epsilon=1e-14)
  model.tib$R2[i] <- rcompanion::nagelkerke(modelList[[i]])$Pseudo.R.squared.for.model.vs.null[3]
  curr.fit <- table(factor(predict(modelList[[i]],dat.d,type='response')>0.5,levels=c(TRUE,FALSE)),
                    factor(dat.d$damage>0,levels=c(TRUE,FALSE)))
  curr.sens <- caret::sensitivity(curr.fit)
  curr.spec <- caret::specificity(curr.fit)
  model.tib$sensitivity[i] <- curr.sens
  model.tib$specificity[i] <- curr.spec
  model.tib$acc[i] <- mean(c(curr.sens,curr.spec))# caret::confusionMatrix(curr.fit)$overall['Accuracy']
  
  print(sprintf('%i out of %i complete.',i,nrow(combs)))
}

# Model comparisons
model.gg <- model.tib %>%
  dplyr::mutate(vars=gsub('\n','+\n',vars),
                vars=gsub('simu_cov','simulated\ncoverage',vars),
                vars=gsub('g_ccontent','GC content',vars),
                vars=gsub('simu_contig_length','simulated\ncontig length',vars),
                R2.mf=1-unlist(lapply(modelList,logLik))/logLik(modelList[[1]])) %>%
  dplyr::slice(1:160) %>%
  ggplot(aes(x=fct_reorder(vars,-R2),y=R2))+
  theme_bw()+
  geom_hline(yintercept=0.5,col='red',linetype='dashed')+
  geom_point(size=3)+
  theme(axis.text.x=element_text(size=6))+
  geom_point(x=1,y=max(model.tib$R2),pch=1,size=5,col='red')+
  xlab('Candidate Variable Set')
model.gg

model.tib %>%
  ggplot(aes(x=sensitivity,y=specificity,col=R2))+
  theme_bw()+
  geom_point()+
  ggrepel::geom_label_repel(aes(label=vars))

model.tib %>%
  rowwise() %>%
  dplyr::mutate(ba=mean(c(sensitivity,specificity))) %>%
  dplyr::arrange(-R2) %>%
  dplyr::filter(specificity>0.15) 
  
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,3,4,5)],applysigns = TRUE)
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,3,5)],applysigns = TRUE)
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,4,5)],applysigns = TRUE)
dat.c %>%rwa::rwa(outcome="sig",predictors=vs[c(1,2,5)],applysigns = TRUE)

dat.glm <- glm(sig~actual_cov+damage+contiglength,
               data=dat.c,#dplyr::slice(dat.c,c(which(dat.c$damage==0),sample(which(dat.c$damage!=0),sum(dat.c$damage==0)))),
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

subN <- 100
sub.Cov <- tibble(damage=mean(dat.c$damage),contiglength=mean(dat.c$contiglength),
       actual_cov=seq(0,max(dat.c$actual_cov,na.rm=T),length.out=subN))
sub.damage <- tibble(actual_cov=mean(dat.c$actual_cov),contiglength=mean(dat.c$contiglength),
                  damage=seq(0,max(dat.c$damage,na.rm=T),length.out=subN))
sub.cl <- tibble(actual_cov=mean(dat.c$actual_cov),damage=mean(dat.c$damage),
                 contiglength=seq(0,max(dat.c$contiglength,na.rm=T),length.out=subN))
gg1 <- sub.Cov %>%
  mutate(pred=predict(dat.glm,sub.Cov,type='response')*100) %>%
  ggplot(aes(x=actual_cov,y=pred))+
  theme_bw()+
  geom_line()+
  ylim(c(0,100))
gg2 <- sub.damage %>%
  mutate(pred=predict(dat.glm,sub.damage,type='response')*100) %>%
  ggplot(aes(x=damage,y=pred))+
  theme_bw()+
  geom_line()+
  ylim(c(0,100))
gg3 <- sub.cl %>%
  mutate(pred=predict(dat.glm,sub.cl,type='response')*100) %>%
  ggplot(aes(x=contiglength,y=pred))+
  theme_bw()+
  geom_line()+
  ylim(c(0,100))

gridExtra::grid.arrange(gg1,gg2,gg3,nrow=1)







