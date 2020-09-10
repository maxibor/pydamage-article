require(ggplot2)
require(tidyverse)
require(cowplot)
require(betareg)
require(lmerTest)
inFolder <- '~/Dropbox/MPI/Pydamage/data/gc/' 
tsvs <- grep('.gz',paste0(inFolder,list.files(inFolder)),invert=T,value=T)
tsvList <- lapply(tsvs,read.delim,sep='\t')
dat <- bind_rows(tsvList) %>%
  janitor::clean_names() %>%
  dplyr::rename(rpt=`repeat`) %>%
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
dat.q <- dat %>%
  dplyr::select(genome,readlength,damage,simu_cov,simu_contig_length,qvalue,g_ccontent) %>%
  dplyr::filter(!is.na(qvalue)) %>%
  dplyr::mutate(sig=qvalue<0.05) %>%
  dplyr::mutate(simu_contig_length=fct_reorder(simu_contig_length,as.numeric(gsub("(-).*",'',simu_contig_length))),
                simu_cov=fct_reorder(simu_cov,as.numeric(gsub("(-).*",'',simu_cov))))
dat.q %>%
  dplyr::filter(damage==0) %>%
  .$sig %>% mean()
dat.q %>%
  dplyr::filter(damage!=0) %>%
  .$sig %>% mean()
dat.obs <- dat.q %>%
  dplyr::group_by(genome,readlength,simu_cov,simu_contig_length,damage) %>%
  dplyr::summarise(obsprop=mean(sig))%>%
  dplyr::mutate(obs=obsprop+1e-3*(obsprop==0)-1e-3*(obsprop==1))
modelList <- list()
vs <- c('readlength','simu_cov','simu_contig_length','damage','g_ccontent')
n <- length(vs)
combs <- as.matrix(expand.grid(as.data.frame(matrix(rep(c(0,1),times=n),nrow=2))))
modelList[[1]] <- glm(sig ~ 1,family='binomial',data=dat.q)
#glmer(sig ~ 1 + (1 | genome),data=dat.q,family='binomial',control=glmerControl(optimizer="bobyqa"),nAGQ=0)
mods <- rep('x',nrow(combs))
for(i in 2:nrow(combs)){
  currF <- formula(paste0('sig~',paste0(vs[which(combs[i,]==1)],collapse='+')))
  #formula(paste0('sig~',paste0(vs[which(combs[i,]==1)],collapse='+'),'+ (1 | genome)'))
  mods[i] <- paste0(vs[which(combs[i,]==1)],collapse='\n')
  modelList[[i]] <- glm(currF,data=dat.q,family='binomial')
  print(sprintf('%i out of %i complete.',i,nrow(combs)))
}
model.tib <- tibble(BIC=unlist(lapply(modelList,BIC)),vars=mods,AIC=unlist(lapply(modelList,AIC)))
model.gg <- model.tib %>%
  dplyr::mutate(vars=gsub('\n','+\n',vars),
                vars=gsub('simu_cov','simulated\ncoverage',vars),
                vars=gsub('g_ccontent','GC content',vars),
                vars=gsub('simu_contig_length','simulated\ncontig length',vars)) %>%
  dplyr::arrange(AIC) %>%
  dplyr::slice(1:10) %>%
  ggplot(aes(x=fct_reorder(vars,AIC),y=AIC))+
  theme_bw()+
  geom_point(size=3)+
  theme(axis.text.x=element_text(size=8))+
  geom_point(x=1,y=min(model.tib$AIC),pch=1,size=5,col='red')+
  xlab('Candidate Variable Set')
model.gg
ggsave("~/Dropbox/MPI/Pydamage/plots/AIC.pdf",model.gg,dp=500,
       height=210,width=297,units='mm')    
dat.glm <- modelList[[which.min(model.tib$BIC)]]
summary(dat.glm)
dat.newtib <- expand.grid(genome=unique(dat.q$genome),
                          readlength=unique(dat.q$readlength),
                          simu_cov=unique(dat.q$simu_cov),
                          simu_contig_length=unique(dat.q$simu_contig_length),
                          damage=unique(dat.q$damage),
                          g_ccontent=seq(min(dat.q$g_ccontent),max(dat.q$g_ccontent),length.out=10)) %>%
  as_tibble()
dat.pred <- dat.newtib  %>%
  dplyr::mutate(minLen=as.numeric(gsub("(-).*",'',simu_contig_length)),
                minCov=as.numeric(gsub("(-).*",'',simu_cov))) %>%
  dplyr::mutate(pred=predict(dat.glm,newdata=dat.newtib,type='response')) %>%
  dplyr::mutate(CF=ifelse(pred>=0.5,'Y','N')) %>%
  inner_join(dat.obs) %>%
  dplyr::filter(damage!=0) %>%
  dplyr::mutate(res.prop=obsprop-pred,CFobs=(obsprop>0.5))  
acc.gg <- ggplot(data=dat.pred,aes(x=fct_reorder(simu_cov,minCov),
                         y=fct_reorder(simu_contig_length,minLen)))+
  theme_bw()+
  geom_tile(aes(fill=pred*100,col=CF))+
  facet_grid(readlength~damage,scale='free')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  xlab('Coverage X')+
  ylab("Contig Length")+
  scale_colour_manual(values=c('red','green'))+
  scale_fill_continuous(limits=c(0,100),low='black',name='Predicted Accuracy')
obs.gg <- ggplot(data=dat.pred,aes(x=fct_reorder(simu_cov,minCov),
                                   y=fct_reorder(simu_contig_length,minLen)))+
  theme_bw()+
  geom_tile(aes(fill=obsprop*100,col=CFobs))+
  facet_grid(readlength~damage,scales='free')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  xlab('Coverage X')+
  ylab("Contig Length")+
  scale_colour_manual(values=c('red','green'),guide=F)+
  scale_fill_continuous(limits=c(0,100),low='black',name='Observed\nAccuracy')
res.gg <- ggplot(data=dat.pred,aes(x=fct_reorder(simu_cov,minCov),
                                   y=fct_reorder(simu_contig_length,minLen)))+
  theme_bw()+
  geom_tile(aes(fill=res.prop,col=CF))+
  facet_grid(genome~damage,scales='free')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  geom_text(aes(label=round(res.prop,2)),col='red',size=2)+
  xlab('Coverage X')+
  ylab("Contig Length")+
  scale_colour_manual(values=c('black','red'))+
  scale_fill_continuous(low='black',name='Residual Accuracy')
# res.gg
# ggsave("~/Dropbox/MPI/Pydamage/plots/Residuals.jpeg",res.gg,dp=500,
#        height=210,width=297,units='mm')  
# plot_grid(obs.gg,acc.gg,nrow=2,labels=c('A','B'))
ggsave("~/Dropbox/MPI/Pydamage/plots/Pred.jpeg",plot_grid(obs.gg,acc.gg,nrow=2,labels=c('A','B')),dp=500,
       height=210,width=297,units='mm') 
ggsave("~/Dropbox/MPI/Pydamage/plots/Obs.jpeg",obs.gg,dp=500,
       height=160,width=297,units='mm') 
acc.cf <- ggplot(data=dat.pred,aes(x=fct_reorder(simu_cov,minCov),
                                   y=fct_reorder(simu_contig_length,minLen)))+
  theme_bw()+
  geom_tile(aes(fill=CF))+
  facet_grid(genome~damage)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  xlab('Coverage X')+
  ylab("Contig Length")
# acc.cf
ggsave("~/Dropbox/MPI/Pydamage/plots/PredAccuracy.jpeg",acc.gg,dp=500,
       height=210,width=297,units='mm')