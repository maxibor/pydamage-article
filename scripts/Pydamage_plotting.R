fit.gg <- model.tib %>%
  dplyr::mutate(vars=gsub('\n','+\n',vars),
                vars=gsub('simu_cov','simulated\ncoverage',vars),
                vars=gsub('g_ccontent','GC content',vars),
                vars=gsub('simu_contig_length','simulated\ncontig length',vars),
                ba=(sensitivity+specificity)/2) %>%
  dplyr::arrange(-R2) %>%
  dplyr::slice(1:4) %>%
  dplyr::select(vars,R2,specificity,sensitivity) %>%
  gather(measure,y,R2:sensitivity) %>%
  dplyr::mutate(numVars=paste0(str_count(vars,"\n")+1,"\nvariables")) %>%
  ggplot(aes(x=vars,y=y,fill=measure))+
  theme_bw()+
  geom_bar(stat='identity',position='dodge',col='black')+
  facet_grid(.~numVars,scales='free')+
  scale_fill_discrete(labels=c(expression(R^2),"Sensitivity","Specificity"),name='')+
  xlab('')+
  ylab('')+
  theme(legend.position='bottom')
facetFix(fit.gg)  
ggsave("~/Dropbox/MPI/Pydamage/plots/ModelFit.png",facetFix(fit.gg),dp=500,
       height=210,width=297,units='mm')

N.marginal <- 20
cl_ac.tib <- as_tibble(expand.grid(contiglength=exp(seq(log(min(dat.c.raw$contiglength)),
                                                        log(50000),#max(dat.c.raw$contiglength),
                                                        length.out=N.marginal)),
                                   actual_cov=seq(min(dat.c.raw$actual_cov),120,
                                                  length.out=N.marginal),
                                   damage=unique(dat.c.raw$damage))) %>%
  dplyr::mutate(median_rl=mean(dat.c.raw$median_rl),
                g_ccontent=mean(dat.c.raw$g_ccontent))


bks <- c(0.1,0.25,0.5,1,2,3,4,5)
cl_ac.gg <- cl_ac.tib %>%
  dplyr::mutate(pred=predict(dat.glm,cl_ac.tib,type='response'),
                sig=(pred>=0.5)) %>%
  ggplot(aes(x=log((contiglength)/1e4),y=(actual_cov),fill=pred))+
  theme_bw()+
  geom_tile(aes(col=sig))+
  facet_wrap(.~damage,nrow=2)+
  coord_cartesian(expand=F)+
  xlab('Contig Length x 10,000')+
  ylab('Coverage')+
  scale_fill_continuous(name='',limits=c(0,1))+
  scale_colour_discrete(guide=F)+
  theme(legend.position='right')+
  scale_x_continuous(breaks=log(bks),labels=bks)
cl_ac.gg  

ggsave("~/Dropbox/MPI/Pydamage/plots/Predicted_Accuracy.png",cl_ac.gg,dp=500,
       height=210,width=297,units='mm')


##### Make observed plot ----
bksobs <- exp(seq(log(min(dat.c.raw$contiglength)),log(50000),length.out=N.marginal+10))/1e4
xlbs <- c(paste0('(',c(0,bksobs[1:(length(bks)-1)]),',',bksobs,']'),paste0('>',bksobs[length(bksobs)]))
rp <- 4
ylb <- c(1,rep('',rp),30,rep('',rp),60,rep('',rp),90,rep('',rp),120)
xlb <- c(rep('',5),0.1,rep('',5),0.25,rep('',4),0.5,rep('',3),1,rep('',2),2,rep('',2),3,rep('',1),4,rep('',1),5)
obs.gg <- dat.c.raw %>%
  dplyr::mutate(contiglength=(contiglength/1e4)) %>%
  dplyr::mutate(actual_cov=cut(actual_cov,c(1e-16,sort(unique(cl_ac.tib$actual_cov)),Inf))) %>%
  dplyr::mutate(contiglength=cut(contiglength,c(-Inf,bksobs,Inf),labels=xlbs)) %>%
  # dplyr::mutate(contiglength=cut(contiglength,c(-Inf,log(bks),Inf),labels=xlbs)) %>%
  dplyr::select(damage,contiglength,actual_cov,sig) %>%
  dplyr::mutate(sig=ifelse(damage==0,!sig,sig)) %>%
  dplyr::group_by(damage,contiglength,actual_cov) %>%
  dplyr::summarise(acc=mean(sig)) %>%
  dplyr::filter(damage!=1000,!is.na(contiglength),!is.na(actual_cov)) %>%
  ggplot(aes(x=(contiglength),y=(actual_cov),fill=acc))+
  theme_bw()+
  geom_tile(aes(col=acc>0.5))+
  facet_wrap(.~damage,nrow=2)+
  coord_cartesian(expand=F)+
  xlab('Contig Length x 10,000')+
  ylab('Coverage')+
  scale_fill_continuous(name='',limits=c(0,1))+
  scale_colour_discrete(guide=F)+
  theme(legend.position='right')+
  scale_y_discrete(labels=ylb)+
  scale_x_discrete(labels=xlb)+
  theme(axis.ticks.y=element_line(size=ifelse(ylb!='',0.5,0)),
        axis.ticks.x=element_line(size=ifelse(xlb!='',0.5,0)))

ggsave("~/Dropbox/MPI/Pydamage/plots/Observed_Accuracy.png",obs.gg,dp=500,
       height=210,width=297,units='mm')
dat.c %>%
  dplyr::filter(damage==0) %>%
  .$sig %>%
  table()