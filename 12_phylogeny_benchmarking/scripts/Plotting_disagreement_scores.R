#Author: Mike Spencer Chapman

library(dplyr)
library(ggplot2)
library(stringr)

load("disagreement_scores")
samples=readLines("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/EM_samples.txt")
names(all_results)<-samples
out_list=lapply(samples,function(sample) {
  list<-all_results[[sample]]
  df=rbind(data.frame(sample=sample,type="data",ds=list$data),data.frame(sample=sample,type="random shuffles",ds=list$rand))
  return(df)
  })

out_df=dplyr::bind_rows(out_list)
out_df<-out_df%>%mutate(sample=str_split(sample,pattern="_",simplify=T)[,1])%>%filter(!grepl("CB",sample)&!grepl("PX",sample))%>%mutate(sample=factor(sample,levels=order_samples))

order_samples=c("KX001", "KX002", "SX001", "AX001", "KX007","KX008", "KX004", "KX003")

breaks <- 10^(-2:2)
minor_breaks <- NULL


plot3<-out_df%>%
  filter(type=="random shuffles")%>%
  ggplot(aes(y=sample,x=ds,col=type))+
  geom_jitter(size=0.15,width = 0.05,alpha=0.5)+
  geom_point(data=out_df%>%filter(type=="data"))+
  scale_x_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  theme_bw()+
  labs(x="Disagreement score",y="Individual",col="Type")#+
  #coord_flip()

plot3<-out_df%>%
  filter(type=="random shuffles")%>%
  ggplot(aes(x=sample,y=ds,col=type))+
  geom_jitter(size=0.15,height = 0.05,alpha=0.5)+
  geom_point(data=out_df%>%filter(type=="data"))+
  scale_y_log10(breaks = breaks,
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                minor_breaks=NULL)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))+
  labs(x="Individual",y="Disagreement score",col="Type")+
  annotation_logticks(sides="l",short=unit(0.05,"cm"),mid=unit(0.05,"cm"),long=unit(0.1,"cm"))

fold_change_list<-lapply(out_list,function(df) {
  sim_ds<-df%>%filter(type=="simulation")%>%pull(ds)
  data_ds<-df%>%filter(type=="data")%>%pull(ds)
  fold_change=sim_ds/data_ds
  return(data.frame(sample=df$sample[1],fold_change=fold_change))
})
fold_change_df<-dplyr::bind_rows(fold_change_list)

plot<-fold_change_df%>%
  mutate(sample=str_split(sample,pattern="_",simplify=T)[,1])%>%
  mutate(ds_change=-log10(1/fold_change))%>%
  filter(!grepl("CB",sample)&!grepl("PX",sample))%>%
  ggplot(aes(x=sample,y=ds_change))+
  geom_jitter(size=0.2,alpha=0.5)+
  #scale_y_log10(limits=c(1,15000))+
  scale_y_continuous(limits=c(0,4.5))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="-log10 Disagreement score of data/ random shuffles",
       x="Sample")

plot2<-fold_change_df%>%
  mutate(sample=str_split(sample,pattern="_",simplify=T)[,1])%>%
  mutate(ds_change=-log10(1/fold_change))%>%
  filter(!grepl("CB",sample)&!grepl("PX",sample))%>%
  group_by(sample)%>%
  summarise(median=median(ds_change))%>%
  ggplot(aes(x=sample,y=median))+
  geom_point(size=1.5)+
  #scale_y_log10(limits=c(1,15000))+
  scale_y_continuous(limits=c(0,4.5))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="-log10 Disagreement score of data/ random shuffles",
       x="Sample")

ggsave(plot,filename = "Disagreement_score_plot.pdf",device="pdf",width=4,height=5)
ggsave(plot2,filename = "Disagreement_score_single_point_plot.pdf",device="pdf",width=4,height=5)
ggsave(plot3,filename = "Disagreement_score.pdf",device="pdf",width=5,height=5)