#CODE ACCOMPANYING DAVIES ET AL (2023)
#J. Davies & M. Kirkland
#This code takes the pre-processed cuckoo migratory timing data and performs the analyses detailed in the manuscript:
#"Spring arrival of the common cuckoo at breeding grounds is strongly determined by environmental conditions in tropical Africa"

#Load packages
require(lubridate) #Various date-related functions needed
require(rgdal)
require(rgeos)
require(geosphere) #For distCosine
require(raster) #For aggregate
require(reshape2) #For melt
require(lme4)
require(nlme)
require(emmeans) #For emmeans
require(MuMIn) #For r.squaredGLMM. R2m = R2 for fixed effects alone; R2c = R2 for the entire model, i.e. fixed + random effects.
require(corrplot) #For corrplot
require(RColorBrewer)
require(brms) #For brm
require(R2jags)
require(mice)
require(VIM)
require(reshape2)
require(fastDummies)
require(ggplot2)

#1. Data preparation
#Load processed cuckoo migratory timing data
dir_path<-getwd()
mig_timing_df<-read.csv(paste0(dir_path,"/mig_timing_df.csv"),header=T)
mig_timing_df_uncertainty<-read.csv(paste0(dir_path,"/mig_timing_df_uncertainty.csv"),header=T)

milestone_names<-c("dep_brgr","dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b","dep_wgr_b","dep_wgr_a","dep_wa","arr_UK","arr_brgr")
threshold<-5 #Threshold number of days' uncertainty arising from measurement error
un_m_ptt<-unique(mig_timing_df$ptt)

#Four different ways of standardising migratory timing:
#First, relative to population median for each milestone
mig_timing_df_std<-mig_timing_df
for(j in 3:16){
  mig_timing_df_std[,j]<-mig_timing_df[,j]-median(mig_timing_df[,j],na.rm=T) 
}
mig_timing_df_std$ptt_num<-as.numeric(as.factor(mig_timing_df_std$ptt)) #For colour
mtd_std_melted<-melt(mig_timing_df_std[,c("ptt","year","mig_dir","br_hab","region",milestone_names)],id=c("ptt","year","mig_dir","br_hab","region"))
mtd_std_melted<-mtd_std_melted[-which(is.na(mtd_std_melted$value) | is.na(mtd_std_melted$mig_dir)),]
mtd_std_melted$value_abs<-abs(mtd_std_melted$value)

#Second, relative to individual mean for each milestone
mig_timing_df_stdind<-mig_timing_df
for(i in 1:length(un_m_ptt)){
  temp_ind_df<-mig_timing_df[mig_timing_df$ptt==un_m_ptt[i],3:16]
  temp_ind_df_std<-temp_ind_df
  for(j in 1:ncol(temp_ind_df)){
    if(any(is.na(temp_ind_df[,j])==F)){
      if(length(which(is.na(temp_ind_df[,j])==F))<2) {temp_ind_df_std[,j]<-NA} else { #If there is only one non-NA value for a given milestone, don't take the mean
        temp_ind_df_std[,j]<-temp_ind_df[,j]-mean(temp_ind_df[,j],na.rm=T)
      }
    } else {
      temp_ind_df_std[,j]<-NA
    }
  }
  mig_timing_df_stdind[mig_timing_df_stdind$ptt==un_m_ptt[i],3:16]<-temp_ind_df_std
}
mig_timing_df_stdind$ptt_num<-as.numeric(as.factor(mig_timing_df_stdind$ptt)) #For colour
mtd_stdind_melted<-melt(mig_timing_df_stdind[,c("ptt","year","mig_dir","outcome","known_age_binary","last_milestone_pre_death",milestone_names)],id=c("ptt","year","mig_dir","known_age_binary","last_milestone_pre_death","outcome"))
mtd_stdind_melted<-mtd_stdind_melted[which(is.na(mtd_stdind_melted$value)==F),]
mtd_stdind_melted$mig_dir<-relevel(as.factor(mtd_stdind_melted$mig_dir),ref = "SW") #Put SW first, so it's on the same side as 'lowland'

#How many birds had multiple mig cycles?
length(unique(mig_timing_df_stdind$ptt[which(apply(is.na(mig_timing_df_stdind[,milestone_names])==F,1,any))]))

#Third, relative to minimum for each milestone
mig_timing_df_stdind_min<-mig_timing_df
for(i in 1:length(un_m_ptt)){
  temp_ind_df<-mig_timing_df[mig_timing_df$ptt==un_m_ptt[i],3:16]
  temp_ind_df_std<-temp_ind_df
  for(j in 1:ncol(temp_ind_df)){
    if(any(is.na(temp_ind_df[,j])==F)){
      if(length(which(is.na(temp_ind_df[,j])==F))<2) {temp_ind_df_std[,j]<-NA} else { #If there is only one non-NA value for a given milestone, don't take the mean
        temp_ind_df_std[,j]<-temp_ind_df[,j]-min(temp_ind_df[,j],na.rm=T)
      }
    } else {
      temp_ind_df_std[,j]<-NA
    }
  }
  mig_timing_df_stdind_min[mig_timing_df_stdind_min$ptt==un_m_ptt[i],3:16]<-temp_ind_df_std
}
mig_timing_df_stdind_min$ptt_num<-as.numeric(as.factor(mig_timing_df_stdind_min$ptt)) #For colour

#Fourth, scale (subtract mean and divide by SD)
mig_timing_df_scaled<-mig_timing_df
for(j in 3:16){
  mig_timing_df_scaled[,j]<-scale(mig_timing_df[,j])[,1]
}

#Fifth, relative to *migratory route* median for each milestone
mig_timing_df_std_bymigdir<-mig_timing_df
mig_timing_df_std_bymigdir[,3:16]<-NA
for(j in 3:16){
  mig_timing_df_std_bymigdir[which(mig_timing_df$mig_dir=="SW"),j]<-
    mig_timing_df[which(mig_timing_df$mig_dir=="SW"),j]-median(mig_timing_df[which(mig_timing_df$mig_dir=="SW"),j],na.rm=T) 
  mig_timing_df_std_bymigdir[which(mig_timing_df$mig_dir=="SE"),j]<-
    mig_timing_df[which(mig_timing_df$mig_dir=="SE"),j]-median(mig_timing_df[which(mig_timing_df$mig_dir=="SE"),j],na.rm=T) 
}

#How many non-NA data per milestone?
data.frame(lapply(apply(is.na(mig_timing_df[,milestone_names])==F,2,which),length))

mtd_melted<-melt(mig_timing_df[,c("ptt","year","mig_dir","mig_dir_3cat","br_hab","last_milestone_pre_death",milestone_names)],id=c("ptt","year","mig_dir","mig_dir_3cat","br_hab","last_milestone_pre_death"))

#Number of events and birds per milestone, for results tables
n_events_per_milestone<-unlist(lapply(apply(is.na(mig_timing_df[,milestone_names])==F,2,which),length))
n_birds_per_milestone<-numeric(length(milestone_names))
for(i in 1:length(milestone_names)){
  n_birds_per_milestone[i]<-length(unique(mig_timing_df$ptt[is.na(mig_timing_df[,milestone_names[i]])==F]))
}

#And again for within-ind anomaly
n_events_per_milestone_withinind<-unlist(lapply(apply(is.na(mig_timing_df_stdind[,milestone_names])==F,2,which),length))
n_birds_per_milestone_withinind<-numeric(length(milestone_names))
for(i in 1:length(milestone_names)){
  n_birds_per_milestone_withinind[i]<-length(unique(mig_timing_df_stdind$ptt[is.na(mig_timing_df_stdind[,milestone_names[i]])==F]))
}

#Fig 1: plot of fixes by migratory stage
#Define layout matrix and plotting order
layout.matrix <- matrix(c(0,0,21,21,0,0,
                          0,0,22,22,0,0,
                          0,0,23,24,0,0,
                          17,17,23,24,1,1,
                          18,18,23,24,2,2,
                          19,20,0,0,3,4,
                          19,20,25,27,3,4,
                          19,20,26,28,3,4,
                          13,13,26,28,5,5,
                          14,14,26,28,6,6,
                          15,16,0,0,7,8,
                          15,16,9,9,7,8,
                          15,16,10,10,7,8,
                          0,0,11,12,0,0,
                          0,0,11,12,0,0,
                          0,0,11,12,0,0),
                        byrow=T, nrow = 16, ncol = 6)
layout(layout.matrix,heights=c(0.5,0.5,0.5,0.5,0.5,
                               0.5,0.5,0.5,0.5,0.5,
                               0.5,0.5,0.5,0.5,0.5,0.5))
for(i in 1:length(un_mig_stage_ordered)){ #Plot fixes by migratory stage
  print(i)
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",ylim=c(0,1))
  text(0,0.1,mig_stage_text_df$title[i],cex=1.4,font=2)
  plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",ylim=c(0,1))
  text(0,0.7,mig_stage_text_df$text[i],cex=1.3)
  par(mar=c(5.1, 4.6, 0.1, 0.1))
  plot(0,0,
       xlim=range(d1_spdf_2_with_mig_stage$location.long),
       ylim=range(d1_spdf_2_with_mig_stage$location.lat),
       col="white",
       xlab="Longitude",
       ylab="Latitude",
       cex.axis=1.5,
       cex.lab=1.5)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightblue1")
  plot(world_shp,add=T,col="light green",border="gray50")
  plot(sah,col="orange",border=NA,add=T)
  plot(ca_rf,col="dark green",border=NA,add=T)
  points(location.lat~location.long,data=d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage==un_mig_stage_ordered[i],],pch=20,cex=0.5)
  if(i==3){
    temp_jul<-d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage==un_mig_stage_ordered[i],]$jul
    temp_jul[temp_jul<150]<-temp_jul[temp_jul<150]+365
    hist(temp_jul,
         breaks=seq(180,550,by=10),
         xlab="Julian day",
         freq=F,
         main="",
         cex.axis=1.4,
         cex.lab=1.4)
    abline(v=365,col="blue")
    abline(v=182,col="red")
  } else {
    hist(d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage==un_mig_stage_ordered[i],]$jul,breaks=seq(0,370,by=10),
         xlab="Julian day",
         freq=F,
         main="",
         cex.axis=1.4,
         cex.lab=1.4)
    abline(v=0,col="blue")
    abline(v=182,col="red")
  }
}
#Northbound migration map
par(mar=c(0.1, 0.1, 0.1, 0.1))
plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",xlim=c(0,1))
text(0.6,0,paste0("Northbound"),font=2,cex=1.4,adj=0.5)
par(mar=c(5.1, 4.6, 0.1, 0.1))
plot(0,0,
     xlim=range(d1_spdf_2_with_mig_stage$location.long),
     ylim=range(d1_spdf_2_with_mig_stage$location.lat),
     col="white",
     xlab="Longitude",
     ylab="Latitude",
     cex.axis=1.4,
     cex.lab=1.4)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightblue1")
plot(world_shp,add=T,col="light green",border="gray50")
plot(sah,col="orange",border=NA,add=T)
plot(ca_rf,col="dark green",border=NA,add=T)
for(j in 1:nrow(mig_timing_df)){
  if(is.na(mig_timing_df$dep_wgr_b_TS[j]) | is.na(mig_timing_df$dep_wgr_a_TS[j]) | is.na(mig_timing_df$mig_dir[j])) next
  temp_post_wgr_b<-mig_timing_df[j,c("dep_wgr_a_TS","dep_wa_TS","arr_UK_TS","arr_brgr_TS")]
  temp_max_avail_post_wgr_b<-as.character(temp_post_wgr_b[max(which(is.na(temp_post_wgr_b)==F))])
  temp_df<-d1[d1$ptt==mig_timing_df$ptt[j] & 
                d1$timestamp>mig_timing_df$dep_wgr_b_TS[j] &
                d1$timestamp<temp_max_avail_post_wgr_b,c("location.long","location.lat")]
  if(nrow(temp_df)<3) next
  temp_df_arrows<-as.matrix(temp_df)
  colnames(temp_df_arrows)<-NULL
  temp_df_arrows<-cbind(temp_df_arrows,rbind(temp_df_arrows[2:nrow(temp_df_arrows),],matrix(0,ncol=2,nrow=1)))
  temp_df_arrows<-temp_df_arrows[-nrow(temp_df_arrows),]
  if(mig_timing_df$mig_dir[j]=="SW") {temp_col<-"yellow"} else {
    temp_col<-"red"
  }
  for(i in 1:nrow(temp_df_arrows)){
    arrows(temp_df_arrows[i,1],temp_df_arrows[i,2],temp_df_arrows[i,3],temp_df_arrows[i,4],length=0,col=temp_col)
  }
}
#Southbound migration map
par(mar=c(0.1, 0.1, 0.1, 0.1))
plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="",xlim=c(0,1))
text(0.6,0,paste0("Southbound"),font=2,cex=1.4)
par(mar=c(5.1, 4.6, 0.1, 0.1))

plot(0,0,
     xlim=range(d1_spdf_2_with_mig_stage$location.long),
     ylim=range(d1_spdf_2_with_mig_stage$location.lat),
     col="white",
     xlab="Longitude",
     ylab="Latitude",
     cex.axis=1.4,
     cex.lab=1.4)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightblue1")
plot(world_shp,add=T,col="light green",border="gray50")
plot(sah,col="orange",border=NA,add=T)
plot(ca_rf,col="dark green",border=NA,add=T)
for(j in 1:nrow(mig_timing_df)){
  if(is.na(mig_timing_df$dep_UK_TS[j]) | is.na(mig_timing_df$mig_dir[j])) next
  temp_pre_wgr_b<-mig_timing_df[j,c("dep_UK_TS","fin_sah_sb_TS","arr_wgr_a_TS","arr_wgr_b_TS")]
  temp_max_avail_pre_wgr_b<-as.character(temp_pre_wgr_b[max(which(is.na(temp_pre_wgr_b)==F))])
  temp_df<-d1[d1$ptt==mig_timing_df$ptt[j] & 
                d1$timestamp>mig_timing_df$dep_brgr_TS[j] &
                d1$timestamp<temp_max_avail_pre_wgr_b,c("location.long","location.lat")]
  if(nrow(temp_df)<3) next
  temp_df_arrows<-as.matrix(temp_df)
  colnames(temp_df_arrows)<-NULL
  temp_df_arrows<-cbind(temp_df_arrows,rbind(temp_df_arrows[2:nrow(temp_df_arrows),],matrix(0,ncol=2,nrow=1)))
  temp_df_arrows<-temp_df_arrows[-nrow(temp_df_arrows),]
  if(mig_timing_df$mig_dir[j]=="SW") {temp_col<-"yellow"} else {
    temp_col<-"red"
  }
  
  for(i in 1:nrow(temp_df_arrows)){
    arrows(temp_df_arrows[i,1],temp_df_arrows[i,2],temp_df_arrows[i,3],temp_df_arrows[i,4],length=0,col=temp_col)
  }
}
par(mfrow=c(1,1)) #To negate unusual layout

################################################################################################################
#2. Variation in migratory timing
#Change in SD between milestones
milestone_name_pairs_2<-data.frame("from"=c("dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr"),
                                   "to"=c("fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr","dep_brgr"))
milestone_name_pairs_2$sd<-apply(mig_timing_df_stdind[,milestone_name_pairs_2$from],2,sd,na.rm=T)
milestone_name_pairs_2$delta_sd<-apply(mig_timing_df_stdind[,milestone_name_pairs_2$to],2,sd,na.rm=T)-milestone_name_pairs_2$sd
milestone_name_pairs_2$var<-apply(mig_timing_df_stdind[,milestone_name_pairs_2$from],2,var,na.rm=T)
milestone_name_pairs_2$delta_var<-apply(mig_timing_df_stdind[,milestone_name_pairs_2$to],2,var,na.rm=T)-milestone_name_pairs_2$var

#F-test: are the variances significantly different between milestones?
milestone_name_pairs_2$var_sig_diff<-NA
for(i in 1:nrow(milestone_name_pairs_2)){
  temp_f_test<-var.test(mig_timing_df_stdind[,milestone_name_pairs_2$from[i]],mig_timing_df_stdind[,milestone_name_pairs_2$to[i]])
  milestone_name_pairs_2$var_sig_diff[i]<-temp_f_test$p.value<0.05
}

mtd_melted_newmilestones<-mtd_melted[which(mtd_melted$variable %in% c("dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr")),]
mtd_melted_newmilestones$variable<-factor(mtd_melted_newmilestones$variable,levels=c("dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr"))
mtd_melted_newmilestones$nice_names<-character(nrow(mtd_melted_newmilestones))
mtd_melted_newmilestones$nice_names[mtd_melted_newmilestones$variable=="dep_brgr"]<-"Depart breeding\ngrounds"
mtd_melted_newmilestones$nice_names[mtd_melted_newmilestones$variable=="fin_sah_sb"]<-"Completion Sahara\ncrossing"
mtd_melted_newmilestones$nice_names[mtd_melted_newmilestones$variable=="arr_wgr_a"]<-"Arrive wintering\ngrounds"
mtd_melted_newmilestones$nice_names[mtd_melted_newmilestones$variable=="dep_wgr_a"]<-"Depart wintering\ngrounds"
mtd_melted_newmilestones$nice_names[mtd_melted_newmilestones$variable=="dep_wa"]<-"Depart West\nAfrica"
mtd_melted_newmilestones$nice_names[mtd_melted_newmilestones$variable=="arr_brgr"]<-"Arrive breeding\ngrounds"
mtd_melted_newmilestones$nice_names<-factor(mtd_melted_newmilestones$nice_names,
                                            levels=c("Depart breeding\ngrounds",
                                                     "Completion Sahara\ncrossing",
                                                     "Arrive wintering\ngrounds",
                                                     "Depart wintering\ngrounds",
                                                     "Depart West\nAfrica",
                                                     "Arrive breeding\ngrounds"))

#Boxplot: raw data, with significance of difference in variance
par(mar=c(8.1, 4.1, 1.1, 1.1))
boxplot(mtd_melted_newmilestones$value~mtd_melted_newmilestones$nice_names,
        las=2,
        frame=F,
        xlab="",
        ylab="Days since start of 1st calendar year of mig cycle",
        ylim=c(min(mtd_melted_newmilestones$value,na.rm=T),max(mtd_melted_newmilestones$value,na.rm=T)+20),
        main="")

#CAUTION HARD CODING! Results of variance comparison test stored in milestone_name_pairs$var_sig_diff.
arrows(1.1,max(mtd_melted_newmilestones$value,na.rm=T)+20,1.9,max(mtd_melted_newmilestones$value,na.rm=T)+20,length=0)
arrows(1.1,max(mtd_melted_newmilestones$value,na.rm=T)+10,1.1,max(mtd_melted_newmilestones$value,na.rm=T)+20,length=0)
arrows(1.9,max(mtd_melted_newmilestones$value,na.rm=T)+10,1.9,max(mtd_melted_newmilestones$value,na.rm=T)+20,length=0)
text(1.5,max(mtd_melted_newmilestones$value,na.rm=T)+30,"*",cex=2)

arrows(4.1,max(mtd_melted_newmilestones$value,na.rm=T)+20,4.9,max(mtd_melted_newmilestones$value,na.rm=T)+20,length=0)
arrows(4.1,max(mtd_melted_newmilestones$value,na.rm=T)+10,4.1,max(mtd_melted_newmilestones$value,na.rm=T)+20,length=0)
arrows(4.9,max(mtd_melted_newmilestones$value,na.rm=T)+10,4.9,max(mtd_melted_newmilestones$value,na.rm=T)+20,length=0)
text(4.5,max(mtd_melted_newmilestones$value,na.rm=T)+30,"*",cex=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))

#Estimation of within- and between-individual variance
#Random effects model (no fixed effects), rather than mixed model: allows us to apportion variance
var_cor_df<-data.frame(matrix(ncol=13,nrow=length(milestone_names)))
names(var_cor_df)<-c("PTT_var","PTT_var_lci","PTT_var_uci","resid_var","resid_var_lci","resid_var_uci","total_var","total_var_lci","total_var_uci","sig_84ci","repeatability","repeatability_lci","repeatability_uci")
row.names(var_cor_df)<-milestone_names

for(i in 1:length(milestone_names)){
  print(i)
  temp_form<-as.formula(paste0(milestone_names[i],"~(1|ptt)"))
  m_all<-brm(temp_form,
             data=mig_timing_df,
             iter=4000,
             control = list(adapt_delta = 0.99999999999,
                            max_treedepth=20));summary(m_all)
  s_m_all<-summary(m_all)
  ptt_var_draws<-as_draws_df(m_all, variable = "sd_ptt__Intercept")$sd_ptt__Intercept^2
  resid_var_draws<-as_draws_df(m_all, variable = "sigma")$sigma^2
  total_var_draws<-ptt_var_draws+resid_var_draws
  
  #Is the difference between between-ind and within-ind variation significant (using overlap of 84% CIs)?
  temp_PTT_lci84<-quantile(ptt_var_draws,0.08)
  temp_PTT_uci84<-quantile(ptt_var_draws,0.92)
  temp_resid_lci84<-quantile(resid_var_draws,0.08)
  temp_resid_uci84<-quantile(resid_var_draws,0.92)
  var_cor_df$sig_84ci[i]<-((temp_resid_lci84>temp_PTT_lci84 & temp_resid_lci84<temp_PTT_uci84) | (temp_resid_uci84>temp_PTT_lci84 & temp_resid_uci84<temp_PTT_uci84))==F
  
  var_cor_df$PTT_var[i]<-mean(ptt_var_draws)
  var_cor_df$PTT_var_lci[i]<-quantile(ptt_var_draws,0.025)
  var_cor_df$PTT_var_uci[i]<-quantile(ptt_var_draws,0.975)
  var_cor_df$resid_var[i]<-mean(resid_var_draws)
  var_cor_df$resid_var_lci[i]<-quantile(resid_var_draws,0.025)
  var_cor_df$resid_var_uci[i]<-quantile(resid_var_draws,0.975)
  var_cor_df$total_var[i]<-mean(total_var_draws)
  var_cor_df$total_var_lci[i]<-quantile(total_var_draws,0.025)
  var_cor_df$total_var_uci[i]<-quantile(total_var_draws,0.975)
  var_cor_df$repeatability[i]<-mean(ptt_var_draws/total_var_draws)
  var_cor_df$repeatability_lci[i]<-quantile(ptt_var_draws/total_var_draws,0.025)
  var_cor_df$repeatability_uci[i]<-quantile(ptt_var_draws/total_var_draws,0.975)
}

#Nice table ('A1') for presentation, with uncertainty included
nice_var_cor_df_tab<-data.frame(matrix(nrow=length(milestone_names),ncol=11))
names(nice_var_cor_df_tab)<-c("Milestone","Between_ind","Between_ind_CI","Within_ind","Within_ind_CI","sig_diff_84","Repeatability","Repeatability_CI","Mean_uncertainty","SD_uncertainty","N")
nice_var_cor_df_tab$Milestone<-milestone_names
nice_var_cor_df_tab$Between_ind<-round(var_cor_df$PTT_var,digits=2)
nice_var_cor_df_tab$Between_ind_CI<-paste0(round(var_cor_df$PTT_var_lci,digits=2),", ",round(var_cor_df$PTT_var_uci,digits=2))
nice_var_cor_df_tab$Within_ind<-round(var_cor_df$resid_var,digits=2)
nice_var_cor_df_tab$Within_ind_CI<-paste0(round(var_cor_df$resid_var_lci,digits=2),", ",round(var_cor_df$resid_var_uci,digits=2))
nice_var_cor_df_tab$sig_diff_84<-var_cor_df$sig_84ci
nice_var_cor_df_tab$Repeatability<-round(var_cor_df$repeatability,digits=2)
nice_var_cor_df_tab$Repeatability_CI<-paste0(round(var_cor_df$repeatability_lci,digits=2),", ",round(var_cor_df$repeatability_uci,digits=2))
mig_timing_df_uncertainty_sub5<-mig_timing_df_uncertainty
mig_timing_df_uncertainty_sub5[,3:16][which(abs(mig_timing_df_uncertainty[,3:16])>threshold,arr.ind=T)]<-NA
nice_var_cor_df_tab$Mean_uncertainty<-round(apply(abs(mig_timing_df_uncertainty_sub5[,milestone_names]),2,mean,na.rm=T),digits=2)
nice_var_cor_df_tab$SD_uncertainty<-round(apply(mig_timing_df_uncertainty_sub5[,milestone_names],2,sd,na.rm=T),digits=2)
nice_var_cor_df_tab$N<-unlist(lapply(apply(is.na(mig_timing_df[,milestone_names])==F,2,which),length))

#Relationship between between-individual variance and total variance
m1<-lm(total_var~PTT_var,
       data=var_cor_df[row.names(var_cor_df) %in% c("dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr"),]);summary(m1);confint(m1)
#Relationship between within-individual variance and total variance
m1<-lm(total_var~resid_var,
       data=var_cor_df[row.names(var_cor_df) %in% c("dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr"),]);summary(m1);confint(m1)


################################################################################################################
#3. Relationships among migratory milestones and with ecological/geographical covariates

# First, get subsequent departure date
dep_next_abs_list <- list()

for(i in 1:length(un_m_ptt)){
  temp_mtd <- mig_timing_df[mig_timing_df$ptt==un_m_ptt[i],]
  if(nrow(temp_mtd)==1) next
  temp_dn_abs_df <- data.frame(matrix(ncol=3,nrow=nrow(temp_mtd)-1))
  names(temp_dn_abs_df) <- c("ptt","year","dep_brgr_next")
  temp_dn_abs_df$ptt <- un_m_ptt[i]
  temp_dn_abs_df$year <- temp_mtd$year[1:(nrow(temp_mtd)-1)]
  temp_dn_abs_df$dep_brgr_next <- temp_mtd$dep_brgr[2:nrow(temp_mtd)]
  dep_next_abs_list[[i]] <- temp_dn_abs_df
}

# Flatten list of dataframes for each bird into single dataframe
dep_next_abs_df <- do.call(rbind,dep_next_abs_list)

# Now repeat for previous arrival date to breeding ground
arr_pre_abs_list <- list()

for(i in 1:length(un_m_ptt)){
  temp_mtd <- mig_timing_df[mig_timing_df$ptt==un_m_ptt[i],]
  if(nrow(temp_mtd)==1) next
  temp_ap_abs_df <- data.frame(matrix(ncol=3,nrow=nrow(temp_mtd)-1))
  names(temp_ap_abs_df) <- c("ptt","year","arr_brgr_pre")
  temp_ap_abs_df$ptt <- un_m_ptt[i]
  temp_ap_abs_df$year <- temp_mtd$year[2:(nrow(temp_mtd))]
  temp_ap_abs_df$arr_brgr_pre <- temp_mtd$arr_brgr[1:nrow(temp_mtd)-1]
  arr_pre_abs_list[[i]] <- temp_ap_abs_df
}

arr_pre_abs_df <- do.call(rbind,arr_pre_abs_list)

# Join with main dataset
mig_timing_df <- merge(mig_timing_df, dep_next_abs_df, all = T, by = c("ptt", "year"))
mig_timing_df <- merge(mig_timing_df, arr_pre_abs_df, all = T, by = c("ptt", "year"))

# What are each of the milestones? ## 6 for simplified causal model
milestone_names <- c("arr_brgr_pre","dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr","dep_brgr_next")

# Name other covariates used in model
cov_names <- c("br_long","br_lat","last_wa_sos_long","lat_last_euro_so_sbound","mig_dir","br_hab")

# Look at patterns of missing data
md.pattern(mig_timing_df[,c(milestone_names,cov_names)], plot = F)
aggr(mig_timing_df[,c(milestone_names,cov_names)], prop = F, numbers = T, sortVars = T, cex.axis = .6)

# Ensure at least one milestone is complete and remove any rows with missing categorical variables and departure from breeding
final_timing_df <- mig_timing_df[!apply(is.na(mig_timing_df[,c("dep_brgr","mig_dir")]), 1, any),] ## No missing data on breeding habitat

# How many individuals?
un_m_ptt <- unique(final_timing_df$ptt) ## 80. Lost 7 individuals

# Standardize migratory timing relative to individual mean for each milestone
## This is the within-individual anomaly
stdind_df <- final_timing_df
for(i in 1:length(un_m_ptt)){
  temp_ind_df <- stdind_df[stdind_df$ptt==un_m_ptt[i],milestone_names] 
  temp_ind_df_std <- temp_ind_df
  for(j in 1:ncol(temp_ind_df)){
    if(any(is.na(temp_ind_df[,j])==F)){
      if(length(which(is.na(temp_ind_df[,j])==F))<2) {temp_ind_df_std[,j] <- NA} else { ## If there is only one non-NA value for a given milestone, don't take the mean
        temp_ind_df_std[,j] <- temp_ind_df[,j]-mean(temp_ind_df[,j],na.rm=T)
      }
    } else {
      temp_ind_df_std[,j] <- NA
    }
  }
  stdind_df[stdind_df$ptt==un_m_ptt[i],milestone_names] <- temp_ind_df_std
}

# Remove cycles with no data for all milestones 
stdind_df <- stdind_df[!apply(is.na(stdind_df[,milestone_names])==T,1,all),] ## No missing data on breeding habitat

# Within-individual anomaly scaled to SD of 1 and center to a mean of 0, then melted 
stdind_scaled <- stdind_df
for(j in milestone_names){
  stdind_scaled[,j] <- as.vector(scale(stdind_scaled[,j], scale = T, center = T))
}

# How many birds had multiple mig cycles? ## 29
length(unique(stdind_scaled$ptt))

# Absolute dates scaled to SD of 1 and center to a mean of 0, then melted 
for(x in milestone_names){
  final_timing_df[,x] <- as.vector(scale(final_timing_df[,x]))
}

# What are the maximum values of scaled longitude of last West African stopover
## This is used in the model to impute missing data
max_long <- max(scale(final_timing_df$last_wa_sos_long), na.rm = T)
min_long <- min(scale(final_timing_df$last_wa_sos_long), na.rm = T) 
max_lat <- max(scale(final_timing_df$lat_last_euro_so_sbound), na.rm = T) 
min_lat <- min(scale(final_timing_df$lat_last_euro_so_sbound), na.rm = T) 

###################################################
#### 1. Model of absolute timing of milestones ####
###################################################

# Formulate model ####
# We want to evaluate direct and indirect effects of several steps of a birds migratory cycle 
{
  sink("sem_tim.jags")		
  cat("
model {

  # Random effects #

  for (j in 1:Nyr) { # Annual effect
    a1[j] ~ dnorm(0, tau_a1) 
    a2[j] ~ dnorm(0, tau_a2) 
    a3[j] ~ dnorm(0, tau_a3) 
    a4[j] ~ dnorm(0, tau_a4) 
    a5[j] ~ dnorm(0, tau_a5) 
    a6[j] ~ dnorm(0, tau_a6) 

  }
  
    tau_a1 <- pow(sd_a1, -2); sd_a1 ~ dunif(0, 10)
    tau_a2 <- pow(sd_a2, -2); sd_a2 ~ dunif(0, 10)
    tau_a3 <- pow(sd_a3, -2); sd_a3 ~ dunif(0, 10)
    tau_a4 <- pow(sd_a4, -2); sd_a4 ~ dunif(0, 10)
    tau_a5 <- pow(sd_a5, -2); sd_a5 ~ dunif(0, 10)
    tau_a6 <- pow(sd_a6, -2); sd_a6 ~ dunif(0, 10)

  for (h in 1:Nreg) { # Regional effect
    c1[h] ~ dnorm(0, tau_c1)
    c2[h] ~ dnorm(0, tau_c2)
    c3[h] ~ dnorm(0, tau_c3)
    c4[h] ~ dnorm(0, tau_c4)
    c5[h] ~ dnorm(0, tau_c5)
    c6[h] ~ dnorm(0, tau_c6)

  }

  for (k in 1:Nind) { # Individual effects nested within region
    d1[k] ~ dnorm(c1[region[ptt[k]]], tau_d1)
    d2[k] ~ dnorm(c2[region[ptt[k]]], tau_d2) 
    d3[k] ~ dnorm(c3[region[ptt[k]]], tau_d3) 
    d4[k] ~ dnorm(c4[region[ptt[k]]], tau_d4) 
    d5[k] ~ dnorm(c5[region[ptt[k]]], tau_d5) 
    d6[k] ~ dnorm(c6[region[ptt[k]]], tau_d6) 

  }
    
    tau_d1 <- pow(sd_d1, -2); sd_d1 ~ dunif(0, 10)
    tau_c1 <-  pow(sd_c1, -2); sd_c1 ~ dunif(0, 10)
    tau_d2 <- pow(sd_d2, -2); sd_d2 ~ dunif(0, 10)
    tau_c2 <-  pow(sd_c2, -2); sd_c2 ~ dunif(0, 10)
    tau_d3 <- pow(sd_d3, -2); sd_d3 ~ dunif(0, 10)
    tau_c3 <-  pow(sd_c3, -2); sd_c3 ~ dunif(0, 10)
    tau_d4 <- pow(sd_d4, -2); sd_d4 ~ dunif(0, 10)
    tau_c4 <-  pow(sd_c4, -2); sd_c4 ~ dunif(0, 10)
    tau_d5 <- pow(sd_d5, -2); sd_d5 ~ dunif(0, 10)
    tau_c5 <-  pow(sd_c5, -2); sd_c5 ~ dunif(0, 10)
    tau_d6 <- pow(sd_d6, -2); sd_d6 ~ dunif(0, 10)
    tau_c6 <-  pow(sd_c6, -2); sd_c6 ~ dunif(0, 10)
    
  # The model #
  
  for (i in 1:N) {
    fin_sah_sb_exp[i] = b1.0 + b1.1*dep_brgr[i] + b1.2*lat_last_euro_so_sbound[i] + b1.3*br_hab[i] + b1.4*mig_dir[i] + a1[year[i]] +
                            d1[ptt[i]]
    fin_sah_sb[i] ~ dnorm(fin_sah_sb_exp[i], fin_sah_sb_tau)
    
    arr_wgr_a_exp[i] =  b2.0 + b2.1*fin_sah_sb[i] + b2.2*br_hab[i] + b2.3*mig_dir[i] + a2[year[i]] + d2[ptt[i]]
    arr_wgr_a[i] ~ dnorm(arr_wgr_a_exp[i], arr_wgr_a_tau)
    
    dep_wgr_a_exp[i] = b3.0 + b3.1*arr_wgr_a[i] + b3.2*br_hab[i] + b3.3*mig_dir[i] + a3[year[i]] + d3[ptt[i]]
    dep_wgr_a[i] ~ dnorm(dep_wgr_a_exp[i], dep_wgr_a_tau)
    
    dep_wa_exp[i] = b4.0 + b4.1*dep_wgr_a[i] + b4.2*last_wa_sos_long[i] + b4.3*br_hab[i] + b4.4*mig_dir[i] +  
                            b4.5*fin_sah_sb[i] + a4[year[i]] + d4[ptt[i]]
    dep_wa[i] ~ dnorm(dep_wa_exp[i], dep_wa_tau)
    
    arr_brgr_exp[i] = b5.0 + b5.1*dep_brgr[i] + b5.2*fin_sah_sb[i] + b5.3*arr_wgr_a[i] + b5.4*dep_wgr_a[i] + b5.5*dep_wa[i] + 
                            b5.6*br_long[i] + b5.7*br_lat[i] + b5.8*br_hab[i] + b5.9*mig_dir[i] + a5[year[i]] + d5[ptt[i]]
    arr_brgr[i] ~ dnorm(arr_brgr_exp[i], arr_brgr_tau)

    dep_brgr_next_exp[i] =  b6.0 + b6.1*arr_brgr[i] + b6.2*br_long[i] + b6.3*br_lat[i] + b6.4*br_hab[i] + b6.5*mig_dir[i] + 
                            a6[year[i]] + d6[ptt[i]]
    dep_brgr_next[i] ~ dnorm(dep_brgr_next_exp[i], dep_brgr_next_tau)  
    
    last_wa_sos_long[i] ~ dnorm(mu_sos_long, tau_sos_long)
    lat_last_euro_so_sbound[i] ~ dnorm(mu_euro_lat, tau_euro_lat)

    }
    
  # Priors #
  
  b1.0 ~ dnorm(0,0.00001); b1.1 ~ dnorm(0,0.00001); b1.2 ~ dnorm(0,0.00001); b1.3 ~ dnorm(0,0.00001); b1.4 ~ dnorm(0,0.00001)
  fin_sah_sb_tau = pow(fin_sah_sb_sigma, -2)
  fin_sah_sb_sigma ~ dunif(0,10)
  
  b2.0 ~ dnorm(0,0.00001); b2.1 ~ dnorm(0,0.00001); b2.2 ~ dnorm(0,0.00001); b2.3 ~ dnorm(0,0.00001)
  arr_wgr_a_tau = pow(arr_wgr_a_sigma, -2)
  arr_wgr_a_sigma  ~ dunif(0,10)
  
  b3.0 ~ dnorm(0,0.00001); b3.1 ~ dnorm(0,0.00001); b3.2 ~ dnorm(0,0.00001); b3.3 ~ dnorm(0,0.00001)
  dep_wgr_a_tau = pow(dep_wgr_a_sigma, -2)
  dep_wgr_a_sigma ~ dunif(0,10)
  
  b4.0 ~ dnorm(0,0.00001); b4.1 ~ dnorm(0,0.00001); b4.2 ~ dnorm(0,0.00001); b4.3 ~ dnorm(0,0.00001); b4.4 ~ dnorm(0,0.00001);
  b4.5 ~ dnorm(0,0.00001); 
  dep_wa_tau = pow(dep_wa_sigma, -2)
  dep_wa_sigma ~ dunif(0,10)
  
  b5.0 ~ dnorm(0,0.00001); b5.1 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001); b5.3 ~ dnorm(0,0.00001); b5.4 ~ dnorm(0,0.00001);
  b5.5 ~ dnorm(0,0.00001); b5.6 ~ dnorm(0,0.00001); b5.7 ~ dnorm(0,0.00001); b5.8 ~ dnorm(0,0.00001); b5.9 ~ dnorm(0,0.00001);
  arr_brgr_tau = pow(arr_brgr_sigma, -2)
  arr_brgr_sigma ~ dunif(0,10)
  
  b6.0 ~ dnorm(0,0.00001); b6.1 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); b6.3 ~ dnorm(0,0.00001); b6.4 ~ dnorm(0,0.00001); 
  b6.5 ~ dnorm(0,0.00001);
  dep_brgr_next_tau = pow(dep_brgr_next_sigma, -2)
  dep_brgr_next_sigma ~ dunif(0,10) 

  mu_sos_long ~ dunif(min_long,max_long)
  tau_sos_long ~ dgamma(.01, .01)
  mu_euro_lat ~ dunif(min_lat,max_lat)
  tau_euro_lat ~ dgamma(.01, .01)


}", 
      
      fill=T)
  
  sink()
}

## Parameters monitored
params_tim <- c('b1.1','b1.2','b1.3','b1.4',
                'b2.1','b2.2','b2.3',
                'b3.1','b3.2','b3.3', 
                'b4.1','b4.2','b4.3','b4.4','b4.5',
                'b5.1','b5.2','b5.3','b5.4','b5.5','b5.6','b5.7','b5.8','b5.9',
                'b6.1','b6.2','b6.3','b6.4','b6.5',
                'last_wa_sos_long','lat_last_euro_so_sbound')

tim_vars <- data.frame(rbind(c('Departure from breeding grounds', 'Southbound Sahara crossing'),
                             c('Last European stopover','Southbound Sahara crossing'),
                             c('Breeding habitat','Southbound Sahara crossing'),
                             c('Migratory direction','Southbound Sahara crossing'),
                             c('Southbound Sahara crossing', 'Arrival to wintering grounds'),
                             c('Breeding habitat', 'Arrival to wintering grounds'),
                             c('Migratory direction', 'Arrival to wintering grounds'),
                             c('Arrival to wintering grounds', 'Departure from wintering grounds'),
                             c('Breeding habitat', 'Departure from wintering grounds'),
                             c('Migratory direction', 'Departure from wintering grounds'),
                             c('Departure from wintering grounds', 'Departure from West Africa'),
                             c('Pre-sahara crossing longitude', 'Departure from West Africa'),
                             c('Breeding habitat', 'Departure from West Africa'),
                             c('Migratory direction', 'Departure from West Africa'),
                             c('Southbound Sahara crossing', 'Departure from West Africa'),
                             c('Departure from breeding grounds', 'Arrival to breeding grounds'),
                             c('Southbound Sahara crossing', 'Arrival to breeding grounds'),
                             c('Arrival to wintering grounds', 'Arrival to breeding grounds'),
                             c('Departure from wintering grounds', 'Arrival to breeding grounds'),
                             c('Departure from West Africa', 'Arrival to breeding grounds'),
                             c('Breeding longitude', 'Arrival to breeding grounds'),
                             c('Breeding latitude', 'Arrival to breeding grounds'),
                             c('Breeding habitat', 'Arrival to breeding grounds'),
                             c('Migratory direction', 'Arrival to breeding grounds'),
                             c('Arrival to breeding grounds', 'Departure from breeding grounds'),
                             c('Breeding longitude', 'Departure from breeding grounds'),
                             c('Breeding latitude', 'Departure from breeding grounds'),
                             c('Breeding habitat', 'Departure from breeding grounds'),
                             c('Migratory direction', 'Departure from breeding grounds')))
names(tim_vars) <- c('Effect of', 'on')

# Absolute timing date ####

# Create data list for jags, which contains anything we are reading into jags as data.
data1 <- list(ptt = as.numeric(as.factor(final_timing_df$ptt)), 
              region = as.numeric(as.factor(final_timing_df$region)),
              year = final_timing_df$year - min(final_timing_df$year) + 1,
              br_hab = as.vector(scale(ifelse(final_timing_df$br_hab == "Lowland",0,1))),
              mig_dir = as.vector(scale(ifelse(final_timing_df$mig_dir == "SE",0,1))),
              br_long = as.vector(scale(final_timing_df$br_long)),
              br_lat = as.vector(scale(final_timing_df$br_lat)), 
              lat_last_euro_so_sbound = as.vector(scale(final_timing_df$lat_last_euro_so_sbound)),
              last_wa_sos_long = as.vector(scale(final_timing_df$last_wa_sos_long)), 
              arr_brgr = final_timing_df$arr_brgr, 
              dep_brgr = final_timing_df$dep_brgr,
              fin_sah_sb = final_timing_df$fin_sah_sb,
              arr_wgr_a = final_timing_df$arr_wgr_a,
              dep_wgr_a = final_timing_df$dep_wgr_a,
              dep_wa = final_timing_df$dep_wa,
              dep_brgr_next = final_timing_df$dep_brgr_next,
              Nyr = length(unique(final_timing_df$year)),
              Nind = length(unique(final_timing_df$ptt)),
              Nreg = length(unique(final_timing_df$region)),
              N = nrow(final_timing_df),
              max_long = max_long,
              min_long = min_long,
              min_lat = min_lat,
              max_lat = max_lat
)

# Run model
(sem_abs <-jags.parallel(data = data1,
                         inits = NULL, 
                         parameters = params_tim, 
                         "sem_tim.jags", 
                         n.chains = 5, 
                         n.thin   = 500, 
                         n.iter   = 1000000,
                         n.burnin = 100000))  
# Check Rhat values, are they close to 1?

# Visually inspect traceplots
traceplot(sem_abs, mfrow = c(2,2))

# Explore model output
param_abs_ind <- which(substr(names(sem_abs$BUGSoutput$sims.list),1,1)=='b')
abs_coefs <- data.frame(
  'mean'=sapply(sem_abs$BUGSoutput$sims.list[param_abs_ind],
                FUN=function(x) mean(x)), 
  'sd'=sapply(sem_abs$BUGSoutput$sims.list[param_abs_ind],
              FUN=function(x) sd(x)), 
  'median'=sapply(sem_abs$BUGSoutput$sims.list[param_abs_ind],
                  FUN=function(x) median(x)),
  '2.5%'=sapply(sem_abs$BUGSoutput$sims.list[param_abs_ind],
                FUN=function(x) quantile(x, 0.025)), 
  '97.5%'=sapply(sem_abs$BUGSoutput$sims.list[param_abs_ind],
                 FUN=function(x) quantile(x, 0.975)), 
  'P(|b|>0)'=sapply(sem_abs$BUGSoutput$sims.list[param_abs_ind],
                    FUN=function(x) max(mean(x>0), mean(x<0))), check.names=FALSE)

# Merge with coefficients
abs_coefs <- cbind(tim_vars, abs_coefs)
write.csv(print(abs_coefs, digits=2), "abs_timing_coeffs.csv", row.names = F)

# Check imputed values
imp_wa_long <- as.data.frame(sem_abs$BUGSoutput$sims.list$last_wa_sos_long)
imp_euro_lat <- as.data.frame(sem_abs$BUGSoutput$sims.list$lat)

# Within individual anomaly ####

# Create data list for jags, which contains anything we are reading into jags as data.
data2 <- list(ptt = as.numeric(as.factor(stdind_scaled$ptt)), 
              region = as.numeric(as.factor(stdind_scaled$region)),
              year = stdind_scaled$year - min(stdind_scaled$year) + 1,
              br_hab = as.vector(scale(ifelse(stdind_scaled$br_hab == "Lowland",0,1))),
              mig_dir = as.vector(scale(ifelse(stdind_scaled$mig_dir == "SE",0,1))),
              br_long = as.vector(scale(stdind_scaled$br_long)),
              br_lat = as.vector(scale(stdind_scaled$br_lat)), 
              lat_last_euro_so_sbound = as.vector(scale(stdind_scaled$lat_last_euro_so_sbound)),
              last_wa_sos_long = as.vector(scale(stdind_scaled$last_wa_sos_long)), 
              arr_brgr = stdind_scaled$arr_brgr, 
              dep_brgr = stdind_scaled$dep_brgr,
              fin_sah_sb = stdind_scaled$fin_sah_sb,
              arr_wgr_a = stdind_scaled$arr_wgr_a,
              dep_wgr_a = stdind_scaled$dep_wgr_a,
              dep_wa = stdind_scaled$dep_wa,
              dep_brgr_next = stdind_scaled$dep_brgr_next,
              Nyr = length(unique(stdind_scaled$year)),
              Nind = length(unique(stdind_scaled$ptt)),
              Nreg = length(unique(stdind_scaled$region)),
              N = nrow(stdind_scaled),
              max_long = max_long,
              min_long = min_long,
              max_lat = max_lat,
              min_lat = min_lat
)

(sem_anom <-jags.parallel(data = data2,
                          inits = NULL, 
                          parameters = params_tim, 
                          "sem_tim.jags", 
                          n.chains = 5, 
                          n.thin   = 500, 
                          n.iter   = 1000000,
                          n.burnin = 100000)) 

# Visually inspect traceplots
traceplot(sem_anom, mfrow = c(2,2))

# Explore model output
param_anom_ind <- which(substr(names(sem_anom$BUGSoutput$sims.list),1,1)=='b')
anom_coefs <- data.frame(
  'mean'=sapply(sem_anom$BUGSoutput$sims.list[param_anom_ind],
                FUN=function(x) mean(x)), 
  'sd'=sapply(sem_anom$BUGSoutput$sims.list[param_anom_ind],
              FUN=function(x) sd(x)), 
  'median'=sapply(sem_anom$BUGSoutput$sims.list[param_anom_ind],
                  FUN=function(x) median(x)),
  '2.5%'=sapply(sem_anom$BUGSoutput$sims.list[param_anom_ind],
                FUN=function(x) quantile(x, 0.025)), 
  '97.5%'=sapply(sem_anom$BUGSoutput$sims.list[param_anom_ind],
                 FUN=function(x) quantile(x, 0.975)), 
  'P(|b|>0)'=sapply(sem_anom$BUGSoutput$sims.list[param_anom_ind],
                    FUN=function(x) max(mean(x>0), mean(x<0))), check.names=FALSE)
anom_coefs <- cbind(tim_vars, anom_coefs)
write.csv(print(anom_coefs, digits=2), "anom_timing_coeffs.csv", row.names = F)

####################################
#### 2. Model of bird mortality ####
####################################

# Melted absolute timing dataset
abs_melted <- melt(final_timing_df[c("ptt","region","year","mig_dir","br_hab","last_milestone_pre_death","outcome",milestone_names)], id = c("ptt","region","year","mig_dir","br_hab","last_milestone_pre_death","outcome"))
abs_melted <- abs_melted[which(is.na(abs_melted$value)==F),]

# Remove repeated rows with previous breeding ground arrival and next departure 
abs_melted <- abs_melted[abs_melted$variable != "arr_brgr_pre" & abs_melted$variable != "dep_brgr_next",]

# If last milestone before death is one of the 'intermediate' milestones not analysed, change to previous milestone
unique(abs_melted$last_milestone_pre_death)
abs_melted$new_milestone_pre_death <- abs_melted$last_milestone_pre_death
abs_melted$new_milestone_pre_death[abs_melted$last_milestone_pre_death == "dep_UK"] <- "dep_brgr"
abs_melted$new_milestone_pre_death[abs_melted$last_milestone_pre_death == "arr_UK"] <- "dep_wa"
abs_melted$new_milestone_pre_death[abs_melted$last_milestone_pre_death == "arr_wgr_b"] <- "arr_wgr_a"
abs_melted$new_milestone_pre_death[abs_melted$last_milestone_pre_death == "dep_wgr_b"] <- "arr_wgr_a"

# Create binary variable indicating whether the bird died or survived that milestone
abs_melted$died <- "Survived"
abs_melted$died[which(abs_melted$new_milestone_pre_death==abs_melted$variable)] <- "Died"
abs_melted$died_bin <- 0 # Create binary variable
abs_melted$died_bin[abs_melted$died=="Died"] <- 1

# Create dummy variable for milestone
abs_melted_dum <- dummy_cols(abs_melted, select_columns = "variable")

# Remove uncertain outcomes
abs_melted_final <- abs_melted_dum[abs_melted_dum$outcome %in% c("UB","UC")==F ,]
table(abs_melted_final$died,abs_melted_final$variable) 
## 56 deaths

# Formulate model ####
# We want to evaluate direct and indirect effects of migratory timing on mortality of cuckoos
{
  sink("sem_mort.jags")		
  cat("
model {

  # Random effects #
  
  for (j in 1:Nyr) { # Annual effect
    a1[j] ~ dnorm(0, tau_a1) 
    a2[j] ~ dnorm(0, tau_a2) 
  }
  
    tau_a1 <- pow(sd_a1, -2)
    sd_a1 ~ dunif(0, 10)
    tau_a2 <- pow(sd_a2, -2)
    sd_a2 ~ dunif(0, 10)

  for (h in 1:Nreg) { # Regional effect
    c1[h] ~ dnorm(0, tau_c1)
    c2[h] ~ dnorm(0, tau_c2)
  }

  for (k in 1:Nind) { # Individual effects nested within region 
    d1[k] ~ dnorm(c1[region[ptt[k]]], tau_d1)
    d2[k] ~ dnorm(c2[region[ptt[k]]], tau_d2) 
  }
    
    tau_d1 <- pow(sd_d1, -2); sd_d1 ~ dunif(0, 10)
    tau_c1 <-  pow(sd_c1, -2); sd_c1 ~ dunif(0, 10)
    tau_d2 <- pow(sd_d2, -2); sd_d2 ~ dunif(0, 10)
    tau_c2 <-  pow(sd_c2, -2); sd_c2 ~ dunif(0, 10)
    
  # The model #
  
  for (i in 1:N) {
    
    # Departure from breeding ground is the intercept
      
    value_exp[i] =  b1.0 + b1.1*mig_dir[i] + b1.2*br_hab[i] + a1[year[i]] + d1[ptt[i]]
    value[i] ~ dnorm(value_exp[i], value_tau)
    
    logit(died_exp[i]) = b2.0 + b2.1*value[i] + b2.2*mig_dir[i] + b2.3*br_hab[i] + b2.4*value[i]*br_hab[i] + a2[year[i]] + d2[ptt[i]]
    died_bin[i] ~ dbern(died_exp[i])

  }
    
  # Priors #
  
  b1.0 ~ dnorm(0,0.00001); b1.1 ~ dnorm(0,0.00001); b1.2 ~ dnorm(0,0.00001); 
  value_tau = pow(value_sigma, -2)
  value_sigma ~ dunif(0,10)
  
  p2.0 ~ dbeta(1, 1); b2.0 <- logit(p2.0);
  b2.1 ~ dnorm(0,0.00001); b2.2 ~ dnorm(0,0.00001); b2.3 ~ dnorm(0,0.00001); b2.4 ~ dnorm(0,0.00001)
  

}", 
      
      fill=T)
  
  sink()
}

## Parameters monitored
params_mort <- c('b1.1','b1.2',
                 'b2.1','b2.2','b2.3','b2.4')

mort_vars <- data.frame(rbind(c('Migratory direction','Timing'),
                              c('Breeding habitat','Timing'),
                              c('Timing', 'Mortality'),
                              c('Migratory direction', 'Mortality'),
                              c('Breeding habitat', 'Mortality'),                             
                              c('Migratory direction x Timing', 'Mortality')))
names(mort_vars) <- c('Effect of', 'on')

# Absolute timing date ####

# Create seperate dataframes for each milestone

dep_brgr_df <- abs_melted_final[abs_melted_final$variable == "dep_brgr",]
fin_sah_df <- abs_melted_final[abs_melted_final$variable == "fin_sah_sb",]
arr_wgr_df <- abs_melted_final[abs_melted_final$variable == "arr_wgr_a",]
dep_wa_df <- abs_melted_final[abs_melted_final$variable == "dep_wa",]
arr_brgr_df <- abs_melted_final[abs_melted_final$variable == "arr_brgr",]

# Create data list
data3 <- list(ptt = as.numeric(as.factor(arr_brgr_df$ptt)), 
              region = as.numeric(as.factor(arr_brgr_df$region)),
              year = arr_brgr_df$year - min(arr_brgr_df$year) + 1,
              br_hab = ifelse(arr_brgr_df$br_hab == "Lowland",0,1),
              mig_dir = ifelse(arr_brgr_df$mig_dir == "SE",1,0),
              died_bin = arr_brgr_df$died_bin,
              value = arr_brgr_df$value,
              Nyr = max(arr_brgr_df$year),
              Nind = length(unique(arr_brgr_df$ptt)),
              Nreg = length(unique(arr_brgr_df$region)),
              N = nrow(arr_brgr_df)
)

# Run model 
(abs_mort <- jags.parallel(data = data3, 
                           inits = NULL, 
                           parameters = params_mort, 
                           "sem_mort.jags", 
                           n.chains = 5, 
                           n.thin   = 500, 
                           n.iter   = 1000000,
                           n.burnin = 100000)) 

# Visually inspect traceplots
traceplot(abs_mort, mfrow = c(2,2))

# Explore model output
param_abs_mort_ind <- which(substr(names(abs_mort$BUGSoutput$sims.list),1,1)=='b')
abs_mort_coefs <- data.frame(
  'mean'=sapply(abs_mort$BUGSoutput$sims.list[param_abs_mort_ind],
                FUN=function(x) mean(x)), 
  'sd'=sapply(abs_mort$BUGSoutput$sims.list[param_abs_mort_ind],
              FUN=function(x) sd(x)), 
  'median'=sapply(abs_mort$BUGSoutput$sims.list[param_abs_mort_ind],
                  FUN=function(x) median(x)),
  '2.5%'=sapply(abs_mort$BUGSoutput$sims.list[param_abs_mort_ind],
                FUN=function(x) quantile(x, 0.025)), 
  '97.5%'=sapply(abs_mort$BUGSoutput$sims.list[param_abs_mort_ind],
                 FUN=function(x) quantile(x, 0.975)), 
  'P(|b|>0)'=sapply(abs_mort$BUGSoutput$sims.list[param_abs_mort_ind],
                    FUN=function(x) max(mean(x>0), mean(x<0))), check.names=FALSE)

abs_mort_coefs <- cbind(mort_vars, abs_mort_coefs)
write.csv(print(abs_mort_coefs, digits=2), "mortality_coeffs.csv", row.names = F)

# Re-run mortality model combining migration stages that positively impacted breeding grounds
{
  sink("sem_mort2.jags")		
  cat("
model {

  # Random effects #
  
  for (j in 1:Nyr) { # Annual effect
    a1[j] ~ dnorm(0, tau_a1) 
    a2[j] ~ dnorm(0, tau_a2) 
  }
  
    tau_a1 <- pow(sd_a1, -2)
    sd_a1 ~ dunif(0, 10)
    tau_a2 <- pow(sd_a2, -2)
    sd_a2 ~ dunif(0, 10)

  for (h in 1:Nreg) { # Regional effect
    c1[h] ~ dnorm(0, tau_c1)
    c2[h] ~ dnorm(0, tau_c2)
  }

  for (k in 1:Nind) { # Individual effects nested within region 
    d1[k] ~ dnorm(c1[region[ptt[k]]], tau_d1)
    d2[k] ~ dnorm(c2[region[ptt[k]]], tau_d2) 
  }
    
    tau_d1 <- pow(sd_d1, -2); sd_d1 ~ dunif(0, 10)
    tau_c1 <-  pow(sd_c1, -2); sd_c1 ~ dunif(0, 10)
    tau_d2 <- pow(sd_d2, -2); sd_d2 ~ dunif(0, 10)
    tau_c2 <-  pow(sd_c2, -2); sd_c2 ~ dunif(0, 10)
    
  # The model #
  
  for (i in 1:N) {
    
    # Departure from breeding ground is the intercept
      
    value_exp[i] =  b1.0 + b1.1*mig_dir[i] + b1.2*br_hab[i] + a1[year[i]] + d1[ptt[i]]
    value[i] ~ dnorm(value_exp[i], value_tau)
    
    logit(died_exp[i]) = b2.0 + b2.1*value[i] + b2.2*mig_dir[i] + b2.3*br_hab[i] + b2.4*value[i]*mig_dir[i] + b2.5*variable_dep_wa[i] + 
    b2.6*variable_arr_brgr[i] + a2[year[i]] + d2[ptt[i]]
    died_bin[i] ~ dbern(died_exp[i])

  }
    
  # Priors #
  
  b1.0 ~ dnorm(0,0.00001); b1.1 ~ dnorm(0,0.00001); b1.2 ~ dnorm(0,0.00001); 
  value_tau = pow(value_sigma, -2)
  value_sigma ~ dunif(0,10)
  
  p2.0 ~ dbeta(1, 1); b2.0 <- logit(p2.0);
  b2.1 ~ dnorm(0,0.00001); b2.2 ~ dnorm(0,0.00001); b2.3 ~ dnorm(0,0.00001); b2.4 ~ dnorm(0,0.00001);  b2.5 ~ dnorm(0,0.00001); 
  b2.6 ~ dnorm(0,0.00001)
  

}", 
      
      fill=T)
  
  sink()
}

## Parameters monitored
params_mort2 <- c('b1.1','b1.2',
                  'b2.1','b2.2','b2.3','b2.4','b2.5','b2.6')

mort_vars2 <- data.frame(rbind(c('Migratory direction','Timing'),
                               c('Breeding habitat','Timing'),
                               c('Timing', 'Mortality'),
                               c('Migratory direction', 'Mortality'),
                               c('Breeding habitat', 'Mortality'),                             
                               c('Breeding habitat x Timing', 'Mortality'),
                               c('Departure West Africa', 'Mortality'),
                               c('Depterture breeding grounds', 'Mortality')))
names(mort_vars2) <- c('Effect of', 'on')

# Create dataframe with relevant milestones
comb_df <- abs_melted_final[abs_melted_final$variable == "fin_sah_sb" | abs_melted_final$variable == "dep_wa" | 
                              abs_melted_final$variable == "arr_brgr",]

# Create data list
data4 <- list(ptt = as.numeric(as.factor(comb_df$ptt)), 
              region = as.numeric(as.factor(comb_df$region)),
              year = comb_df$year - min(comb_df$year) + 1,
              br_hab = ifelse(comb_df$br_hab == "Lowland",0,1),
              mig_dir = ifelse(comb_df$mig_dir == "SE",1,0),
              died_bin = comb_df$died_bin,
              value = comb_df$value,
              variable_dep_wa = comb_df$variable_dep_wa,
              variable_arr_brgr = comb_df$variable_arr_brgr,
              Nyr = max(comb_df$year),
              Nind = length(unique(comb_df$ptt)),
              Nreg = length(unique(comb_df$region)),
              N = nrow(comb_df))

# Run model 
(abs_mort2 <- jags.parallel(data = data4, 
                            inits = NULL, 
                            parameters = params_mort2, 
                            "sem_mort2.jags", 
                            n.chains = 5, 
                            n.thin   = 500, 
                            n.iter   = 1000000,
                            n.burnin = 100000)) 

# Visually inspect traceplots
traceplot(abs_mort2, mfrow = c(2,2))

# Explore model output
param_abs_mort_ind2 <- which(substr(names(abs_mort2$BUGSoutput$sims.list),1,1)=='b')
abs_mort_coefs2 <- data.frame(
  'mean'=sapply(abs_mort2$BUGSoutput$sims.list[param_abs_mort_ind2],
                FUN=function(x) mean(x)), 
  'sd'=sapply(abs_mort2$BUGSoutput$sims.list[param_abs_mort_ind2],
              FUN=function(x) sd(x)), 
  'median'=sapply(abs_mort2$BUGSoutput$sims.list[param_abs_mort_ind2],
                  FUN=function(x) median(x)),
  '2.5%'=sapply(abs_mort2$BUGSoutput$sims.list[param_abs_mort_ind2],
                FUN=function(x) quantile(x, 0.025)), 
  '97.5%'=sapply(abs_mort2$BUGSoutput$sims.list[param_abs_mort_ind2],
                 FUN=function(x) quantile(x, 0.975)), 
  'P(|b|>0)'=sapply(abs_mort2$BUGSoutput$sims.list[param_abs_mort_ind2],
                    FUN=function(x) max(mean(x>0), mean(x<0))), check.names=FALSE)

abs_mort_coefs2 <- cbind(mort_vars2, abs_mort_coefs2)
write.csv(print(abs_mort_coefs2, digits=2), "mortality_coeffs2c.csv", row.names = F)


################################################################################################################
#4. Effects of age on migratory timing

age_coeff_df<-data.frame(matrix(ncol=8,nrow=length(milestone_names)))
row.names(age_coeff_df)<-milestone_names
names(age_coeff_df)<-c("first_ad_cyc_est","first_ad_cyc_lci","first_ad_cyc_uci",
                       "later_cyc_est","later_cyc_lci","later_cyc_uci","sig.05","sig.1")

for(i in 1:nrow(age_coeff_df)){
  print(i)
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(known_age_binary)+(1|ptt)"))
  m1<-brm(temp_form_age,
          silent=1,
          data=mig_timing_df_stdind,
          iter=4000,
          control = list(adapt_delta = 0.99999999999, #Close to 1 = short step length to avoid divergent transitions
                         max_treedepth=20));summary(m1) #Check Rhat<1.1 for all parameters
  sm1<-summary(m1)
  age_coeff_df$first_ad_cyc_est[i]<-sm1$fixed[row.names(sm1$fixed)=="as.factorknown_age_binary0",]$Estimate
  age_coeff_df$first_ad_cyc_lci[i]<-sm1$fixed[row.names(sm1$fixed)=="as.factorknown_age_binary0",]$'l-95% CI'
  age_coeff_df$first_ad_cyc_uci[i]<-sm1$fixed[row.names(sm1$fixed)=="as.factorknown_age_binary0",]$'u-95% CI'
  age_coeff_df$later_cyc_est[i]<-sm1$fixed[row.names(sm1$fixed)=="as.factorknown_age_binary1",]$Estimate
  age_coeff_df$later_cyc_lci[i]<-sm1$fixed[row.names(sm1$fixed)=="as.factorknown_age_binary1",]$'l-95% CI'
  age_coeff_df$later_cyc_uci[i]<-sm1$fixed[row.names(sm1$fixed)=="as.factorknown_age_binary1",]$'u-95% CI'
  
  temp_form_age_sign<-as.formula(paste0(milestone_names[i],"~as.factor(known_age_binary)+(1|ptt)")) #Remove intercept to test significance of difference
  m2<-brm(temp_form_age_sign,
          silent=1,
          data=mig_timing_df_stdind,
          iter=4000,
          control = list(adapt_delta = 0.99999999999, #Close to 1 = short step length to avoid divergent transitions
                         max_treedepth=20));summary(m1) #Check Rhat<1.1 for all parameters
  later_draws<-as_draws_df(m2,variable="b_as.factorknown_age_binary1")$b_as.factorknown_age_binary1
  temp_lci_95<-quantile(later_draws,0.025)
  temp_uci_95<-quantile(later_draws,0.975)
  temp_lci_90<-quantile(later_draws,0.05)
  temp_uci_90<-quantile(later_draws,0.95)
  age_coeff_df$sig.05[i]<-(temp_lci_95<0 & temp_uci_95>0)==F
  age_coeff_df$sig.1[i]<-(temp_lci_90<0 & temp_uci_90>0)==F
}

age_coeff_df_newmilestones<-age_coeff_df[row.names(age_coeff_df) %in% c("dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr"),]
age_coeff_df_newmilestones$nice_names<-character(nrow(age_coeff_df_newmilestones))
age_coeff_df_newmilestones$nice_names[row.names(age_coeff_df_newmilestones)=="dep_brgr"]<-"Depart breeding\ngrounds"
age_coeff_df_newmilestones$nice_names[row.names(age_coeff_df_newmilestones)=="fin_sah_sb"]<-"Completion Sahara\ncrossing"
age_coeff_df_newmilestones$nice_names[row.names(age_coeff_df_newmilestones)=="arr_wgr_a"]<-"Arrive wintering\ngrounds"
age_coeff_df_newmilestones$nice_names[row.names(age_coeff_df_newmilestones)=="dep_wgr_a"]<-"Depart wintering\ngrounds"
age_coeff_df_newmilestones$nice_names[row.names(age_coeff_df_newmilestones)=="dep_wa"]<-"Depart West\nAfrica"
age_coeff_df_newmilestones$nice_names[row.names(age_coeff_df_newmilestones)=="arr_brgr"]<-"Arrive breeding\ngrounds"
age_coeff_df_newmilestones$nice_names<-factor(age_coeff_df_newmilestones$nice_names,
                                              levels=c("Depart breeding\ngrounds",
                                                       "Completion Sahara\ncrossing",
                                                       "Arrive wintering\ngrounds",
                                                       "Depart wintering\ngrounds",
                                                       "Depart West\nAfrica",
                                                       "Arrive breeding\ngrounds"))
milestone_names_new<-row.names(age_coeff_df_newmilestones)

par(mar=c(8.6, 4.1, 1.1, 1.1))
plot(1:nrow(age_coeff_df_newmilestones),
     age_coeff_df_newmilestones$first_ad_cyc_est,
     xlim=c(0.5,nrow(age_coeff_df_newmilestones)+0.5),
     ylim=c(min(c(age_coeff_df_newmilestones$first_ad_cyc_lci,age_coeff_df_newmilestones$first_ad_cyc_lci)),
            max(c(age_coeff_df_newmilestones$first_ad_cyc_uci,age_coeff_df_newmilestones$first_ad_cyc_uci))),
     col="white",xaxt="n",xlab="",pch=19,ylab="Within-individual anomaly",main="")
points((1:length(milestone_names_new))-0.1,age_coeff_df_newmilestones$first_ad_cyc_est,pch=19,col="cornflowerblue")
arrows((1:nrow(age_coeff_df_newmilestones))-0.1,age_coeff_df_newmilestones$first_ad_cyc_lci,(1:nrow(age_coeff_df_newmilestones))-0.1,age_coeff_df_newmilestones$first_ad_cyc_uci,length=0,col="cornflowerblue")
points((1:length(milestone_names_new))+0.1,age_coeff_df_newmilestones$later_cyc_est,pch=19,col="chocolate")
arrows((1:nrow(age_coeff_df_newmilestones))+0.1,age_coeff_df_newmilestones$later_cyc_lci,(1:nrow(age_coeff_df_newmilestones))+0.1,age_coeff_df_newmilestones$later_cyc_uci,length=0,col="chocolate")
age_sig1_ind<-which(age_coeff_df_newmilestones$sig.1==T & age_coeff_df_newmilestones$sig.05==F)
arrows(age_sig1_ind-0.3,20,age_sig1_ind+0.3,20,length=0)
arrows(age_sig1_ind-0.3,20,age_sig1_ind-0.3,19,length=0)
arrows(age_sig1_ind+0.3,20,age_sig1_ind+0.3,19,length=0)
text(age_sig1_ind,22,"*",cex=2)
axis(1, at=1:length(milestone_names_new),labels=age_coeff_df_newmilestones$nice_names,las=2)
legend("topleft",c("Cycle after first breeding season","All subsequent cycles"),col=c("cornflowerblue","chocolate"),pch=c(19,19),bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))

#Cycles since fledging
par(mfrow=c(2,3))
for(i in 1:length(milestone_names_new)){
  print(i)
  if(milestone_names_new[i]=="dep_wgr_a") next #Not enough data for this milestone: trying to estimate an individual random effect and four parameters from 6 data
  temp_form_age<-as.formula(paste0(milestone_names_new[i],"~0+as.factor(known_age)+(1|ptt)"))
  m1<-brm(temp_form_age,
          silent=1,
          data=mig_timing_df_stdind[which(is.na(mig_timing_df_stdind[,milestone_names_new[i]])==F),],
          iter=4000,
          control = list(adapt_delta = 0.99999999999, #Close to 1 = short step length to avoid divergent transitions
                         max_treedepth=20));summary(m1) #Check Rhat<1.1 for all parameters
  sm1_f<-summary(m1)$fixed
  #Calculate 84% CIs for each estimate
  sm1_f$sig_84<-sm1_f$uci84<-sm1_f$lci84<-NA
  for(j in 1:nrow(sm1_f)){
    temp_factor_name<-paste0("b_",row.names(sm1_f[j,]))
    temp_draws<-data.frame(as_draws_df(m1,temp_factor_name))[,temp_factor_name]
    sm1_f$lci84[j]<-quantile(temp_draws,0.08)
    sm1_f$uci84[j]<-quantile(temp_draws,0.92)
  }
  #Does the estimate overlap with the 84% CIs for any other parameter?
  for(j in 1:nrow(sm1_f)){
    sm1_f$sig_84[j]<-any(sm1_f$Estimate[j]>sm1_f$lci84[-j] & sm1_f$Estimate[j]<sm1_f$uci84[-j])==F
  }
  
  temp_lci<-sm1_f$`l-95% CI`
  temp_uci<-sm1_f$`u-95% CI`
  temp_ages<-as.numeric(substr(row.names(sm1_f),nchar(row.names(sm1_f)),nchar(row.names(sm1_f))))
  plot(temp_ages,sm1_f$Estimate,
       main=age_coeff_df_newmilestones$nice_names[i],
       col="blue",
       xlab="Cycles since first breeding season",
       xlim=c(1,6),
       ylab="Within-individual anomaly",
       ylim=range(c(temp_lci,temp_uci)))
  arrows(temp_ages,temp_lci,temp_ages,temp_uci,length=0,col="blue")
  abline(h=0,lty=3)
  points(temp_ages[sm1_f$sig_84],sm1_f$Estimate[sm1_f$sig_84],pch=19,col="blue")
}
par(mfrow=c(1,1))

#Fig A9 but for non-known age birds
par(mfrow=c(2,3))
for(i in 1:length(milestone_names_new)){
  print(i)
  if(milestone_names_new[i]=="dep_wgr_a") next #Not enough data for this milestone: trying to estimate an individual random effect and four parameters from 6 data
  temp_form_age<-as.formula(paste0(milestone_names_new[i],"~0+as.factor(cycles_since_tagging_4s)+(1|ptt)"))
  m1<-brm(temp_form_age,
          silent=1,
          data=mig_timing_df_stdind[which(is.na(mig_timing_df_stdind[,milestone_names_new[i]])==F),],
          iter=4000,
          control = list(adapt_delta = 0.99999999999, #Close to 1 = short step length to avoid divergent transitions
                         max_treedepth=20));summary(m1) #Check Rhat<1.1 for all parameters
  sm1_f<-summary(m1)$fixed
  #Calculate 84% CIs for each estimate
  sm1_f$sig_84<-sm1_f$uci84<-sm1_f$lci84<-NA
  for(j in 1:nrow(sm1_f)){
    temp_factor_name<-paste0("b_",row.names(sm1_f[j,]))
    temp_draws<-data.frame(as_draws_df(m1,temp_factor_name))[,temp_factor_name]
    sm1_f$lci84[j]<-quantile(temp_draws,0.08)
    sm1_f$uci84[j]<-quantile(temp_draws,0.92)
  }
  #Does the estimate overlap with the 84% CIs for any other parameter?
  for(j in 1:nrow(sm1_f)){
    sm1_f$sig_84[j]<-any(sm1_f$Estimate[j]>sm1_f$lci84[-j] & sm1_f$Estimate[j]<sm1_f$uci84[-j])==F
  }
  
  temp_lci<-sm1_f$`l-95% CI`
  temp_uci<-sm1_f$`u-95% CI`
  temp_ages<-as.numeric(substr(row.names(sm1_f),nchar(row.names(sm1_f)),nchar(row.names(sm1_f))))
  plot(temp_ages,sm1_f$Estimate,
       main=age_coeff_df_newmilestones$nice_names[i],
       col="blue",
       xlab="Cycles since tagging",
       xlim=c(0,4),
       ylab="Within-individual anomaly",
       ylim=range(c(temp_lci,temp_uci)))
  arrows(temp_ages,temp_lci,temp_ages,temp_uci,length=0,col="blue")
  abline(h=0,lty=3)
  points(temp_ages[sm1_f$sig_84],sm1_f$Estimate[sm1_f$sig_84],pch=19,col="blue")
}
par(mfrow=c(1,1))

################################################################################################################
#5. Effect of outcome status on migratory timing
mig_timing_df_scaled$last_milestone_pre_death_newmilestones<-mig_timing_df_scaled$last_milestone_pre_death
mig_timing_df_scaled$last_milestone_pre_death_newmilestones[which(mig_timing_df_scaled$last_milestone_pre_death=="dep_UK")]<-"dep_brgr"
mig_timing_df_scaled$last_milestone_pre_death_newmilestones[which(mig_timing_df_scaled$last_milestone_pre_death=="arr_wgr_b")]<-"arr_wgr_a"
mig_timing_df_scaled$last_milestone_pre_death_newmilestones[which(mig_timing_df_scaled$last_milestone_pre_death=="dep_wgr_b")]<-"arr_wgr_a"
mig_timing_df_scaled$last_milestone_pre_death_newmilestones[which(mig_timing_df_scaled$last_milestone_pre_death=="arr_UK")]<-"dep_wa"

mtd_melted_scaled_newmilestones<-melt(mig_timing_df_scaled[,c("ptt","year","mig_dir","br_hab","region","known_age_binary","last_milestone_pre_death_newmilestones","outcome",milestone_names)],id=c("ptt","year","mig_dir","br_hab","region","known_age_binary","last_milestone_pre_death_newmilestones","outcome"))
mtd_melted_scaled_newmilestones<-mtd_melted_scaled_newmilestones[which(is.na(mtd_melted_scaled_newmilestones$value)==F),]
mtd_melted_scaled_newmilestones$died<-"Survived"
mtd_melted_scaled_newmilestones$died[which(mtd_melted_scaled_newmilestones$last_milestone_pre_death_newmilestones==mtd_melted_scaled_newmilestones$variable)]<-"Died"
mtd_melted_scaled_newmilestones$died<-relevel(as.factor(mtd_melted_scaled_newmilestones$died),ref = "Survived") #Put Survived first - more data for intercept
mtd_melted_scaled_newmilestones$died_bin<-0
mtd_melted_scaled_newmilestones$died_bin[mtd_melted_scaled_newmilestones$died=="Died"]<-1
mtd_melted_scaled_newmilestones$season<-NA
mtd_melted_scaled_newmilestones[mtd_melted_scaled_newmilestones$variable %in% c("dep_brgr","fin_sah_sb","arr_wgr_a"),]$season<-"Autumn"
mtd_melted_scaled_newmilestones[mtd_melted_scaled_newmilestones$variable %in% c("dep_wgr_a","dep_wa","arr_brgr"),]$season<-"Spring"

#Plot by outcome
core_milestone_names<-c("dep_brgr","fin_sah_sb","arr_wgr_a","dep_wgr_a","dep_wa","arr_brgr")
lcmn<-length(core_milestone_names)
died_m_x_seq<-seq(1,6*lcmn,by=6)
died_ua_x_seq<-died_m_x_seq+1
died_ub_x_seq<-died_m_x_seq+2
died_uc_x_seq<-died_m_x_seq+3
surv_x_seq<-died_m_x_seq+4

par(mar=c(8.6, 4.1, 1.1, 1.1))
plot(rep(0,lcmn*6),ylim=c(min(mtd_melted_scaled_newmilestones$value),max(mtd_melted_scaled_newmilestones$value)+2),col="white",xaxt="n",xlab="",ylab="Date scaled to mean 0 and SD 1")
for(i in 1:lcmn){
  temp_data<-mtd_melted_scaled_newmilestones[mtd_melted_scaled_newmilestones$variable==core_milestone_names[i],]
  points(rep(died_m_x_seq[i],nrow(temp_data[temp_data$died=="Died" & temp_data$outcome=="M",])),
         temp_data[temp_data$died=="Died" & temp_data$outcome=="M",]$value)
  points(rep(died_ua_x_seq[i],nrow(temp_data[temp_data$died=="Died" & temp_data$outcome=="UA",])),
         temp_data[temp_data$died=="Died" & temp_data$outcome=="UA",]$value,col="blue")
  points(rep(died_ub_x_seq[i],nrow(temp_data[temp_data$died=="Died" & temp_data$outcome=="UB",])),
         temp_data[temp_data$died=="Died" & temp_data$outcome=="UB",]$value,col="orange")
  points(rep(died_uc_x_seq[i],nrow(temp_data[temp_data$died=="Died" & temp_data$outcome=="UC",])),
         temp_data[temp_data$died=="Died" & temp_data$outcome=="UC",]$value,col="red")
  points(rep(surv_x_seq[i],nrow(temp_data[temp_data$died=="Survived",])),temp_data[temp_data$died=="Survived",]$value,col="green")
}
axis(1, at=died_ub_x_seq,labels=c("Depart breeding\ngrounds",
                                  "Completion Sahara\ncrossing",
                                  "Arrive wintering\ngrounds",
                                  "Depart wintering\ngrounds",
                                  "Depart West\nAfrica",
                                  "Arrive breeding\ngrounds"),las=2)
abline(h=0,lty=3)
legend("topright",col=c("black",
                        "blue",
                        "orange",
                        "red",
                        "green"),pch=c(1,1),legend=c("Died after milestone",
                                                     "Unknown fate after milestone: A",
                                                     "Unknown fate after milestone: B",
                                                     "Unknown fate after milestone: C",
                                                     "Survived to next milestone"))
par(mar=c(5.1, 4.1, 4.1, 2.1))
