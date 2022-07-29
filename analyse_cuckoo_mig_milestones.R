#ANALYSE CUCKOO DATA

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


mig_timing_df<-read.csv("H:/JGD H drive backup 160320/Cuckoos/Output data/mig_timing_df.csv",header=T)
mig_timing_df_uncertainty<-read.csv("H:/JGD H drive backup 160320/Cuckoos/Output data/mig_timing_df_uncertainty.csv",header=T)

milestone_names<-c("dep_brgr","dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b","dep_wgr_b","dep_wgr_a","dep_wa","arr_UK","arr_brgr")

threshold<-5
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

#4A: NATURE OF VARIATION IN MILESTONES
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

#Nice table for presentation, with uncertainty included
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
# write.csv(nice_var_cor_df_tab,"H:/JGD H drive backup 160320/Cuckoos/Output data/nice_var_cor_df_tab.csv",row.names=F)

#Relationship between between-individual variance and total variance
m1<-lm(var_cor_df$total_var~var_cor_df$PTT_var);summary(m1);confint(m1)
#Relationship between within-individual variance and total variance
m1<-lm(var_cor_df$total_var~var_cor_df$resid_var);summary(m1);confint(m1)
#Relationship between between-individual variance and within-individual variance
m1<-lm(PTT_var~resid_var,data=var_cor_df);summary(m1)
#Relationship between between-individual variance and within-individual variance, without arr_wgr_b
m1<-lm(PTT_var~resid_var,data=var_cor_df[row.names(var_cor_df)!="arr_wgr_b",]);summary(m1)

#Double check - how well a job has the model done of characterising total variance?
plot(apply(mig_timing_df[,milestone_names],2,sd,na.rm=T),sqrt(var_cor_df$total_var));abline(0,1) #A very good job.

#Barplots
var_abs_df<-data.frame(rbind(var_cor_df$PTT_var,var_cor_df$resid_var))
repeatability_df<-data.frame(rbind(var_cor_df$PTT_var/(var_cor_df$PTT_var+var_cor_df$resid_var),
                                   var_cor_df$resid_var/(var_cor_df$PTT_var+var_cor_df$resid_var)))
names(var_abs_df)<-names(repeatability_df)<-milestone_names
row.names(var_abs_df)<-row.names(repeatability_df)<-c("Between-individual","Within-individual")

range(repeatability_df["Between-individual",])
mean(as.numeric(repeatability_df["Between-individual",]))

coul<-brewer.pal(nrow(var_abs_df),"Set1") 


########################################################################################################################

#4B: CAUSES OF VARIATION IN MILESTONES
#Migratory direction, breeding habitat, year, age
mtd_melted$mig_dir<-relevel(as.factor(mtd_melted$mig_dir),ref = "SW") #Put SW first, so it's on the same side as 'lowland'

##UPDATED BAYESIAN MIXED MODEL INCLUDING PTT NESTED WITHIN REGION, AND YEAR AS A CONTINUOUS FIXED EFFECT.
#Scale year to avoid numerical issues
mig_timing_df$year_s<-scale(mig_timing_df$year)[,1]
mig_timing_df$year_f<-as.factor(mig_timing_df$year)

#Create data frame to populate with results
brhab_migdir_res_df<-data.frame(matrix(ncol=18,nrow=length(milestone_names)))
names(brhab_migdir_res_df)<-c("bh_CI_overlap_0","br_CI_overlap_0_90","md_CI_overlap_0","md_CI_overlap_0_90",
                              "R2c","R2m",
                              "bh_Upland_est_marg","bh_Upland_lci_marg","bh_Upland_uci_marg","bh_Lowland_est_marg",
                              "bh_Lowland_lci_marg","bh_Lowland_uci_marg","md_SE_est_marg","md_SE_lci_marg",
                              "md_SE_uci_marg","md_SW_est_marg","md_SW_lci_marg","md_SW_uci_marg")
row.names(brhab_migdir_res_df)<-milestone_names

for(i in 1:length(milestone_names)){
  #Have to do two full models to output coefficients for each category: one with breeding habitat first, the second with migratory direction first.
  
  print(i)
  if(milestone_names[i]=="dep_wgr_a"){ #Can't handle having year_f in there as well for low sample size at dep_wgr_a
    temp_form_br_hab<-as.formula(paste0(milestone_names[i],"~br_hab+mig_dir+(1|region/ptt)"))
  } else {
    temp_form_br_hab<-as.formula(paste0(milestone_names[i],"~br_hab+mig_dir+(1|year_f)+(1|region/ptt)"))
  }
  
  m_all_br_hab<-NULL
  m_all_br_hab<-brm(temp_form_br_hab,
                    silent=1,
                    data=mig_timing_df,
                    iter=4000,
                    control = list(adapt_delta = 0.99999999999, #Close to 1 = short step length to avoid divergent transitions
                                   max_treedepth=20));summary(m_all_br_hab) #Check Rhat<1.1 for all parameters
  #Calculate R2 - just do this once, it's effectively the same for both models
  brhab_migdir_res_df$R2c[i]<-round(bayes_R2(m_all_br_hab,re.form=NULL)[1],digits=3) #Including random effects. R2c / conditional R2.
  brhab_migdir_res_df$R2m[i]<-round(bayes_R2(m_all_br_hab,re.form=NA)[1],digits=3) #Not including random effects. R2m / marginal R2.
  
  s_m_all_br_hab<-summary(m_all_br_hab)
  s_m_all_br_hab_f<-s_m_all_br_hab$fixed
  brhab_migdir_res_df$md_CI_overlap_0[i]<-s_m_all_br_hab_f[row.names(s_m_all_br_hab_f)=="mig_dirSW",]$'l-95% CI'<0 & 
    s_m_all_br_hab_f[row.names(s_m_all_br_hab_f)=="mig_dirSW",]$'u-95% CI'>0
  
  #Get 90% CI
  mig_dirSW_draws<-as_draws_df(m_all_br_hab, variable = "b_mig_dirSW")$b_mig_dirSW
  temp_mig_dirSW_lci90<-quantile(mig_dirSW_draws,0.05)
  temp_mig_dirSW_uci90<-quantile(mig_dirSW_draws,0.95)
  brhab_migdir_res_df$md_CI_overlap_0_90[i]<-temp_mig_dirSW_lci90<0 & temp_mig_dirSW_uci90>0
  
  if(milestone_names[i]=="dep_wgr_a"){ #Can't handle having year_f in there as well for low sample size at dep_wgr_a
    temp_form_mig_dir<-as.formula(paste0(milestone_names[i],"~mig_dir+br_hab+(1|region/ptt)"))
  } else {
    temp_form_mig_dir<-as.formula(paste0(milestone_names[i],"~mig_dir+br_hab+(1|year_f)+(1|region/ptt)"))
  }
  m_all_mig_dir<-NULL
  m_all_mig_dir<-brm(temp_form_mig_dir,
                     silent=1,
                     data=mig_timing_df,
                     iter=4000,
                     control = list(adapt_delta = 0.99999999999,
                                    max_treedepth=20));summary(m_all_mig_dir) #Check Rhat<1.1 for all parameters
  
  s_m_all_mig_dir<-summary(m_all_mig_dir)
  s_m_all_mig_dir_f<-s_m_all_mig_dir$fixed
  brhab_migdir_res_df$bh_CI_overlap_0[i]<-s_m_all_mig_dir_f[row.names(s_m_all_mig_dir_f)=="br_habUpland",]$'l-95% CI'<0 & 
    s_m_all_mig_dir_f[row.names(s_m_all_mig_dir_f)=="br_habUpland",]$'u-95% CI'>0
  
  #Get 90% CI
  br_habUpland_draws<-as_draws_df(m_all_mig_dir, variable = "b_br_habUpland")$b_br_habUpland
  temp_br_habUpland_lci90<-quantile(br_habUpland_draws,0.05)
  temp_br_habUpland_uci90<-quantile(br_habUpland_draws,0.95)
  brhab_migdir_res_df$br_CI_overlap_0_90[i]<-temp_br_habUpland_lci90<0 & temp_br_habUpland_uci90>0
  
  #Marginal means - taking into account of other variable
  em_bh<-data.frame(emmeans(m_all_br_hab,specs="br_hab"))
  brhab_migdir_res_df$bh_Lowland_est_marg[i]<-round(em_bh$emmean[em_bh$br_hab=="Lowland"],digits=3)
  brhab_migdir_res_df$bh_Lowland_lci_marg[i]<-round(em_bh$lower.HPD[em_bh$br_hab=="Lowland"],digits=3)
  brhab_migdir_res_df$bh_Lowland_uci_marg[i]<-round(em_bh$upper.HPD[em_bh$br_hab=="Lowland"],digits=3)
  brhab_migdir_res_df$bh_Upland_est_marg[i]<-round(em_bh$emmean[em_bh$br_hab=="Upland"],digits=3)
  brhab_migdir_res_df$bh_Upland_lci_marg[i]<-round(em_bh$lower.HPD[em_bh$br_hab=="Upland"],digits=3)
  brhab_migdir_res_df$bh_Upland_uci_marg[i]<-round(em_bh$upper.HPD[em_bh$br_hab=="Upland"],digits=3)
  
  em_md<-data.frame(emmeans(m_all_br_hab,specs="mig_dir"))
  brhab_migdir_res_df$md_SE_est_marg[i]<-round(em_md$emmean[em_md$mig_dir=="SE"],digits=3)
  brhab_migdir_res_df$md_SE_lci_marg[i]<-round(em_md$lower.HPD[em_md$mig_dir=="SE"],digits=3)
  brhab_migdir_res_df$md_SE_uci_marg[i]<-round(em_md$upper.HPD[em_md$mig_dir=="SE"],digits=3)
  brhab_migdir_res_df$md_SW_est_marg[i]<-round(em_md$emmean[em_md$mig_dir=="SW"],digits=3)
  brhab_migdir_res_df$md_SW_lci_marg[i]<-round(em_md$lower.HPD[em_md$mig_dir=="SW"],digits=3)
  brhab_migdir_res_df$md_SW_uci_marg[i]<-round(em_md$upper.HPD[em_md$mig_dir=="SW"],digits=3)
}

par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names)-0.1,brhab_migdir_res_df$bh_Lowland_est_marg,pch=19,ylab="Julian day",bty="n",xaxt="n",xlab="",col="light green")
points(1:length(milestone_names)+0.1,brhab_migdir_res_df$bh_Upland_est_marg,pch=19,col="brown")
arrows(1:length(milestone_names)-0.1,brhab_migdir_res_df$bh_Lowland_lci_marg,1:length(milestone_names)-0.1,brhab_migdir_res_df$bh_Lowland_uci_marg,length=0,col="light green")
arrows(1:length(milestone_names)+0.1,brhab_migdir_res_df$bh_Upland_lci_marg,1:length(milestone_names)+0.1,brhab_migdir_res_df$bh_Upland_uci_marg,length=0,col="brown")
legend("topleft",legend=c("Lowland","Upland"),pch=c(19,19),bty="n",col=c("light green","brown"))
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))

par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names)-0.1,brhab_migdir_res_df$md_SW_est_marg,col="blue",pch=19,ylab="Julian day",bty="n",xaxt="n",xlab="")
points(1:length(milestone_names)+0.1,brhab_migdir_res_df$md_SE_est_marg,pch=19,col="darkorange1")
arrows(1:length(milestone_names)-0.1,brhab_migdir_res_df$md_SW_lci_marg,1:length(milestone_names)-0.1,brhab_migdir_res_df$md_SW_uci_marg,length=0,col="blue")
arrows(1:length(milestone_names)+0.1,brhab_migdir_res_df$md_SE_lci_marg,1:length(milestone_names)+0.1,brhab_migdir_res_df$md_SE_uci_marg,length=0,col="darkorange1")
legend("topleft",legend=c("SW","SE"),pch=c(19,19),bty="n",col=c("blue","darkorange1"))
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))

#Make prettier for table
brhab_migdir_res_df$Upland<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$bh_Upland_est_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$Upland_lci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$bh_Upland_lci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$Upland_uci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$bh_Upland_uci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$Lowland<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$bh_Lowland_est_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$Lowland_lci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$bh_Lowland_lci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$Lowland_uci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$bh_Lowland_uci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$SE<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$md_SE_est_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$SE_lci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$md_SE_lci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$SE_uci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$md_SE_uci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$SW<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$md_SW_est_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$SW_lci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$md_SW_lci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$SW_uci<-paste0("x",substr(as.Date(round(brhab_migdir_res_df$md_SW_uci_marg),origin="1990-01-01"),6,10))
brhab_migdir_res_df$n_events<-n_events_per_milestone
brhab_migdir_res_df$n_birds<-n_birds_per_milestone
brhab_migdir_res_df$Milestone<-milestone_names

#Relate timing to lat and long
#Scale lat and long
mig_timing_df_std$br_lat_s<-scale(mig_timing_df_std$br_lat)[,1]
mig_timing_df_std$br_long_s<-scale(mig_timing_df_std$br_long)[,1]
mig_timing_df_std$year_f<-as.factor(mig_timing_df_std$year)

#Create data frame to populate with results
latlong_res_df_noregion<-latlong_res_df<-data.frame(matrix(ncol=11,nrow=length(milestone_names)))
names(latlong_res_df_noregion)<-names(latlong_res_df)<-c("br_lat_est","br_lat_lci","br_lat_uci","br_long_est","br_long_lci","br_long_uci","br_latlongint_est","br_latlongint_lci","br_latlongint_uci","R2c","R2m")
row.names(latlong_res_df_noregion)<-row.names(latlong_res_df)<-milestone_names

#With region, first
for(i in 1:length(milestone_names)){
  print(i)
  temp_form_latlong<-as.formula(paste0(milestone_names[i],"~0+br_lat_s*br_long_s+(1|year_f)+(1|region/ptt)"))
  m_all_br_latlong<-NULL
  m_all_br_latlong<-brm(temp_form_latlong,
                        silent=1,
                        data=mig_timing_df_std,
                        iter=3000,
                        control = list(adapt_delta = 0.99999999999, #Close to 1 = short step length to avoid divergent transitions
                                       max_treedepth=20));summary(m_all_br_latlong) #Check Rhat<1.1 for all parameters
  
  latlong_res_df$R2c[i]<-round(bayes_R2(m_all_br_latlong,re.form=NULL)[1],digits=3) #Including random effects. R2c / conditional R2.
  latlong_res_df$R2m[i]<-round(bayes_R2(m_all_br_latlong,re.form=NA)[1],digits=3) #Not including random effects. R2m / marginal R2.
  
  s_m_all_br_latlong<-summary(m_all_br_latlong)
  s_m_all_br_latlong_f<-s_m_all_br_latlong$fixed
  latlong_res_df$br_lat_est[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_lat_s",]$Estimate,digits=3)
  latlong_res_df$br_lat_lci[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_lat_s",]$'l-95% CI',digits=3)
  latlong_res_df$br_lat_uci[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_lat_s",]$'u-95% CI',digits=3)
  latlong_res_df$br_long_est[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_long_s",]$Estimate,digits=3)
  latlong_res_df$br_long_lci[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_long_s",]$'l-95% CI',digits=3)
  latlong_res_df$br_long_uci[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_long_s",]$'u-95% CI',digits=3)
  latlong_res_df$br_latlongint_est[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_lat_s:br_long_s",]$Estimate,digits=3)
  latlong_res_df$br_latlongint_lci[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_lat_s:br_long_s",]$'l-95% CI',digits=3)
  latlong_res_df$br_latlongint_uci[i]<-round(s_m_all_br_latlong_f[row.names(s_m_all_br_latlong_f)=="br_lat_s:br_long_s",]$'u-95% CI',digits=3)
}

latlong_res_df_range<-range(latlong_res_df[,(names(latlong_res_df) %in% c("R2c","R2m"))==F])
coul3<-brewer.pal(3,"Set1") 

par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names)-0.2,latlong_res_df$br_lat_est,pch=19,ylab="Coefficient",bty="n",xaxt="n",xlab="",col=coul3[1],ylim=latlong_res_df_range,xlim=c(0.7,10.5))
points(1:length(milestone_names),latlong_res_df$br_long_est,pch=19,col=coul3[2])
points(1:length(milestone_names)+0.2,latlong_res_df$br_latlongint_est,pch=19,col=coul3[3])
arrows(1:length(milestone_names)-0.2,latlong_res_df$br_lat_lci,1:length(milestone_names)-0.2,latlong_res_df$br_lat_uci,length=0,col=coul3[1])
arrows(1:length(milestone_names)-0,latlong_res_df$br_long_lci,1:length(milestone_names)-0,latlong_res_df$br_long_uci,length=0,col=coul3[2])
arrows(1:length(milestone_names)+0.2,latlong_res_df$br_latlongint_lci,1:length(milestone_names)+0.2,latlong_res_df$br_latlongint_uci,length=0,col=coul3[3])
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
legend("topright",c("Latitude","Longitude","Latitude:Longitude"),pch=c(19,19,19),col=coul3,bty="n")
abline(h=0,lty=3)
par(mar=c(5.1, 4.1, 4.1, 2.1))

#Now without region
for(i in 1:length(milestone_names)){
  print(i)
  temp_form_latlong_noregion<-as.formula(paste0(milestone_names[i],"~0+br_lat_s*br_long_s+(1|year_f)+(1|ptt)"))
  m_all_br_latlong_noregion<-NULL
  m_all_br_latlong_noregion<-brm(temp_form_latlong_noregion,
                                 silent=1,
                                 data=mig_timing_df_std,
                                 iter=3000,
                                 control = list(adapt_delta = 0.99999999999, #Close to 1 = short step length to avoid divergent transitions
                                                max_treedepth=20));summary(m_all_br_latlong_noregion) #Check Rhat<1.1 for all parameters
  
  latlong_res_df_noregion$R2c[i]<-round(bayes_R2(m_all_br_latlong_noregion,re.form=NULL)[1],digits=3) #Including random effects. R2c / conditional R2.
  latlong_res_df_noregion$R2m[i]<-round(bayes_R2(m_all_br_latlong_noregion,re.form=NA)[1],digits=3) #Not including random effects. R2m / marginal R2.
  
  s_m_all_br_latlong_noregion<-summary(m_all_br_latlong_noregion)
  s_m_all_br_latlong_noregion_f<-s_m_all_br_latlong_noregion$fixed
  latlong_res_df_noregion$br_lat_est[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_lat_s",]$Estimate,digits=3)
  latlong_res_df_noregion$br_lat_lci[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_lat_s",]$'l-95% CI',digits=3)
  latlong_res_df_noregion$br_lat_uci[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_lat_s",]$'u-95% CI',digits=3)
  latlong_res_df_noregion$br_long_est[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_long_s",]$Estimate,digits=3)
  latlong_res_df_noregion$br_long_lci[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_long_s",]$'l-95% CI',digits=3)
  latlong_res_df_noregion$br_long_uci[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_long_s",]$'u-95% CI',digits=3)
  latlong_res_df_noregion$br_latlongint_est[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_lat_s:br_long_s",]$Estimate,digits=3)
  latlong_res_df_noregion$br_latlongint_lci[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_lat_s:br_long_s",]$'l-95% CI',digits=3)
  latlong_res_df_noregion$br_latlongint_uci[i]<-round(s_m_all_br_latlong_noregion_f[row.names(s_m_all_br_latlong_noregion_f)=="br_lat_s:br_long_s",]$'u-95% CI',digits=3)
}
#Now do equivalent plots

par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names)-0.2,latlong_res_df_noregion$br_lat_est,pch=19,ylab="Coefficient",bty="n",xaxt="n",xlab="",col=coul3[1],ylim=latlong_res_df_range,xlim=c(0.7,10.5))
points(1:length(milestone_names),latlong_res_df_noregion$br_long_est,pch=19,col=coul3[2])
points(1:length(milestone_names)+0.2,latlong_res_df_noregion$br_latlongint_est,pch=19,col=coul3[3])
arrows(1:length(milestone_names)-0.2,latlong_res_df_noregion$br_lat_lci,1:length(milestone_names)-0.2,latlong_res_df_noregion$br_lat_uci,length=0,col=coul3[1])
arrows(1:length(milestone_names)-0,latlong_res_df_noregion$br_long_lci,1:length(milestone_names)-0,latlong_res_df_noregion$br_long_uci,length=0,col=coul3[2])
arrows(1:length(milestone_names)+0.2,latlong_res_df_noregion$br_latlongint_lci,1:length(milestone_names)+0.2,latlong_res_df_noregion$br_latlongint_uci,length=0,col=coul3[3])
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
legend("topright",c("Latitude","Longitude","Latitude:Longitude"),pch=c(19,19,19),col=coul3,bty="n")
abline(h=0,lty=3)
par(mar=c(5.1, 4.1, 4.1, 2.1))


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

#Write out param estimates for SI results tables
age_coeff_df_nicer<-age_coeff_df
age_coeff_df_nicer$first_ad_cyc_est<-round(age_coeff_df$first_ad_cyc_est,digits=3)
age_coeff_df_nicer$first_ad_cyc_lci<-round(age_coeff_df$first_ad_cyc_lci,digits=3)
age_coeff_df_nicer$first_ad_cyc_uci<-round(age_coeff_df$first_ad_cyc_uci,digits=3)
age_coeff_df_nicer$later_cyc_est<-round(age_coeff_df$later_cyc_est,digits=3)
age_coeff_df_nicer$later_cyc_lci<-round(age_coeff_df$later_cyc_lci,digits=3)
age_coeff_df_nicer$later_cyc_uci<-round(age_coeff_df$later_cyc_uci,digits=3)

par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:nrow(age_coeff_df),
     age_coeff_df$first_ad_cyc_est,
     xlim=c(0.5,nrow(age_coeff_df)+0.5),
     ylim=c(min(c(age_coeff_df$first_ad_cyc_lci,age_coeff_df$first_ad_cyc_lci)),
            max(c(age_coeff_df$first_ad_cyc_uci,age_coeff_df$first_ad_cyc_uci))),
     col="white",xaxt="n",xlab="",pch=19,ylab="Within-individual anomaly",main="")
points((1:length(milestone_names))-0.1,age_coeff_df$first_ad_cyc_est,pch=19,col="cornflowerblue")
arrows((1:nrow(age_coeff_df))-0.1,age_coeff_df$first_ad_cyc_lci,(1:nrow(age_coeff_df))-0.1,age_coeff_df$first_ad_cyc_uci,length=0,col="cornflowerblue")
points((1:length(milestone_names))+0.1,age_coeff_df$later_cyc_est,pch=19,col="chocolate")
arrows((1:nrow(age_coeff_df))+0.1,age_coeff_df$later_cyc_lci,(1:nrow(age_coeff_df))+0.1,age_coeff_df$later_cyc_uci,length=0,col="chocolate")
age_sig05_ind<-which(age_coeff_df$sig.05)
age_sig1_ind<-which(age_coeff_df$sig.1==T & age_coeff_df$sig.05==F)
arrows(age_sig05_ind-0.3,20,age_sig05_ind+0.3,20,length=0) #NB hard coding!
arrows(age_sig05_ind-0.3,20,age_sig05_ind-0.3,19,length=0) #NB hard coding!
arrows(age_sig05_ind+0.3,20,age_sig05_ind+0.3,19,length=0) #NB hard coding!
text(age_sig05_ind,21,"**",cex=2) #NB hard coding!
arrows(age_sig1_ind-0.3,20,age_sig1_ind+0.3,20,length=0) #NB hard coding!
arrows(age_sig1_ind-0.3,20,age_sig1_ind-0.3,19,length=0) #NB hard coding!
arrows(age_sig1_ind+0.3,20,age_sig1_ind+0.3,19,length=0) #NB hard coding!
text(age_sig1_ind,21,"*",cex=2) #NB hard coding!
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
legend("topleft",c("Cycle after first breeding season","All subsequent cycles"),col=c("cornflowerblue","chocolate"),pch=c(19,19),bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))

#Cycles since fledging
ord_age_param_est_list<-list()
par(mfrow=c(3,3))
for(i in 1:length(milestone_names)){
  print(i)
  if(milestone_names[i]=="dep_wgr_a") next #Not enough data for this milestone: trying to estimate an individual random effect and four parameters from 6 data
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(known_age)+(1|ptt)"))
  m1<-brm(temp_form_age,
          silent=1,
          data=mig_timing_df_stdind[which(is.na(mig_timing_df_stdind[,milestone_names[i]])==F),],
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
       main=milestone_names[i],
       col="blue",
       xlab="Cycles since first breeding season",
       xlim=c(1,5),
       ylab="Within-individual anomaly",
       ylim=range(c(temp_lci,temp_uci)))
  arrows(temp_ages,temp_lci,temp_ages,temp_uci,length=0,col="blue")
  abline(h=0,lty=3)
  points(temp_ages[sm1_f$sig_84],sm1_f$Estimate[sm1_f$sig_84],pch=19,col="blue")
  
  temp_param_est_df<-sm1_f[,c("Estimate","l-95% CI","u-95% CI")]
  row.names(temp_param_est_df)<-paste0(milestone_names[i],"_age",temp_ages)
  ord_age_param_est_list[[i]]<-temp_param_est_df[temp_ages<6,]
}
par(mfrow=c(1,1))
ord_age_param_est_df<-do.call(rbind,ord_age_param_est_list)

ord_age_param_est_df_nicer<-ord_age_param_est_df
ord_age_param_est_df_nicer$Estimate<-round(ord_age_param_est_df$Estimate,digits=3)
ord_age_param_est_df_nicer$`l-95% CI`<-round(ord_age_param_est_df$`l-95% CI`,digits=3)
ord_age_param_est_df_nicer$`u-95% CI`<-round(ord_age_param_est_df$`u-95% CI`,digits=3)

#Fig A9 but for non-known age birds
par(mfrow=c(3,3))
for(i in 1:length(milestone_names)){
  print(i)
  if(milestone_names[i]=="dep_wgr_a") next #Not enough data for this milestone: trying to estimate an individual random effect and four parameters from 6 data
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(cycles_since_tagging_4s)+(1|ptt)"))
  m1<-brm(temp_form_age,
          silent=1,
          data=mig_timing_df_stdind[which(is.na(mig_timing_df_stdind[,milestone_names[i]])==F),],
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
       main=milestone_names[i],
       col="blue",
       xlab="Cycles since tagging",
       xlim=c(0,4),
       ylab="Within-individual anomaly",
       ylim=range(c(temp_lci,temp_uci)))
  arrows(temp_ages,temp_lci,temp_ages,temp_uci,length=0,col="blue")
  abline(h=0,lty=3)
  points(temp_ages[sm1_f$sig_84],sm1_f$Estimate[sm1_f$sig_84],pch=19,col="blue")
  
  temp_param_est_df<-sm1_f[,c("Estimate","l-95% CI","u-95% CI")]
  row.names(temp_param_est_df)<-paste0(milestone_names[i],"_age",temp_ages)
}
par(mfrow=c(1,1))


########################################################################################################################
#4C: RELATIONSHIPS BETWEEN MILESTONES
#What is the relationship between subsequent anomalies?
milestone_name_pairs<-data.frame("from"=c("dep_brgr","dep_UK","fin_sah_sb","fin_sah_sb","arr_wgr_a","arr_wgr_b",
                                          "dep_wgr_a","dep_wgr_b","dep_wa","arr_UK","arr_brgr"),
                                 "to"=c("dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b",
                                        "dep_wa","dep_wa","arr_UK","arr_brgr","dep_brgr"))
milestone_name_pairs$R2m<-milestone_name_pairs$slope<-milestone_name_pairs$slope_lci<-milestone_name_pairs$slope_uci<-milestone_name_pairs$sign<-NA

#Make data frame of arrival to the breeding grounds and the following departure from the breeding grounds
un_mtdsi_ptt<-unique(mig_timing_df_stdind$ptt)
arrdepbrgr_df_list<-list()

for(i in 1:length(un_mtdsi_ptt)){
  temp_mtdsi<-mig_timing_df_stdind[mig_timing_df_stdind$ptt==un_mtdsi_ptt[i],]
  if(nrow(temp_mtdsi)==1) next
  temp_arrdepbrgr_df<-data.frame(matrix(ncol=6,nrow=nrow(temp_mtdsi)-1))
  names(temp_arrdepbrgr_df)<-c("ptt","year_1st_cycle","arr_brgr","dep_brgr","mig_dir_1st_cycle","region")
  temp_arrdepbrgr_df$ptt<-un_mtdsi_ptt[i]
  temp_arrdepbrgr_df$year_1st_cycle<-temp_mtdsi$year[1:(nrow(temp_mtdsi)-1)]
  temp_arrdepbrgr_df$arr_brgr<-temp_mtdsi$arr_brgr[1:(nrow(temp_mtdsi)-1)]
  temp_arrdepbrgr_df$dep_brgr<-temp_mtdsi$dep_brgr[2:nrow(temp_mtdsi)]
  temp_arrdepbrgr_df$mig_dir_1st_cycle<-temp_mtdsi$mig_dir[1:(nrow(temp_mtdsi)-1)]
  temp_arrdepbrgr_df$region<-temp_mtdsi$region[1]
  arrdepbrgr_df_list[[i]]<-temp_arrdepbrgr_df
}
arrdepbrgr_df<-do.call(rbind,arrdepbrgr_df_list)

par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){ #NB ONLY GO AS FAR AS ARR_UK - THE ARR_BRGR:DEP_BRGR STEP HAS TO BE DONE MANUALLY
  print(i)
  temp_data_df<-data.frame("t"=mig_timing_df_stdind[,milestone_name_pairs$from[i]],
                           "tplus1"=mig_timing_df_stdind[,milestone_name_pairs$to[i]],
                           "ptt"=mig_timing_df_stdind$ptt,
                           "region"=mig_timing_df_stdind$region)
  
  max_abs_vals<-max(abs(temp_data_df[,c("t","tplus1")]),na.rm=T)
  plot(temp_data_df$t,temp_data_df$tplus1,
       xlim=c(-max_abs_vals,max_abs_vals),
       ylim=c(-max_abs_vals,max_abs_vals),
       bty="n",
       pch=21,
       bg="seagreen3",
       xlab=milestone_name_pairs$from[i],
       ylab=milestone_name_pairs$to[i])
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  m2<-brm(tplus1~t+(1|region/ptt),
          data=temp_data_df,
          control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
  sm2<-summary(m2)
  R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
  milestone_name_pairs$R2m[i]<-R2m
  milestone_name_pairs$slope[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$Estimate
  milestone_name_pairs$slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'
  milestone_name_pairs$slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'
  milestone_name_pairs$sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'>0)==F
  newx<-seq(-max_abs_vals,max_abs_vals)
  intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
  t_draws<-as_draws_df(m2,variable="b_t")$b_t
  predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
  for(i in 1:nrow(predicted_lines_df)){
    predicted_lines_df[i,]<-(newx*t_draws[i])+intercept_draws[i]
  }
  lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
  lines(newx,apply(predicted_lines_df,2,mean))
  lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
  text(-max_abs_vals,max_abs_vals,paste0("R2m = ",R2m),pos=4)
}

#Add-on for arr_brgr vs dep_brgr
i<-nrow(milestone_name_pairs)
max_abs_vals<-max(abs(arrdepbrgr_df[,c("arr_brgr","dep_brgr")]),na.rm=T)
plot(arrdepbrgr_df$arr_brgr,arrdepbrgr_df$dep_brgr,
     xlim=c(-max_abs_vals,max_abs_vals),
     ylim=c(-max_abs_vals,max_abs_vals),
     bty="n",
     pch=21,
     bg="seagreen3",
     xlab="arr_brgr",
     ylab="dep_brgr")
abline(h=0,lty=3)
abline(v=0,lty=3)
m2<-brm(dep_brgr~arr_brgr+(1|region/ptt),
        data=arrdepbrgr_df,
        control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
sm2<-summary(m2)
R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
milestone_name_pairs$R2m[i]<-R2m
milestone_name_pairs$slope[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$Estimate
milestone_name_pairs$slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'l-95% CI'
milestone_name_pairs$slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'u-95% CI'
milestone_name_pairs$sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'u-95% CI'>0)==F
newx<-seq(-max_abs_vals,max_abs_vals)
intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
arr_brgr_draws<-as_draws_df(m2,variable="b_arr_brgr")$b_arr_brgr
predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
for(i in 1:nrow(predicted_lines_df)){
  predicted_lines_df[i,]<-(newx*arr_brgr_draws[i])+intercept_draws[i]
}
lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
lines(newx,apply(predicted_lines_df,2,mean))
lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
text(-max_abs_vals,max_abs_vals,paste0("R2m = ",R2m),pos=4)

par(mfrow=c(1,1))


milestone_name_pairs$sd<-apply(mig_timing_df_stdind[,milestone_name_pairs$from],2,sd,na.rm=T)
milestone_name_pairs$delta_sd<-apply(mig_timing_df_stdind[,milestone_name_pairs$to],2,sd,na.rm=T)-milestone_name_pairs$sd
milestone_name_pairs$var<-apply(mig_timing_df_stdind[,milestone_name_pairs$from],2,var,na.rm=T)
milestone_name_pairs$delta_var<-apply(mig_timing_df_stdind[,milestone_name_pairs$to],2,var,na.rm=T)-milestone_name_pairs$var

#F-test: are the variances significantly different between milestones?
milestone_name_pairs$var_sig_diff<-NA
for(i in 1:nrow(milestone_name_pairs)){
  temp_f_test<-var.test(mig_timing_df_stdind[,milestone_name_pairs$from[i]],mig_timing_df_stdind[,milestone_name_pairs$to[i]])
  milestone_name_pairs$var_sig_diff[i]<-temp_f_test$p.value<0.05
}

var_sig_diff_df<-data.frame(matrix(nrow=length(milestone_names),ncol=length(milestone_names)))
row.names(var_sig_diff_df)<-names(var_sig_diff_df)<-milestone_names
for(i in 1:length(milestone_names)){
  for(j in 1:length(milestone_names)){
    temp_var_test<-var.test(mig_timing_df_stdind[,row.names(var_sig_diff_df)[i]],mig_timing_df_stdind[,names(var_sig_diff_df)[j]])
    var_sig_diff_df[i,j]<-temp_var_test$p.value<0.05
  }
}

#Boxplot: raw data, with significance of difference in variance
boxplot(mtd_melted$value~mtd_melted$variable,
        las=2,
        frame=F,
        xlab="",
        ylab="Days since start of 1st calendar year of mig cycle",
        ylim=c(min(mtd_melted$value,na.rm=T),max(mtd_melted$value,na.rm=T)+50),
        main="")
#CAUTION HARD CODING!
arrows(2.1,max(mtd_melted$value,na.rm=T)+20,2.9,max(mtd_melted$value,na.rm=T)+20,length=0)
arrows(2.1,max(mtd_melted$value,na.rm=T)+10,2.1,max(mtd_melted$value,na.rm=T)+20,length=0)
arrows(2.9,max(mtd_melted$value,na.rm=T)+10,2.9,max(mtd_melted$value,na.rm=T)+20,length=0)
text(2.5,max(mtd_melted$value,na.rm=T)+30,"*",cex=2)

arrows(3.1,max(mtd_melted$value,na.rm=T)+20,4.9,max(mtd_melted$value,na.rm=T)+20,length=0)
arrows(3.1,max(mtd_melted$value,na.rm=T)+10,3.1,max(mtd_melted$value,na.rm=T)+20,length=0)
arrows(4.9,max(mtd_melted$value,na.rm=T)+10,4.9,max(mtd_melted$value,na.rm=T)+20,length=0)
text(4,max(mtd_melted$value,na.rm=T)+30,"*",cex=2)

arrows(6.1,max(mtd_melted$value,na.rm=T)+20+30,7.9,max(mtd_melted$value,na.rm=T)+20+30,length=0)
arrows(6.1,max(mtd_melted$value,na.rm=T)+10+30,6.1,max(mtd_melted$value,na.rm=T)+20+30,length=0)
arrows(7.9,max(mtd_melted$value,na.rm=T)+10+30,7.9,max(mtd_melted$value,na.rm=T)+20+30,length=0)
text(7,max(mtd_melted$value,na.rm=T)+30+30,"*",cex=2)

arrows(7.1,max(mtd_melted$value,na.rm=T)+20,7.9,max(mtd_melted$value,na.rm=T)+20,length=0)
arrows(7.1,max(mtd_melted$value,na.rm=T)+10,7.1,max(mtd_melted$value,na.rm=T)+20,length=0)
arrows(7.9,max(mtd_melted$value,na.rm=T)+10,7.9,max(mtd_melted$value,na.rm=T)+20,length=0)
text(7.5,max(mtd_melted$value,na.rm=T)+30,"*",cex=2)

#Relationship between anomalies and arr_brgr anomaly
milestone_df_arrbrgr<-data.frame("from"=milestone_names[milestone_names!="arr_brgr"])
milestone_df_arrbrgr$R2m<-milestone_df_arrbrgr$slope<-milestone_df_arrbrgr$slope_lci<-milestone_df_arrbrgr$slope_uci<-milestone_df_arrbrgr$sign<-NA

par(mfrow=c(3,3))
for(i in 1:nrow(milestone_df_arrbrgr)){
  print(i)
  temp_data_df<-data.frame("t"=mig_timing_df_stdind[,milestone_df_arrbrgr$from[i]],
                           "arr_brgr"=mig_timing_df_stdind$arr_brgr,
                           "ptt"=mig_timing_df_stdind$ptt,
                           "region"=mig_timing_df_stdind$region)
  
  max_abs_vals<-max(abs(temp_data_df[,c("t","arr_brgr")]),na.rm=T)
  plot(temp_data_df$t,temp_data_df$arr_brgr,
       xlim=c(-max_abs_vals,max_abs_vals),
       ylim=c(-max_abs_vals,max_abs_vals),
       # main="Relationship between anomalies\nfor subsequent milestones in same year",
       bty="n",
       pch=21,
       bg="seagreen3",
       xlab=milestone_df_arrbrgr$from[i],
       ylab="arr_brgr")
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  m2<-brm(arr_brgr~t+(1|region/ptt),
          data=temp_data_df,
          control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
  sm2<-summary(m2)
  R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
  milestone_df_arrbrgr$R2m[i]<-R2m
  milestone_df_arrbrgr$slope[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$Estimate
  milestone_df_arrbrgr$slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'
  milestone_df_arrbrgr$slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'
  milestone_df_arrbrgr$sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'>0)==F
  newx<-seq(-max_abs_vals,max_abs_vals)
  intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
  t_draws<-as_draws_df(m2,variable="b_t")$b_t
  predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
  for(i in 1:nrow(predicted_lines_df)){
    predicted_lines_df[i,]<-(newx*t_draws[i])+intercept_draws[i]
  }
  lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
  lines(newx,apply(predicted_lines_df,2,mean))
  lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
  
  text(-max_abs_vals,max_abs_vals,paste0("R2m = ",R2m),pos=4)
}
par(mfrow=c(1,1))

ms_table<-milestone_name_pairs[,c("from","var","to","delta_var","var_sig_diff","R2m","sign")]
ms_table$var<-round(milestone_name_pairs$var,digits=2)
ms_table$delta_var<-round(milestone_name_pairs$delta_var,digits=2)
ms_table$R2m<-round(milestone_name_pairs$R2m,digits=2)
ms_table$R2m_arrbrgr_sign<-ms_table$R2m_arrbrgr<-NA
for(i in 1:nrow(milestone_df_arrbrgr)){
  ms_table[ms_table$from==milestone_df_arrbrgr[i,]$from,]$R2m_arrbrgr<-round(milestone_df_arrbrgr[i,]$R2m,digits=2)
  ms_table[ms_table$from==milestone_df_arrbrgr[i,]$from,]$R2m_arrbrgr_sign<-milestone_df_arrbrgr[i,]$sign
}

#Make pretty table for paper
sd_r2_tab_within_ind<-milestone_name_pairs[,c("from","sd","to","delta_sd","R2m")]
sd_r2_tab_within_ind$R2m_arrbrgr<-NA
for(i in 1:nrow(milestone_df_arrbrgr)){
  sd_r2_tab_within_ind[sd_r2_tab_within_ind$from==milestone_df_arrbrgr[i,]$from,]$R2m_arrbrgr<-milestone_df_arrbrgr[i,]$R2m
}
sd_r2_tab_within_ind[,c("sd","delta_sd","R2m","R2m_arrbrgr")]<-round(sd_r2_tab_within_ind[,c("sd","delta_sd","R2m","R2m_arrbrgr")],digits=3)


#What is the relationship between subsequent milestones (absolute values)?
milestone_name_pairs_abs<-data.frame("from"=c("dep_brgr","dep_UK","fin_sah_sb","fin_sah_sb","arr_wgr_a","arr_wgr_b",
                                              "dep_wgr_a","dep_wgr_b","dep_wa","arr_UK","arr_brgr"),
                                     "to"=c("dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b",
                                            "dep_wa_a","dep_wa","arr_UK","arr_brgr","dep_brgr"))
milestone_name_pairs_abs$R2m<-milestone_name_pairs_abs$slope<-milestone_name_pairs_abs$slope_lci<-milestone_name_pairs_abs$slope_uci<-milestone_name_pairs_abs$sign<-NA

#Make data frame of arrival to the breeding grounds and the following departure from the breeding grounds
un_mtd_ptt<-unique(mig_timing_df$ptt)
arrdepbrgr_abs_df_list<-list()

for(i in 1:length(un_mtd_ptt)){
  temp_mtd<-mig_timing_df[mig_timing_df$ptt==un_mtd_ptt[i],]
  if(nrow(temp_mtd)==1) next
  temp_arrdepbrgr_abs_df<-data.frame(matrix(ncol=6,nrow=nrow(temp_mtd)-1))
  names(temp_arrdepbrgr_abs_df)<-c("ptt","year_1st_cycle","arr_brgr","dep_brgr","mig_dir_1st_cycle","region")
  temp_arrdepbrgr_abs_df$ptt<-un_mtd_ptt[i]
  temp_arrdepbrgr_abs_df$year_1st_cycle<-temp_mtd$year[1:(nrow(temp_mtd)-1)]
  temp_arrdepbrgr_abs_df$arr_brgr<-temp_mtd$arr_brgr[1:(nrow(temp_mtd)-1)]
  temp_arrdepbrgr_abs_df$dep_brgr<-temp_mtd$dep_brgr[2:nrow(temp_mtd)]
  temp_arrdepbrgr_abs_df$mig_dir_1st_cycle<-temp_mtd$mig_dir[1:(nrow(temp_mtd)-1)]
  temp_arrdepbrgr_abs_df$region<-temp_mtd$region[1]
  arrdepbrgr_abs_df_list[[i]]<-temp_arrdepbrgr_abs_df
}
arrdepbrgr_abs_df<-do.call(rbind,arrdepbrgr_abs_df_list)

par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){ #NB ONLY GO AS FAR AS ARR_UK - THE ARR_BRGR:DEP_BRGR STEP HAS TO BE DONE MANUALLY
  print(i)
  temp_data_df<-data.frame("t"=mig_timing_df[,milestone_name_pairs_abs$from[i]],
                           "tplus1"=mig_timing_df[,milestone_name_pairs_abs$to[i]],
                           "ptt"=mig_timing_df$ptt,
                           "region"=mig_timing_df$region)
  
  max_abs_vals<-max(abs(temp_data_df[,c("t","tplus1")]),na.rm=T)
  plot(temp_data_df$t,temp_data_df$tplus1,
       xlim=range(temp_data_df$t,na.rm=T),
       ylim=range(temp_data_df$tplus1,na.rm=T),
       # main="Relationship between dates\nfor subsequent milestones in same year",
       bty="n",
       pch=21,
       bg="seagreen3",
       xlab=milestone_name_pairs_abs$from[i],
       ylab=milestone_name_pairs_abs$to[i])
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  m2<-brm(tplus1~t+(1|region/ptt),
          iter=4000,
          data=temp_data_df,
          control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
  sm2<-summary(m2)
  R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
  milestone_name_pairs_abs$R2m[i]<-R2m
  milestone_name_pairs_abs$slope[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$Estimate
  milestone_name_pairs_abs$slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'
  milestone_name_pairs_abs$slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'
  milestone_name_pairs_abs$sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'>0)==F
  newx<-seq(min(temp_data_df$t,na.rm=T),max(temp_data_df$t,na.rm=T))
  intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
  t_draws<-as_draws_df(m2,variable="b_t")$b_t
  predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
  for(i in 1:nrow(predicted_lines_df)){
    predicted_lines_df[i,]<-(newx*t_draws[i])+intercept_draws[i]
  }
  lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
  lines(newx,apply(predicted_lines_df,2,mean))
  lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
  
  text(min(temp_data_df$t,na.rm=T),max(temp_data_df$tplus1,na.rm=T),paste0("R2m = ",R2m),pos=4)
}

#Add-on for arr_brgr vs dep_brgr
i<-nrow(milestone_name_pairs_abs)
max_abs_vals<-max(abs(arrdepbrgr_abs_df[,c("arr_brgr","dep_brgr")]),na.rm=T)
plot(arrdepbrgr_abs_df$arr_brgr,arrdepbrgr_abs_df$dep_brgr,
     xlim=range(arrdepbrgr_abs_df$arr_brgr,na.rm=T),
     ylim=range(arrdepbrgr_abs_df$dep_brgr,na.rm=T),
     # main="Relationship between anomalies\nfor subsequent milestones in same year",
     bty="n",
     pch=21,
     bg="seagreen3",
     xlab="arr_brgr",
     ylab="dep_brgr")
abline(h=0,lty=3)
abline(v=0,lty=3)
m2<-brm(dep_brgr~arr_brgr+(1|region/ptt),
        data=arrdepbrgr_abs_df,
        control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
sm2<-summary(m2)
R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
milestone_name_pairs_abs$R2m[i]<-R2m
milestone_name_pairs_abs$slope[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$Estimate
milestone_name_pairs_abs$slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'l-95% CI'
milestone_name_pairs_abs$slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'u-95% CI'
milestone_name_pairs_abs$sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="arr_brgr",]$'u-95% CI'>0)==F
newx<-seq(min(arrdepbrgr_abs_df$arr_brgr,na.rm=T),max(arrdepbrgr_abs_df$arr_brgr,na.rm=T))
intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
arr_brgr_draws<-as_draws_df(m2,variable="b_arr_brgr")$b_arr_brgr
predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
for(i in 1:nrow(predicted_lines_df)){
  predicted_lines_df[i,]<-(newx*arr_brgr_draws[i])+intercept_draws[i]
}
lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
lines(newx,apply(predicted_lines_df,2,mean))
lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
text(min(arrdepbrgr_abs_df$arr_brgr,na.rm=T),max(arrdepbrgr_abs_df$dep_brgr,na.rm=T),paste0("R2m = ",R2m),pos=4)
par(mfrow=c(1,1))


png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig7_2.png",6,6,res=800,units="in")
par(mar=c(10.1, 4.1, 1.1, 1.1))
plot(1:nrow(milestone_name_pairs_abs),milestone_name_pairs_abs$slope,ylim=c(min(milestone_name_pairs_abs$slope_lci),max(milestone_name_pairs_abs$slope_uci)),col="white",xaxt="n",xlab="",pch=19,ylab="Slope coefficient",main="")
points(which(milestone_name_pairs_abs$sign==F),milestone_name_pairs_abs$slope[which(milestone_name_pairs_abs$sign==F)])
points(which(milestone_name_pairs_abs$sign),milestone_name_pairs_abs$slope[which(milestone_name_pairs_abs$sign)],pch=19)
arrows(1:nrow(milestone_name_pairs_abs),milestone_name_pairs_abs$slope_lci,1:nrow(milestone_name_pairs_abs),milestone_name_pairs_abs$slope_uci,length=0)
abline(h=0,lty=3)
axis(1, at=1:nrow(milestone_name_pairs_abs),labels=paste0(milestone_name_pairs_abs$from," : ",milestone_name_pairs_abs$to),las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Combine slope coeff plots for both within-ind and absolute
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig2a.png",6,6,res=800,units="in")
par(mar=c(10.1, 4.1, 1.1, 1.1))
plot(1:nrow(milestone_name_pairs),
     milestone_name_pairs$slope,
     xlim=c(0.5,nrow(milestone_name_pairs)+0.5),
     ylim=c(min(c(milestone_name_pairs$slope_lci,milestone_name_pairs_abs$slope_lci)),
            max(c(milestone_name_pairs$slope_uci,milestone_name_pairs_abs$slope_uci))),
     col="white",xaxt="n",xlab="",pch=19,ylab="Slope coefficient",main="")
points(which(milestone_name_pairs$sign==F)-0.1,milestone_name_pairs$slope[which(milestone_name_pairs$sign==F)],col="cornflowerblue")
points(which(milestone_name_pairs$sign)-0.1,milestone_name_pairs$slope[which(milestone_name_pairs$sign)],pch=19,col="cornflowerblue")
arrows((1:nrow(milestone_name_pairs))-0.1,milestone_name_pairs$slope_lci,(1:nrow(milestone_name_pairs))-0.1,milestone_name_pairs$slope_uci,length=0,col="cornflowerblue")
points(which(milestone_name_pairs_abs$sign==F)+0.1,milestone_name_pairs_abs$slope[which(milestone_name_pairs_abs$sign==F)],col="chocolate")
points(which(milestone_name_pairs_abs$sign)+0.1,milestone_name_pairs_abs$slope[which(milestone_name_pairs_abs$sign)],pch=19,col="chocolate")
arrows((1:nrow(milestone_name_pairs_abs))+0.1,milestone_name_pairs_abs$slope_lci,(1:nrow(milestone_name_pairs_abs))+0.1,milestone_name_pairs_abs$slope_uci,length=0,col="chocolate")
abline(h=0,lty=3)
axis(1, at=1:nrow(milestone_name_pairs),labels=paste0(milestone_name_pairs$from," : ",milestone_name_pairs$to),las=2)
legend("topright",c("Within-individual anomaly","Absolute date"),col=c("cornflowerblue","chocolate"),pch=c(19,19),bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#What is the relationship between milestones and arrival to the breeding grounds (absolute values)?
milestone_df_arrbrgr_abs<-data.frame("from"=milestone_names[milestone_names!="arr_brgr"])
milestone_df_arrbrgr_abs$R2m<-milestone_df_arrbrgr_abs$slope<-milestone_df_arrbrgr_abs$slope_lci<-milestone_df_arrbrgr_abs$slope_uci<-milestone_df_arrbrgr_abs$sign<-NA

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA4.png",10,10,res=800,units="in")
par(mfrow=c(3,3))
for(i in 1:nrow(milestone_df_arrbrgr_abs)){
  print(i)
  temp_data_df<-data.frame("t"=mig_timing_df[,milestone_df_arrbrgr_abs$from[i]],
                           "arr_brgr"=mig_timing_df$arr_brgr,
                           "ptt"=mig_timing_df$ptt,
                           "region"=mig_timing_df$region)
  max_abs_vals<-max(abs(temp_data_df[,c("t","arr_brgr")]),na.rm=T)
  plot(temp_data_df$t,temp_data_df$arr_brgr,
       xlim=range(temp_data_df$t,na.rm=T),
       ylim=range(temp_data_df$arr_brgr,na.rm=T),
       # main="Relationship between anomalies\nfor subsequent milestones in same year",
       bty="n",
       pch=21,
       bg="seagreen3",
       xlab=milestone_df_arrbrgr_abs$from[i],
       ylab="arr_brgr")
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  m2<-brm(arr_brgr~t+(1|region/ptt),
          iter=4000,
          data=temp_data_df,
          control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
  sm2<-summary(m2)
  R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
  milestone_df_arrbrgr_abs$R2m[i]<-R2m
  milestone_df_arrbrgr_abs$slope[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$Estimate
  milestone_df_arrbrgr_abs$slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'
  milestone_df_arrbrgr_abs$slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'
  milestone_df_arrbrgr_abs$sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="t",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="t",]$'u-95% CI'>0)==F
  newx<-seq(-max_abs_vals,max_abs_vals)
  intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
  t_draws<-as_draws_df(m2,variable="b_t")$b_t
  predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
  for(i in 1:nrow(predicted_lines_df)){
    predicted_lines_df[i,]<-(newx*t_draws[i])+intercept_draws[i]
  }
  lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
  lines(newx,apply(predicted_lines_df,2,mean))
  lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
  
  text(min(temp_data_df$t,na.rm=T),max(temp_data_df$arr_brgr,na.rm=T),paste0("R2m = ",R2m),pos=4)
}
par(mfrow=c(1,1))
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig8_2.png",6,6,res=800,units="in")
par(mar=c(10.1, 4.1, 1.1, 1.1))
plot(1:nrow(milestone_df_arrbrgr_abs),milestone_df_arrbrgr_abs$slope,ylim=c(min(milestone_df_arrbrgr_abs$slope_lci),max(milestone_df_arrbrgr_abs$slope_uci)),col="white",xaxt="n",xlab="",pch=19,ylab="Slope coefficient",main="")
points(which(milestone_df_arrbrgr_abs$sign==F),milestone_df_arrbrgr_abs$slope[which(milestone_df_arrbrgr_abs$sign==F)])
points(which(milestone_df_arrbrgr_abs$sign),milestone_df_arrbrgr_abs$slope[which(milestone_df_arrbrgr_abs$sign)],pch=19)
arrows(1:nrow(milestone_df_arrbrgr_abs),milestone_df_arrbrgr_abs$slope_lci,1:nrow(milestone_df_arrbrgr_abs),milestone_df_arrbrgr_abs$slope_uci,length=0)
abline(h=0,lty=3)
axis(1, at=1:nrow(milestone_df_arrbrgr_abs),labels=milestone_df_arrbrgr_abs$from,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Combine slope coeff plots for both within-ind and absolute
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig2b.png",6,6,res=800,units="in")
par(mar=c(10.1, 4.1, 1.1, 1.1))
plot(1:nrow(milestone_df_arrbrgr),
     milestone_df_arrbrgr$slope,
     xlim=c(0.5,nrow(milestone_df_arrbrgr)+0.5),
     ylim=c(min(c(milestone_df_arrbrgr$slope_lci,milestone_df_arrbrgr_abs$slope_lci)),
            max(c(milestone_df_arrbrgr$slope_uci,milestone_df_arrbrgr_abs$slope_uci))),
     col="white",xaxt="n",xlab="",pch=19,ylab="Slope coefficient",main="")
points(which(milestone_df_arrbrgr$sign==F)-0.1,milestone_df_arrbrgr$slope[which(milestone_df_arrbrgr$sign==F)],col="cornflowerblue")
points(which(milestone_df_arrbrgr$sign)-0.1,milestone_df_arrbrgr$slope[which(milestone_df_arrbrgr$sign)],col="cornflowerblue",pch=19)
arrows((1:nrow(milestone_df_arrbrgr))-0.1,milestone_df_arrbrgr$slope_lci,(1:nrow(milestone_df_arrbrgr))-0.1,milestone_df_arrbrgr$slope_uci,length=0,col="cornflowerblue")
points(which(milestone_df_arrbrgr_abs$sign==F)+0.1,milestone_df_arrbrgr_abs$slope[which(milestone_df_arrbrgr_abs$sign==F)],col="chocolate")
points(which(milestone_df_arrbrgr_abs$sign)+0.1,milestone_df_arrbrgr_abs$slope[which(milestone_df_arrbrgr_abs$sign)],pch=19,col="chocolate")
arrows((1:nrow(milestone_df_arrbrgr_abs))+0.1,milestone_df_arrbrgr_abs$slope_lci,(1:nrow(milestone_df_arrbrgr_abs))+0.1,milestone_df_arrbrgr_abs$slope_uci,length=0,col="chocolate")
abline(h=0,lty=3)
axis(1, at=1:nrow(milestone_df_arrbrgr),labels=milestone_df_arrbrgr$from,las=2)
legend("topleft",c("Within-individual anomaly","Absolute date"),col=c("cornflowerblue","chocolate"),pch=c(19,19),bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Make pretty table for paper
milestone_name_pairs_abs$sd<-apply(mig_timing_df[,milestone_name_pairs_abs$from],2,sd,na.rm=T)
milestone_name_pairs_abs$delta_sd<-apply(mig_timing_df[,milestone_name_pairs_abs$to],2,sd,na.rm=T)-milestone_name_pairs_abs$sd
milestone_name_pairs_abs$var<-apply(mig_timing_df[,milestone_name_pairs_abs$from],2,var,na.rm=T)
milestone_name_pairs_abs$delta_var<-apply(mig_timing_df[,milestone_name_pairs_abs$to],2,var,na.rm=T)-milestone_name_pairs_abs$var

#F-test: are the variances significantly different between milestones?
milestone_name_pairs_abs$var_sig_diff<-NA
for(i in 1:nrow(milestone_name_pairs_abs)){
  temp_f_test<-var.test(mig_timing_df_stdind[,milestone_name_pairs_abs$from[i]],mig_timing_df_stdind[,milestone_name_pairs_abs$to[i]])
  milestone_name_pairs_abs$var_sig_diff[i]<-temp_f_test$p.value<0.05
}

var_sig_diff_df_abs<-data.frame(matrix(nrow=length(milestone_names),ncol=length(milestone_names)))
row.names(var_sig_diff_df_abs)<-names(var_sig_diff_df_abs)<-milestone_names
for(i in 1:length(milestone_names)){
  for(j in 1:length(milestone_names)){
    temp_var_test<-var.test(mig_timing_df[,row.names(var_sig_diff_df_abs)[i]],mig_timing_df[,names(var_sig_diff_df_abs)[j]])
    var_sig_diff_df_abs[i,j]<-temp_var_test$p.value<0.05
  }
}

ms_table_abs<-milestone_name_pairs_abs[,c("from","var","to","delta_var","var_sig_diff","R2m","sign")]
ms_table_abs$var<-round(milestone_name_pairs_abs$var,digits=2)
ms_table_abs$delta_var<-round(milestone_name_pairs_abs$delta_var,digits=2)
ms_table_abs$R2m<-round(milestone_name_pairs_abs$R2m,digits=2)
ms_table_abs$R2m_arrbrgr_sign<-ms_table_abs$R2m_arrbrgr<-NA
for(i in 1:nrow(milestone_df_arrbrgr_abs)){
  ms_table_abs[ms_table_abs$from==milestone_df_arrbrgr_abs[i,]$from,]$R2m_arrbrgr<-round(milestone_df_arrbrgr_abs[i,]$R2m,digits=2)
  ms_table_abs[ms_table_abs$from==milestone_df_arrbrgr_abs[i,]$from,]$R2m_arrbrgr_sign<-milestone_df_arrbrgr_abs[i,]$sign
}
write.csv(ms_table_abs,"H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/ms_table_abs.csv",row.names=F)


#Relationship between anomalies - by migratory direction - SW first
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig4_21.png",10,10,res=800,units="in")
par(mfrow=c(4,3))
for(i in 1:nrow(milestone_name_pairs)){
  print(i)
  temp_data_df<-data.frame("t"=mig_timing_df_stdind[,milestone_name_pairs$from[i]],
                           "tplus1"=mig_timing_df_stdind[,milestone_name_pairs$to[i]],
                           "ptt"=mig_timing_df_stdind$ptt,
                           "region"=mig_timing_df_stdind$region,
                           "mig_dir"=mig_timing_df_stdind$mig_dir)
  temp_data_df<-temp_data_df[temp_data_df$mig_dir=="SW",]
  max_abs_vals<-max(abs(temp_data_df[,c("t","tplus1")]),na.rm=T)
  plot(temp_data_df$t,temp_data_df$tplus1,
       xlim=c(-max_abs_vals,max_abs_vals),
       ylim=c(-max_abs_vals,max_abs_vals),
       # main="Relationship between anomalies\nfor subsequent milestones in same year",
       bty="n",
       pch=21,
       bg="seagreen3",
       xlab=milestone_name_pairs$from[i],
       ylab=milestone_name_pairs$to[i])
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  if(length(which(apply(is.na(temp_data_df[,c("t","tplus1")]),1,any)==F))<10) next
  m2<-brm(tplus1~t+(1|region),
          data=temp_data_df,
          control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
  newx<-seq(-max_abs_vals,max_abs_vals)
  intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
  t_draws<-as_draws_df(m2,variable="b_t")$b_t
  predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
  for(i in 1:nrow(predicted_lines_df)){
    predicted_lines_df[i,]<-(newx*t_draws[i])+intercept_draws[i]
  }
  lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
  lines(newx,apply(predicted_lines_df,2,mean))
  lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
  
  R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
  text(-max_abs_vals,max_abs_vals,paste0("R2m = ",R2m),pos=4)
}

#Add-on for arr_brgr vs dep_brgr
i<-nrow(milestone_name_pairs)
arrdepbrgr_df_SW<-arrdepbrgr_df[arrdepbrgr_df$mig_dir_1st_cycle=="SW",]
max_abs_vals<-max(abs(arrdepbrgr_df_SW[,c("arr_brgr","dep_brgr")]),na.rm=T)
plot(arrdepbrgr_df_SW$arr_brgr,arrdepbrgr_df_SW$dep_brgr,
     xlim=c(-max_abs_vals,max_abs_vals),
     ylim=c(-max_abs_vals,max_abs_vals),
     # main="Relationship between anomalies\nfor subsequent milestones in same year",
     bty="n",
     pch=21,
     bg="seagreen3",
     xlab="arr_brgr",
     ylab="dep_brgr")
abline(h=0,lty=3)
abline(v=0,lty=3)
m2<-brm(dep_brgr~arr_brgr+(1|region/ptt),
        data=arrdepbrgr_df_SW,
        control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
sm2<-summary(m2)
R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
newx<-seq(-max_abs_vals,max_abs_vals)
intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
arr_brgr_draws<-as_draws_df(m2,variable="b_arr_brgr")$b_arr_brgr
predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
for(i in 1:nrow(predicted_lines_df)){
  predicted_lines_df[i,]<-(newx*arr_brgr_draws[i])+intercept_draws[i]
}
lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
lines(newx,apply(predicted_lines_df,2,mean))
lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
text(-max_abs_vals,max_abs_vals,paste0("R2m = ",R2m),pos=4)

par(mfrow=c(1,1))
dev.off()

#Now SE
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig4_22.png",10,10,res=800,units="in")
par(mfrow=c(4,3))
for(i in 1:nrow(milestone_name_pairs)){
  print(i)
  temp_data_df<-data.frame("t"=mig_timing_df_stdind[,milestone_name_pairs$from[i]],
                           "tplus1"=mig_timing_df_stdind[,milestone_name_pairs$to[i]],
                           "ptt"=mig_timing_df_stdind$ptt,
                           "region"=mig_timing_df_stdind$region,
                           "mig_dir"=mig_timing_df_stdind$mig_dir)
  temp_data_df<-temp_data_df[temp_data_df$mig_dir=="SE",]
  max_abs_vals<-max(abs(temp_data_df[,c("t","tplus1")]),na.rm=T)
  plot(temp_data_df$t,temp_data_df$tplus1,
       xlim=c(-max_abs_vals,max_abs_vals),
       ylim=c(-max_abs_vals,max_abs_vals),
       # main="Relationship between anomalies\nfor subsequent milestones in same year",
       bty="n",
       pch=21,
       bg="seagreen3",
       xlab=milestone_name_pairs$from[i],
       ylab=milestone_name_pairs$to[i])
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  if(length(which(apply(is.na(temp_data_df[,c("t","tplus1")]),1,any)==F))<10) next
  m2<-brm(tplus1~t+(1|region),
          data=temp_data_df,
          control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
  newx<-seq(-max_abs_vals,max_abs_vals)
  intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
  t_draws<-as_draws_df(m2,variable="b_t")$b_t
  predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
  for(i in 1:nrow(predicted_lines_df)){
    predicted_lines_df[i,]<-(newx*t_draws[i])+intercept_draws[i]
  }
  lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
  lines(newx,apply(predicted_lines_df,2,mean))
  lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
  
  R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
  text(-max_abs_vals,max_abs_vals,paste0("R2m = ",R2m),pos=4)
}
#Add-on for arr_brgr vs dep_brgr
i<-nrow(milestone_name_pairs)
arrdepbrgr_df_SE<-arrdepbrgr_df[arrdepbrgr_df$mig_dir_1st_cycle=="SE",]
max_abs_vals<-max(abs(arrdepbrgr_df_SE[,c("arr_brgr","dep_brgr")]),na.rm=T)
plot(arrdepbrgr_df_SE$arr_brgr,arrdepbrgr_df_SE$dep_brgr,
     xlim=c(-max_abs_vals,max_abs_vals),
     ylim=c(-max_abs_vals,max_abs_vals),
     # main="Relationship between anomalies\nfor subsequent milestones in same year",
     bty="n",
     pch=21,
     bg="seagreen3",
     xlab="arr_brgr",
     ylab="dep_brgr")
abline(h=0,lty=3)
abline(v=0,lty=3)
m2<-brm(dep_brgr~arr_brgr+(1|region/ptt),
        data=arrdepbrgr_df_SE,
        control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
sm2<-summary(m2)
R2m<-round(bayes_R2(m2,re.form=NA)[1],digits=3)
newx<-seq(-max_abs_vals,max_abs_vals)
intercept_draws<-as_draws_df(m2,variable="b_Intercept")$b_Intercept
arr_brgr_draws<-as_draws_df(m2,variable="b_arr_brgr")$b_arr_brgr
predicted_lines_df<-data.frame(matrix(ncol=length(newx),nrow=length(intercept_draws)))
for(i in 1:nrow(predicted_lines_df)){
  predicted_lines_df[i,]<-(newx*arr_brgr_draws[i])+intercept_draws[i]
}
lines(newx,apply(predicted_lines_df,2,quantile,0.025),lty=2)
lines(newx,apply(predicted_lines_df,2,mean))
lines(newx,apply(predicted_lines_df,2,quantile,0.975),lty=2)
text(-max_abs_vals,max_abs_vals,paste0("R2m = ",R2m),pos=4)
par(mfrow=c(1,1))
dev.off()

#Check to see whether the interaction between migratory direction and last anomaly is ever significant, for any milestone.
milestone_name_pairs_migdirint<-data.frame("from"=c("dep_brgr","dep_UK","fin_sah_sb","fin_sah_sb","arr_wgr_a","arr_wgr_b",
                                                    "dep_wgr_a","dep_wgr_b","dep_wa","arr_UK","arr_brgr"),
                                           "to"=c("dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b",
                                                  "dep_wa","dep_wa","arr_UK","arr_brgr","dep_brgr"))
milestone_name_pairs_migdirint$int_slope<-milestone_name_pairs_migdirint$int_slope_lci<-milestone_name_pairs_migdirint$int_slope_uci<-milestone_name_pairs_migdirint$int_sign<-NA

for(i in 1:length(milestone_names)){ #NB ONLY GO AS FAR AS ARR_UK - THE ARR_BRGR:DEP_BRGR STEP HAS TO BE DONE MANUALLY
  print(i)
  temp_data_df<-data.frame("t"=mig_timing_df_stdind[,milestone_name_pairs_migdirint$from[i]],
                           "tplus1"=mig_timing_df_stdind[,milestone_name_pairs_migdirint$to[i]],
                           "mig_dir"=mig_timing_df_stdind$mig_dir,
                           "ptt"=mig_timing_df_stdind$ptt,
                           "region"=mig_timing_df_stdind$region)
  m2<-brm(tplus1~t*mig_dir+(1|region/ptt),
          data=temp_data_df,
          control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
  sm2<-summary(m2)
  milestone_name_pairs_migdirint$int_slope[i]<-sm2$fixed[row.names(sm2$fixed)=="t:mig_dirSW",]$Estimate
  milestone_name_pairs_migdirint$int_slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="t:mig_dirSW",]$'l-95% CI'
  milestone_name_pairs_migdirint$int_slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="t:mig_dirSW",]$'u-95% CI'
  milestone_name_pairs_migdirint$int_sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="t:mig_dirSW",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="t:mig_dirSW",]$'u-95% CI'>0)==F
}

#Add-on for arr_brgr vs dep_brgr
i<-nrow(milestone_name_pairs_migdirint)
m2<-brm(dep_brgr~arr_brgr*mig_dir_1st_cycle+(1|region/ptt),
        data=arrdepbrgr_df,
        control = list(adapt_delta = 0.99999999999,max_treedepth=20));summary(m2)
sm2<-summary(m2)
milestone_name_pairs_migdirint$int_slope[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr:mig_dir_1st_cycleSW",]$Estimate
milestone_name_pairs_migdirint$int_slope_lci[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr:mig_dir_1st_cycleSW",]$'l-95% CI'
milestone_name_pairs_migdirint$int_slope_uci[i]<-sm2$fixed[row.names(sm2$fixed)=="arr_brgr:mig_dir_1st_cycleSW",]$'u-95% CI'
milestone_name_pairs_migdirint$int_sign[i]<-(sm2$fixed[row.names(sm2$fixed)=="arr_brgr:mig_dir_1st_cycleSW",]$'l-95% CI'<0 & sm2$fixed[row.names(sm2$fixed)=="arr_brgr:mig_dir_1st_cycleSW",]$'u-95% CI'>0)==F


#################################################
##Effect of location of last stopover N of Sahara on timing of completion of Sahara crossing
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/Fig4a.png",6,6,res=800,units="in")
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(mig_timing_df$lat_last_euro_so_sbound,mig_timing_df$fin_sah_sb,xlab="Latitude of last stopover N of Sahara",ylab="fin_sah_sb")
abline(v=46.2,lty=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#Is there a relationship between latitude of the last stopover in Europe and the date of completion of the Sahara crossing southbound?
m1<-lm(mig_timing_df$fin_sah_sb~mig_timing_df$lat_last_euro_so_sbound);summary(m1);confint(m1)
m2<-brm(fin_sah_sb~lat_last_euro_so_sbound+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=mig_timing_df);summary(m2)

#Do birds either side of 46.2N differ in their timing?
m1<-lm(mig_timing_df$fin_sah_sb~(mig_timing_df$lat_last_euro_so_sbound>46.2));summary(m1);confint(m1)
mig_timing_df$lat_last_euro_so_sbound_462<-mig_timing_df$lat_last_euro_so_sbound>46.2
m2<-brm(fin_sah_sb~lat_last_euro_so_sbound_462+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=mig_timing_df);summary(m2)

#Can this be explained by individual flexibility - do we see the same relationship within individuals?
#Which birds *ever* stopped north of 46.2N?
n_stopping_ptts<-mig_timing_df_stdind[which(mig_timing_df_stdind$lat_last_euro_so_sbound>46.2 & is.na(mig_timing_df_stdind$fin_sah_sb)==F),]$ptt
n_stopping_ptts_df<-mig_timing_df_stdind[mig_timing_df_stdind$ptt %in% n_stopping_ptts,c("ptt","year","fin_sah_sb","lat_last_euro_so_sbound")]
n_stopping_ptts_df<-n_stopping_ptts_df[is.na(n_stopping_ptts_df$fin_sah_sb)==F,]
n_stopping_ptts_df$lat_last_euro_so_sbound_462<-n_stopping_ptts_df$lat_last_euro_so_sbound>46.2

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/Fig4b.png",6,6,res=800,units="in")
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(n_stopping_ptts_df$lat_last_euro_so_sbound,n_stopping_ptts_df$fin_sah_sb,col=n_stopping_ptts_df$ptt,
     pch=20,cex=2,xlab="Latitude of last stopover N of Sahara",ylab="fin_sah_sb within-individual anomaly")
abline(v=46.2,lty=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

m1<-lm(fin_sah_sb~lat_last_euro_so_sbound_462,data=n_stopping_ptts_df);summary(m1);confint(m1)
m2<-brm(fin_sah_sb~lat_last_euro_so_sbound_462+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=n_stopping_ptts_df);summary(m2) #CI of coefficient for lat_last_euro_so_sbound_462 doesn't overlap zero.

#What is the migratory direction/breeding habitat of birds that made their first stopover further north?
n_stopping_ptts_df$fin_sah_sb_abs<-n_stopping_ptts_df$mig_dir<-n_stopping_ptts_df$br_hab<-NA
for(i in 1:nrow(n_stopping_ptts_df)){
  n_stopping_ptts_df$mig_dir[i]<-mig_timing_df[mig_timing_df$ptt==n_stopping_ptts_df[i,]$ptt & mig_timing_df$year==n_stopping_ptts_df[i,]$year,]$mig_dir
  n_stopping_ptts_df$br_hab[i]<-mig_timing_df[mig_timing_df$ptt==n_stopping_ptts_df[i,]$ptt & mig_timing_df$year==n_stopping_ptts_df[i,]$year,]$br_hab
  n_stopping_ptts_df$fin_sah_sb_abs[i]<-mig_timing_df[mig_timing_df$ptt==n_stopping_ptts_df[i,]$ptt & mig_timing_df$year==n_stopping_ptts_df[i,]$year,]$fin_sah_sb
}

#Now look at absolute timing for these individuals
plot(n_stopping_ptts_df$lat_last_euro_so_sbound,n_stopping_ptts_df$fin_sah_sb_abs,col=n_stopping_ptts_df$ptt,
     pch=20,cex=2,xlab="Latitude of last stopover N of Sahara",ylab="fin_sah_sb within-individual anomaly")
abline(v=46.2,lty=2)

m1<-lm(fin_sah_sb_abs~lat_last_euro_so_sbound_462,data=n_stopping_ptts_df);summary(m1);confint(m1)
m2<-brm(fin_sah_sb_abs~lat_last_euro_so_sbound_462+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=n_stopping_ptts_df);summary(m2) #CI of coefficient for lat_last_euro_so_sbound_462 doesn't overlap zero.


#################################################
#Effect of location of last fix/stopover in West Africa on timing of departure from West Africa

#Fix first
plot(mig_timing_df$last_wa_fix_long,mig_timing_df$dep_wa,ylab="Departure from West Africa",xlab="Longitude of last fix in W Africa")
m1<-lm(mig_timing_df$dep_wa~mig_timing_df$last_wa_fix_long);summary(m1)
m2<-brm(dep_wa~last_wa_fix_long+last_wa_fix_lat+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=mig_timing_df);summary(m2)
plot(mig_timing_df$last_wa_fix_lat,mig_timing_df$dep_wa,ylab="Departure from West Africa",xlab="Latitude of last fix in W Africa")
m1<-lm(mig_timing_df$dep_wa~mig_timing_df$last_wa_fix_lat);summary(m1)
m2<-brm(dep_wa~last_wa_fix_lat+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=mig_timing_df);summary(m2)

plot(mig_timing_df_stdind$last_wa_fix_long,mig_timing_df_stdind$dep_wa,ylab="Departure from West Africa (within-individual anomaly)",xlab="Longitude of last fix in W Africa")
m1<-lm(mig_timing_df_stdind$dep_wa~mig_timing_df_stdind$last_wa_fix_long);summary(m1)
plot(mig_timing_df_stdind$last_wa_fix_lat,mig_timing_df_stdind$dep_wa,ylab="Departure from West Africa (within-individual anomaly)",xlab="Latitude of last fix in W Africa")
m1<-lm(mig_timing_df_stdind$dep_wa~mig_timing_df_stdind$last_wa_fix_lat);summary(m1)
m1<-lm(mig_timing_df_stdind$dep_wa~mig_timing_df_stdind$last_wa_fix_long+mig_timing_df_stdind$last_wa_fix_lat);summary(m1)

temp_std_depwac<-mig_timing_df$dep_wa-min(mig_timing_df$dep_wa,na.rm=T)
temp_std_depwac_submax<-temp_std_depwac/max(temp_std_depwac,na.rm=T)
par(mar=c(5.1, 4.1, 0.1, 0.1))
plot(last_wa_fix_lat~last_wa_fix_long,
     xlab="",
     ylab="",
     data=mig_timing_df[is.na(mig_timing_df$dep_wa)==F,],
     col=rgb(temp_std_depwac_submax[is.na(mig_timing_df$dep_wa)==F],0,1-temp_std_depwac_submax[is.na(mig_timing_df$dep_wa)==F],1),pch=20,cex=2)
plot(world_shp,add=T)
plot(sah_spdf,add=T,col="yellow")

#Now last stopover
plot(mig_timing_df$last_wa_sos_long,mig_timing_df$dep_wa,ylab="Departure from West Africa",xlab="Longitude of last stopover in W Africa")
m1<-lm(mig_timing_df$dep_wa~mig_timing_df$last_wa_sos_long);summary(m1)
m2<-brm(dep_wa~last_wa_sos_long+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=mig_timing_df);summary(m2)
plot(mig_timing_df$last_wa_sos_lat,mig_timing_df$dep_wa,ylab="Departure from West Africa",xlab="Latitude of last stopover in W Africa")
m1<-lm(mig_timing_df$dep_wa~mig_timing_df$last_wa_sos_lat);summary(m1)
m2<-brm(dep_wa~last_wa_sos_lat+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=mig_timing_df);summary(m2)
m2<-brm(dep_wa~last_wa_sos_long+last_wa_sos_lat+(1|ptt),
        control = list(adapt_delta = 0.99999999999,max_treedepth=20),
        data=mig_timing_df);summary(m2)

plot(dep_wa~last_wa_sos_long,data=mig_timing_df_stdind[is.na(mig_timing_df_stdind$dep_wa)==F,],ylab="Departure from West Africa (within-individual anomaly)",xlab="Longitude of last stopover in W Africa")
m1<-lm(mig_timing_df_stdind$dep_wa~mig_timing_df_stdind$last_wa_sos_long);summary(m1);confint(m1)
plot(dep_wa~last_wa_sos_lat,data=mig_timing_df_stdind[is.na(mig_timing_df_stdind$dep_wa)==F,],ylab="Departure from West Africa (within-individual anomaly)",xlab="Latitude of last stopover in W Africa")
m1<-lm(mig_timing_df_stdind$dep_wa~mig_timing_df_stdind$last_wa_sos_lat);summary(m1);confint(m1)
m1<-lm(mig_timing_df_stdind$dep_wa~mig_timing_df_stdind$last_wa_sos_long+mig_timing_df_stdind$last_wa_sos_lat);summary(m1);confint(m1)

temp_std_depwac<-mig_timing_df$dep_wa-min(mig_timing_df$dep_wa,na.rm=T)
temp_std_depwac_submax<-temp_std_depwac/max(temp_std_depwac,na.rm=T)
par(mar=c(5.1, 4.1, 0.1, 0.1))
plot(last_wa_sos_lat~last_wa_sos_long,
     xlab="",
     ylab="",
     data=mig_timing_df[is.na(mig_timing_df$dep_wa)==F,],
     col=rgb(temp_std_depwac_submax[is.na(mig_timing_df$dep_wa)==F],0,1-temp_std_depwac_submax[is.na(mig_timing_df$dep_wa)==F],1),pch=20,cex=2)
plot(world_shp,add=T)
plot(sah_spdf,add=T,col="yellow")


########################################################################################################################
#4D: CONSEQUENCES OF VARIATION IN MILESTONES
#Were birds that left the breeding grounds later more likely to die?
boxplot(dep_brgr~died_this_cycle,data=mig_timing_df,ylab="dep_brgr, absolute date")
boxplot(dep_brgr~died_this_cycle,data=mig_timing_df_stdind,ylab="dep_brgr, within-individual anomaly")
boxplot(dep_brgr~died_this_cycle,data=mig_timing_df_std,ylab="dep_brgr, relative to population median")

boxplot(dep_brgr~season_died,data=mig_timing_df,ylab="Departure from the breeding grounds (absolute)")
boxplot(dep_brgr~season_died,data=mig_timing_df_stdind,ylab="Departure from the breeding grounds (within-ind anomaly)")
boxplot(dep_brgr~season_died,data=mig_timing_df_std,ylab="Departure from the breeding grounds (relative to pop median)")


#Relationship between anomaly and death
mtd_stdind_melted$died<-"Survived"
mtd_stdind_melted$died[which(mtd_stdind_melted$last_milestone_pre_death==mtd_stdind_melted$variable)]<-"Died"
mtd_stdind_melted$died_bin<-0
mtd_stdind_melted$died_bin[mtd_stdind_melted$died=="Died"]<-1
mtd_stdind_melted$season<-NA
mtd_stdind_melted[mtd_stdind_melted$variable %in% c("dep_brgr","dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b"),]$season<-"Autumn"
mtd_stdind_melted[mtd_stdind_melted$variable %in% c("dep_wgr_b","dep_wgr_a","dep_wa","arr_UK","arr_brgr"),]$season<-"Spring"

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig6_1.png",7,7,res=800,units="in")
par(mar=c(10.1, 4.1, 4.1, 2.1))
boxplot(value~died*variable,data=mtd_stdind_melted,col=c("brown","green"),las=2,xlab="",
        ylab="Within-individual anomaly")
legend("topright",fill=c("brown","green"),legend=c("Died after milestone","Survived to next milestone"))
abline(h=0,lty=3)
dev.off()

boxplot(value~died*variable,main="SE",data=mtd_stdind_melted[mtd_stdind_melted$mig_dir=="SE",],col=c("brown","green"),las=2,xlab="")
abline(h=0,lty=3)
boxplot(value~died*variable,main="SW",data=mtd_stdind_melted[mtd_stdind_melted$mig_dir=="SW",],col=c("brown","green"),las=2,xlab="")
abline(h=0,lty=3)
par(mar=c(5.1, 4.1, 4.1, 2.1))

m1<-glm(died_bin~value,family="binomial",data=mtd_stdind_melted);summary(m1)
m0<-glm(died_bin~1,family="binomial",data=mtd_stdind_melted);summary(m0)
m1<-glm(died_bin~value*variable,family="binomial",data=mtd_stdind_melted);summary(m1)


par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(as.numeric(table(mtd_stdind_melted$variable[mtd_stdind_melted$died=="Died" & mtd_stdind_melted$mig_dir=="SE"])/
                  table(mtd_stdind_melted$variable[mtd_stdind_melted$died=="Survived" & mtd_stdind_melted$mig_dir=="SE"])),
     col="darkorange1",
     xaxt="n",
     xlab="",
     bty="n",
     ylab="Proportion of birds dying after milestone",
     pch=19,
     ylim=c(0,0.3))
points(as.numeric(table(mtd_stdind_melted$variable[mtd_stdind_melted$died=="Died" & mtd_stdind_melted$mig_dir=="SW"])/
                    table(mtd_stdind_melted$variable[mtd_stdind_melted$died=="Survived" & mtd_stdind_melted$mig_dir=="SW"])),
       col="blue",pch=19)
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
legend("topright",legend=c("SE","SW"),pch=c(19,19),col=c("darkorange1","blue"))


#Repeat death analyses, but using across-population anomaly rather than within-individual anomaly.
mtd_melted_stdmil<-melt(mig_timing_df_std[,c("ptt","year","mig_dir","known_age_binary","last_milestone_pre_death","outcome",milestone_names)],id=c("ptt","year","mig_dir","known_age_binary","last_milestone_pre_death","outcome"))
mtd_melted_stdmil<-mtd_melted_stdmil[which(is.na(mtd_melted_stdmil$value)==F),]
mtd_melted_stdmil$died<-"Survived"
mtd_melted_stdmil$died[which(mtd_melted_stdmil$last_milestone_pre_death==mtd_melted_stdmil$variable)]<-"Died"
mtd_melted_stdmil$died_bin<-0
mtd_melted_stdmil$died_bin[mtd_melted_stdmil$died=="Died"]<-1
mtd_melted_stdmil$season<-NA
mtd_melted_stdmil[mtd_melted_stdmil$variable %in% c("dep_brgr","dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b"),]$season<-"Autumn"
mtd_melted_stdmil[mtd_melted_stdmil$variable %in% c("dep_wgr_b","dep_wgr_a","dep_wa","arr_UK","arr_brgr"),]$season<-"Spring"


png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig6_2.png",7,7,res=800,units="in")
par(mar=c(10.1, 4.1, 4.1, 2.1))
boxplot(value~died*variable,data=mtd_melted_stdmil,col=c("brown","green"),las=2,xlab="",
        ylab="Within-population anomaly")
legend("topright",fill=c("brown","green"),legend=c("Died after milestone","Survived to next milestone"))
abline(h=0,lty=3)
dev.off()

m1<-glm(died_bin~value,family="binomial",data=mtd_melted_stdmil);summary(m1) #No overall effect of population anomaly on death probability
m0<-glm(died_bin~1,family="binomial",data=mtd_melted_stdmil);summary(m0)
m1<-glm(died_bin~value*variable,family="binomial",data=mtd_melted_stdmil);summary(m1) #AIC is lower than m0, but no effect of anomaly on death probability for individual milestones

#Effect of timing on mortality conditional on season
#Within-individual anomaly first
#Plot raw data. Boxplot no good here because only three birds w/repeated measures died in spring.
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig6_3.png",7,7,res=800,units="in")
plot(rep(1,length(which(mtd_stdind_melted$season=="Autumn" & mtd_stdind_melted$died=="Died"))),
     mtd_stdind_melted[mtd_stdind_melted$season=="Autumn" & mtd_stdind_melted$died=="Died",]$value,
     ylim=range(mtd_stdind_melted$value),xlim=c(0.5,4.5),
     ylab="Within-individual anomaly",xaxt="n",xlab="Season")
points(rep(2,length(which(mtd_stdind_melted$season=="Autumn" & mtd_stdind_melted$died=="Survived"))),
       mtd_stdind_melted[mtd_stdind_melted$season=="Autumn" & mtd_stdind_melted$died=="Survived",]$value,col="green")
points(rep(3,length(which(mtd_stdind_melted$season=="Spring" & mtd_stdind_melted$died=="Died"))),
       mtd_stdind_melted[mtd_stdind_melted$season=="Spring" & mtd_stdind_melted$died=="Died",]$value)
points(rep(4,length(which(mtd_stdind_melted$season=="Spring" & mtd_stdind_melted$died=="Survived"))),
       mtd_stdind_melted[mtd_stdind_melted$season=="Spring" & mtd_stdind_melted$died=="Survived",]$value,col="green")
axis(1, at=c(1.5,3.5),labels=c("Autumn","Spring"))
abline(h=0,lty=3)
legend("bottomright",col=c("black","green"),pch=c(1,1),legend=c("Died after milestone","Survived to next milestone"))
dev.off()

#Model effect of within-individual anomaly on mortality probability
m1<-glm(died_bin~value,family="binomial",data=mtd_stdind_melted[mtd_stdind_melted$season=="Spring",]);summary(m1)
m1_xseq<-seq(floor(min(mtd_stdind_melted$value)),ceiling(max(mtd_stdind_melted$value)))
pm1<-predict(m1,newdata=data.frame(value=m1_xseq),type="link",se.fit=T)
critval <- 1.96 ## approx 95% CI
upr <- pm1$fit + (critval * pm1$se.fit)
lwr <- pm1$fit - (critval * pm1$se.fit)
fit <- pm1$fit
upr_probs<-1/(1+exp(-upr))
lwr_probs<-1/(1+exp(-lwr))
fit_probs<-1/(1+exp(-fit))

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig6_5.png",7,7,res=800,units="in")
plot(died_bin~value,data=mtd_stdind_melted[mtd_stdind_melted$season=="Spring",],main="Spring",
     xlab="Within-individual anomaly",
     ylab="Probability of dying before next milestone")
lines(m1_xseq,fit_probs)
lines(m1_xseq,upr_probs,lty=2)
lines(m1_xseq,lwr_probs,lty=2)
dev.off()

#Now population anomaly
#Plot raw data
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig6_4.png",7,7,res=800,units="in")
plot(rep(1,length(which(mtd_melted_stdmil$season=="Autumn" & mtd_melted_stdmil$died=="Died"))),
     mtd_melted_stdmil[mtd_melted_stdmil$season=="Autumn" & mtd_melted_stdmil$died=="Died",]$value,
     ylim=range(mtd_melted_stdmil$value),xlim=c(0.5,4.5),
     ylab="Anomaly relative to population median",xaxt="n",xlab="Season")
points(rep(2,length(which(mtd_melted_stdmil$season=="Autumn" & mtd_melted_stdmil$died=="Survived"))),
       mtd_melted_stdmil[mtd_melted_stdmil$season=="Autumn" & mtd_melted_stdmil$died=="Survived",]$value,col="green")
points(rep(3,length(which(mtd_melted_stdmil$season=="Spring" & mtd_melted_stdmil$died=="Died"))),
       mtd_melted_stdmil[mtd_melted_stdmil$season=="Spring" & mtd_melted_stdmil$died=="Died",]$value)
points(rep(4,length(which(mtd_melted_stdmil$season=="Spring" & mtd_melted_stdmil$died=="Survived"))),
       mtd_melted_stdmil[mtd_melted_stdmil$season=="Spring" & mtd_melted_stdmil$died=="Survived",]$value,col="green")
axis(1, at=c(1.5,3.5),labels=c("Autumn","Spring"))
abline(h=0,lty=3)
legend("topright",col=c("black","green"),pch=c(1,1),legend=c("Died after milestone","Survived to next milestone"))
dev.off()

#Now plot by milestone
lmn<-length(milestone_names)
died_x_seq<-seq(1,3*lmn,by=3)
surv_x_seq<-died_x_seq+1
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig6_7.png",8,6,res=800,units="in")
par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(rep(0,lmn*3),ylim=range(mtd_melted_stdmil$value),col="white",xaxt="n",xlab="",ylab="Anomaly relative to population median")
for(i in 1:lmn){
  temp_data<-mtd_melted_stdmil[mtd_melted_stdmil$variable==milestone_names[i] &
                                 (mtd_melted_stdmil$outcome %in% c("UB","UC"))==F,]
  points(rep(died_x_seq[i],nrow(temp_data[temp_data$died=="Died",])),temp_data[temp_data$died=="Died",]$value)
  points(rep(surv_x_seq[i],nrow(temp_data[temp_data$died=="Survived",])),temp_data[temp_data$died=="Survived",]$value,col="green")
}
axis(1, at=died_x_seq+0.5,labels=milestone_names,las=2)
abline(h=0,lty=3)
legend("topright",col=c("black","green"),pch=c(1,1),legend=c("Died after milestone","Survived to next milestone"))
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Model effect of population anomaly on mortality probability
m1<-glm(died_bin~value,family="binomial",data=mtd_melted_stdmil[mtd_melted_stdmil$season=="Spring",]);summary(m1)
m1_xseq<-seq(floor(min(mtd_melted_stdmil$value)),ceiling(max(mtd_melted_stdmil$value)))
pm1<-predict(m1,newdata=data.frame(value=m1_xseq),type="link",se.fit=T)
critval <- 1.96 ## approx 95% CI
upr <- pm1$fit + (critval * pm1$se.fit)
lwr <- pm1$fit - (critval * pm1$se.fit)
fit <- pm1$fit
upr_probs<-1/(1+exp(-upr))
lwr_probs<-1/(1+exp(-lwr))
fit_probs<-1/(1+exp(-fit))

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig6_6.png",7,7,res=800,units="in")
plot(died_bin~value,data=mtd_melted_stdmil[mtd_melted_stdmil$season=="Spring",],main="Spring",
     xlab="Date anomaly relative to population median",
     ylab="Probability of dying before next milestone")
lines(m1_xseq,fit_probs)
lines(m1_xseq,upr_probs,lty=2)
lines(m1_xseq,lwr_probs,lty=2)
dev.off()


m1<-glm(died_bin~value,data=mtd_melted_stdmil[mtd_melted_stdmil$season=="Spring" &
                                                (mtd_melted_stdmil$outcome %in% c("UA","UB","UC"))==F ,],family="binomial");summary(m1)
m1<-glmer(died_bin~value+(1|variable)+(1|year)+(1|ptt),data=mtd_melted_stdmil[mtd_melted_stdmil$season=="Spring",],family="binomial");summary(m1)
m1<-brm(died_bin~value+(1|ptt),
        control = list(adapt_delta = 0.99999999999,
                       max_treedepth=20),
        data=mtd_melted_stdmil[mtd_melted_stdmil$season=="Spring" &
                                 (mtd_melted_stdmil$outcome %in% c("UA","UB","UC"))==F ,],family="bernoulli");summary(m1)
#Standardise by SD as well to take into account of variable variance between milestones

#Now for scaled dataset
mtd_melted_scaled<-melt(mig_timing_df_scaled[,c("ptt","year","mig_dir","br_hab","region","known_age_binary","last_milestone_pre_death","outcome",milestone_names)],id=c("ptt","year","mig_dir","br_hab","region","known_age_binary","last_milestone_pre_death","outcome"))
mtd_melted_scaled<-mtd_melted_scaled[which(is.na(mtd_melted_scaled$value)==F),]
mtd_melted_scaled$died<-"Survived"
mtd_melted_scaled$died[which(mtd_melted_scaled$last_milestone_pre_death==mtd_melted_scaled$variable)]<-"Died"
mtd_melted_scaled$died<-relevel(as.factor(mtd_melted_scaled$died),ref = "Survived") #Put Survived first - more data for intercept
mtd_melted_scaled$died_bin<-0
mtd_melted_scaled$died_bin[mtd_melted_scaled$died=="Died"]<-1
mtd_melted_scaled$season<-NA
mtd_melted_scaled[mtd_melted_scaled$variable %in% c("dep_brgr","dep_UK","fin_sah_sb","arr_wgr_a","arr_wgr_b"),]$season<-"Autumn"
mtd_melted_scaled[mtd_melted_scaled$variable %in% c("arr_wgr_a","arr_wgr_b","dep_wa","arr_UK","arr_brgr"),]$season<-"Spring"

m1<-lmer(value~died+(1|year)+(1|region/ptt),
         data=mtd_melted_scaled[mtd_melted_scaled$season=="Spring" & 
                                  (mtd_melted_scaled$outcome %in% c("UB","UC"))==F ,]);summary(m1);confint(m1)
m2<-lmer(value~died+(1|year)+(1|region/ptt),
         data=mtd_melted_scaled[mtd_melted_scaled$season=="Spring" & 
                                  (mtd_melted_scaled$outcome %in% c("UA","UB","UC"))==F ,]);summary(m2);confint(m2)
