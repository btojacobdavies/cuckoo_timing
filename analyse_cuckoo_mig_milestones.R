#CLEAN, COMPILE AND ANALYSE CUCKOO DATA
#1. Load and clean data; 2. Define stopovers; 3. Define migratory milestones; 4. Analysis
#This file uses code from: 'cuckoo_compile_clean_movebank_30Jan.R', 'function to add stopover duration data for analysis.R', 'cuckoo_stopovers.R' and 'paper_analyses.R'

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

world_shp<-readOGR("H:/JGD H drive backup 160320/Raw data/GADM world/TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp")

#########################################################################################################
### 1. LOAD AND CLEAN DATA ###
#Lots of code here from cuckoo_compile_clean_movebank_30Jan.R
#Load raw data from Movebank. First, all fixes from loggers - files from Chris H
all_fixes_1<-read.csv("H:/JGD H drive backup 160320/Raw data/Cuckoos/Derived CSVs/BTO Cuckoo migration study (1).csv")
all_fixes_2<-read.csv("H:/JGD H drive backup 160320/Raw data/Cuckoos/Derived CSVs/BTO Cuckoos - Argos webservice feed (3).csv")
all_fixes_3<-read.csv("H:/JGD H drive backup 160320/Raw data/Cuckoos/Derived CSVs/BTO Cuckoos - Argos webservice feed_from 010120_for Jacob.csv")

#Remove all fixes that overlap between data files - otherwise can end up with multiple fixes per day
#Find the fixes in all_fixes_1 that have algorithm.marked.outlier==NA that would get dropped by this method
all_fixes_1_ptt_day<-paste0(all_fixes_1$tag.local.identifier," ",substr(all_fixes_1$timestamp,1,10))
all_fixes_2_ptt_day<-paste0(all_fixes_2$tag.local.identifier," ",substr(all_fixes_2$timestamp,1,10))
af2pd<-unique(all_fixes_1_ptt_day[is.na(all_fixes_1$algorithm.marked.outlier)==F])
af3pd<-unique(all_fixes_2_ptt_day[is.na(all_fixes_3$algorithm.marked.outlier)==F])

af1_af2_valid_ids<-all_fixes_1$event.id[all_fixes_1$timestamp>min(all_fixes_2$timestamp) & ((all_fixes_1$ptt_day %in% af2pd)==F)]
af2_af3_valid_ids<-all_fixes_2$event.id[all_fixes_2$timestamp>min(all_fixes_3$timestamp) & ((all_fixes_2$ptt_day %in% af3pd)==F)]

all_fixes_1<-all_fixes_1[which((all_fixes_1$timestamp<min(all_fixes_2$timestamp)) | (all_fixes_1$event.id %in% af1_af2_valid_ids)),]
all_fixes_2<-all_fixes_2[which((all_fixes_2$timestamp<min(all_fixes_3$timestamp)) | (all_fixes_2$event.id %in% af2_af3_valid_ids)),]

afn1<-names(all_fixes_1)
afn2<-names(all_fixes_2)

#Column in all_fixes_1 and _2 differ from each other - trim down to shared columns
af_col_names_df<-data.frame(table(c(names(all_fixes_1),names(all_fixes_2))))
af_col_names_df_shared<-as.character(af_col_names_df$Var1[af_col_names_df$Freq==2])

#Combine all_fixes_1 and _2 and _3
all_fixes<-rbind(all_fixes_1[,af_col_names_df_shared],
                 all_fixes_2[,af_col_names_df_shared],
                 all_fixes_3[,af_col_names_df_shared])

max(all_fixes$timestamp)

#Remove fixes for Idemili. This bird was tagged for a short time and then the tag was removed. 
all_fixes<-all_fixes[which(all_fixes$individual.local.identifier!="Idemili"),]

# #Tidy up naming convention for Disco Tony
# all_fixes[which(all_fixes$individual.local.identifier=="Disco_Tony"),]$individual.local.identifier<-"Disco Tony"

d1<-all_fixes
d1$timestamp <- as.POSIXct(strptime(substr(d1$timestamp,1,19), "%Y-%m-%d %H:%M:%S"), "UTC")
d1$year <- year(d1$timestamp)
d1$month <- month(d1$timestamp)
d1$day <- mday(d1$timestamp)
d1$julian <- yday(d1$timestamp)

names(d1)[names(d1)=="tag.local.identifier"]<-"ptt"

#Remove argos flagged erroneous points to give best-of-day
d1_nbod<-d1 #First keep hold of full dataset, will be useful later
d1<-d1[which(is.na(d1$algorithm.marked.outlier)),] #Switched

#Remove duplicates
#NB this must be done AFTER removing outliers
d1$dupID<-paste(d1$ptt, d1$timestamp)
d1$dupLC<-duplicated(d1$dupID)
d1<-d1[which(d1$dupLC==F),]

#Remove fixes with no location data
d1<-d1[which(is.na(d1$location.long)==F),]
d1_nbod<-d1_nbod[which(is.na(d1_nbod$location.long)==F),]

#REMOVE SOME PTTS
#Remove woodcock and Mongolian/Beijing cuckoos or tags which seem to never have been fitted to a bird. **CHECK ALL THESE WITH CHRIS**
#51897 is potentially not a UK cuckoo - just points in S & E Africa
#53388-92 are woodcocks
#161311-7 & 170436-8 are cuckoos in China
#115587 & 179110 just have one point in UK (near Cambridge) - ever fitted to bird?
#162704 has no location info associated with it
non_cuckoo_ptts<-c(51897,53388:53392,115587,161311:161317,170436:170438,179110)
d1<-d1[which((d1$ptt %in% non_cuckoo_ptts)==F),]
#Also remove Idemili. Idemili never left the UK, and her tag 115590 was passed on to Waller in 2013.
d1<-d1[which((d1$ptt==115590 & d1$year==2012)==F),]
#Also remove PTTs 170429, 170432, 170433, 170434, 170435 (2G tags with greater timing uncertainty)
ptts_2g<-c(170429,170432,170433,170434,170435)
d1<-d1[which((d1$ptt %in% ptts_2g)==F),]
#And PTT 128301 (unclear migratory route / just one migratory milestone)
d1<-d1[which((d1$ptt==128301)==F),]
#Omit cuckoos tagged in 2021
d1<-d1[which((d1$ptt %in% 213801:213810)==F),]

#Not all of the fixes in d1 have names, even though they are associated with cuckoos that do have names. 
#Fill these in from the fixes that do have names. This is a bit unnecessary but does help with error checking.
un_d1_ptt<-unique(d1$ptt)
un_d1_tab<-data.frame(matrix(nrow=length(un_d1_ptt),ncol=2))
names(un_d1_tab)<-c("ptt","name")
un_d1_tab$ptt<-un_d1_ptt
for(i in 1:nrow(un_d1_tab)){ #Create conversion table for ptt to name
  temp_names<-unique(d1[which(d1$ptt==un_d1_tab$ptt[i]),]$individual.local.identifier)
  if(any(nchar(temp_names)>0)) {un_d1_tab$name[i]<-temp_names[which(nchar(temp_names)>0)]} else {
    un_d1_tab$name[i]<-""} #Some cuckoos were never named
}

#Need to know tagging date in order to remove fixes from before tagging (from testing of the tag).
#Use tagging info file from Chris H. Fill in un_d1_tab with date of tagging, and use that to fill in d1.
tag_details<-read.csv("H:/JGD H drive backup 160320/Raw data/Cuckoos/Tagged Cuckoo stats - Sheet1.csv")
names(tag_details)[names(tag_details)=="X"]<-"ptt"
#Remove Idemili from tag_details so that there is only one bird per tag, to allow matching.
tag_details<-tag_details[-which(tag_details$Name=="Idemili"),]
un_d1_tab$tagging_date<-NA
for(i in 1:nrow(un_d1_tab)){
  temp_date<-tag_details[which(tag_details$ptt==un_d1_tab$ptt[i]),]$Date.caught
  temp_year<-strsplit(temp_date,"/")[[1]][3]
  if(length(temp_date)==0) {next}else{
    if(nchar(temp_year)==2) un_d1_tab$tagging_date[i]<-substr(as.character(strptime(temp_date,format="%d/%m/%y")),1,10) #Two different date formats used
    if(nchar(temp_year)==4) un_d1_tab$tagging_date[i]<-substr(as.character(strptime(temp_date,format="%d/%m/%Y")),1,10)
  }
}
d1$tagging_date<-character(nrow(d1))
for(i in 1:nrow(un_d1_tab)){ #Fill in main database
  d1[which(d1$ptt==un_d1_tab$ptt[i]),]$individual.local.identifier<-un_d1_tab$name[i]
  d1[which(d1$ptt==un_d1_tab$ptt[i]),]$tagging_date<-un_d1_tab$tagging_date[i]
}
#Remove tag test locations (where fix was before the date the bird was tagged)
d1$date_tagged_posix <- as.POSIXct(strptime(d1$tagging_date,format="%Y-%m-%d"))
d1<-d1[which((d1$timestamp<d1$date_tagged_posix)==F),]

#Remove locations outside Afro-Palearctic (one of the fixes for Chance is outside the Microwave Telemetry office in Columbia, MD)
d1<-d1[which(d1$location.long>-30 & d1$location.long<45),]
plot(d1$location.long,d1$location.lat,col=rgb(1-(d1$julian/(max(d1$julian))),d1$julian/(max(d1$julian)),1,1))

#Remove points from after birds had died
#I put this list of known deaths together. Dead birds were identified by lack of movement for long periods at uncharacteristic times.
#Outcome code: M, definitely died; U, unknown - UA least likely to be tag failure, UC most likely to be tag failure.
death_tab<-data.frame(do.call("rbind",list(c(50017,"2018-10-05","UB"),
                                           c(50023,"2019-03-08","M"), #Updated >1WK EARLIER
                                           c(50024,"2019-04-11","M"), #Checked 
                                           c(50026,"2021-02-26","UB"), #MUCH LATER 
                                           c(50031,"2018-08-19","M"), #Checked 
                                           c(50032,"2019-07-24","UC"),
                                           c(50042,"2019-06-19","M"), #Checked 
                                           c(50046,"2019-07-12","UA"),
                                           c(50051,"2019-08-03","M"), #Checked 
                                           c(50053,"2019-08-14","M"), #Checked (how do we know this one died?)
                                           c(62518,"2012-04-11","UA"),
                                           c(62520,"2012-02-20","M"), #Checked 
                                           c(62602,"2012-04-04","M"), #Checked 
                                           c(62608,"2015-08-08","M"), #Checked >1WK EARLIER
                                           c(62688,"2012-08-05","M"), #Checked 
                                           c(115586,"2015-08-19","UA"),
                                           c(115588,"2012-09-15","UC"),
                                           c(115589,"2012-08-05","M"), #Checked
                                           c(115590,"2014-09-29","M"), #Checked >1WK EARLIER
                                           c(115591,"2014-11-23","UA"),
                                           c(115592,"2012-06-02","M"), #Checked >1WK EARLIER
                                           c(115593,"2012-10-07","UA"),
                                           c(115594,"2017-01-12","UA"),
                                           c(115595,"2012-07-24","M"), #Later than spreadsheet - definitely still moving on 24th July
                                           c(115596,"2012-08-05","M"), #Checked >1WK LATER
                                           c(115597,"2013-04-25","M"), #Checked >1WK EARLIER
                                           c(115598,"2012-09-17","M"), #Checked 
                                           c(115599,"2014-03-14","M"), #Checked >1WK EARLIER
                                           c(115600,"2012-08-02","UA"),
                                           c(115602,"2013-12-05","UB"),
                                           c(121791,"2013-09-10","M"), #Checked >1WK EARLIER
                                           c(121792,"2016-04-20","M"), #Checked >1WK EARLIER
                                           c(128295,"2013-06-11","M"), #Checked 
                                           c(128296,"2014-04-04","M"), #Checked >1WK EARLIER
                                           c(128297,"2015-01-23","UA"),
                                           c(128298,"2013-07-28","UB"),
                                           c(128299,"2013-08-20","M"), #Checked >1WK EARLIER
                                           c(128300,"2015-03-17","UA"),
                                           c(128301,"2013-07-29","UC"),
                                           c(128302,"2014-11-07","UB"),
                                           c(128303,"2014-04-21","M"), #Checked 
                                           c(128304,"2013-08-15","M"), #Checked
                                           c(134950,"2014-07-15","M"), #Checked 1WK EARLIER (given date of death in ss is later than last fix - how?)
                                           c(134951,"2015-05-01","M"), #Checked - no way of telling when this bird died during breeding season
                                           c(134952,"2015-07-10","M"), #Checked >1WK EARLIER
                                           c(134953,"2014-07-14","M"), #Checked >1WK EARLIER (given date of death in ss is later than last fix - how?)
                                           c(134954,"2014-05-27","M"), #Checked 
                                           c(134955,"2015-06-12","UC"),
                                           c(134956,"2015-08-14","UA"),
                                           c(134957,"2015-06-24","M"), #Checked >1WK LATER
                                           c(134958,"2015-04-20","UB"),
                                           c(134959,"2014-10-16","UC"),
                                           c(134960,"2014-08-10","UB"),
                                           c(134961,"2014-07-17","M"), #Checked 
                                           c(134962,"2014-05-29","M"), #Checked - may have died before this potentially but doesn't matter as never left UK
                                           c(134963,"2015-07-01","UA"),
                                           c(134964,"2014-07-12","M"), #Checked
                                           c(146753,"2015-06-29","M"), #Checked >1WK EARLIER
                                           c(146754,"2016-05-04","UA"),
                                           c(146755,"2015-05-25","M"), #Checked - may have died before this potentially but doesn't matter as never left UK
                                           c(146756,"2015-07-26","M"), #Checked 
                                           c(146757,"2016-04-23","M"), #Checked
                                           c(146758,"2018-04-03","UC"),
                                           c(146759,"2019-07-19","M"), #Checked 1WK EARLIER
                                           c(146760,"2016-06-22","M"), #Checked >1WK EARLIER
                                           c(146761,"2015-09-15","M"), #Checked 
                                           c(146762,"2015-05-21","M"), #Checked
                                           c(161319,"2016-08-16","UA"),
                                           c(161320,"2016-10-06","M"), #Checked
                                           c(161321,"2019-04-03","M"), #Checked
                                           c(161322,"2016-08-07","UA"),
                                           c(161323,"2017-04-05","M"), #Checked >1WK EARLIER
                                           c(161324,"2019-04-10","M"), #Checked
                                           c(161325,"2016-12-22","M"), #Checked >1WK EARLIER
                                           c(179106,"2020-05-03","M"), 
                                           c(179107,"2021-01-31","M"),
                                           c(179108,"2020-04-06","UB"), #From Chris - uncertain outcome, not necessarily death
                                           c(179109,"2019-08-09","M") #Checked 1WK EARLIER
)))
names(death_tab)<-c("ptt","last_known_alive_fix","outcome_code")
# write.csv(death_tab,"H:/JGD H drive backup 160320/Cuckoos/Output data/mortality_tab.csv",row.names=F)

#For each bird that we know to have died, remove fixes that came after the last fix when it was known to be alive
temp_dead_fix_ids_list<-list()
for(i in 1:nrow(death_tab)){
  # if(death_tab[i,]$outcome_code != "M") next
  temp_d1<-d1[which(d1$ptt==death_tab[i,]$ptt),]
  temp_alive_ts<-temp_d1$timestamp[substr(temp_d1$timestamp,1,10)<=death_tab[i,]$last_known_alive_fix]
  if(length(temp_alive_ts)==0) next
  temp_max_alive_ts<-max(temp_alive_ts)
  temp_dead_fix_ids_list[[i]]<-temp_d1$event.id[which(temp_d1$timestamp>temp_max_alive_ts)]
} #Returns warning if there are no fixes that come after the last fix when it was known to be alive
temp_dead_fix_ids<-unlist(temp_dead_fix_ids_list)
d1<-d1[which((d1$event.id %in% temp_dead_fix_ids)==F),]


# #Check for outliers visually and remove outliers manually
# #Only do this once - 
# d6<-d1[d1$year>=2020,]
# d6$ptt_yearID<-paste(d6$ptt, d6$year, sep="_")
# d6$outlier_flagged<-""
# world_shp_simp<-gSimplify(world_shp,0.2)
# out<-NULL
# for(i in unique(d6$ptt_yearID)){ #Use this code to look for outliers
#   print(i)
#   print(table(d6[d6$ptt==strsplit(i, "_")[[1]][1],]$year))
#   temp<- d6[d6$ptt_yearID==i,]
#   plot(location.lat~location.long, temp, col="red",pch=20)
#   lines(location.lat~location.long, temp, col=3)
#   plot(world_shp_simp,add=T,border="gray50")
#   a1<-"N"
#   a1<-readline("are there outliers? (y?)")
#   if(a1=="y")
#   {id<-identify(x=temp$location.long, y=temp$location.lat)
#   temp[id,]$outlier_flagged<-"Y"
#   points(location.lat~location.long,
#          temp[temp$outlier_flagged=="Y",], col=4)}
#   out<-rbind(out, temp)
# }
# out[out$outlier_flagged=="Y",]

# 
# #Remove outliers identified above
# #Event ID 192215226 (Waller) is erroneous - a fix over the North Sea on the day of tagging in the Highlands
# #Event ID 3532431591 (Boris) seems erroneous - a fix in Namibia two days before fixes in Senegal
# d1<-d1[which((d1$event.id %in% c(192215226,3532431591))==F),]
# 
# #Hard to see UK outliers at scale above - repeat but zoomed in
# d7<-out
# out2<-NULL
# for(i in unique(d7$ptt_yearID)[41:167]){ #Use this code to look for outliers
#   print(i)
#   print(table(d7[d7$ptt==strsplit(i, "_")[[1]][1],]$year))
#   temp<- d7[d7$ptt_yearID==i,]
#   plot(location.lat~location.long, temp, col="red",pch=20,ylim=c(45,62),xlim=c(-10,5))
#   lines(location.lat~location.long, temp, col=3)
#   plot(world_shp_simp,add=T,border="gray50")
#   a1<-"N"
#   a1<-readline("are there outliers? (y?)")
#   if(a1=="y")
#   {id<-identify(x=temp$location.long, y=temp$location.lat)
#   temp[id,]$outlier_flagged<-"Y"
#   points(location.lat~location.long,
#          temp[temp$outlier_flagged=="Y",], col=4)}
#   out2<-rbind(out2, temp)
# }
# 
# #Save file of outliers to be read back in, so only have to do this once
# outliers<-out2$event.id[which(out2$outlier_flagged=="Y")]
# # save(outliers,file="H:/JGD H drive backup 160320/Raw data/Cuckoos/outliers.RData")
# #262 outliers identified

load("H:/JGD H drive backup 160320/Raw data/Cuckoos/outliers.RData")
additional_outliers<-c(10257499576,10273573510,10274634045,14158249882,15194553858,325476489,14714124246,10666314632,3819981836,194874779,728272626,1512652704,1546375775,14769737138,3819966333) #Various extra position errors subsequently noticed
d1<-d1[which((d1$event.id %in% c(outliers,additional_outliers))==F),]
d1_nbod<-d1_nbod[which((d1_nbod$event.id %in% c(outliers,additional_outliers))==F),] #Same for full dataset

#Remove birds that didn't make it out of the UK
#Mark does this using 'migratory.strategy' from strategy.dat; I automate it here - which birds don't have any points outside GB?
#Which points fall outside of GB? - using a polygon that allows some extra space for position error
gb_poly<-readOGR("H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation/GB+ polygon/GB+ polygon.shp")
d1_spdf<-d1
coordinates(d1_spdf)<-~location.long+location.lat
proj4string(d1_spdf)<-CRS("+init=epsg:4326")
d1$in_gb<-numeric(nrow(d1))
for(i in 1:nrow(d1)){ #Takes a couple of minutes
  d1$in_gb[i]<-over(d1_spdf[i,],gb_poly)$id
}
ever_left_gb_ptts<-unique(d1[which(is.na(d1$in_gb)),]$ptt)
d1<-d1[which(d1$ptt %in% ever_left_gb_ptts),] #Limit dataset to birds that ever left GB
un_d1_ptt<-unique(d1$ptt) #Updated

#Make spatial version of d1_nbod too
d1_nbod_spdf<-d1_nbod
coordinates(d1_nbod_spdf)<-~location.long+location.lat
proj4string(d1_nbod_spdf)<-CRS("+init=epsg:4326")

#Remove unnecessary variables
d1$visible<-NULL
d1$argos.best.level<-NULL
d1$argos.calcul.freq<-NULL
d1$argos.gdop<-NULL
d1$argos.iq<-NULL
d1$argos.lat1<-NULL
d1$argos.lat2<-NULL
d1$argos.lon1<-NULL
d1$argos.lon2<-NULL
d1$argos.nb.mes<-NULL
d1$argos.nb.mes.120<-NULL
d1$argos.nopc<-NULL
d1$argos.pass.duration<-NULL
d1$argos.sat.id<-NULL
d1$argos.sensor.1<-NULL
d1$argos.sensor.2<-NULL
d1$argos.sensor.3<-NULL
d1$argos.sensor.4<-NULL
d1$argos.transmission.timestamp<-NULL
d1$argos.valid.location.algorithm<-NULL
d1$sensor.type<-NULL
d1$individual.taxon.canonical.name<-NULL
d1$study.name<-NULL
d1$dupID<-NULL
d1$dupLC<-NULL
d1$test.location<-NULL

## Define breeding area buffers - 50km buffer around the mean location during breeding season
#Need to define breeding locations and put 50km buffer around them
#Limit to birds that ever returned to the UK
# *NB - mention to Chris* Due to using the non-best-of-day data, there is a lot of error in the breeding area location points. This shouldn't affect the mean location though, which is all we need these for.
# d1_spdf$breeding_binary<-0

#Write out shapefiles for each bird. 
# for(i in 1:length(un_d1_ptt)){
#   temp_spdf<-d1_spdf[d1_spdf$ptt==un_d1_ptt[i],]
#   writeOGR(temp_spdf,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation/Shapefiles by ptt",layer=un_d1_ptt[i],driver="ESRI Shapefile")
# }
#Then manually (in QGIS) limit the points down to those that are in a clear cluster denoting the breeding area. 
#This was straightforward to identify for all birds, but there was clearly wide position error around the mean.
#  - this is not a problem though as we're just taking the mean.
#Read back in and fill in table of mean coordinates, giving the centre of the cluster
br_cent_coords<-data.frame(matrix(ncol=3,nrow=length(un_d1_ptt)))
names(br_cent_coords)<-c("ptt","long","lat")
br_cent_coords$ptt<-un_d1_ptt
for(i in 1:length(un_d1_ptt)){
  temp_spdf<-readOGR(paste0("H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation/Shapefiles by ptt/",un_d1_ptt[i],".shp"))
  br_cent_coords[i,2:3]<-as.numeric(apply(coordinates(temp_spdf),2,mean))
}

br_cent_coords_spdf<-br_cent_coords
coordinates(br_cent_coords_spdf)<-~long+lat
proj4string(br_cent_coords_spdf)<-CRS("+proj=longlat +datum=WGS84")
br_cent_coords_spdf_30utm<-spTransform(br_cent_coords_spdf,CRS("+init=epsg:32630"))
plot(br_cent_coords_spdf_30utm)
br_buffs<-gBuffer(br_cent_coords_spdf_30utm,width=50000,byid=T) #50km buffer around centre of breeding points cluster
br_buffs_ll<-spTransform(br_buffs,CRS("+init=epsg:4326"))
plot(br_buffs)
points(br_cent_coords_spdf_30utm,pch=20)
plot(br_buffs[br_buffs$ptt==50026,],col="green",add=T)
points(br_cent_coords_spdf_30utm[br_cent_coords_spdf_30utm$ptt==50026,],pch=20,col="red")
#Next, use this to work out when birds made it back to their breeding area (not just the UK)

# writeOGR(br_buffs_ll,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation",layer="br_buffs_ll",driver="ESRI Shapefile")

#Cuckoo satellite tag dataset now cleaned


#########################################################################################################

### 2. CALCULATE STOPOVERS ###

##2A. Calculate migratory status, on the basis of time and distance between points

d1<-d1[order(d1$ptt, d1$timestamp),] #Order by ptt and timestamp: timing calculations require fixes being in time sequence

#Give each fix a transmission cycle ID. Fixes get unique TC IDs if they are more than 10hrs apart
#Add transmission cycle info
tcyclefunc <- function(x) {
  ##function diff calculates the difference between consecutive values in a vector
  ##so- is the difference between the consecutive times greater than 10/24
  tcycle <- as.numeric(diff(x$timestamp), units="days") > 10/24
  ##if yes gets 1 else gets zero
  ##and add a one to start with
  tcycle <- c(1,ifelse(tcycle==TRUE,1,0))
  ##summuative sum as you go along the vector 
  tcycle <- cumsum(tcycle)
  ##new dataframe with data, tcycle and column showing the number of days between successive timestamps
  ##add an NA for the first record for each bird as can't have travelled any distance for that one
  newdataframe <- data.frame(x,tcycle,days.btwn.trans=c(NA,as.numeric(diff(x$timestamp), units="days")))
  ##return dataframe
  return(newdataframe)
}
d2 <- do.call(rbind,by(d1, list(d1$ptt), tcyclefunc))
rownames(d2) <- d2$id <- 1:nrow(d2) #Change rownames, otherwise hard to view whole dataset in console. Also add id column.

names(d2)[names(d2)=="location.long"]<-"long"
names(d2)[names(d2)=="location.lat"]<-"lat"

#Add movement data
addmovedata <- function(input) { # distance, bearing, etc function
  # for each individual, do the following:
  ### --- DISTANCE & BEARING CALCULATION --- ###
  # use long & lat to calculate distances based on original measured location point
  dist.cuckoo <- input[,c("long","lat")]
  coordinates(dist.cuckoo) <- c("long", "lat")
  proj4string(dist.cuckoo) <- CRS("+proj=longlat +datum=WGS84")
  d1 <- dist.cuckoo@coords
  d2 <- dist.cuckoo@coords[-1,] # create new matrix minus the first row so that d2 starts at the second observation
  d2 <- rbind(d2, c(0,0)) # create a placeholder last row in the second distance matrix
  dist <- distCosine(d1,d2)/1000 # distance between points, in km
  bear <- bearing(d1,d2) # bearing between points
  dist <- c(NA,dist[-length(dist)])
  bear <- c(NA,bear[-length(bear)])
  distbear <- data.frame(input, distance=dist, bearing=bear)
  ### --- MOVEMENT GROUPS --- ###
  ##classifying as migratory if at least 30 km (now 50 km) from previous location otherwise is stationary
  mgroup <- rep(NA,nrow(distbear))
  for (n in 1:nrow(distbear)){
    if (is.na(distbear$distance[n])) {mgroup[n] <- 1 #Start with 1
    ##if under 50km gets same group as previous point # 
    ## Previously used 30, changed to 50 as per nature paper
    } else if (distbear$distance[n] <= 50) {mgroup[n] <- mgroup[n-1]
    ##if greater than 50km gets new point
    } else {mgroup[n] <- mgroup[n-1] + 1}
  }
  ### --- MOVEMENT TYPES --- ###
  mtype <- c("C", rep(NA, nrow(distbear)-1))
  for (n in 2:(nrow(distbear)-1)){
    if (mgroup[n] != mgroup[n-1] & mgroup[n] != mgroup[n+1]) {
      mtype[n] <- "M"
    } else {
      mtype[n] <- "S"
    }
  }
  completedata <- data.frame(distbear,mgroup,mtype)
  return(completedata)
}
withmovedata <- by(d2, list(d2$ptt), addmovedata)
withmovedata.all <- do.call(rbind, withmovedata)

#Some visualisation of stationary and migratory points
plot(world_shp,border="gray50",xlim=range(withmovedata.all$long),ylim=range(withmovedata.all$lat))
points(withmovedata.all$long[withmovedata.all$mtype=="S"],withmovedata.all$lat[withmovedata.all$mtype=="S"],col="red",pch=20)
points(withmovedata.all$long[withmovedata.all$mtype=="M"],withmovedata.all$lat[withmovedata.all$mtype=="M"],col="blue",pch=20)

#Write out for data exploration
# withmovedata.all_S<-withmovedata.all[which(withmovedata.all$mtype=="S"),]
# withmovedata.all_M<-withmovedata.all[which(withmovedata.all$mtype=="M"),]
# coordinates(withmovedata.all_S)<-coordinates(withmovedata.all_M)<-~long+lat
# proj4string(withmovedata.all_S) <- proj4string(withmovedata.all_M) <- CRS("+proj=longlat +datum=WGS84")
# writeOGR(withmovedata.all_S,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation",layer="withmovedata.all_S",driver="ESRI Shapefile")
# writeOGR(withmovedata.all_M,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation",layer="withmovedata.all_M",driver="ESRI Shapefile")

#More visualisation of stationary and migratory points, by month
months<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
plot(world_shp,border="gray50",xlim=range(withmovedata.all$long),ylim=range(withmovedata.all$lat),main=months[i])
points(withmovedata.all[withmovedata.all$month==i,]$long,withmovedata.all[withmovedata.all$month==i,]$lat,pch=20,col="red")


##2B. Calculate stopover durations, on the basis of migratory status
#This looks like a roundabout way of summing days.btwn.trans within mgroup!
calculate.LOS <- function(data.list.byname) {
  fulldata <- data.list.byname[order(data.list.byname$timestamp),]
  timediffs <- tapply(fulldata$timestamp, fulldata$mgroup, function(x) {
    x <- rev(x)
    difftime(x[1:(length(x)-1)] , x[2:length(x)], units="days")
  })  
  LOS <- do.call(rbind, lapply(timediffs, function(x) ifelse(length(x)==0, NA, sum(x))))
  LOStable <- data.frame(mgroup=levels(as.factor(fulldata$mgroup)), LOS)
  newdat.LOS <- merge(fulldata, LOStable, by.x="mgroup", by.y="mgroup")
  return(newdat.LOS)
}
withLOSdata <- lapply(withmovedata, calculate.LOS)
withLOSdata.all <- do.call(rbind, withLOSdata)

# need to extend LOS if next stopover > 1 duty cycle of time away.
# If so, we assume the bird stayed at the current stopover until
# the final duty cycle (2 days) prior to the next stopover.
# Note: we need to consider 'M' classed points as limiters to stopover duration
withLOSdata.all$LOS.recalc<-withLOSdata.all$LOS
for(j in 1:nrow(withLOSdata.all)){
  if(withLOSdata.all[j,]$mgroup==withLOSdata.all[j+1,]$mgroup | is.na(withLOSdata.all[j,]$days.btwn.trans+withLOSdata.all[j+1,]$days.btwn.trans) | withLOSdata.all[j+1,]$days.btwn.trans<3) {next} else {
    withLOSdata.all[which(withLOSdata.all$ptt==withLOSdata.all[j,]$ptt & withLOSdata.all$mgroup==withLOSdata.all[j,]$mgroup),]$LOS.recalc<-withLOSdata.all[which(withLOSdata.all$ptt==withLOSdata.all[j,]$ptt & withLOSdata.all$mgroup==withLOSdata.all[j,]$mgroup),]$LOS+(withLOSdata.all[j+1,]$days.btwn.trans-2);
  }
}

wld_spdf<-withLOSdata.all
coordinates(wld_spdf)<-~long+lat
proj4string(wld_spdf)<-CRS("+init=epsg:4326")

#The following is lines 172-309 from cuckoo_stopovers.R (Mark's GitHub version)
dat<-withLOSdata.all
birds<-unique(dat$ptt)
stopovers_tab<-NULL
for(i in birds){
  # for each bird give the mgroups which are not migratory ones
  bird_out<-NULL
  mgz<-na.omit(unique(dat[dat$ptt==i & dat$mtype=="S",]$mgroup))
  # removes 'M' migratory and 'C' deployment point
  for (j in mgz){
    ptt_mgroup<-dat[dat$ptt==i & dat$mgroup==j,]
    ###########################################################
    ########## Allows a stopover to be defined based on time #
    # Here we use 1 hour - we're doing this on non duty-cycle data # THIS ACTUALLY ELIMINATES STOPOVERS OF LESS THAN 1 DAY - CORRECT?
    # as it gives more accurate timings 
    # if(unique(ptt_mgroup$LOS)<1){next} #Using LOS rather than LOS.recalc, for now
    if(unique(ptt_mgroup$LOS.recalc)<1){next} 
    ###########################################################
    # new stopover end datetime from LOS.recalc
    endTT_LOS<-as.double(min(ptt_mgroup$timestamp))+(unique(ptt_mgroup$LOS)*(24*60*60))
    endtimestamp_LOS<-paste(as.Date(as.POSIXlt(endTT_LOS, origin="1970-01-01", "UTC")), 
                            format((as.POSIXlt(endTT_LOS, origin="1970-01-01", "UTC")), "%H:%M:%S"))
    endTT<-as.double(min(ptt_mgroup$timestamp))+(unique(ptt_mgroup$LOS.recalc)*(24*60*60))
    endtimestamp<-paste(as.Date(as.POSIXlt(endTT, origin="1970-01-01", "UTC")), 
                        format((as.POSIXlt(endTT, origin="1970-01-01", "UTC")), "%H:%M:%S"))
    ###########################################
    ptt_mgroup_out<-data.frame(ptt=i, mgroup=j, #Got rid of code asking for ptt_mgroup$name or $comments - no such column in my data
                               SO_start=min(ptt_mgroup$timestamp), #NB There is no need for SO_start_LOS, because the definition of the start of the stopover doesn't change using LOS.recalc
                               SO_end=endtimestamp,
                               SO_end_LOS=endtimestamp_LOS,
                               SO_days=unique(ptt_mgroup$LOS.recalc),
                               SO_days_LOS=unique(ptt_mgroup$LOS),
                               SO_median_long=median(ptt_mgroup$long),
                               SO_median_lat=median(ptt_mgroup$lat),
                               SO_month=round(median(ptt_mgroup$month)),
                               SO_year=round(median(ptt_mgroup$year)))
    bird_out<-rbind(bird_out, ptt_mgroup_out) 
  }
  stopovers_tab<-rbind(stopovers_tab, bird_out)
}
stopovers_tab$SO_start<-as.POSIXct(stopovers_tab$SO_start,tz="UTC") #Need timezone in there - in Mark's code it was treating summer UTC as BST!
stopovers_tab$SO_end<-as.POSIXct(stopovers_tab$SO_end,tz="UTC")
stopovers_tab$SO_end_LOS<-as.POSIXct(stopovers_tab$SO_end_LOS,tz="UTC")

#Identify whether there was a missing duty cycle preceding or following each stopover
#Specifically, were there  - a) any fixes in the 48hrs preceding the start of the stopover
#                          - b) any fixes in the 48hrs following the end of the stopover
stopovers_tab$miss_foll_dcs<-stopovers_tab$miss_prec_dcs<-logical(nrow(stopovers_tab))
for(i in 1:nrow(stopovers_tab)){
  temp_prec_48hr_start<-stopovers_tab[i,]$SO_start-(48*60*60)
  temp_foll_48hr_end<-stopovers_tab[i,]$SO_end_LOS+(48*60*60)
  stopovers_tab$miss_prec_dcs[i]<-any(d1$ptt==stopovers_tab$ptt[i] & d1$timestamp>temp_prec_48hr_start & d1$timestamp<stopovers_tab[i,]$SO_start)==F 
  stopovers_tab$miss_foll_dcs[i]<-any(d1$ptt==stopovers_tab$ptt[i] & d1$timestamp>stopovers_tab[i,]$SO_end_LOS & d1$timestamp<temp_foll_48hr_end)==F
}
stopovers_tab$miss_prec_dcs[which(stopovers_tab$mgroup==1)]<-NA #If mgroup is 1, it is by definition the first stopover: make miss_prec_dcs NA
#stopovers_tab prepared.

#Frequencies of missing duty cycle preceding or following stopover
table(stopovers_tab$miss_prec_dcs) #Preceding
table(stopovers_tab$miss_foll_dcs) #Following
#I.e. most stopovers are followed or preceded by missing duty cycles

#Visualise how stopovers are defined
i<-4
temp_ptt<-un_d1_ptt[i]
temp_d1<-withmovedata.all[withmovedata.all$ptt==temp_ptt,]
plot(temp_d1$timestamp[temp_d1$mtype=="S"],temp_d1$lat[temp_d1$mtype=="S"],xlab="Date",ylab="Latitude",main=temp_ptt)
points(temp_d1$timestamp[temp_d1$mtype=="M"],temp_d1$lat[temp_d1$mtype=="M"],col="blue")
temp_st<-stopovers_tab[stopovers_tab$ptt==temp_ptt,]
for(i in 1:nrow(temp_st)){
  arrows(temp_st[i,]$SO_start,temp_st$SO_median_lat[i],temp_st[i,]$SO_end,temp_st$SO_median_lat[i],col="orange",lwd=2,length=0)
  points(temp_st[i,]$SO_start,temp_st$SO_median_lat[i],col="dark green",pch=4,cex=1.2)
  points(temp_st[i,]$SO_end,temp_st$SO_median_lat[i],col="red",pch=4,cex=1.2)
  # arrows(temp_st[i,]$SO_start,temp_st$SO_median_lat[i],temp_st[i,]$SO_end_LOS,temp_st$SO_median_lat[i],col="orange",lwd=2,length=0)
  # points(temp_st[i,]$SO_start,temp_st$SO_median_lat[i],col="dark green",pch=4,cex=1.2)
  # points(temp_st[i,]$SO_end_LOS,temp_st$SO_median_lat[i],col="red",pch=4,cex=1.2)
}

#Not all fixes where a bird is clearly moving are correctly assigned movement type "M".
#Briefly tried assigning movement type based on bearing, turning angle, and velocity.
#There were too many false positives with this approach. Probably a HMM approach is needed.
# temp_d1$bearing_360<-temp_d1$bearing
# temp_d1$bearing_360[which(temp_d1$bearing_360<0)]<-temp_d1$bearing_360[which(temp_d1$bearing_360<0)]+360
# temp_d1$turning_angle<-numeric(nrow(temp_d1))
# temp_d1$turning_angle[2:nrow(temp_d1)]<-temp_d1$bearing_360[2:nrow(temp_d1)]-temp_d1$bearing_360[1:(nrow(temp_d1)-1)]
# temp_d1$time_elapsed<-numeric(nrow(temp_d1))
# temp_d1$time_elapsed[2:nrow(temp_d1)]<-(temp_d1$timestamp[2:nrow(temp_d1)]-temp_d1$timestamp[1:(nrow(temp_d1)-1)])/3600
# temp_d1$velocity<-temp_d1$distance/temp_d1$time_elapsed
# 
# temp_d1$probably_moving<-temp_d1$turning_angle>-60 & temp_d1$turning_angle<60 & temp_d1$velocity>15 & temp_d1$velocity<150
# temp_d1$probably_moving_mult<-logical(nrow(temp_d1))
# for(i in 1:(nrow(temp_d1)-1)){
#   temp_d1$probably_moving_mult[i]<-temp_d1$probably_moving[i] & temp_d1$probably_moving[i+1]
# }
# temp_d1[,c("timestamp","distance","velocity","bearing_360","turning_angle","mtype","probably_moving","probably_moving_mult")][400:500,]

# Add country and biome data to stopovers
countries<-readOGR("H:/JGD H drive backup 160320/Raw data/Cuckoos/cuckoo_tracking/sourced_data/country_borders/TM_WORLD_BORDERS-0.3.shp")
biomes<-readOGR("H:/JGD H drive backup 160320/Raw data/Cuckoos/cuckoo_tracking/sourced_data/biomes_TNC/tnc_terr_ecoregions.shp")

#Use TNC terrestrial ecoregions shapefiles to define a) the Sahara and b) the wintering grounds
ecoregs<-readOGR("H:/JGD H drive backup 160320/Raw data/Cuckoos/cuckoo_tracking/sourced_data/biomes_TNC/tnc_terr_ecoregions.shp")
#The Sahara
#10674 Atlantic coastal desert
#10609 Saharan halophytics
#10691 North Saharan steppe and woodlands
#10698 Sahara desert
#10703 West Saharan montane xeric woodlands
#10702 Tibesti-Jebel Uweinat montane xeric woodlands
#10700 South Saharan steppe and woodlands
#NB omitting Nile delta, which extends some way inland - no cuckoo points fall in it though
sah_polys<-ecoregs[ecoregs$ECO_ID_U %in% c(10674,10609,10691,10698,10703,10702,10700),]
sah<-aggregate(sah_polys) #Ignore geometry is invalid warning
sah<-SpatialPolygons(list(Polygons(list(sah@polygons[[1]]@Polygons[[8]]),ID=1))) #For some reason it retains holes - get outer polygon only



#The wintering grounds - Central African rainforest, and countries to the south
#10091 #Cross-Sanaga-Bioko coastal forests #REMOVED
#10086 #Atlantic Equatorial coastal forests
#10108 #Northeastern Congolian lowland forests
#10110 #Northwestern Congolian lowland forests
#10113 #Western Congolian swamp forests
#10088 #Central Congolian lowland forests
#10140 #Western Congolian forest-savanna mosaic
#10135 #Southern Congolian forest-savanna mosaic
#10156 #Angolan Scarp savanna and woodlands
#10118 #Angolan Miombo woodlands
ca_rf_polys<-ecoregs[ecoregs$ECO_ID_U %in% c(10086,10108,10110,10113,10088,10140,10135,10156,10118),]
#Add in Central African mangroves south of 3N
cam<-ecoregs[ecoregs$ECO_ID_U==10193,]
cam_polys<-disaggregate(cam)
cam_centroids<-gCentroid(cam_polys,byid=T)
cam_polys_sub3<-cam_polys[as.numeric(which(coordinates(cam_centroids)[,2]<3)),]
ca_rf_polys<-rbind(ca_rf_polys,cam_polys_sub3)
ca_rf<-aggregate(ca_rf_polys)
ca_rf<-SpatialPolygons(list(Polygons(list(ca_rf@polygons[[1]]@Polygons[[1]]),ID=1))) #For some reason it retains holes - get outer polygon only

sah_spdf<-SpatialPolygonsDataFrame(sah,data.frame(reg="Sahara"))
ca_rf_spdf<-SpatialPolygonsDataFrame(ca_rf,data.frame(reg="CAfr rainforest"))
regs_spdf<-rbind(sah_spdf,ca_rf_spdf)
crs(sah_spdf)<-crs("+init=epsg:4326")
crs(ca_rf_spdf)<-crs("+init=epsg:4326")
crs(regs_spdf)<-crs("+init=epsg:4326")

#New code 091221 - buffer north edge of wintering grounds inwards by 55km. NB this eventually overwrites ca_rf_spdf
ca_rf_spdf_t<-spTransform(ca_rf_spdf, CRS("+init=epsg:32633"))
ca_rf_spdf_t_b<-gBuffer(ca_rf_spdf_t,width=-55000)

# #Write out projected wintering grounds and internal buffer as points shapefiles, for identifying vertices at which polygons should be joined
ca_rf_pts<-data.frame(ca_rf_spdf_t@polygons[[1]]@Polygons[[1]]@coords)
ca_rf_b_pts<-data.frame(ca_rf_spdf_t_b@polygons[[1]]@Polygons[[1]]@coords)
names(ca_rf_pts)<-names(ca_rf_b_pts)<-c("x","y")
coordinates(ca_rf_pts)<-coordinates(ca_rf_b_pts)<-c("x","y")
crs(ca_rf_pts)<-crs(ca_rf_b_pts)<-crs(ca_rf_spdf_t)
ca_rf_pts$order<-1:length(ca_rf_pts)
ca_rf_b_pts$order<-1:length(ca_rf_b_pts)

# writeOGR(ca_rf_pts,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation/Wintering grounds",layer="ca_rf_pts2",driver="ESRI Shapefile",overwrite_layer=T)
# writeOGR(ca_rf_b_pts,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation/Wintering grounds",layer="ca_rf_b_pts2",driver="ESRI Shapefile",overwrite_layer=T)

#Define new wintering grounds polygon based on identified vertices
# ca_rf_new_pts<-rbind(ca_rf_pts[1:1442,],ca_rf_b_pts[86:462,],ca_rf_pts[2206:4269,]) #Previously, when Central African Mangroves were not included
ca_rf_new_pts<-rbind(ca_rf_pts[1:2533,],ca_rf_pts[2566:2862,],ca_rf_b_pts[161:537,],ca_rf_pts[3626:4926,]) #Skip 2534:2565 because they delineate an estuary in which a fix falls, causing exit of the forest zone to be estimated too early
ca_rf_new_poly<-SpatialPolygons(list(Polygons(list(Polygon(ca_rf_new_pts)),1)))
crs(ca_rf_new_poly)<-crs(ca_rf_spdf_t)
ca_rf_new_poly_unt<-spTransform(ca_rf_new_poly,crs(sah_spdf))
ca_rf_spdf<-SpatialPolygonsDataFrame(ca_rf_new_poly_unt,data.frame(reg="CAfr rainforest"))
plot(ca_rf_spdf)

dat_spdf<-stopovers_tab
coordinates(dat_spdf)<-~SO_median_long+SO_median_lat
proj4string(dat_spdf)<-CRS("+init=epsg:4326")
dat_spdf$region<-character(nrow(dat_spdf))
# dat_spdf$biome2<-dat_spdf$biome1<-dat_spdf$country<-character(nrow(dat_spdf))

for(i in 1:nrow(dat_spdf)){
  # dat_spdf$country[i]<-over(dat_spdf[i,],countries)$NAME
  # dat_spdf[i,c("biome1","biome2")]<-over(dat_spdf[i,],biomes)[,c("WWF_MHTNUM","ECO_NAME")]
  dat_spdf$region[i]<-over(dat_spdf[i,],regs_spdf)$reg
}

# #Some points lie just off the coast, so need to use the *nearest* country to the point, rather than the point the country lies in
# #Do this for biome/regions NAs too if important later on
# country_nas<-which(is.na(dat_spdf$country))
# closest_countries<-character(length(country_nas))
# for(i in 1:length(country_nas)){
#   temp_dists_to_countries<-numeric(nrow(world_shp))
#   for(j in 1:nrow(world_shp)){
#     temp_dists_to_countries[j]<-gDistance(dat_spdf[country_nas[i],],world_shp[j,])
#   }
#   closest_countries[i]<-world_shp$NAME[which.min(temp_dists_to_countries)]
# } #Ignore warnings - it will complain about coordinates not being planar, but distance calcs will work for this
# dat_spdf$country[country_nas]<-closest_countries

# stopovers_tab$country<-dat_spdf$country
# stopovers_tab$biome1<-dat_spdf$biome1
# stopovers_tab$biome2<-dat_spdf$biome2
stopovers_tab$region<-dat_spdf$region

#Plot showing Sahara and C African rainforest (s. lato) biomes, and stopover locations within each
png("H:/JGD H drive backup 160320/Cuckoos/Plots/sah_ca_rf_map.png",7,7,res=800,units="in")
plot(stopovers_tab$SO_median_long,stopovers_tab$SO_median_lat,pch=20,col="white",xlab="Longitude",ylab="Latitude")
plot(world_shp,add=T,border="gray80")
plot(sah,col="orange",border=NA,add=T)
plot(ca_rf,col="dark green",border=NA,add=T)
points(stopovers_tab$SO_median_long,stopovers_tab$SO_median_lat,pch=20,col="royalblue2")
dev.off()

stopovers_tab_spdf<-stopovers_tab
coordinates(stopovers_tab_spdf)<-~SO_median_long+SO_median_lat
proj4string(stopovers_tab_spdf)<-CRS("+init=epsg:4326")
stopovers_tab$country<-stopovers_tab_spdf$country<-over(stopovers_tab_spdf,world_shp)$NAME #So much quicker than looping 'over' over all d1_spdf!

# #Some points lie just off the coast, so need to use the *nearest* country to the point, rather than the point the country lies in
country_nas<-which(is.na(stopovers_tab_spdf$country))
closest_countries<-character(length(country_nas))
for(i in 1:length(country_nas)){
  temp_dists_to_countries<-numeric(nrow(world_shp))
  for(j in 1:nrow(world_shp)){
    temp_dists_to_countries[j]<-gDistance(stopovers_tab_spdf[country_nas[i],],world_shp[j,])
  }
  closest_countries[i]<-world_shp$NAME[which.min(temp_dists_to_countries)]
} #Ignore warnings - it will complain about coordinates not being planar, but distance calcs will work for this
stopovers_tab$country[country_nas]<-stopovers_tab_spdf$country[country_nas]<-closest_countries

stopovers_tab$id<-stopovers_tab_spdf$id<-1:nrow(stopovers_tab) #Give unique ID to each stopover, to allow picking out the last stopover

#Stopover data ready for contributing to migratory milestones data 

# write.csv(stopovers_tab,"H:/JGD H drive backup 160320/Cuckoos/Output data/stopovers_tab.csv",row.names=F)


#########################################################################################################
### 3. DEFINE MIGRATORY MILESTONES ###
#Two tables:
# - A: all migratory milestones for all birds
# - B: all migratory milestones for complete migratory cycles
# - C: all migratory milestones for complete migratory cycles for which we have high confidence (i.e. preceded/followed by minimal missing duty cycles)

#Migratory milestones and definitions
# *Need a polygon of Africa south of the Sahara

#Define for all migratory cycles, complete *or incomplete* (a migratory cycle being a GB-Africa-GB movement)
#Want a table with one row for each *started* migratory cycle

d1$country<-character(nrow(d1))
d1_spdf<-d1 #Needs redefining because d1 has changed in length
coordinates(d1_spdf)<-~location.long+location.lat
proj4string(d1_spdf)<-CRS("+init=epsg:4326")
d1_spdf$country<-d1$country<-over(d1_spdf,world_shp)$NAME #So much quicker than looping 'over' over all d1_spdf!

#For each unique ptt-year combination, were there any fixes in countries other than the United Kingdom after the fixes in the United Kingdom?
d1$ptt_year<-paste0(d1$ptt,"_",d1$year)
un_ptt_year<-unique(d1$ptt_year)
ever_left_uk_ptt_year_log<-logical(length(un_ptt_year))

for(i in 1:length(un_ptt_year)){
  temp_d1<-d1[d1$ptt_year==un_ptt_year[i],]
  if(nrow(temp_d1[which(temp_d1$country=="United Kingdom"),])==0) next #Skip (default FALSE) if never in UK
  temp_max_uk_date<-max(temp_d1[temp_d1$country=="United Kingdom",]$julian,na.rm=T)
  ever_left_uk_ptt_year_log[i]<-length(na.omit(temp_d1[temp_d1$julian>temp_max_uk_date,]$country))>0
}
ever_left_uk_ptt_year<-un_ptt_year[ever_left_uk_ptt_year_log]

#MIGRATORY MILESTONES
plot(location.lat~timestamp,data=d1)
#NB these migratory milestones use either first fix/stopover or last fix/stopover criteria. 
# - If using the 'last' criteria, I need an indicator for whether there were any subsequent fixes/stopovers
# - If not (because the tag failed or the bird died), I'm not interested in using them as criteria.
#One row for each ptt-year combination in ever_left_uk_ptt_year

mig_timing_df<-data.frame(do.call(rbind,strsplit(ever_left_uk_ptt_year,"_")))
names(mig_timing_df)<-c("ptt","year")
mig_timing_df$year<-as.numeric(mig_timing_df$year)
mig_timing_df[,c("dep_brgr","dep_UK","dep_euna_sb","fin_sah_sb","arr_wgr_a","arr_wgr_b","dep_wgr_b","dep_wgr_a","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr","br_hab","mig_dir","dep_brgr_TS","min_lat_stopover","any_so_s_europe_sbound","lat_last_euro_so_sbound","last_wa_fix_long","last_wa_fix_lat","last_wa_sos_long","last_wa_sos_lat","br_lat","br_long","first_fix_TS","dep_UK_TS","fin_sah_sb_TS","arr_wgr_a_TS","arr_wgr_b_TS","dep_wgr_b_TS","min_lat_stopover_id","dep_wgr_a_TS","dep_wgr_a_TS_ford1plot","dep_wa_TS","arr_UK_TS","arr_brgr_TS","last_fix_TS")]<-NA
mig_timing_df_uncertainty<-mig_timing_df

for(i in 1:nrow(mig_timing_df)){
  temp_yr<-mig_timing_df[i,]
  temp_d1_spdf<-d1_spdf[d1_spdf$ptt==mig_timing_df$ptt[i],]
  temp_d1_nbod_spdf<-d1_nbod_spdf[d1_nbod_spdf$ptt==mig_timing_df$ptt[i],]
  temp_so_tab_spdf<-stopovers_tab_spdf[stopovers_tab_spdf$ptt==mig_timing_df$ptt[i],]
  # writeOGR(temp_d1_spdf,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation",layer="temp_d1_spdf",driver="ESRI Shapefile")
  # writeOGR(temp_so_tab_spdf,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation",layer="temp_so_tab_spdf",driver="ESRI Shapefile")
  
  mig_timing_df$br_lat[i]<-br_cent_coords[br_cent_coords$ptt==mig_timing_df$ptt[i],]$lat #For comparison with UK milestone timing
  mig_timing_df$br_long[i]<-br_cent_coords[br_cent_coords$ptt==mig_timing_df$ptt[i],]$long #For comparison with UK milestone timing
  
  #Which is the final stopover and fix for this bird?
  temp_last_so_id<-temp_so_tab_spdf[which.max(temp_so_tab_spdf$SO_end),]$id #Previously $SO_end_hybrid
  temp_last_fix_id<-temp_d1_spdf[which.max(temp_d1_spdf$timestamp),]$event.id
  mig_timing_df$first_fix_TS[i]<-as.character(min(temp_d1_spdf$timestamp))
  mig_timing_df$last_fix_TS[i]<-as.character(max(temp_d1_spdf$timestamp))
  
  #M1. 'Departure from breeding grounds'
  # First fix outside breeding grounds buffers [* this differs from M8, which is first fix *inside* breeding grounds buffer]
  # BUT the position error in the non-best-of-day data takes genuine breeding season points outside the 50km buffer
  # To avoid problems with position error in non-best-of-day data, I use the LAST fix from breeding grounds buffer in the first calendar year of migratory cycle.
  # - Alternatively I could use best-of-day data 
  # - LAST-TYPE MILESTONE - but all birds included here left the UK, so don't need to check whether this is the final fix.
  temp_br_buff<-br_buffs_ll[br_buffs_ll$ptt==mig_timing_df$ptt[i],]
  # #Previous approach - last point outside buffer
  # temp_inbr_log<-is.na(over(temp_d1_spdf,temp_br_buff))==F #Which points are within the breeding grounds?
  # temp_dep_brgr_ts<-max(temp_d1_spdf[which(temp_d1_spdf$year==mig_timing_df$year[i] & temp_inbr_log),]$timestamp)
  # mig_timing_df$dep_brgr[i]<-yday(temp_dep_brgr_ts)
  # # UNCERTAINTY: how long until next fix?
  # mig_timing_df_uncertainty$dep_brgr[i]<-as.numeric(difftime(min(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp>temp_dep_brgr_ts]),temp_dep_brgr_ts,units="days"))
  
  #NEW APPROACH: the day before the first 50km movement between successive BOD fixes in the UK (unless that is before the day of the last fix from the breeding grounds).
  #NB this is after the longest period without 50km movements, otherwise it will pick up birds e.g. 50026 making a UK stopover before the breeding grounds
  #NB have to add that this period has to be within the breeding area, otherwise it will pick up birds like 62608 which spent longer on stopover in GB than on the breeding grounds
  temp_wld_spdf<-wld_spdf[which(wld_spdf$ptt==mig_timing_df$ptt[i] & wld_spdf$year==mig_timing_df$year[i] & wld_spdf$in_gb==1 & is.na(over(wld_spdf,temp_br_buff))==F),]
  if(any(is.na(unique(temp_wld_spdf$LOS))==F)==F){temp_end_longest_period_sub50km<-temp_wld_spdf$timestamp[1]} else { #For PTT 161323, there are no periods without 50km movements: the bird immediately leaves
    temp_longest_period_sub50km<-max(unique(temp_wld_spdf$LOS),na.rm=T) #The longest period without 50km movements. Leverage fact that migratory groups are delineated by >50km movements.
    temp_end_longest_period_sub50km<-max(temp_wld_spdf[which(temp_wld_spdf$LOS==temp_longest_period_sub50km),]$timestamp)
  }
  #Which is the next fix?
  temp_subs_fix<-temp_d1_spdf[which(temp_d1_spdf$timestamp==temp_end_longest_period_sub50km)+1,]$timestamp #NB use d1 as wld is subsetted to GB and fix might not be in GB
  
  #Check against full dataset: were there any fixes from outside the breeding area on the same day as temp_end_longest_period_sub50km?
  #Which are the fixes from the full dataset which fell on the same day as and after temp_end_longest_period_sub50km?
  temp_nbod_end_longest_period_sub50km<-temp_d1_nbod_spdf[which(year(temp_d1_nbod_spdf$timestamp)==year(temp_end_longest_period_sub50km) & 
                                                                  yday(temp_d1_nbod_spdf$timestamp)==yday(temp_end_longest_period_sub50km) &
                                                                  temp_d1_nbod_spdf$timestamp>temp_end_longest_period_sub50km),]
  if(length(temp_nbod_end_longest_period_sub50km)==0){temp_dep_brgr_ts<-max(temp_end_longest_period_sub50km,temp_subs_fix-(24*3600))} else {
    if(any(is.na(over(temp_nbod_end_longest_period_sub50km,temp_br_buff)))==F){temp_dep_brgr_ts<-max(temp_end_longest_period_sub50km,temp_subs_fix-(24*3600))} else {
      temp_dep_brgr_ts<-temp_end_longest_period_sub50km
    }
  }
  mig_timing_df$dep_brgr[i]<-yday(temp_dep_brgr_ts)
  mig_timing_df$dep_brgr_TS[i]<-as.character(temp_dep_brgr_ts) #For definition of migratory cycles for attribution of death timing
  
  # UNCERTAINTY: how long until next fix?
  mig_timing_df_uncertainty$dep_brgr[i]<-as.numeric(difftime(temp_subs_fix,temp_end_longest_period_sub50km,units="days"))
  
  
  #M1B 'Departure from the UK'
  #Day before first fix outside of UK (and over land), unless that is before the day of the last fix in the UK
  temp_uk_fixes<-temp_d1_spdf[which(temp_d1_spdf$year==mig_timing_df$year[i] & temp_d1_spdf$in_gb==1),]
  temp_last_uk_fix<-max(temp_uk_fixes$timestamp)
  temp_nonuk_fixes<-temp_d1_spdf[temp_d1_spdf$year==mig_timing_df$year[i] & temp_d1_spdf$timestamp>temp_last_uk_fix & is.na(over(temp_d1_spdf,world_shp)$NAME)==F,]
  temp_first_nonuk_fix<-min(temp_nonuk_fixes$timestamp)
  mig_timing_df$dep_UK[i]<-ifelse(yday(temp_first_nonuk_fix)>yday(temp_last_uk_fix),yday(temp_first_nonuk_fix)-1,yday(temp_first_nonuk_fix))
  mig_timing_df$dep_UK_TS[i]<-as.character(temp_first_nonuk_fix)
  
  # UNCERTAINTY: how long between last fix in UK and first fix outside UK?
  mig_timing_df_uncertainty$dep_UK[i]<-difftime(temp_last_uk_fix,temp_first_nonuk_fix,units="days")
  
  
  #M2. 'Southbound departure from Europe / N Africa'
  #The day before the first fix south of, >50km from and after the last major stopover north of the Sahara (or same day, if fix was in the evening of departure).
  #NB THIS CAN BE QUITE CONSERVATIVE - e.g. ptt 50017 left Europe/N Africa on the 10th Aug, but this method declares it the 9th
  #First, find the stopovers north of the Sahara in the first calendar year of the migratory cycle
  # - LAST-TYPE MILESTONE -
  so_n_of_sah<-temp_so_tab_spdf[as.logical(temp_so_tab_spdf$SO_median_lat>22 & temp_so_tab_spdf$SO_end>temp_dep_brgr_ts & temp_so_tab_spdf$SO_year==mig_timing_df$year[i] & is.na(over(temp_so_tab_spdf,sah_spdf))),]
  #If there are no stopovers between leaving the UK and the Sahara (e.g. for 134957) then skip
  if(nrow(so_n_of_sah)==0) {mig_timing_df$dep_euna_sb[i]<-NA} else {
    mig_timing_df$lat_last_euro_so_sbound[i]<-so_n_of_sah$SO_median_lat[nrow(so_n_of_sah)]
    mig_timing_df$any_so_s_europe_sbound[i]<-any(so_n_of_sah$SO_median_lat<47) #For interest - did the bird make any stopovers in S Europe (defined as <47N, i.e. S of north border of Switzerland)
    #Is the last stopover north of the Sahara the last stopover of all for this bird? If so skip.
    if(so_n_of_sah[which.max(so_n_of_sah$SO_end),]$id==temp_last_so_id) {mig_timing_df$dep_euna_sb[i]<-NA} else {
      #What is the earliest timestamp south of, >50km from and after the last stopover north of the Sahara?
      #NB no need to create buffer to identify points >50km from stopover. Stopover is by definition delineated by >50km movement
      fixes_aft_last_so_n_of_sah<-temp_d1_spdf$timestamp[as.numeric(which(temp_d1_spdf$timestamp>max(so_n_of_sah$SO_end) & temp_d1_spdf$location.lat<so_n_of_sah[which.max(so_n_of_sah$SO_end),]$SO_median_lat))]
      if(length(fixes_aft_last_so_n_of_sah)==0) {mig_timing_df$dep_euna_sb[i]<-NA} else {
        earl_fix_aft_last_so_n_of_sah<-yday(min(fixes_aft_last_so_n_of_sah))
        #UNCERTAINTY: how long since previous fix?
        mig_timing_df_uncertainty$dep_euna_sb[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<min(fixes_aft_last_so_n_of_sah)]),min(fixes_aft_last_so_n_of_sah),units="days"))
        #If fix is same day as departure from stopover, make this milestone on that day; otherwise make this milestone the day preceding the fix
        if(earl_fix_aft_last_so_n_of_sah==yday(max(so_n_of_sah$SO_end))) {mig_timing_df$dep_euna_sb[i]<-earl_fix_aft_last_so_n_of_sah} else {
          mig_timing_df$dep_euna_sb[i]<-earl_fix_aft_last_so_n_of_sah-1
        }
      }
    }
  }
  
  #M3. 'Completion of Sahara crossing southbound'
  # The day of the first location in the first major stopover south of the Sahara.
  # NB restricting to stopovers in the first calendar year of the migratory cycle. 
  so_s_of_sah<-temp_so_tab_spdf[as.logical(temp_so_tab_spdf$SO_median_lat<22 & temp_so_tab_spdf$SO_end>temp_dep_brgr_ts & temp_so_tab_spdf$SO_year==mig_timing_df$year[i] & is.na(over(temp_so_tab_spdf,sah_spdf))),]
  if(nrow(so_s_of_sah)==0) {mig_timing_df$fin_sah_sb[i]<-NA} else {
    mig_timing_df$fin_sah_sb_TS[i]<-as.character(min(so_s_of_sah$SO_start))
    mig_timing_df$fin_sah_sb[i]<-yday(min(so_s_of_sah$SO_start)) #NB using sensu stricto version of stopover start
    #UNCERTAINTY: how long since previous fix?
    mig_timing_df_uncertainty$fin_sah_sb[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<min(so_s_of_sah$SO_start)]),min(so_s_of_sah$SO_start),units="days"))
  }
  
  
  #M4. 'Arrival in wintering grounds'
  #M4a. The day after the final fix north of the forest zone (due to problems with charging and transmission within the forest zone)
  #Which was the first fix in the forest zone?
  forest_zone_tss<-temp_d1_spdf[which(over(temp_d1_spdf,ca_rf_spdf)=="CAfr rainforest" & temp_d1_spdf$timestamp>temp_dep_brgr_ts & temp_d1_spdf$timestamp<(temp_dep_brgr_ts+(365*60*60*24))),]$timestamp
  if(length(forest_zone_tss)==0) {mig_timing_df$arr_wgr_a[i]<-NA} else {
    mig_timing_df$arr_wgr_a_TS[i]<-as.character(min(forest_zone_tss))
    #Then find the preceding fix
    #Take the day after that unless that would put it on the day after the first fix in the forest zone.
    #UNCERTAINTY: how long until the next fix?
    mig_timing_df_uncertainty$arr_wgr_a[i]<-as.numeric(difftime(min(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp>max(temp_d1_spdf[temp_d1_spdf$timestamp<min(forest_zone_tss),]$timestamp)]),max(temp_d1_spdf[temp_d1_spdf$timestamp<min(forest_zone_tss),]$timestamp),units="days"))
    #Refine using NBOD data - there are a couple of birds (121792 & 146759, both in 2015) for which the first fix in the forest zone is the same day as the final fix north of the forest zone on the NBOD data, but not in the BOD data
    forest_zone_nbod_tss<-temp_d1_nbod_spdf[which(over(temp_d1_nbod_spdf,ca_rf_spdf)=="CAfr rainforest" & 
                                                    temp_d1_nbod_spdf$timestamp>temp_dep_brgr_ts & 
                                                    temp_d1_nbod_spdf$timestamp<(temp_dep_brgr_ts+(365*60*60*24))),]$timestamp    
    if(yday(max(temp_d1_spdf[temp_d1_spdf$timestamp<min(forest_zone_tss),]$timestamp))==yday(min(forest_zone_tss))) {mig_timing_df$arr_wgr_a[i]<-yday(min(forest_zone_tss)) } else {
      if(yday(max(temp_d1_spdf[temp_d1_spdf$timestamp<min(forest_zone_tss),]$timestamp))==yday(min(forest_zone_nbod_tss))) {mig_timing_df$arr_wgr_a[i]<-yday(min(forest_zone_nbod_tss)) } else {
        mig_timing_df$arr_wgr_a[i]<-yday(max(temp_d1_spdf[temp_d1_spdf$timestamp<min(forest_zone_tss),]$timestamp))+1
      }
    } 
  }
  
  #M4b. The first day at the stopover with minimum latitude
  #Which are the stopovers within the forest zone (NB forest zone defined as rainforest and any area further south)?
  #NB these are typically quite a lot later than M4a.
  forest_zone_sos<-matrix(nrow=0,ncol=1)
  if(length(forest_zone_tss)==0) {mig_timing_df$arr_wgr_b[i]<-NA} else{
    mig_timing_df$dep_wgr_a_TS_ford1plot[i]<-as.character(max(forest_zone_tss)) #For plotting only
    forest_zone_sos<-temp_so_tab_spdf[temp_so_tab_spdf$SO_start>=min(forest_zone_tss) & temp_so_tab_spdf$SO_start<max(forest_zone_tss),]
    if(nrow(forest_zone_sos)==0) {mig_timing_df$arr_wgr_b[i]<-NA} else{ #Need this clause in because there might be fixes but no stopovers in the forest zone
      #Of these stopovers in the forest zone, which is the furthest south?
      mig_timing_df$arr_wgr_b[i]<-yday(forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_start)
      mig_timing_df$arr_wgr_b_TS[i]<-as.character(forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_start)
      mig_timing_df$dep_wgr_b_TS[i]<-as.character(forest_zone_sos[which.max(forest_zone_sos$SO_median_lat),]$SO_start) #For visualisation (Fig 1) only
      mig_timing_df$min_lat_stopover[i]<-min(forest_zone_sos$SO_median_lat)
      mig_timing_df$min_lat_stopover_id[i]<-forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$id
      #UNCERTAINTY: how long since previous fix?
      mig_timing_df_uncertainty$arr_wgr_b[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_start]),forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_start,units="days"))
    }
  }
  
  #M5. 'Departure from wintering grounds'
  #M5a. Day before the first day outside of the forest zone.
  #NB this doesn't quite take into account uncertainty in the same way as M4a
  if(length(forest_zone_tss)==0){mig_timing_df$dep_wgr_a[i]<-NA} else {
    tss_after_forest_zone<-temp_d1_spdf[temp_d1_spdf$timestamp>max(forest_zone_tss) &
                                          is.na(over(temp_d1_spdf,world_shp)$NAME)==F,]$timestamp #Need the clause about being over land, otherwise point marginally over water will count
    if(length(tss_after_forest_zone)==0) {mig_timing_df$dep_wgr_a[i]<-NA} else {
      mig_timing_df$dep_wgr_a[i]<-yday(min(tss_after_forest_zone))-1
      mig_timing_df$dep_wgr_a_TS[i]<-as.character(min(tss_after_forest_zone))
      #UNCERTAINTY: how long since the previous fix?
      mig_timing_df_uncertainty$dep_wgr_a[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<min(tss_after_forest_zone)]),min(tss_after_forest_zone),units="days"))
    }
  }
  
  
  #M5b. Last day at stopover with minimum latitude (interpolate if missing duty cycles)
  #NB these typically fall before M5a
  # - LAST-TYPE MILESTONE - 
  if(nrow(forest_zone_sos)==0) {mig_timing_df$dep_wgr_b[i]<-NA} else { #Need this clause in because there might be fixes but no stopovers in the forest zone
    #Of these stopovers in the forest zone, which is the furthest south?
    #If this stopover is the final stopover altogether for the bird, and if there are no fixes after the forest zone, skip
    if(forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$id==temp_last_so_id &
       length(tss_after_forest_zone)==0) {mig_timing_df$dep_wgr_b[i]<-NA} else {
         #UNCERTAINTY: how long until the next fix?
         mig_timing_df_uncertainty$dep_wgr_b[i]<-as.numeric(difftime(min(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp>forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_end_LOS]),forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_end_LOS,units="days")) #NB using $SO_end_LOS here, because otherwise (using $SO_end) by definition the end of the stopover would be dependent on the next fix
         if(mig_timing_df_uncertainty$dep_wgr_b[i]<2){mig_timing_df$dep_wgr_b[i]<-yday(forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_end_LOS)} else { #No need to interpolate if no missing duty cycles
           mig_timing_df$dep_wgr_b[i]<-yday(min(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp>forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_end_LOS]))-2
           # mig_timing_df$dep_wgr_b[i]<-yday(forest_zone_sos[which.min(forest_zone_sos$SO_median_lat),]$SO_end_LOS+((mig_timing_df_uncertainty$dep_wgr_b[i]/2)*24*3600)) #Interpolate: half way between sensu-stricto end of stopover and the next fix
         }
       }
  }
  
  #M6. 'Departure from West Africa'
  #M6a. Last day at final stopover before Sahara crossing (interpolate if missing duty cycles)
  # - LAST-TYPE MILESTONE - 
  #Which is the last stopover south of the Sahara (and not in the forest zone) in the 2nd calendar year of the migratory cycle?
  wa_sos<-temp_so_tab_spdf[which(temp_so_tab_spdf$SO_year==(mig_timing_df$year[i]+1) & temp_so_tab_spdf$SO_end<(temp_dep_brgr_ts+(365*60*60*24)) & temp_so_tab_spdf$SO_median_lat<22 & is.na(over(temp_so_tab_spdf,sah_spdf)) & is.na(over(temp_so_tab_spdf,ca_rf_spdf))),]
  if(nrow(wa_sos)>0){
    mig_timing_df$last_wa_sos_long[i]<-wa_sos[which.max(wa_sos$SO_end),]$SO_median_long
    mig_timing_df$last_wa_sos_lat[i]<-wa_sos[which.max(wa_sos$SO_end),]$SO_median_lat
  }
  if(nrow(wa_sos)==0) {mig_timing_df$dep_wa_a[i]<-NA} else {
    #Are there any fixes north of West Africa?
    tss_after_wa<-temp_d1_spdf[which(temp_d1_spdf$timestamp>max(wa_sos$SO_end) & temp_d1_spdf$timestamp<(temp_dep_brgr_ts+(365*60*60*24)) & (temp_d1_spdf$location.lat>22 | over(temp_d1_spdf,sah_spdf)=="Sahara")),]
    #If the last stopover in West Africa is the very last for this bird AND there are no fixes after West Africa, skip
    if(wa_sos[which.max(wa_sos$SO_end),]$id==temp_last_so_id & nrow(tss_after_wa)==0) {mig_timing_df$dep_wa_a[i]<-NA} else {
      #UNCERTAINTY: how long until the next fix?
      mig_timing_df_uncertainty$dep_wa_a[i]<-as.numeric(difftime(min(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp>max(wa_sos$SO_end_LOS)]),max(wa_sos$SO_end_LOS),units="days"))  #NB using $SO_end_LOS here, because otherwise (using $SO_end) by definition the end of the stopover would be dependent on the next fix
      if(mig_timing_df_uncertainty$dep_wa_a[i]<2){mig_timing_df$dep_wa_a[i]<-yday(max(wa_sos$SO_end))} else { #No need to interpolate if no missing duty cycles
        mig_timing_df$dep_wa_a[i]<-yday(max(wa_sos$SO_end_LOS)+((mig_timing_df_uncertainty$dep_wa_a[i]/2)*24*3600)) #Interpolate: half way between sensu-stricto end of stopover and the next fix
      }
    }
  }
  
  #M6b. First fix during desert crossing - check bird is moving though, rather than at short stopover in Sahara
  # - *Use fixes with movement type 'M'? Or any fix that was not during a stopover?
  #This will mostly be blanks
  #Which fixes are in the Sahara?
  fixes_in_sah_nb<-temp_d1_spdf[which(temp_d1_spdf$year==(mig_timing_df$year[i]+1) & temp_d1_spdf$timestamp<(temp_dep_brgr_ts+(365*60*60*24)) & is.na(over(temp_d1_spdf,sah_spdf))==F),]
  if(nrow(fixes_in_sah_nb)==0) {mig_timing_df$dep_wa_b[i]<-NA} else {
    #Which of these fixes are given movement type 'M' - i.e. geographically distant from proceeding and following fixes
    moving_tss_in_sah_nb<-withLOSdata.all[(withLOSdata.all$event.id %in% fixes_in_sah_nb$event.id) & withLOSdata.all$mtype=="M",]$timestamp
    if(length(moving_tss_in_sah_nb)==0) {mig_timing_df$dep_wa_b[i]<-NA} else {
      mig_timing_df$dep_wa_b[i]<-yday(min(moving_tss_in_sah_nb))
      #UNCERTAINTY: how long since last fix?
      mig_timing_df_uncertainty$dep_wa_b[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<min(moving_tss_in_sah_nb)]),min(moving_tss_in_sah_nb),units="days"))
    }
  }
  
  #M6c. New hybrid milestone for departure from West Africa - using Saharan fixes preferentially, otherwise WA fixes.
  wa_fixes<-temp_d1_spdf[which(temp_d1_spdf$year==(mig_timing_df$year[i]+1) & temp_d1_spdf$timestamp<(temp_dep_brgr_ts+(365*60*60*24)) & temp_d1_spdf$location.lat<22 & is.na(over(temp_d1_spdf,sah_spdf)) & is.na(over(temp_d1_spdf,ca_rf_spdf))),]
  #Add in lat/long of last fix in West Africa
  if(nrow(wa_fixes)>0){
    mig_timing_df$last_wa_fix_long[i]<-wa_fixes[which.max(wa_fixes$timestamp),]$location.long
    mig_timing_df$last_wa_fix_lat[i]<-wa_fixes[which.max(wa_fixes$timestamp),]$location.lat
    mig_timing_df$dep_wa_TS[i]<-as.character(max(wa_fixes$timestamp))
  }
  moving_fixes_in_sah_nb<-withLOSdata.all[which((withLOSdata.all$event.id %in% fixes_in_sah_nb$event.id) & withLOSdata.all$distance>50),] #Use distance>50 rather than mtype="M", because the last fix of all (e.g. for 161321) will get a blank mtype
  if(nrow(wa_fixes)==0 & nrow(moving_fixes_in_sah_nb)==0){ mig_timing_df$dep_wa[i]<-NA 
  #If no fixes in Sahara, dep_wa is day of last fix in West Africa
  } else if(nrow(moving_fixes_in_sah_nb)==0) {
    if(wa_fixes[which.max(wa_fixes$timestamp),]$event.id==temp_last_fix_id) {mig_timing_df$dep_wa[i]<-NA} else {
      mig_timing_df$dep_wa[i]<-yday(max(wa_fixes$timestamp))
      #UNCERTAINTY: how long until the next fix?
      mig_timing_df_uncertainty$dep_wa[i]<-as.numeric(difftime(min(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp>max(wa_fixes$timestamp)]),max(wa_fixes$timestamp),units="days"))
    }
  } else {
    #UNCERTAINTY: if we have fixes in the Sahara, we are relatively confident that they are max 2 days after leaving W Africa.
    #Uniquely (among milestones), we give uncertainty a fixed value here, regardless of how long since the last fix.
    mig_timing_df_uncertainty$dep_wa[i]<-2
    first_moving_sah_nb_tss<-moving_fixes_in_sah_nb[1,]$timestamp
    if(moving_fixes_in_sah_nb[1,]$lat<20){
      #If there are fixes in Sahara, then dep_wa is 1 or 2 days before first Sahara fix: if N of 20deg then bird left 2 days earlier, otherwise 1 day.
      mig_timing_df$dep_wa[i]<-max(yday(first_moving_sah_nb_tss)-1,max(yday(wa_fixes$timestamp)))
    } else {
      mig_timing_df$dep_wa[i]<-max(yday(first_moving_sah_nb_tss)-2,max(yday(wa_fixes$timestamp)))
    }
  }
  
  
  #M7. 'Completion of Sahara crossing northbound'
  # The day of the first location in the first major stopover north of Sahara (and not in the UK, and not after returning to the UK)
  #So what is the first return date to the UK breeding grounds?
  first_return_br_tss<-temp_d1_spdf[as.logical(temp_d1_spdf$year==(mig_timing_df$year[i]+1) & is.na(over(temp_d1_spdf,temp_br_buff))==F),]$timestamp
  if(length(first_return_br_tss)==0){ #If the bird never returns to the UK, then ignore the clause about this milestone having to be before the bird returns to the UK
    so_n_of_sah_nb<-temp_so_tab_spdf[as.logical(temp_so_tab_spdf$SO_median_lat>22 & temp_so_tab_spdf$SO_year==(mig_timing_df$year[i]+1) & temp_so_tab_spdf$SO_start<(temp_dep_brgr_ts+(365*60*60*24)) & temp_so_tab_spdf$country!="United Kingdom" & is.na(over(temp_so_tab_spdf,sah_spdf))),]
    if(nrow(so_n_of_sah_nb)==0) {mig_timing_df$fin_sah_nb[i]<-NA} else {
      # so_n_of_sah_nb_list[[i]]<-so_n_of_sah_nb[which.min(so_n_of_sah_nb$SO_start),]
      mig_timing_df$fin_sah_nb[i]<-yday(min(so_n_of_sah_nb$SO_start))
      #UNCERTAINTY: how long since previous fix?
      mig_timing_df_uncertainty$fin_sah_nb[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<min(so_n_of_sah_nb$SO_start)]),min(so_n_of_sah_nb$SO_start),units="days"))
    }
  } else { #If the bird does return to the UK, then apply the clause that this milestone must be before the bird returns to the UK
    so_n_of_sah_nb<-temp_so_tab_spdf[as.logical(temp_so_tab_spdf$SO_median_lat>22 & temp_so_tab_spdf$SO_year==(mig_timing_df$year[i]+1) & temp_so_tab_spdf$SO_start<min(first_return_br_tss) & temp_so_tab_spdf$country!="United Kingdom" & is.na(over(temp_so_tab_spdf,sah_spdf))),]
    if(nrow(so_n_of_sah_nb)==0) {mig_timing_df$fin_sah_nb[i]<-NA} else {
      # so_n_of_sah_nb_list[[i]]<-so_n_of_sah_nb[which.min(so_n_of_sah_nb$SO_start),]
      mig_timing_df$fin_sah_nb[i]<-yday(min(so_n_of_sah_nb$SO_start))
      #UNCERTAINTY: how long since previous fix?
      mig_timing_df_uncertainty$fin_sah_nb[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<min(so_n_of_sah_nb$SO_start)]),min(so_n_of_sah_nb$SO_start),units="days"))
    }
  }
  
  #M8A. 'Arrival to the UK'
  #Reinstating a previously used milestone. The first fix in the UK in the second calendar year of the migratory cycle.
  # temp_UK_return_fixes<-temp_d1_spdf[which(temp_d1_spdf$year==mig_timing_df$year[i]+1 & temp_d1_spdf$in_gb==1),]$timestamp #Previous
  temp_UK_return_fixes<-temp_d1_spdf[which(temp_d1_spdf$year==mig_timing_df$year[i]+1 & temp_d1_spdf$country=="United Kingdom"),]$timestamp
  
  if(length(temp_UK_return_fixes)==0) {mig_timing_df$arr_UK[i]<-NA} else {
    temp_arr_UK_ts<-min(temp_UK_return_fixes)
    mig_timing_df$arr_UK[i]<-yday(temp_arr_UK_ts)
    mig_timing_df$arr_UK_TS[i]<-as.character(temp_arr_UK_ts)
    
    # UNCERTAINTY: how long since the last fix?
    mig_timing_df_uncertainty$arr_UK[i]<-difftime(temp_arr_UK_ts,max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<temp_arr_UK_ts]),units="days")
  }
  
  #M8. 'Arrival in breeding grounds'
  # Previously first fix within breeding grounds buffer in second calendar year of migratory cycle.
  #Now start of the longest period in GB (and in breeding grounds) free of >50km movements in 2nd calendar year of mig cycle
  
  if(length(first_return_br_tss)==0) {mig_timing_df$arr_brgr[i]<-NA} else {
    mig_timing_df$arr_brgr_TS[i]<-as.character(min(first_return_br_tss))
    # mig_timing_df$arr_brgr[i]<-yday(min(first_return_br_tss))
    # #UNCERTAINTY: how long since previous fix?
    # mig_timing_df_uncertainty$arr_brgr[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<min(first_return_br_tss)]),min(first_return_br_tss),units="days"))
    # 
    #New approach: start of the longest period in GB (and in breeding grounds) free of >50km movements in 2nd calendar year of mig cycle
    temp_wld_spdf_yr2<-wld_spdf[which(wld_spdf$ptt==mig_timing_df$ptt[i] & wld_spdf$year==mig_timing_df$year[i]+1 & wld_spdf$in_gb==1 & is.na(over(wld_spdf,temp_br_buff))==F),]
    temp_longest_period_yr2_sub50km<-max(unique(temp_wld_spdf_yr2$LOS),na.rm=T) #The longest period without 50km movements. Leverage fact that migratory groups are delineated by >50km movements.
    temp_start_longest_period_yr2_sub50km<-min(temp_wld_spdf_yr2[which(temp_wld_spdf_yr2$LOS==temp_longest_period_yr2_sub50km),]$timestamp)
    mig_timing_df$arr_brgr[i]<-yday(temp_start_longest_period_yr2_sub50km)
    
    #UNCERTAINTY: how long since previous fix?
    mig_timing_df_uncertainty$arr_brgr[i]<-as.numeric(difftime(max(temp_d1_spdf$timestamp[temp_d1_spdf$timestamp<temp_start_longest_period_yr2_sub50km]),temp_start_longest_period_yr2_sub50km,units="days"))
  }
  
  #Write out location of different milestones
}

#We know Vigilamus made it back to the UK
mig_timing_df[mig_timing_df$ptt=="146757",]$arr_UK<-114
mig_timing_df_uncertainty[mig_timing_df_uncertainty$ptt=="146757",]$arr_UK<-3

#Which milestones are we interested in?
milestone_names<-c("dep_brgr","dep_UK","dep_euna_sb","fin_sah_sb","arr_wgr_a","arr_wgr_b","dep_wgr_b","dep_wgr_a","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr")

#Omit milestones we're no longer using
milestone_names<-milestone_names[(milestone_names %in% c("dep_euna_sb",
                                                         # "arr_wgr_b", #Maybe keep in
                                                         # "dep_wgr_b", #Maybe keep in
                                                         "dep_wa_a",
                                                         "dep_wa_b",
                                                         "fin_sah_nb"))==F]

milestone_names_inds<-which(names(mig_timing_df)[3:16] %in% milestone_names)
#After which milestone did the bird die?
death_tab$milestone_prec_death_name<-death_tab$milestone_prec_death<-death_tab$last_mig_cycle<-NA
for(i in 1:nrow(death_tab)){
  temp_ptt_dt<-death_tab[i,]$ptt
  if(death_tab[i,]$ptt %in% mig_timing_df$ptt==F) next
  temp_mtd<-mig_timing_df[mig_timing_df$ptt==temp_ptt_dt,]
  temp_last_year<-temp_mtd[nrow(temp_mtd),]
  if(any(is.na(temp_last_year[3:16])==F)==F) next #If all milestones are NA then skip
  if(is.na(temp_last_year$arr_brgr)==F) {
    death_tab$milestone_prec_death[i]<-13
    death_tab$milestone_prec_death_name[i]<-"arr_brgr"
    death_tab$last_mig_cycle[i]<-temp_mtd[nrow(temp_mtd),]$year
  } else { #i.e. if there are values for all milestones, arr_brgr is the last one before the bird's death
    death_tab$milestone_prec_death[i]<-max(which(is.na(temp_last_year[3:16])==F))
    #We omit some milestones from the analysis. *Of the milestones we're interested in*, which is the last milestone before death?
    death_tab$milestone_prec_death_name[i]<-names(mig_timing_df)[3:16][max(milestone_names_inds[milestone_names_inds<=death_tab$milestone_prec_death[i]])]
    death_tab$last_mig_cycle[i]<-temp_last_year$year
  }
}

# so_n_of_sah_nb_list_ul<-unlist(so_n_of_sah_nb_list)
# so_n_of_sah_nb_all<-do.call(rbind,so_n_of_sah_nb_list_ul)
# writeOGR(so_n_of_sah_nb_all,dsn="H:/JGD H drive backup 160320/Raw data/Cuckoos/Temporary shapefiles for visualisation",layer="so_n_of_sah_nb_all",driver="ESRI Shapefile")

#UNCERTAINTY
#How does uncertainty vary across the migratory milestones?
colMeans(mig_timing_df_uncertainty[,milestone_names],na.rm=T)
#Remove uncertain milestones: where the time since the previous fix or until the next fix is more than a threshold number of days

#Number of events and birds per milestone BEFORE REMOVING UNCERTAIN MILESTONES, for results tables
n_events_per_milestone_incuncertain<-unlist(lapply(apply(is.na(mig_timing_df[,milestone_names])==F,2,which),length))
n_birds_per_milestone_incuncertain<-numeric(length(milestone_names))
for(i in 1:length(milestone_names)){
  n_birds_per_milestone_incuncertain[i]<-length(unique(mig_timing_df$ptt[is.na(mig_timing_df[,milestone_names[i]])==F]))
}

threshold<-5
mig_timing_df[,3:16][which(abs(mig_timing_df_uncertainty[,3:16])>threshold,arr.ind=T)]<-NA #Remove uncertain milestones

#Add breeding habitat (lowland/upland)
mig_timing_df[which(mig_timing_df$ptt %in% c("50017","50023","50024","50026","50031","50032","50042","50046","50051","50053","62518","62520","62602","62608","62688","115586","115589","115595","121791","121792","128296","128299","128302","134950","134951","134952","134953","134955","134956","134957","134959","134960","134963","134964","146761","161318","161320","161321","161322","161323","161324","161325","179106","179107","179108","179109")),]$br_hab<-"Lowland"
mig_timing_df[which(mig_timing_df$ptt %in% c("115588","115590","115591","115593","115594","115596","115597","115598","115599","115600","115602","128297","128298","128300","128303","128304","134958","134961","146753","146754","146756","146757","146758","146759","146760","161319")),]$br_hab<-"Upland"

#Add breeding region
mig_timing_df$region<-character(nrow(mig_timing_df))
mig_timing_df[which(mig_timing_df$ptt %in% c("50026",
                                             "62518",
                                             "62602",
                                             "62688",
                                             "115586",
                                             "115589",
                                             "115595",
                                             "128296",
                                             "128301",
                                             "128302",
                                             "134964",
                                             "161322")),]$region<-"Norfolk Broads"
mig_timing_df[which(mig_timing_df$ptt %in% c("50031",
                                             "50032",
                                             "62520",
                                             "62608",
                                             "121792",
                                             "128299",
                                             "134950",
                                             "134960",
                                             "161318",
                                             "161324",
                                             "161325",
                                             "170429",
                                             "179106",
                                             "179107",
                                             "179108")),]$region<-"Thetford Forest"
mig_timing_df[which(mig_timing_df$ptt %in% c("50046",
                                             "50051",
                                             "50053",
                                             "121791",
                                             "134951",
                                             "134963",
                                             "170432")),]$region<-"South Downs" #Anything on south coast that is clearly east of the New Forest - fair?
mig_timing_df[which(mig_timing_df$ptt %in% c("50017",
                                             "50024",
                                             "134953",
                                             "134956",
                                             "134959",
                                             "161320",
                                             "161321",
                                             "170434",
                                             "179109")),]$region<-"New Forest"
mig_timing_df[which(mig_timing_df$ptt %in% c("115599",
                                             "128297",
                                             "128298",
                                             "128304",
                                             "134958",
                                             "134961",
                                             "161319")),]$region<-"Dartmoor"
mig_timing_df[which(mig_timing_df$ptt %in% c("115594",
                                             "115596",
                                             "115597",
                                             "115598",
                                             "146760")),]$region<-"Wales"
mig_timing_df[which(mig_timing_df$ptt %in% c("50023",
                                             "50042",
                                             "134952",
                                             "134955",
                                             "134957",
                                             "146761",
                                             "161323",
                                             "170433",
                                             "170435")),]$region<-"Sherwood Forest"
mig_timing_df[which(mig_timing_df$ptt %in% c("115590",
                                             "128300",
                                             "128303")),]$region<-"Skye & Lochalsh"
mig_timing_df[which(mig_timing_df$ptt %in% c("115588",
                                             "115591",
                                             "115593",
                                             "115600",
                                             "115602")),]$region<-"Trossachs"
mig_timing_df[which(mig_timing_df$ptt %in% c("146754",
                                             "146756",
                                             "146757")),]$region<-"North York Moors"
mig_timing_df[which(mig_timing_df$ptt %in% c("146753",
                                             "146758",
                                             "146759")),]$region<-"Forest of Bowland"
mig_timing_df$ptt[which(mig_timing_df$region=="")]


##Assign migratory direction southbound out of Europe
#Visualise migration tracks for each bird
world_shp_simp<-gSimplify(world_shp,0.2)
# un_ptt<-unique(mig_timing_df$ptt)
# i<-70
# mig_timing_df[which(mig_timing_df$ptt==un_ptt[i]),]
# temp_df<-d1[d1$ptt==un_ptt[i],c("location.long","location.lat")]
# temp_df_arrows<-as.matrix(temp_df)
# colnames(temp_df_arrows)<-NULL
# temp_df_arrows<-cbind(temp_df_arrows,rbind(temp_df_arrows[2:nrow(temp_df_arrows),],matrix(0,ncol=2,nrow=1)))
# temp_df_arrows<-temp_df_arrows[-nrow(temp_df_arrows),]
# plot(temp_df,main=un_ptt[i])
# plot(world_shp,add=T,border="gray70",col="pale green")
# for(i in 1:nrow(temp_df_arrows)){
#   arrows(temp_df_arrows[i,1],temp_df_arrows[i,2],temp_df_arrows[i,3],temp_df_arrows[i,4],length=0.1,col="red")
# }

#Based on region from which they made their final crossing to Africa
#NB FINAL crossing (161322 and 161323 made initial crossing from France to Algeria or near Algeria, before crossing back to Spain and from there to Morocco)
#If incomplete southward journey but following a complete southward journey in previous cycle, assume same direction
mig_timing_df[which(mig_timing_df$ptt %in% c("50017","50023","50026","50031","50032",
                                             "50046","50051","50053","62520",
                                             "62688","115586","115589","115595",
                                             "128296","128298","128299","128302","128304",
                                             "134951","134952","134953",
                                             "134955","134957","134960","134961",
                                             "134964","161321","161322",
                                             "161323","161324","179106",
                                             "179107","179108","179109")),]$mig_dir<-"SW"
mig_timing_df[which(mig_timing_df$ptt %in% c("50042","62518","62602","62608","115588",
                                             "115590","115591","115593",
                                             "115594","115596","115597",
                                             "115598","115599","115600",
                                             "115602","121791","121792","128300",
                                             "128303","134956","134958","134959",
                                             "134963","146753","146754","146756",
                                             "146757","146758","146759",
                                             "146760","146761","161319",
                                             "161320")),]$mig_dir<-"SE"
mig_timing_df[which(mig_timing_df$ptt %in% c("50024","134950","161325")),]$mig_dir<-NA #Unclear route
mig_timing_df[mig_timing_df$ptt=="128297" & mig_timing_df$year==2013,]$mig_dir<-"SW" #Started as if to go SE
mig_timing_df[mig_timing_df$ptt=="128297" & mig_timing_df$year==2014,]$mig_dir<-"SE"
mig_timing_df[mig_timing_df$ptt=="161318" & mig_timing_df$year==2016,]$mig_dir<-"SE"
mig_timing_df[mig_timing_df$ptt=="161318" & mig_timing_df$year==2017,]$mig_dir<-"SE"
mig_timing_df[mig_timing_df$ptt=="161318" & mig_timing_df$year==2018,]$mig_dir<-"SW"
mig_timing_df[mig_timing_df$ptt=="161318" & mig_timing_df$year==2019,]$mig_dir<-"SW"
mig_timing_df[mig_timing_df$ptt=="161318" & mig_timing_df$year==2020,]$mig_dir<-"SE"
mig_timing_df[mig_timing_df$ptt=="161318" & mig_timing_df$year==2021,]$mig_dir<-"SW" 

#Now add a 3-category migratory direction classification: anything that doesn't have a fix in mainland Iberia or mainland Italy/Balkans is 'Central'
mig_timing_df$mig_dir_3cat<-NA
mig_timing_df[which(mig_timing_df$ptt %in% c("50023","50026","50031","50032",
                                             "50046","50051","50053","62520",
                                             "62688","115586","115589","115595","128296",
                                             "128299","128302","128304","134951",
                                             "134952","134953","134955","134957","134960",
                                             "134961","134964","161321","161322","161323",
                                             "161324","179106","179107","179109")),]$mig_dir_3cat<-"SW"
mig_timing_df[which(mig_timing_df$ptt %in% c("50017","50024","115599","121791",
                                             "128298","134958","161325",
                                             "179108")),]$mig_dir_3cat<-"C"
mig_timing_df[which(mig_timing_df$ptt %in% c("50042","62518","62602","62608",
                                             "115588","115590","115591",
                                             "115593","115594","115596","115597",
                                             "115598","115600","115602","128300",
                                             "128303","134956","134959",
                                             "134963","146753","146754","146756",
                                             "146757","146758","146759","146760",
                                             "146761","161319","161320")),]$mig_dir_3cat<-"SE"
mig_timing_df[which(mig_timing_df$ptt=="128297" & mig_timing_df$year==2013),]$mig_dir_3cat<-"SW"
mig_timing_df[which(mig_timing_df$ptt=="128297" & mig_timing_df$year==2014),]$mig_dir_3cat<-"SE"
mig_timing_df[which(mig_timing_df$ptt=="121792" & mig_timing_df$year==2014),]$mig_dir_3cat<-"SE"
mig_timing_df[which(mig_timing_df$ptt=="121792" & mig_timing_df$year==2015),]$mig_dir_3cat<-"C"
mig_timing_df[which(mig_timing_df$ptt=="161318" & mig_timing_df$year %in% c(2016,2017)),]$mig_dir_3cat<-"C"
mig_timing_df[which(mig_timing_df$ptt=="161318" & mig_timing_df$year %in% c(2018,2019,2021)),]$mig_dir_3cat<-"SW"
mig_timing_df[which(mig_timing_df$ptt=="161318" & mig_timing_df$year==2020),]$mig_dir_3cat<-"SE"
mig_timing_df[which(mig_timing_df$ptt=="134950"),]$mig_dir_3cat<-NA

mig_dir_3cat_col<-character(nrow(mig_timing_df))
mig_dir_3cat_col[mig_timing_df$mig_dir_3cat=="SW"]<-"orange"
mig_dir_3cat_col[mig_timing_df$mig_dir_3cat=="C"]<-"blue"
mig_dir_3cat_col[mig_timing_df$mig_dir_3cat=="SE"]<-"red"
mig_dir_3cat_col[is.na(mig_timing_df$mig_dir_3cat)]<-"black"


# # #Plot different broad migratory routes
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/SE.png",10,7,res=800,units="in")
# plot(world_shp,border="gray70",col="pale green",xlim=range(d1$location.long),ylim=range(d1$location.lat),main="SE (3 category)")
# for(j in which(mig_timing_df$mig_dir_3cat=="SE")){
#   temp_df<-d1[d1$ptt==mig_timing_df$ptt[j] & d1$year==mig_timing_df$year[j] & d1$julian>150,c("location.long","location.lat")]
#   temp_df_arrows<-as.matrix(temp_df)
#   colnames(temp_df_arrows)<-NULL
#   temp_df_arrows<-cbind(temp_df_arrows,rbind(temp_df_arrows[2:nrow(temp_df_arrows),],matrix(0,ncol=2,nrow=1)))
#   temp_df_arrows<-temp_df_arrows[-nrow(temp_df_arrows),]
#   for(i in 1:nrow(temp_df_arrows)){
#     arrows(temp_df_arrows[i,1],temp_df_arrows[i,2],temp_df_arrows[i,3],temp_df_arrows[i,4],length=0,col=mig_dir_3cat_col[j])
#   }
# }
# dev.off()


# #Plot out migration routes for checking migratory route classifications
# for(j in 1:nrow(mig_timing_df)){
#   png(paste0("H:/JGD H drive backup 160320/Cuckoos/Plots/Mig dirs by ptt year/",mig_timing_df$ptt[j],"_",mig_timing_df$year[j],".png"),10,7,res=800,units="in")
#   plot(world_shp,border="gray70",col="pale green",xlim=range(d1$location.long),ylim=range(d1$location.lat),
#        main=paste0("PTT ",mig_timing_df$ptt[j],", ",mig_timing_df$year[j],"\n",mig_timing_df$mig_dir[j]," (2 cat), ",mig_timing_df$mig_dir_3cat[j]," (3 cat)"))
#   temp_df<-d1[d1$ptt==mig_timing_df$ptt[j] & d1$year==mig_timing_df$year[j] & d1$julian>150,c("location.long","location.lat")]
#   temp_df_arrows<-as.matrix(temp_df)
#   colnames(temp_df_arrows)<-NULL
#   temp_df_arrows<-cbind(temp_df_arrows,rbind(temp_df_arrows[2:nrow(temp_df_arrows),],matrix(0,ncol=2,nrow=1)))
#   temp_df_arrows<-temp_df_arrows[-nrow(temp_df_arrows),]
#   for(i in 1:nrow(temp_df_arrows)){
#     arrows(temp_df_arrows[i,1],temp_df_arrows[i,2],temp_df_arrows[i,3],temp_df_arrows[i,4],length=0,col="black")
#   }
#   dev.off()
# }


#Sanity-checking against Mark's equivalent table
mm_mig_timing_df<-read.csv("H:/JGD H drive backup 160320/Raw data/Cuckoos/Additional data from MM/stopover_bestofday_2018_1daymin_recalc_ANNUAL_mig_summary_dead_final.csv")

mm_mig_timing_df$arr_brgr_mine<-mm_mig_timing_df$dep_wgr_b_mine<-NA
for(i in 1:nrow(mm_mig_timing_df)){
  mm_mig_timing_df$dep_brgr_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$dep_brgr
  mm_mig_timing_df$dep_UK_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$dep_UK
  mm_mig_timing_df$arr_UK_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$arr_UK
  mm_mig_timing_df$arr_brgr_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$arr_brgr
  mm_mig_timing_df$arr_wgr_a_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$arr_wgr_a
  mm_mig_timing_df$dep_wgr_a_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$dep_wgr_a
  mm_mig_timing_df$dep_wgr_b_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$dep_wgr_b
  mm_mig_timing_df$dep_wa_a_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$dep_wa_a
  mm_mig_timing_df$fin_sah_sb_mine[i]<-mig_timing_df[mig_timing_df$year==mm_mig_timing_df[i,]$year-1 & mig_timing_df$ptt==mm_mig_timing_df[i,]$ptt,]$fin_sah_sb
}
mm_mig_timing_df$depart_winterSO365<-ifelse(mm_mig_timing_df$depart_winterSO<1,mm_mig_timing_df$depart_winterSO+365,mm_mig_timing_df$depart_winterSO)

# #6 milestones in common
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/dep_UK.png",6,6,res=800,units="in")
# plot(mm_mig_timing_df$dep_UK_mine,mm_mig_timing_df$DEPuk,xlab="Mine",ylab="Mark's",
#      main="Departure from the UK");abline(0,1)
# dev.off()
# # plot(mm_mig_timing_df$arr_wgr_a_mine,mm_mig_timing_df$DEPsahel,xlab="Mine",ylab="Mark's", #Not quite equivalent - Mark's is departure from the Sahel
# #      main="Arrival to the wintering grounds");abline(0,1)
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/dep_wint.png",6,6,res=800,units="in")
# plot(mm_mig_timing_df$dep_wgr_a_mine,mm_mig_timing_df$DEPcentralAF,xlab="Mine",ylab="Mark's",
#      main="Departure from the wintering grounds");abline(0,1)
# dev.off()
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/dep_WA.png",6,6,res=800,units="in")
# plot(mm_mig_timing_df$dep_wa_a_mine,mm_mig_timing_df$DEPwestAF,xlab="Mine",ylab="Mark's",
#      main="Departure from West Africa");abline(0,1)
# dev.off()
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/arr_UK.png",6,6,res=800,units="in")
# plot(mm_mig_timing_df$arr_UK_mine,mm_mig_timing_df$arrive_uk,xlab="Mine",ylab="Mark's",
#      main="Arrival to the UK");abline(0,1)
# dev.off()
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/arr_br.png",6,6,res=800,units="in")
# plot(mm_mig_timing_df$arr_brgr_mine,mm_mig_timing_df$arrive_breeding,xlab="Mine",ylab="Mark's",
#      main="Arrival to the breeding grounds");abline(0,1)
# dev.off()

# #Map departures from breeding grounds to compare with Mark's
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/dep_wint_map_mine.png",8,8,res=800,units="in")
# plot(world_shp,xlim=c(-10,30),ylim=c(-12,10),col="gray90")
# plot(ca_rf_spdf,col="green",add=T)
# for(i in 1:nrow(mm_mig_timing_df)){
#   if(is.na(mm_mig_timing_df$dep_wgr_a_mine[i])) next
#   temp_ds<-d1_spdf[d1_spdf$ptt==mm_mig_timing_df$ptt[i] & d1_spdf$year==mm_mig_timing_df$year[i],]
#   pre_dep_ind<-max(which(temp_ds$julian<mm_mig_timing_df$dep_wgr_a_mine[i]))
#   post_dep_ind<-min(which(temp_ds$julian>mm_mig_timing_df$dep_wgr_a_mine[i]))
#   points(temp_ds[pre_dep_ind,],col="yellow",pch=20)
#   points(temp_ds[post_dep_ind,],col="red",pch=20)
#   arrow_coords<-c(as.numeric(coordinates(temp_ds[pre_dep_ind,])),as.numeric(coordinates(temp_ds[post_dep_ind,])))
#   arrows(arrow_coords[1],arrow_coords[2],arrow_coords[3],arrow_coords[4],col="blue",length=0)
# }
# legend("topright",legend=c("Last fix before day of milestone","First fix after day of milestone"),pch=c(20,20),col=c("yellow","red"))
# dev.off()
# 
# png("H:/JGD H drive backup 160320/Cuckoos/Plots/dep_wint_map_marks.png",8,8,res=800,units="in")
# plot(world_shp,xlim=c(-12,30),ylim=c(-12,22),col="gray90")
# plot(ca_rf_spdf,col="green",add=T)
# for(i in 1:nrow(mm_mig_timing_df)){
#   if(is.na(mm_mig_timing_df$DEPcentralAF[i])) next
#   temp_ds<-d1_spdf[d1_spdf$ptt==mm_mig_timing_df$ptt[i] & d1_spdf$year==mm_mig_timing_df$year[i],]
#   pre_dep_ind<-max(which(temp_ds$julian<mm_mig_timing_df$DEPcentralAF[i]))
#   post_dep_ind<-min(which(temp_ds$julian>mm_mig_timing_df$DEPcentralAF[i]))
#   points(temp_ds[pre_dep_ind,],col="yellow",pch=20)
#   points(temp_ds[post_dep_ind,],col="red",pch=20)
#   arrow_coords<-c(as.numeric(coordinates(temp_ds[pre_dep_ind,])),as.numeric(coordinates(temp_ds[post_dep_ind,])))
#   arrows(arrow_coords[1],arrow_coords[2],arrow_coords[3],arrow_coords[4],col="blue",length=0)
#   print(i)
#   print(arrow_coords[3])
# }
# legend("topright",legend=c("Last fix before day of milestone","First fix after day of milestone"),pch=c(20,20),col=c("yellow","red"))
# dev.off()

mig_timing_df_365<-mig_timing_df

#Convert 2nd calendar year dates to days after the start of the 1st calendar year (differs for leap year)
sub_200_ind_ly<-which(mig_timing_df[leap_year(mig_timing_df$year),c("arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr")]<200,arr.ind=T)
sub_200_ind_nly<-which(mig_timing_df[leap_year(mig_timing_df$year)==F,c("arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr")]<200,arr.ind=T)
mig_timing_df[leap_year(mig_timing_df$year),c("arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr")][sub_200_ind_ly]<-mig_timing_df[leap_year(mig_timing_df$year),c("arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr")][sub_200_ind_ly]+366
mig_timing_df[leap_year(mig_timing_df$year)==F,c("arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr")][sub_200_ind_nly]<-mig_timing_df[leap_year(mig_timing_df$year)==F,c("arr_wgr_a","arr_wgr_b","dep_wgr_a","dep_wgr_b","dep_wa_a","dep_wa_b","dep_wa","fin_sah_nb","arr_UK","arr_brgr")][sub_200_ind_nly]+365

# mig_timing_df$time_brgr_2_sah<-mig_timing_df$dep_euna_sb-mig_timing_df$dep_brgr
# mig_timing_df$time_in_sah_sb<-mig_timing_df$fin_sah_sb-mig_timing_df$dep_euna_sb
# mig_timing_df$time_sah_2_wgr<-mig_timing_df$arr_wgr_a-mig_timing_df$fin_sah_sb
# mig_timing_df$time_in_wgr<-mig_timing_df$dep_wgr_a-mig_timing_df$arr_wgr_a
# mig_timing_df$time_wgr_2_sah<-mig_timing_df$dep_wa_a-mig_timing_df$dep_wgr_a
# mig_timing_df$time_in_sah_nb<-mig_timing_df$fin_sah_nb-mig_timing_df$dep_wa_a
# mig_timing_df$time_sah_2_brgr<-mig_timing_df$arr_brgr-mig_timing_df$fin_sah_nb

#Birds taking SE route spend longer going from the UK to the Sahara
#Birds taking SW route spend longer in Sahel on way south before wintering grounds, and less time in wintering grounds
#Birds taking central route spend nearly twice as long in W Africa after leaving wintering grounds
#Birds from all routes spend the same amount of time in both Sahara crossings, and travelling from the Sahara to the UK

#After removing uncertain milestones, how many migratory cycles do we have?
length(which(is.na(as.numeric(as.matrix(mig_timing_df[,milestone_names]))==F)))
length(which(apply(is.na(mig_timing_df[,milestone_names])==F,1,any))) #118 cycles
length(unique(mig_timing_df$ptt[which(apply(is.na(mig_timing_df[,milestone_names])==F,1,any))])) #72 birds
mtd_ptt_freq_df<-data.frame(table(mig_timing_df$ptt[which(apply(is.na(mig_timing_df[,milestone_names])==F,1,any))]))
length(which(mtd_ptt_freq_df$Freq>1))/nrow(mtd_ptt_freq_df)

plot(as.numeric(mig_timing_df[i,3:16]),type="l",ylim=c(153,509))
for(i in 1:nrow(mig_timing_df)){
  lines(as.numeric(mig_timing_df[i,3:16]))
}
#When do birds die?
mig_timing_df$last_milestone_pre_death<-mig_timing_df$outcome<-NA
for(i in 1:nrow(death_tab)){
  if(is.na(death_tab$last_mig_cycle[i])) next
  mig_timing_df[mig_timing_df$ptt==death_tab[i,]$ptt & 
                  mig_timing_df$year==death_tab[i,]$last_mig_cycle,]$last_milestone_pre_death<-death_tab[i,]$milestone_prec_death_name
  mig_timing_df[mig_timing_df$ptt==death_tab[i,]$ptt & 
                  mig_timing_df$year==death_tab[i,]$last_mig_cycle,]$outcome<-death_tab[i,]$outcome_code
}

mig_timing_df$died_this_cycle<-FALSE
mig_timing_df$died_this_cycle[which(is.na(mig_timing_df$last_milestone_pre_death)==F)]<-TRUE

mig_timing_df$season_died<-NA
mig_timing_df$season_died[mig_timing_df$last_milestone_pre_death %in% milestone_names[1:5]]<-"Autumn"
mig_timing_df$season_died[mig_timing_df$last_milestone_pre_death %in% milestone_names[6:10]]<-"Spring"


for(i in 1:nrow(death_tab)){
  points(death_tab$milestone_prec_death[i],mig_timing_df[which(mig_timing_df$ptt==death_tab[i,]$ptt & mig_timing_df$year==death_tab[i,]$last_mig_cycle),death_tab[i,]$milestone_prec_death+2],col="red",pch=20,cex=2)
}
#BIRDS THAT LEAVE THE UK LATE ARE MUCH MORE LIKELY TO DIE
temp_app<-apply(is.na(mig_timing_df[,3:16])==F,2,which)
milestones_surv_df<-data.frame("milestones_made"=as.numeric(unlist(lapply(temp_app,length))))
milestones_surv_df$ind<-1:nrow(milestones_surv_df)
death_freq_df<-data.frame(table(death_tab$milestone_prec_death))

milestones_surv_df$death_freq<-0
for(i in 1:nrow(death_freq_df)){
  milestones_surv_df[as.numeric(as.character(death_freq_df[i,]$Var1)),]$death_freq<-death_freq_df[i,]$Freq
}
milestones_surv_df$death_prob<-milestones_surv_df$death_freq/milestones_surv_df$milestones_made
milestones_surv_df$milestone<-names(mig_timing_df[,3:16])

died_after_ms2_tab<-death_tab[which(death_tab$milestone_prec_death==2),]
mig_timing_df$died_after_ms2<-logical(nrow(mig_timing_df))

for(i in 1:nrow(died_after_ms2_tab)){
  mig_timing_df[which(mig_timing_df$ptt==died_after_ms2_tab[i,]$ptt & mig_timing_df$year==died_after_ms2_tab[i,]$last_mig_cycle),]$died_after_ms2<-TRUE
}
boxplot(mig_timing_df$dep_UK~mig_timing_df$died_after_ms2)

#Add in age at ringing. NB not possible to age birds as 6 - all birds in spreadsheet as '6' given '4' here.
mig_timing_df$age_at_ring<-NA
mig_timing_df[which(mig_timing_df$ptt %in% c("50017","50024","50026","50031",
                                             "50032","50042","50046","50051",
                                             "50053","170429","170435","179108",
                                             "115598","161320","161321","121792",
                                             "146753","146754","146755","146756",
                                             "146757","146758","146759","146761",
                                             "146762","115599","128296","128297",
                                             "128298","128301","134951","134963",
                                             "170432","170433","170434","179106",
                                             "179107","179109","115590","62602",
                                             "115586","115588","115592","115593",
                                             "115594","115595","115600","121791",
                                             "128300","128302","128303","128304",
                                             "134952","134957","170430","170431")),]$age_at_ring<-4
mig_timing_df[which(mig_timing_df$ptt %in% c("50023","161319","161323","62518",
                                             "62520","115589","115590","115591",
                                             "115596","115597","115602","128295",
                                             "128299","134960","62608","62688",
                                             "134953","146760","161318","161324",
                                             "161325","134964")),]$age_at_ring<-5

#Add known age (migratory cycles since fledging)
mig_timing_df$known_age<-numeric(nrow(mig_timing_df))
mig_timing_df$known_age<-NA
ever_known_5s<-unique(mig_timing_df[which(mig_timing_df$age_at_ring==5),]$ptt)
for(i in 1:length(ever_known_5s)){
  temp_5s_df<-mig_timing_df[mig_timing_df$ptt==ever_known_5s[i],]
  mig_timing_df[mig_timing_df$ptt==ever_known_5s[i],]$known_age<-temp_5s_df$year-(min(temp_5s_df$year)-1)
}

#Add known age (0 for mig cycle definitely started year after fledging, 1 for mig cycle definitely started more than one year after fledging, NA for unknown)
mig_timing_df$known_age_binary<-numeric(nrow(mig_timing_df))
mig_timing_df$known_age_binary<-NA
un_m_ptt<-unique(mig_timing_df$ptt)
for(i in 1:length(un_m_ptt)){
  temp_agebin_df<-mig_timing_df[mig_timing_df$ptt==un_m_ptt[i],]
  if(nrow(temp_agebin_df)==1) next
  mig_timing_df[mig_timing_df$ptt==un_m_ptt[i],][which(temp_agebin_df$year>min(temp_agebin_df$year)),]$known_age_binary<-1
}
mig_timing_df[which(mig_timing_df$known_age==1),]$known_age_binary<-0

#Add cycles since tagging for birds ringed as 4s, to compare with effects of known_age, to tease apart effects of age and time since taggin
mig_timing_df$cycles_since_tagging_4s<-numeric(nrow(mig_timing_df))
mig_timing_df$cycles_since_tagging_4s<-NA
ringed_as_4s<-unique(mig_timing_df[which(mig_timing_df$age_at_ring==4),]$ptt)
for(i in 1:length(ringed_as_4s)){
  temp_4s_df<-mig_timing_df[mig_timing_df$ptt==ringed_as_4s[i],]
  mig_timing_df[mig_timing_df$ptt==ringed_as_4s[i],]$cycles_since_tagging_4s<-temp_4s_df$year-(min(temp_4s_df$year))
}

#Now add cycles since tagging for all birds
mig_timing_df$cycles_since_tagging<-numeric(nrow(mig_timing_df))
mig_timing_df$cycles_since_tagging<-NA
un_ptt<-unique(mig_timing_df$ptt)
for(i in 1:length(un_ptt)){
  temp_df<-mig_timing_df[mig_timing_df$ptt==un_ptt[i],]
  mig_timing_df[mig_timing_df$ptt==un_ptt[i],]$cycles_since_tagging<-temp_df$year-(min(temp_df$year))
}


#Is there a relationship between latitude of minimum latitude stopover, and date of arrival/departure to/from that stopover?
plot(mig_timing_df$min_lat_stopover,
     mig_timing_df$arr_wgr_b,
     xlab="Latitude of min lat stopover",
     ylab="Date of arrival to min lat stopover",
     main="Arrival")
plot(mig_timing_df$min_lat_stopover,
     mig_timing_df$dep_wgr_b,
     xlab="Latitude of min lat stopover",
     ylab="Date of departure from to min lat stopover",
     main="Departure")


#Milestone timing vs breeding latitude
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/lat vs milestone.png",10,10,res=800,units="in")
par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){
  plot(mig_timing_df$br_lat,mig_timing_df[,milestone_names[i]],xlab="Latitude of breeding site",ylab=milestone_names[i])
}
par(mfrow=c(1,1))
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/lat vs dep_brgr to dep_UK.png",7,7,res=800,units="in")
plot(mig_timing_df$br_lat,mig_timing_df$dep_UK-mig_timing_df$dep_brgr,xlab="Latitude of breeding site",ylab="Days between departing breeding site and departing UK")
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/lat vs arr_UK to arr_brgr.png",7,7,res=800,units="in")
plot(mig_timing_df$br_lat,mig_timing_df$arr_brgr-mig_timing_df$arr_UK,xlab="Latitude of breeding site",ylab="Days between arriving to UK and arriving to breeding site")
dev.off()

#Milestone timing vs year
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/year vs milestone.png",10,10,res=800,units="in")
par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){
  plot(mig_timing_df$year,mig_timing_df[,milestone_names[i]],xlab="Year",ylab=milestone_names[i])
}
par(mfrow=c(1,1))
dev.off()

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
# cbind(n_events_per_milestone_withinind,n_birds_per_milestone_withinind)


#Migratory milestone and migration stage length data now ready for analysis

# temp_mtd_forchecking<-mig_timing_df[,c("ptt","br_hab","region")]
# temp_mtd_forchecking_unique<-aggregate(temp_mtd_forchecking,by=list(temp_mtd_forchecking$ptt),FUN=unique)
# temp_mtd_forchecking_unique<-temp_mtd_forchecking_unique[order(as.numeric(temp_mtd_forchecking_unique$ptt)),]
# write.csv(temp_mtd_forchecking_unique[,c("ptt","br_hab","region")],"H:/JGD H drive backup 160320/Cuckoos/Output data/Cuckoo br_hab and regions for checking.csv",row.names=F)

# write.csv(mig_timing_df,"H:/JGD H drive backup 160320/Cuckoos/Output data/mig_timing_df.csv",row.names=F)


########################################################################################################################
######################################### 4. ANALYSIS OF MIGRATORY TIMING ################################################
#4A: NATURE OF VARIATION IN MILESTONES
#Histograms of variation within each milestone
par(mfrow=c(length(milestone_names),1))
par(mar=c(2.1, 4.1, 0.1, 2.1))
for(i in 1:length(milestone_names)){
  h<-hist(mig_timing_df_std[,milestone_names[i]],
          breaks=seq(range(mig_timing_df_std[,milestone_names],na.rm=T)[1]-1,
                     range(mig_timing_df_std[,milestone_names],na.rm=T)[2]+1,by=2),
          main="",xlab="")
  text(x=range(mig_timing_df_std[,milestone_names],na.rm=T)[2]-20,
       y=max(h$counts)/2,
       labels=milestone_names[i])
}
par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,1))

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

plot(var_cor_df$PTT_var,var_cor_df$total_var)
plot(var_cor_df$resid_var,var_cor_df$total_var)

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
write.csv(nice_var_cor_df_tab,"H:/JGD H drive backup 160320/Cuckoos/Output data/nice_var_cor_df_tab.csv",row.names=F)

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

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig3a.png",6,6,res=800,units="in")
par(mar=c(5.1, 4.1, 1.1, 1.1))
barplot(as.matrix(var_abs_df),
        col=coul,
        las=2,
        ylab="Variance",
        main="")
legend("topright",legend=row.names(var_abs_df),fill=coul)
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig3c.png",6,6,res=800,units="in")
par(mar=c(5.1, 4.1, 1.1, 1.1))
barplot(as.matrix(repeatability_df),
        ylim=c(0,1.1),
        col=coul,
        las=2,
        ylab="Proportion of variance",
        main="")
legend("topright",legend=row.names(repeatability_df),fill=coul)
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig3_5.png",7,7,res=800,units="in")
par(mar=c(10.1, 4.1, 4.1, 2.1))
plot(apply(var_abs_df,2,sum)[2:ncol(var_abs_df)]-apply(var_abs_df,2,sum)[1:(ncol(var_abs_df)-1)],ylab="Change in total variance of migratory timing",xlab="",xaxt="n")
abline(h=0,lty=3)
axis(1, at=1:(length(milestone_names)-1),
     labels=paste0(milestone_names[1:(length(milestone_names)-1)]," : ",milestone_names[2:length(milestone_names)]),las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Between- vs within-individual variation
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA2.png",7,7,res=800,units="in")
par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names)-0.1,var_cor_df$PTT_var,ylim=c(0,max(c(var_cor_df$PTT_var_uci,var_cor_df$resid_var_uci))),pch=19,bty="n",xlab="",ylab="Variance",xaxt="n",col=coul[1])
points(1:length(milestone_names)+0.1,var_cor_df$resid_var,pch=19,col=coul[2])
arrows(1:length(milestone_names)-0.1,var_cor_df$PTT_var_lci,1:length(milestone_names)-0.1,var_cor_df$PTT_var_uci,length=0,col=coul[1])
arrows(1:length(milestone_names)+0.1,var_cor_df$resid_var_lci,1:length(milestone_names)+0.1,var_cor_df$resid_var_uci,length=0,col=coul[2])
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
legend("topright",legend=c("Between-individual","Within-individual"),pch=19,col=coul,bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig3b.png",6,6,res=800,units="in")
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(var_cor_df$PTT_var,var_cor_df$total_var,
     pch=21,
     bg=coul[1],
     bty="n",
     xlim=c(min(var_cor_df$PTT_var_lci),max(var_cor_df$PTT_var_uci)),
     ylim=c(min(var_cor_df$total_var_lci),max(var_cor_df$total_var_uci)),
     xlab="Between-individual variance",ylab="Total variance")
arrows(var_cor_df$PTT_var_lci,var_cor_df$total_var,var_cor_df$PTT_var_uci,var_cor_df$total_var,length=0,col="gray50")
arrows(var_cor_df$PTT_var,var_cor_df$total_var_lci,var_cor_df$PTT_var,var_cor_df$total_var_uci,length=0,col="gray50")
m1<-lm(total_var~PTT_var,data=var_cor_df)
newx<-seq(min(var_cor_df$PTT_var_lci),max(var_cor_df$PTT_var_uci),by=0.1)
ps<-predict(m1, newdata=data.frame(PTT_var=newx), interval="confidence", level = 0.95)
lines(newx,ps[,1],lty=1)
lines(newx,ps[,2],lty=2)
lines(newx,ps[,3],lty=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig3d.png",6,6,res=800,units="in")
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(var_cor_df$resid_var,var_cor_df$total_var,
     pch=21,
     bg=coul[2],
     bty="n",
     xlim=c(min(var_cor_df$resid_var_lci),max(var_cor_df$resid_var_uci)),
     ylim=c(min(var_cor_df$total_var_lci),max(var_cor_df$total_var_uci)),
     xlab="Within-individual variance",ylab="Total variance")
arrows(var_cor_df$resid_var_lci,var_cor_df$total_var,var_cor_df$resid_var_uci,var_cor_df$total_var,length=0,col="gray50")
arrows(var_cor_df$resid_var,var_cor_df$total_var_lci,var_cor_df$resid_var,var_cor_df$total_var_uci,length=0,col="gray50")
m1<-lm(total_var~resid_var,data=var_cor_df)
newx<-seq(min(var_cor_df$resid_var_lci),max(var_cor_df$resid_var_uci),by=0.1)
ps<-predict(m1, newdata=data.frame(resid_var=newx), interval="confidence", level = 0.95)
lines(newx,ps[,1],lty=1)
lines(newx,ps[,2],lty=2)
lines(newx,ps[,3],lty=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#What is the distribution of *individual* anomalies for each milestone?
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_1.png",7,5,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
boxplot(mtd_stdind_melted$value~mtd_stdind_melted$variable,ylab="Within-individual anomaly",xlab="",xaxt="n",frame=F)
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Visualising individual birds' tracks through the year
mig_timing_df_stdind$ptt_num<-as.numeric(as.factor(mig_timing_df_stdind$ptt)) #For colour

#First: all birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig1d.png",6,6,res=800,units="in")
par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly",xlab="",xaxt="n",main="")
abline(h=0,lty=3)
for(i in 1:nrow(mig_timing_df_stdind)){
  temp_mtdsi<-mig_timing_df_stdind[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=mig_timing_df_stdind$ptt_num[i],lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),
          col=mig_timing_df_stdind$ptt_num[i],lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
# for(i in 1:nrow(mig_timing_df_stdind)){
#   if(is.na(mig_timing_df_stdind[i,]$last_milestone_pre_death)) next
#   if(mig_timing_df_stdind[i,]$outcome=="M"){temp_col<-"black"} else 
#     if(mig_timing_df_stdind[i,]$outcome=="UA"){temp_col<-"blue"} else {
#       temp_col<-"red"
#     }
#   points(which(milestone_names==mig_timing_df_stdind[i,]$last_milestone_pre_death),
#          mig_timing_df_stdind[i,mig_timing_df_stdind[i,]$last_milestone_pre_death],pch=20,cex=2,col=temp_col)
# }
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Just SW birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_4a.png",7,7,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly",xlab="",xaxt="n",main="SW")
abline(h=0,lty=3)
for(i in which(mig_timing_df_stdind$mig_dir=="SW")){
  temp_mtdsi<-mig_timing_df_stdind[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=mig_timing_df_stdind$ptt_num[i],lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),col=mig_timing_df_stdind$ptt_num[i],lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
for(i in which(mig_timing_df_stdind$mig_dir=="SW")){
  if(is.na(mig_timing_df_stdind[i,]$last_milestone_pre_death)) next
  points(which(milestone_names==mig_timing_df_stdind[i,]$last_milestone_pre_death),
         mig_timing_df_stdind[i,mig_timing_df_stdind[i,]$last_milestone_pre_death],pch=20,cex=2)
}
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Just SE birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_4b.png",7,7,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly",xlab="",xaxt="n",main="SE")
abline(h=0,lty=3)
for(i in which(mig_timing_df_stdind$mig_dir=="SE")){
  temp_mtdsi<-mig_timing_df_stdind[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=mig_timing_df_stdind$ptt_num[i],lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),col=mig_timing_df_stdind$ptt_num[i],lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
for(i in which(mig_timing_df_stdind$mig_dir=="SE")){
  if(is.na(mig_timing_df_stdind[i,]$last_milestone_pre_death)) next
  points(which(milestone_names==mig_timing_df_stdind[i,]$last_milestone_pre_death),
         mig_timing_df_stdind[i,mig_timing_df_stdind[i,]$last_milestone_pre_death],pch=20,cex=2)
}
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Colour by the order in which birds left the breeding grounds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_5.png",7,7,res=800,units="in")
mig_timing_df_stdind_noNA_depbrgr<-mig_timing_df_stdind[which(is.na(mig_timing_df_stdind$dep_brgr)==F),]
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind_noNA_depbrgr[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly",xlab="",xaxt="n")
abline(h=0,lty=3)
for(i in 1:nrow(mig_timing_df_stdind_noNA_depbrgr)){
  temp_mtdsi<-mig_timing_df_stdind_noNA_depbrgr[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=rgb(rank(mig_timing_df_stdind_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_noNA_depbrgr),
                0,1-(rank(mig_timing_df_stdind_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_noNA_depbrgr))),
        lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),
          col=rgb(rank(mig_timing_df_stdind_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_noNA_depbrgr),
                  0,1-(rank(mig_timing_df_stdind_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_noNA_depbrgr))),
          lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
for(i in 1:nrow(mig_timing_df_stdind)){
  if(is.na(mig_timing_df_stdind[i,]$last_milestone_pre_death)) next
  points(which(milestone_names==mig_timing_df_stdind[i,]$last_milestone_pre_death),
         mig_timing_df_stdind[i,mig_timing_df_stdind[i,]$last_milestone_pre_death],pch=20,cex=2)
}
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#Now do again but relative to population median
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig1b.png",6,6,res=800,units="in")
par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_std[,milestone_names],na.rm=T),
     ylab="Anomaly relative to population median",xlab="",xaxt="n",main="")
abline(h=0,lty=3)
for(i in 1:nrow(mig_timing_df_std)){
  temp_mtdsi<-mig_timing_df_std[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=mig_timing_df_std$ptt_num[i],lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),
          col=mig_timing_df_std$ptt_num[i],lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
# for(i in 1:nrow(mig_timing_df_std)){
#   if(is.na(mig_timing_df_std[i,]$last_milestone_pre_death)) next
#   if(mig_timing_df_stdind[i,]$outcome=="M"){temp_col<-"black"} else 
#     if(mig_timing_df_stdind[i,]$outcome=="UA"){temp_col<-"blue"} else {
#     temp_col<-"red"
#   }
#   points(which(milestone_names==mig_timing_df_std[i,]$last_milestone_pre_death),
#          mig_timing_df_std[i,mig_timing_df_std[i,]$last_milestone_pre_death],pch=20,cex=2,col=temp_col)
# }
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Just SW birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_20.png",7,7,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_std[,milestone_names],na.rm=T),
     ylab="Anomaly relative to population median",xlab="",xaxt="n",main="SW")
abline(h=0,lty=3)
for(i in which(mig_timing_df_std$mig_dir=="SW")){
  temp_mtdsi<-mig_timing_df_std[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=mig_timing_df_std$ptt_num[i],lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),col=mig_timing_df_std$ptt_num[i],lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
for(i in which(mig_timing_df_std$mig_dir=="SW")){
  if(is.na(mig_timing_df_std[i,]$last_milestone_pre_death)) next
  points(which(milestone_names==mig_timing_df_std[i,]$last_milestone_pre_death),
         mig_timing_df_std[i,mig_timing_df_std[i,]$last_milestone_pre_death],pch=20,cex=2)
}
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Just SE birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA10.png",7,7,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_std[,milestone_names],na.rm=T),
     ylab="Anomaly relative to population median",xlab="",xaxt="n",main="SE")
abline(h=0,lty=3)
for(i in which(mig_timing_df_std$mig_dir=="SE")){
  temp_mtdsi<-mig_timing_df_std[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=mig_timing_df_std$ptt_num[i],lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),col=mig_timing_df_std$ptt_num[i],lwd=2,lty=2)
  }
}
for(i in which(mig_timing_df_std$mig_dir=="SE")){
  if(is.na(mig_timing_df_std[i,]$last_milestone_pre_death)) next
  points(which(milestone_names==mig_timing_df_std[i,]$last_milestone_pre_death),
         mig_timing_df_std[i,mig_timing_df_std[i,]$last_milestone_pre_death],pch=20,cex=2)
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Colour by the order in which birds left the breeding grounds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA11.png",7,7,res=800,units="in")
mig_timing_df_std_noNA_depbrgr<-mig_timing_df_std[which(is.na(mig_timing_df_std$dep_brgr)==F),]
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_std_noNA_depbrgr[,milestone_names],na.rm=T),
     ylab="Anomaly relative to population median",xlab="",xaxt="n")
abline(h=0,lty=3)
for(i in 1:nrow(mig_timing_df_std_noNA_depbrgr)){
  temp_mtdsi<-mig_timing_df_std_noNA_depbrgr[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=rgb(rank(mig_timing_df_std_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_std_noNA_depbrgr),
                0,1-(rank(mig_timing_df_std_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_std_noNA_depbrgr))),
        lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),
          col=rgb(rank(mig_timing_df_std_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_std_noNA_depbrgr),
                  0,1-(rank(mig_timing_df_std_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_std_noNA_depbrgr))),
          lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
for(i in 1:nrow(mig_timing_df_std)){
  if(is.na(mig_timing_df_std[i,]$last_milestone_pre_death)) next
  points(which(milestone_names==mig_timing_df_std[i,]$last_milestone_pre_death),
         mig_timing_df_std[i,mig_timing_df_std[i,]$last_milestone_pre_death],pch=20,cex=2)
}
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#Now do again but relative to earliest
#First: all birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig1c.png",6,6,res=800,units="in")
par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind_min[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly relative to earliest date",xlab="",xaxt="n",main="")
for(i in 1:nrow(mig_timing_df_stdind_min)){
  temp_mtdsi<-mig_timing_df_stdind_min[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=mig_timing_df_stdind_min$ptt_num[i],
        lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),
          col=mig_timing_df_stdind_min$ptt_num[i],
          lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Just SW birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_16.png",7,7,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind_min[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly relative to earliest date",xlab="",xaxt="n",main="SW")
for(i in which(mig_timing_df_stdind_min$mig_dir=="SW")){
  temp_mtdsi<-mig_timing_df_stdind_min[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=i,lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),col=i,lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Just SE birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_17.png",7,7,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind_min[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly relative to earliest date",xlab="",xaxt="n",main="SE")
abline(h=0,lty=3)
for(i in which(mig_timing_df_stdind_min$mig_dir=="SE")){
  temp_mtdsi<-mig_timing_df_stdind_min[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=i,lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),col=i,lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Colour by the order in which birds left the breeding grounds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_18.png",7,7,res=800,units="in")
mig_timing_df_stdind_min_noNA_depbrgr<-mig_timing_df_stdind_min[which(is.na(mig_timing_df_stdind_min$dep_brgr)==F),]
par(mar=c(7.1, 4.1, 4.1, 2.1))
plot(1:length(milestone_names),rep(0,length(milestone_names)),col="white",
     ylim=range(mig_timing_df_stdind_min_noNA_depbrgr[,milestone_names],na.rm=T),
     ylab="Within-individual anomaly relative to earliest date",xlab="",xaxt="n")
abline(h=0,lty=3)
for(i in 1:nrow(mig_timing_df_stdind_min_noNA_depbrgr)){
  temp_mtdsi<-mig_timing_df_stdind_min_noNA_depbrgr[i,milestone_names]
  if(any(is.na(temp_mtdsi)==F)==F) next
  lines(1:length(milestone_names),as.numeric(temp_mtdsi),
        col=rgb(rank(mig_timing_df_stdind_min_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_min_noNA_depbrgr),
                0,1-(rank(mig_timing_df_stdind_min_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_min_noNA_depbrgr))),
        lwd=2)
  #If there is a gap, draw a dashed line
  temp_na_inds<-which(is.na(temp_mtdsi))
  if(length(temp_na_inds)==0) next
  temp_split<-split(temp_na_inds, cumsum(c(1, diff(temp_na_inds) != 1)))
  for(j in 1:length(temp_split)){
    temp_start<-temp_split[[j]][1]-1
    temp_end<-temp_split[[j]][length(temp_split[[j]])]+1
    if(temp_start<1 | temp_end>10) next
    lines(c(temp_start,temp_end),c(temp_mtdsi[temp_start],temp_mtdsi[temp_end]),
          col=rgb(rank(mig_timing_df_stdind_min_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_min_noNA_depbrgr),
                  0,1-(rank(mig_timing_df_stdind_min_noNA_depbrgr$dep_brgr)[i]/nrow(mig_timing_df_stdind_min_noNA_depbrgr))),
          lwd=2,lty=2)
  }
}
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

########################################################################################################################

#4B: CAUSES OF VARIATION IN MILESTONES
#Migratory direction, breeding habitat, year, age
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig2_2.png",7,5,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
boxplot(mtd_melted$value~mtd_melted$br_hab*mtd_melted$variable,
        at=(1:(3*length(milestone_names)))[-3*(1:length(milestone_names))], 
        col=c("light green","brown"),
        frame=F,
        xlab="",
        ylab="Days since start of 1st calendar year of mig cycle",
        main="Variation in migratory milestone date, by breeding habitat",
        xaxt="n")
axis(1,at=seq(1.5,3*length(milestone_names),3),labels=milestone_names,las=2)
legend("topleft",c("Lowland","Upland"),fill=c("light green","brown"),title="Breeding habitat",bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

mtd_melted$mig_dir<-relevel(as.factor(mtd_melted$mig_dir),ref = "SW") #Put SW first, so it's on the same side as 'lowland'
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig2_3.png",7,5,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
boxplot(mtd_melted$value~mtd_melted$mig_dir*mtd_melted$variable,
        at=(1:(3*length(milestone_names)))[-3*(1:length(milestone_names))], 
        col=c("blue","darkorange1"),
        frame=F,
        xlab="",
        ylab="Days since start of 1st calendar year of mig cycle",
        main="Variation in migratory milestone date, by migratory direction (2 category)",
        xaxt="n")
axis(1,at=seq(1.5,3*length(milestone_names),3),labels=milestone_names,las=2)
legend("topleft",c("SW","SE"),fill=c("blue","darkorange1"),title="Migratory direction",bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

mtd_melted$mig_dir_3cat<-relevel(as.factor(mtd_melted$mig_dir_3cat),ref = "SW") #Put SW first, so it's on the same side as 'lowland'

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig2_4.png",7,5,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
boxplot(mtd_melted$value~mtd_melted$mig_dir_3cat*mtd_melted$variable,
        at=(1:(4*length(milestone_names)))[-4*(1:length(milestone_names))], 
        col=c("blue","green","darkorange1"),
        frame=F,
        xlab="",
        ylab="Days since start of 1st calendar year of mig cycle",
        main="Variation in migratory milestone date, by migratory direction (3 category)",
        xaxt="n")
axis(1,at=seq(2,4*length(milestone_names),4),labels=milestone_names,las=2)
legend("topleft",c("SW","C","SE"),fill=c("blue","green","darkorange1"),title="Migratory direction",bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

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

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig5a.png",6,6,res=800,units="in")
par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names)-0.1,brhab_migdir_res_df$bh_Lowland_est_marg,pch=19,ylab="Julian day",bty="n",xaxt="n",xlab="",col="light green")
points(1:length(milestone_names)+0.1,brhab_migdir_res_df$bh_Upland_est_marg,pch=19,col="brown")
arrows(1:length(milestone_names)-0.1,brhab_migdir_res_df$bh_Lowland_lci_marg,1:length(milestone_names)-0.1,brhab_migdir_res_df$bh_Lowland_uci_marg,length=0,col="light green")
arrows(1:length(milestone_names)+0.1,brhab_migdir_res_df$bh_Upland_lci_marg,1:length(milestone_names)+0.1,brhab_migdir_res_df$bh_Upland_uci_marg,length=0,col="brown")
legend("topleft",legend=c("Lowland","Upland"),pch=c(19,19),bty="n",col=c("light green","brown"))
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig5b.png",6,6,res=800,units="in")
par(mar=c(7.1, 4.1, 1.1, 1.1))
plot(1:length(milestone_names)-0.1,brhab_migdir_res_df$md_SW_est_marg,col="blue",pch=19,ylab="Julian day",bty="n",xaxt="n",xlab="")
points(1:length(milestone_names)+0.1,brhab_migdir_res_df$md_SE_est_marg,pch=19,col="darkorange1")
arrows(1:length(milestone_names)-0.1,brhab_migdir_res_df$md_SW_lci_marg,1:length(milestone_names)-0.1,brhab_migdir_res_df$md_SW_uci_marg,length=0,col="blue")
arrows(1:length(milestone_names)+0.1,brhab_migdir_res_df$md_SE_lci_marg,1:length(milestone_names)+0.1,brhab_migdir_res_df$md_SE_uci_marg,length=0,col="darkorange1")
legend("topleft",legend=c("SW","SE"),pch=c(19,19),bty="n",col=c("blue","darkorange1"))
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

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

write.csv(brhab_migdir_res_df[,c("Milestone","n_events","n_birds",
                                 "Upland","Upland_lci","Upland_uci",
                                 "Lowland","Lowland_lci","Lowland_uci",
                                 "SE","SE_lci","SE_uci",
                                 "SW","SW_lci","SW_uci",
                                 "R2c","R2m")],file="H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/brhab_migdir_res_df.csv",row.names=F)

#How does the distribution of *individual* anomalies for each milestone vary between migratory pathways?
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 090322/fig4_2.png",7,5,res=800,units="in")
par(mar=c(7.1, 4.1, 4.1, 2.1))
boxplot(mtd_stdind_melted$value~mtd_stdind_melted$mig_dir*mtd_stdind_melted$variable,
        at=(1:(3*length(milestone_names)))[-3*(1:length(milestone_names))], 
        col=c("blue","darkorange1"),
        ylab="Within-individual anomaly",
        xlab="",
        xaxt="n",
        frame=F)
axis(1,at=seq(1.5,3*length(milestone_names),3),labels=milestone_names,las=2)
legend("topleft",c("SW","SE"),fill=c("blue","darkorange1"),title="Migratory direction",bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

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

#Write out param estimates for results tables
write.csv(latlong_res_df,"H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/latlong_res_df.csv")

latlong_res_df_range<-range(latlong_res_df[,(names(latlong_res_df) %in% c("R2c","R2m"))==F])
coul3<-brewer.pal(3,"Set1") 

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig6.png",7,7,res=800,units="in")
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
dev.off()

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

#Write out param estimates for SI results tables
write.csv(latlong_res_df_noregion,"H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/latlong_res_df_noregion.csv")

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig7.png",7,7,res=800,units="in")
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
dev.off()


#Age effects
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig5_1.png",10,10,res=800,units="in")
par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(known_age)"))
  boxplot(temp_form_age,data=mig_timing_df_stdind,
          main=milestone_names[i],
          xlab="Cycles since first year",
          ylab="Within-individual anomaly")
}
par(mfrow=c(1,1))
dev.off()

#Cycles since fledging
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig5_3.png",10,10,res=800,units="in")
par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(known_age)"))
  # boxplot(temp_form_age,data=mig_timing_df_stdind,
  #         main=milestone_names[i],
  #         xlab="Cycles since fledging",
  #         ylab="Within-individual anomaly")
  m1<-lm(temp_form_age,data=mig_timing_df_stdind);summary(m1)
  csm1<-data.frame(coefficients(summary(m1)))
  ci_m1<-data.frame(confint(m1))
  temp_ages<-as.numeric(do.call(rbind,strsplit(row.names(csm1),")"))[,2])
  plot(temp_ages,csm1$Estimate,
       main=milestone_names[i],
       col="blue",
       xlab="Cycles since first year",
       xlim=c(1,6),
       ylab="Within-individual anomaly",
       ylim=range(ci_m1))
  arrows(temp_ages,ci_m1$X2.5..,temp_ages,ci_m1$X97.5..,length=0,col="blue")
  abline(h=0,lty=3)
  points(temp_ages[csm1$Pr...t..<0.05],csm1$Estimate[csm1$Pr...t..<0.05],pch=19,col="blue")
}
par(mfrow=c(1,1))
dev.off()

#Check that the lateness of 2cy birds is a juv effect, rather than a year-of-tagging effect
par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){
  temp_form_cycles<-as.formula(paste0(milestone_names[i],"~0+as.factor(cycles_since_tagging_4s)"))
  m1<-lm(temp_form_cycles,data=mig_timing_df_stdind);summary(m1)
  csm1<-data.frame(coefficients(summary(m1)))
  temp_ages<-as.numeric(do.call(rbind,strsplit(row.names(csm1),")"))[,2])
  ci_m1<-data.frame(confint(m1))
  plot(temp_ages,csm1$Estimate,
       main=milestone_names[i],
       col="blue",
       xlab="Cycles since tagging",
       ylab="Within-individual anomaly",
       xlim=c(0,4),
       ylim=range(ci_m1))
  arrows(temp_ages,ci_m1$X2.5..,temp_ages,ci_m1$X97.5..,length=0,col="blue")
  abline(h=0,lty=3)
  points(temp_ages[csm1$Pr...t..<0.05],csm1$Estimate[csm1$Pr...t..<0.05],pch=19,col="blue")
}
par(mfrow=c(1,1))


#Cycle after fledging
par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(known_age_binary)"))
  # boxplot(temp_form_age,data=mig_timing_df_stdind,
  #         main=milestone_names[i],
  #         xlab="Cycle after fledging",
  #         ylab="Within-individual anomaly")
  m1<-lm(temp_form_age,data=mig_timing_df_stdind);summary(m1)
  csm1<-data.frame(coefficients(summary(m1)))
  ci_m1<-data.frame(confint(m1))
  temp_ages<-as.numeric(do.call(rbind,strsplit(row.names(csm1),")"))[,2])
  plot(temp_ages,csm1$Estimate,
       main=milestone_names[i],
       col="blue",
       xlab="Cycles since fledging",
       xlim=c(0,1),
       ylab="Within-individual anomaly",
       ylim=range(ci_m1))
  arrows(temp_ages,ci_m1$X2.5..,temp_ages,ci_m1$X97.5..,length=0,col="blue")
  abline(h=0,lty=3)
  points(temp_ages[csm1$Pr...t..<0.05],csm1$Estimate[csm1$Pr...t..<0.05],pch=19,col="blue")
}
par(mfrow=c(1,1))

#Print p-value for difference
for(i in 1:length(milestone_names)){
  print(milestone_names[i])
  temp_form_age<-as.formula(paste0(milestone_names[i],"~as.factor(known_age_binary)"))
  m1<-lm(temp_form_age,data=mig_timing_df_stdind);summary(m1)
  csm1<-data.frame(coefficients(summary(m1)))
  print(csm1$Pr...t..[2])
}

#Significant difference for arr_wgr_b...is this to do with where young birds overwinter?
m1<-lm(arr_wgr_b~as.factor(known_age_binary)+min_lat_stopover,data=mig_timing_df_stdind);summary(m1) #No
m1<-lm(min_lat_stopover~as.factor(known_age_binary),data=mig_timing_df_stdind);summary(m1) #NO


png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig5_2.png",10,10,res=800,units="in")
par(mfrow=c(4,3))
for(i in 1:length(milestone_names)){
  temp_form_age<-as.formula(paste0(milestone_names[i],"~as.factor(known_age_binary)"))
  boxplot(temp_form_age,data=mig_timing_df_stdind,
          main=milestone_names[i],
          xlab="Cycle after first year",
          ylab="Within-individual anomaly")
}
par(mfrow=c(1,1))
dev.off()

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

write.csv(age_coeff_df_nicer,"H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/age_coeff_df.csv")

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA12.png",7,7,res=800,units="in")
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
arrows(age_sig05_ind-0.3,20,age_sig05_ind+0.3,20,length=0)
arrows(age_sig05_ind-0.3,20,age_sig05_ind-0.3,19,length=0)
arrows(age_sig05_ind+0.3,20,age_sig05_ind+0.3,19,length=0)
text(age_sig05_ind,21,"**",cex=2)
arrows(age_sig1_ind-0.3,20,age_sig1_ind+0.3,20,length=0)
arrows(age_sig1_ind-0.3,20,age_sig1_ind-0.3,19,length=0)
arrows(age_sig1_ind+0.3,20,age_sig1_ind+0.3,19,length=0)
text(age_sig1_ind,21,"*",cex=2)
axis(1, at=1:length(milestone_names),labels=milestone_names,las=2)
legend("topleft",c("Cycle after first breeding season","All subsequent cycles"),col=c("cornflowerblue","chocolate"),pch=c(19,19),bty="n")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

#Cycles since fledging
ord_age_param_est_list<-list()
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA9.png",10,10,res=800,units="in")
par(mfrow=c(3,3))
for(i in 1:length(milestone_names)){
  print(i)
  if(milestone_names[i]=="dep_wgr_a") next #Not enough data for this milestone: trying to estimate an individual random effect and four parameters from 6 data
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(known_age)+(1|ptt)"))
  m1<-brm(temp_form_age,
          silent=1,
          data=mig_timing_df_stdind[which(is.na(mig_timing_df_stdind[,milestone_names[i]])==F),],
          # data=mig_timing_df_stdind[which(mig_timing_df_stdind$known_age<=5 & is.na(mig_timing_df_stdind[,milestone_names[i]])==F),],
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
dev.off()
ord_age_param_est_df<-do.call(rbind,ord_age_param_est_list)

ord_age_param_est_df_nicer<-ord_age_param_est_df
ord_age_param_est_df_nicer$Estimate<-round(ord_age_param_est_df$Estimate,digits=3)
ord_age_param_est_df_nicer$`l-95% CI`<-round(ord_age_param_est_df$`l-95% CI`,digits=3)
ord_age_param_est_df_nicer$`u-95% CI`<-round(ord_age_param_est_df$`u-95% CI`,digits=3)


#Write out param estimates for SI results tables
write.csv(ord_age_param_est_df_nicer,"H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/ord_age_param_est_df.csv")


#Fig A9 but for non-known age birds
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA9_for_unknown_age_birds_.png",10,10,res=800,units="in")
par(mfrow=c(3,3))
for(i in 1:length(milestone_names)){
  print(i)
  if(milestone_names[i]=="dep_wgr_a") next #Not enough data for this milestone: trying to estimate an individual random effect and four parameters from 6 data
  temp_form_age<-as.formula(paste0(milestone_names[i],"~0+as.factor(cycles_since_tagging_4s)+(1|ptt)"))
  m1<-brm(temp_form_age,
          silent=1,
          data=mig_timing_df_stdind[which(is.na(mig_timing_df_stdind[,milestone_names[i]])==F),],
          # data=mig_timing_df_stdind[which(mig_timing_df_stdind$known_age<=5 & is.na(mig_timing_df_stdind[,milestone_names[i]])==F),],
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
dev.off()


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

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA5.png",10,10,res=800,units="in")
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
       # main="Relationship between anomalies\nfor subsequent milestones in same year",
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
     # main="Relationship between anomalies\nfor subsequent milestones in same year",
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
dev.off()

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
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig1_a.png",6,6,res=800,units="in")
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
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig4_8.png",6,6,res=800,units="in")
par(mar=c(10.1, 4.1, 1.1, 1.1))
plot(1:nrow(milestone_name_pairs),milestone_name_pairs$slope,ylim=c(min(milestone_name_pairs$slope_lci),max(milestone_name_pairs$slope_uci)),col="white",xaxt="n",xlab="",pch=19,ylab="Slope coefficient",main="")
points(which(milestone_name_pairs$sign==F),milestone_name_pairs$slope[which(milestone_name_pairs$sign==F)])
points(which(milestone_name_pairs$sign),milestone_name_pairs$slope[which(milestone_name_pairs$sign)],pch=19)
arrows(1:nrow(milestone_name_pairs),milestone_name_pairs$slope_lci,1:nrow(milestone_name_pairs),milestone_name_pairs$slope_uci,length=0)
abline(h=0,lty=3)
axis(1, at=1:nrow(milestone_name_pairs),labels=paste0(milestone_name_pairs$from," : ",milestone_name_pairs$to),las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


#Relationship between anomalies and arr_brgr anomaly
milestone_df_arrbrgr<-data.frame("from"=milestone_names[milestone_names!="arr_brgr"])
milestone_df_arrbrgr$R2m<-milestone_df_arrbrgr$slope<-milestone_df_arrbrgr$slope_lci<-milestone_df_arrbrgr$slope_uci<-milestone_df_arrbrgr$sign<-NA

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA6.png",10,10,res=800,units="in")
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
dev.off()

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/fig4_12.png",6,6,res=800,units="in")
par(mar=c(10.1, 4.1, 1.1, 1.1))
plot(1:nrow(milestone_df_arrbrgr),milestone_df_arrbrgr$slope,ylim=c(min(milestone_df_arrbrgr$slope_lci),max(milestone_df_arrbrgr$slope_uci)),col="white",xaxt="n",xlab="",pch=19,ylab="Slope coefficient",main="")
points(which(milestone_df_arrbrgr$sign==F),milestone_df_arrbrgr$slope[which(milestone_df_arrbrgr$sign==F)])
points(which(milestone_df_arrbrgr$sign),milestone_df_arrbrgr$slope[which(milestone_df_arrbrgr$sign)],pch=19)
arrows(1:nrow(milestone_df_arrbrgr),milestone_df_arrbrgr$slope_lci,1:nrow(milestone_df_arrbrgr),milestone_df_arrbrgr$slope_uci,length=0)
abline(h=0,lty=3)
axis(1, at=1:nrow(milestone_df_arrbrgr),labels=milestone_df_arrbrgr$from,las=2)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

ms_table<-milestone_name_pairs[,c("from","var","to","delta_var","var_sig_diff","R2m","sign")]
ms_table$var<-round(milestone_name_pairs$var,digits=2)
ms_table$delta_var<-round(milestone_name_pairs$delta_var,digits=2)
ms_table$R2m<-round(milestone_name_pairs$R2m,digits=2)
ms_table$R2m_arrbrgr_sign<-ms_table$R2m_arrbrgr<-NA
for(i in 1:nrow(milestone_df_arrbrgr)){
  ms_table[ms_table$from==milestone_df_arrbrgr[i,]$from,]$R2m_arrbrgr<-round(milestone_df_arrbrgr[i,]$R2m,digits=2)
  ms_table[ms_table$from==milestone_df_arrbrgr[i,]$from,]$R2m_arrbrgr_sign<-milestone_df_arrbrgr[i,]$sign
}
write.csv(ms_table,"H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/ms_table.csv",row.names=F)


#Make pretty table for paper
sd_r2_tab_within_ind<-milestone_name_pairs[,c("from","sd","to","delta_sd","R2m")]
sd_r2_tab_within_ind$R2m_arrbrgr<-NA
for(i in 1:nrow(milestone_df_arrbrgr)){
  sd_r2_tab_within_ind[sd_r2_tab_within_ind$from==milestone_df_arrbrgr[i,]$from,]$R2m_arrbrgr<-milestone_df_arrbrgr[i,]$R2m
}
sd_r2_tab_within_ind[,c("sd","delta_sd","R2m","R2m_arrbrgr")]<-round(sd_r2_tab_within_ind[,c("sd","delta_sd","R2m","R2m_arrbrgr")],digits=3)
write.csv(sd_r2_tab_within_ind,"H:/JGD H drive backup 160320/Cuckoos/Plots/Results 250322/sd_r2_tab_within_ind.csv")



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

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/figA3.png",10,10,res=800,units="in")
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
dev.off()

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

#Fix first (doesn't make sense to do this, as last fix can be anywhere)
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


#Were birds that died 

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

m1<-glm(died_bin~value,family="binomial",data=mtd_stdind_melted);summary(m1) #No overall effect of anomaly on death probability
m0<-glm(died_bin~1,family="binomial",data=mtd_stdind_melted);summary(m0)
m1<-glm(died_bin~value*variable,family="binomial",data=mtd_stdind_melted);summary(m1) #No effect of anomaly on death probability for individual milestones
#SO: apart from a non-significant hint that birds leaving the UK late, or departing the wintering grounds early (in fact all birds that died after leaving the wintering grounds left *early*), have a higher probability of death...
# ...there is no reason to suspect that the effect of migratory timing on fitness is passed through survival

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
#SW birds die more on way out, SE birds die more on way back 
#Do we know this already - check Hewson et al 2016


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
                                  # (mtd_melted_scaled$ptt %in% c("134951","128303")==F) & 
                                  (mtd_melted_scaled$outcome %in% c("UB","UC"))==F ,]);summary(m1);confint(m1)
m2<-lmer(value~died+(1|year)+(1|region/ptt),
         data=mtd_melted_scaled[mtd_melted_scaled$season=="Spring" & 
                                  # (mtd_melted_scaled$ptt %in% c("134951","128303")==F) & 
                                  (mtd_melted_scaled$outcome %in% c("UA","UB","UC"))==F ,]);summary(m2);confint(m2)


##############################################################################################
#Map stopovers by migration stage
#Categorise stopovers by migratory stage
stopovers_tab
stopovers_tab$SO_mean<-NA
for(i in 1:nrow(stopovers_tab)){stopovers_tab$SO_mean[i]<-as.character(mean(c(stopovers_tab$SO_start[i],stopovers_tab$SO_end_LOS[i])))}
stopovers_tab$mig_stage<-NA

#NB this currently leaves out 4 stopovers for birds whose last milestone was dep_wa_a, but for which there were no fixes in WA
stopovers_tab_with_mig_stage_list<-list()
for(i in 1:nrow(mig_timing_df)){
  temp_SO_tab<-stopovers_tab[stopovers_tab$ptt==mig_timing_df[i,]$ptt,]
  if(i==which(mig_timing_df$ptt==mig_timing_df[i,]$ptt)[1]){ #i.e. for first year only
    pre_dep_brgr_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$first_fix_TS & 
                              temp_SO_tab$SO_mean<mig_timing_df[i,]$dep_brgr_TS)
    if(length(pre_dep_brgr_ind)>0){temp_SO_tab[pre_dep_brgr_ind,]$mig_stage<-"pre_dep_brgr"}
  }
  dep_brgr_2_dep_UK_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_brgr_TS & 
                                 temp_SO_tab$SO_mean<mig_timing_df[i,]$dep_UK_TS)
  if(length(dep_brgr_2_dep_UK_ind)>0){temp_SO_tab[dep_brgr_2_dep_UK_ind,]$mig_stage<-"dep_brgr_2_dep_UK"}
  dep_UK_2_fin_sah_sb_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_UK_TS & 
                                   temp_SO_tab$SO_mean<mig_timing_df[i,]$fin_sah_sb_TS)
  if(length(dep_UK_2_fin_sah_sb_ind)>0){temp_SO_tab[dep_UK_2_fin_sah_sb_ind,]$mig_stage<-"dep_UK_2_fin_sah_sb"}
  fin_sah_sb_2_arr_wgr_a_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$fin_sah_sb_TS & 
                                      temp_SO_tab$SO_mean<mig_timing_df[i,]$arr_wgr_a_TS)
  if(length(fin_sah_sb_2_arr_wgr_a_ind)>0){temp_SO_tab[fin_sah_sb_2_arr_wgr_a_ind,]$mig_stage<-"fin_sah_sb_2_arr_wgr_a"}
  arr_wgr_a_2_arr_wgr_b_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_wgr_a_TS & 
                                     temp_SO_tab$id<mig_timing_df[i,]$min_lat_stopover_id)
  if(length(arr_wgr_a_2_arr_wgr_b_ind)>0){temp_SO_tab[arr_wgr_a_2_arr_wgr_b_ind,]$mig_stage<-"arr_wgr_a_2_arr_wgr_b"}
  arr_wgr_b_2_dep_wgr_b_ind<-which(temp_SO_tab$id==mig_timing_df[i,]$min_lat_stopover_id)
  if(length(arr_wgr_b_2_dep_wgr_b_ind)>0){temp_SO_tab[arr_wgr_b_2_dep_wgr_b_ind,]$mig_stage<-"arr_wgr_b_2_dep_wgr_b"}
  dep_wgr_b_2_dep_wgr_a_ind<-which(temp_SO_tab$id>mig_timing_df[i,]$min_lat_stopover_id & 
                                     temp_SO_tab$SO_mean<mig_timing_df[i,]$dep_wgr_a_TS)
  if(length(dep_wgr_b_2_dep_wgr_a_ind)>0){temp_SO_tab[dep_wgr_b_2_dep_wgr_a_ind,]$mig_stage<-"dep_wgr_b_2_dep_wgr_a"}
  dep_wgr_a_2_dep_wa_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_wgr_a_TS & 
                                  temp_SO_tab$SO_mean<mig_timing_df[i,]$dep_wa_TS)
  if(length(dep_wgr_a_2_dep_wa_ind)>0){temp_SO_tab[dep_wgr_a_2_dep_wa_ind,]$mig_stage<-"dep_wgr_a_2_dep_wa"}
  dep_wa_2_arr_UK_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_wa_TS & 
                               temp_SO_tab$SO_mean<mig_timing_df[i,]$arr_UK_TS)
  if(length(dep_wa_2_arr_UK_ind)>0){temp_SO_tab[dep_wa_2_arr_UK_ind,]$mig_stage<-"dep_wa_2_arr_UK"}
  arr_UK_2_arr_brgr_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_UK_TS & 
                                 temp_SO_tab$SO_mean<mig_timing_df[i,]$arr_brgr_TS)
  if(length(arr_UK_2_arr_brgr_ind)>0){temp_SO_tab[arr_UK_2_arr_brgr_ind,]$mig_stage<-"arr_UK_2_arr_brgr"}
  if(i!=nrow(mig_timing_df)){ #In which case there would be no mig_timing_df$ptt[i+1]
    if(mig_timing_df$ptt[i]==mig_timing_df$ptt[i+1]){
      arr_brgr_2_dep_brgr_ind<-which(temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_brgr_TS & 
                                       temp_SO_tab$SO_mean<mig_timing_df[i+1,]$dep_brgr_TS)
      if(length(arr_brgr_2_dep_brgr_ind)>0){temp_SO_tab[arr_brgr_2_dep_brgr_ind,]$mig_stage<-"arr_brgr_2_dep_brgr"}
    }
  }
  
  #Add clause for died after milestone
  lmpd<-mig_timing_df$last_milestone_pre_death[i]
  if(mig_timing_df$ptt[i]=="146757"){lmpd<-"dep_wa_a"} #Because it had been manually adjusted to arr_UK for Vigilamus, which wouldn't work with code below
  
  if(is.na(lmpd)==F){
    if(lmpd=="dep_UK" & 
       any(temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_UK_TS)){
      temp_SO_tab[temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_UK_TS,]$mig_stage<-"dep_UK_2_fin_sah_sb"
    }
    if(lmpd=="fin_sah_sb" & 
       any(temp_SO_tab$SO_mean>mig_timing_df[i,]$fin_sah_sb_TS)){
      temp_SO_tab[temp_SO_tab$SO_mean>mig_timing_df[i,]$fin_sah_sb_TS,]$mig_stage<-"fin_sah_sb_2_arr_wgr_a"
    }
    if(lmpd=="arr_wgr_a" & 
       any(temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_wgr_a_TS)){
      temp_SO_tab[temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_wgr_a_TS,]$mig_stage<-"arr_wgr_a_2_arr_wgr_b"
    }
    if(lmpd=="dep_wgr_b" & 
       any(temp_SO_tab$id>mig_timing_df[i,]$min_lat_stopover_id)){
      temp_SO_tab[temp_SO_tab$id>mig_timing_df[i,]$min_lat_stopover_id,]$mig_stage<-"dep_wgr_b_2_dep_wa"
    }
    if(lmpd=="dep_wgr_a" & 
       any(temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_wgr_a_TS)){
      temp_SO_tab[temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_wgr_a_TS,]$mig_stage<-"dep_wgr_a_2_dep_wa"
    }
    if(lmpd=="dep_wa_a" & 
       any(temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_wa_a_TS)){
      temp_SO_tab[temp_SO_tab$SO_mean>mig_timing_df[i,]$dep_wa_a_TS,]$mig_stage<-"dep_wa_2_arr_UK"
    }
    if(lmpd=="arr_UK" & 
       any(temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_UK_TS)){
      temp_SO_tab[temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_UK_TS,]$mig_stage<-"arr_UK_2_arr_brgr"
    }
    if(lmpd=="arr_brgr" & 
       any(temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_brgr_TS)){
      temp_SO_tab[temp_SO_tab$SO_mean>mig_timing_df[i,]$arr_brgr_TS,]$mig_stage<-"arr_brgr_2_dep_brgr"
    }
  }
  
  stopovers_tab_with_mig_stage_list[[i]]<-temp_SO_tab
}
stopovers_tab_with_mig_stage<-do.call(rbind,stopovers_tab_with_mig_stage_list)
stopovers_tab_with_mig_stage<-stopovers_tab_with_mig_stage[is.na(stopovers_tab_with_mig_stage$mig_stage)==F,]

stopovers_tab_with_mig_stage<-stopovers_tab_with_mig_stage[stopovers_tab_with_mig_stage$mig_stage!="pre_dep_brgr",]

#Order migratory stages for plotting - NB this misses out pre_dep_brgr
un_mig_stage_ordered<-c("dep_brgr_2_dep_UK",
                        "dep_UK_2_fin_sah_sb",
                        "fin_sah_sb_2_arr_wgr_a",
                        "arr_wgr_a_2_arr_wgr_b",
                        "arr_wgr_b_2_dep_wgr_b",
                        "dep_wgr_b_2_dep_wgr_a",
                        "dep_wgr_a_2_dep_wa",
                        "dep_wa_2_arr_UK",
                        "arr_UK_2_arr_brgr",
                        "arr_brgr_2_dep_brgr")

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig_A1.png",9,9,res=800,units="in")
layout.matrix <- matrix(c(1,2,3,4,
                          10,0,0,5,
                          9,8,7,6),byrow=T, nrow = 3, ncol = 4)
layout(layout.matrix)
for(i in 1:length(un_mig_stage_ordered)){
  temp_ss<-strsplit(un_mig_stage_ordered[i],"_2_")
  plot(0,0,
       xlim=range(stopovers_tab_with_mig_stage$SO_median_long),
       ylim=range(stopovers_tab_with_mig_stage$SO_median_lat),
       col="white",
       xlab="Longitude",
       ylab="Latitude",
       main=paste0(temp_ss[[1]][1]," to ",temp_ss[[1]][2]))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightblue1")
  plot(world_shp,add=T,col="light green",border="gray50")
  plot(sah,col="orange",border=NA,add=T)
  plot(ca_rf,col="dark green",border=NA,add=T)
  points(SO_median_lat~SO_median_long,data=stopovers_tab_with_mig_stage[stopovers_tab_with_mig_stage$mig_stage==un_mig_stage_ordered[i],],pch=20)
}
dev.off()
par(mfrow=c(1,1)) #To negate unusual layout

#Equivalent histograms
stopovers_tab_with_mig_stage$jul<-yday(stopovers_tab_with_mig_stage$SO_mean)

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig_A13.png",9,9,res=800,units="in")
layout.matrix <- matrix(c(1,2,3,4,
                          10,0,0,5,
                          9,8,7,6),byrow=T, nrow = 3, ncol = 4)
layout(layout.matrix)
for(i in 1:length(un_mig_stage_ordered)){
  temp_ss<-strsplit(un_mig_stage_ordered[i],"_2_")
  hist(stopovers_tab_with_mig_stage[stopovers_tab_with_mig_stage$mig_stage==un_mig_stage_ordered[i],]$jul,breaks=seq(0,370,by=10),
       xlab="Julian day",
       freq=F,
       main=paste0(temp_ss[[1]][1]," to ",temp_ss[[1]][2]))
}
dev.off()
par(mfrow=c(1,1)) #To negate unusual layout


#Both plots together
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig_Ax1.png",18,12,res=800,units="in")
layout.matrix <- matrix(c(1,1,4,4,7,7,10,10,
                          2,3,5,6,8,9,11,12,
                          28,28,0,0,0,0,13,13,
                          29,30,0,0,0,0,14,15,
                          25,25,22,22,19,19,16,16,
                          26,27,23,24,20,21,17,18),
                        byrow=T, nrow = 6, ncol = 8)
layout(layout.matrix,heights=c(1,3,1,3,1,3))
for(i in 1:length(un_mig_stage_ordered)){
  print(i)
  temp_ss<-strsplit(un_mig_stage_ordered[i],"_2_")
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="")
  text(0,0,paste0(temp_ss[[1]][1]," to ",temp_ss[[1]][2]),cex=2,font=2)
  par(mar=c(5.1, 3.1, 0.1, 1.1))
  plot(0,0,
       xlim=range(stopovers_tab_with_mig_stage$SO_median_long),
       ylim=range(stopovers_tab_with_mig_stage$SO_median_lat),
       col="white",
       xlab="Longitude",
       ylab="Latitude",
       main="")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightblue1")
  plot(world_shp,add=T,col="light green",border="gray50")
  plot(sah,col="orange",border=NA,add=T)
  plot(ca_rf,col="dark green",border=NA,add=T)
  points(SO_median_lat~SO_median_long,data=stopovers_tab_with_mig_stage[stopovers_tab_with_mig_stage$mig_stage==un_mig_stage_ordered[i],],pch=20)
  hist(stopovers_tab_with_mig_stage[stopovers_tab_with_mig_stage$mig_stage==un_mig_stage_ordered[i],]$jul,breaks=seq(0,370,by=10),
       xlab="Julian day",
       freq=F,
       main="")
}
dev.off()
par(mfrow=c(1,1)) #To negate unusual layout

##############################################################################################
#Map fixes by migration stage
#Categorise fixes by migratory stage
d1_spdf_2<-data.frame(d1_spdf)
d1_spdf_2$mig_stage<-NA
d1_spdf_2<-d1_spdf_2[which((d1_spdf_2$event.id==3819970477)==F),] #Point over water

d1_spdf_2_with_mig_stage_list<-list()
for(i in 1:nrow(mig_timing_df)){
  temp_d1_spdf_2<-d1_spdf_2[d1_spdf_2$ptt==mig_timing_df[i,]$ptt,]
  if(i==which(mig_timing_df$ptt==mig_timing_df[i,]$ptt)[1]){ #i.e. for first year only
    pre_dep_brgr_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$first_fix_TS & 
                              temp_d1_spdf_2$timestamp<mig_timing_df[i,]$dep_brgr_TS)
    if(length(pre_dep_brgr_ind)>0){temp_d1_spdf_2$mig_stage[pre_dep_brgr_ind]<-"pre_dep_brgr"}
  }
  dep_brgr_2_dep_UK_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_brgr_TS & 
                                 temp_d1_spdf_2$timestamp<mig_timing_df[i,]$dep_UK_TS)
  if(length(dep_brgr_2_dep_UK_ind)>0){temp_d1_spdf_2$mig_stage[dep_brgr_2_dep_UK_ind]<-"dep_brgr_2_dep_UK"}
  dep_UK_2_fin_sah_sb_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_UK_TS & 
                                   temp_d1_spdf_2$timestamp<mig_timing_df[i,]$fin_sah_sb_TS)
  if(length(dep_UK_2_fin_sah_sb_ind)>0){temp_d1_spdf_2$mig_stage[dep_UK_2_fin_sah_sb_ind]<-"dep_UK_2_fin_sah_sb"}
  fin_sah_sb_2_arr_wgr_a_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$fin_sah_sb_TS & 
                                      temp_d1_spdf_2$timestamp<mig_timing_df[i,]$arr_wgr_a_TS)
  if(length(fin_sah_sb_2_arr_wgr_a_ind)>0){temp_d1_spdf_2$mig_stage[fin_sah_sb_2_arr_wgr_a_ind]<-"fin_sah_sb_2_arr_wgr_a"}
  arr_wgr_a_2_arr_wgr_b_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_wgr_a_TS & 
                                     temp_d1_spdf_2$timestamp<mig_timing_df[i,]$arr_wgr_b_TS)
  if(length(arr_wgr_a_2_arr_wgr_b_ind)>0){temp_d1_spdf_2$mig_stage[arr_wgr_a_2_arr_wgr_b_ind]<-"arr_wgr_a_2_arr_wgr_b"}
  arr_wgr_b_2_dep_wgr_b_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_wgr_b_TS & 
                                     temp_d1_spdf_2$timestamp<mig_timing_df[i,]$dep_wgr_b_TS)
  if(length(arr_wgr_b_2_dep_wgr_b_ind)>0){temp_d1_spdf_2$mig_stage[arr_wgr_b_2_dep_wgr_b_ind]<-"arr_wgr_b_2_dep_wgr_b"}
  dep_wgr_b_2_dep_wgr_a_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wgr_b_TS & 
                                     temp_d1_spdf_2$timestamp<mig_timing_df[i,]$dep_wgr_a_TS_ford1plot)
  if(length(dep_wgr_b_2_dep_wgr_a_ind)>0){temp_d1_spdf_2[dep_wgr_b_2_dep_wgr_a_ind,]$mig_stage<-"dep_wgr_b_2_dep_wgr_a"}
  dep_wgr_a_2_dep_wa_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wgr_a_TS_ford1plot & 
                                  temp_d1_spdf_2$timestamp<mig_timing_df[i,]$dep_wa_TS)
  if(length(dep_wgr_a_2_dep_wa_ind)>0){temp_d1_spdf_2[dep_wgr_a_2_dep_wa_ind,]$mig_stage<-"dep_wgr_a_2_dep_wa"}
  dep_wa_2_arr_UK_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wa_TS & 
                               temp_d1_spdf_2$timestamp<mig_timing_df[i,]$arr_UK_TS)
  if(length(dep_wa_2_arr_UK_ind)>0){temp_d1_spdf_2[dep_wa_2_arr_UK_ind,]$mig_stage<-"dep_wa_2_arr_UK"}
  arr_UK_2_arr_brgr_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_UK_TS & 
                                 temp_d1_spdf_2$timestamp<mig_timing_df[i,]$arr_brgr_TS)
  if(length(arr_UK_2_arr_brgr_ind)>0){temp_d1_spdf_2[arr_UK_2_arr_brgr_ind,]$mig_stage<-"arr_UK_2_arr_brgr"}
  if(i!=nrow(mig_timing_df)){ #In which case there would be no mig_timing_df$ptt[i+1]
    if(mig_timing_df$ptt[i]==mig_timing_df$ptt[i+1]){
      arr_brgr_2_dep_brgr_ind<-which(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_brgr_TS & 
                                       temp_d1_spdf_2$timestamp<mig_timing_df[i+1,]$dep_brgr_TS)
      if(length(arr_brgr_2_dep_brgr_ind)>0){temp_d1_spdf_2[arr_brgr_2_dep_brgr_ind,]$mig_stage<-"arr_brgr_2_dep_brgr"}
    }
  }
  
  #Add clause for died after milestone
  lmpd<-mig_timing_df$last_milestone_pre_death[i]
  if(mig_timing_df$ptt[i]=="146757"){lmpd<-"dep_wa_a"} #Because it had been manually adjusted to arr_UK for Vigilamus, which wouldn't work with code below
  
  if(is.na(lmpd)==F){
    if(lmpd=="dep_UK" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_UK_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_UK_TS,]$mig_stage<-"dep_UK_2_fin_sah_sb"
    }
    if(lmpd=="fin_sah_sb" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$fin_sah_sb_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$fin_sah_sb_TS,]$mig_stage<-"fin_sah_sb_2_arr_wgr_a"
    }
    if(lmpd=="arr_wgr_a" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_wgr_a_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_wgr_a_TS,]$mig_stage<-"arr_wgr_a_2_arr_wgr_b"
    }
    if(lmpd=="dep_wgr_b" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wgr_b_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wgr_b_TS,]$mig_stage<-"dep_wgr_b_2_dep_wgr_a"
    }
    if(lmpd=="dep_wgr_a" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wgr_a_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wgr_a_TS,]$mig_stage<-"dep_wgr_a_2_dep_wa"
    }
    if(lmpd=="dep_wa_a" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wa_a_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$dep_wa_a_TS,]$mig_stage<-"dep_wa_2_arr_UK"
    }
    if(lmpd=="arr_UK" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_UK_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_UK_TS,]$mig_stage<-"arr_UK_2_arr_brgr"
    }
    if(lmpd=="arr_brgr" & 
       any(temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_brgr_TS)){
      temp_d1_spdf_2[temp_d1_spdf_2$timestamp>mig_timing_df[i,]$arr_brgr_TS,]$mig_stage<-"arr_brgr_2_dep_brgr"
    }
  }
  
  d1_spdf_2_with_mig_stage_list[[i]]<-temp_d1_spdf_2
}
d1_spdf_2_with_mig_stage<-do.call(rbind,d1_spdf_2_with_mig_stage_list)
d1_spdf_2_with_mig_stage<-d1_spdf_2_with_mig_stage[is.na(d1_spdf_2_with_mig_stage$mig_stage)==F,]
d1_spdf_2_with_mig_stage<-d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage!="pre_dep_brgr",]

#Order migratory stages for plotting - NB this misses out pre_dep_brgr
un_mig_stage_ordered<-c("dep_brgr_2_dep_UK",
                        "dep_UK_2_fin_sah_sb",
                        "fin_sah_sb_2_arr_wgr_a",
                        "arr_wgr_a_2_arr_wgr_b",
                        "arr_wgr_b_2_dep_wgr_b",
                        "dep_wgr_b_2_dep_wgr_a",
                        "dep_wgr_a_2_dep_wa",
                        "dep_wa_2_arr_UK",
                        "arr_UK_2_arr_brgr",
                        "arr_brgr_2_dep_brgr")

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig_A1_fixes.png",9,9,res=800,units="in")
layout.matrix <- matrix(c(1,2,3,4,
                          10,0,0,5,
                          9,8,7,6),byrow=T, nrow = 3, ncol = 4)
layout(layout.matrix)
for(i in 1:length(un_mig_stage_ordered)){
  temp_ss<-strsplit(un_mig_stage_ordered[i],"_2_")
  plot(0,0,
       xlim=range(d1_spdf_2_with_mig_stage$location.long),
       ylim=range(d1_spdf_2_with_mig_stage$location.lat),
       col="white",
       xlab="Longitude",
       ylab="Latitude",
       main=paste0(temp_ss[[1]][1]," to ",temp_ss[[1]][2]))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightblue1")
  plot(world_shp,add=T,col="light green",border="gray50")
  plot(sah,col="orange",border=NA,add=T)
  plot(ca_rf,col="dark green",border=NA,add=T)
  points(location.lat~location.long,data=d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage==un_mig_stage_ordered[i],],pch=20)
}
dev.off()
par(mfrow=c(1,1)) #To negate unusual layout

#Equivalent histograms
d1_spdf_2_with_mig_stage$jul<-yday(d1_spdf_2_with_mig_stage$timestamp)

png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig_A13_fixes.png",9,9,res=800,units="in")
layout.matrix <- matrix(c(1,2,3,4,
                          10,0,0,5,
                          9,8,7,6),byrow=T, nrow = 3, ncol = 4)
layout(layout.matrix)
for(i in 1:length(un_mig_stage_ordered)){
  temp_ss<-strsplit(un_mig_stage_ordered[i],"_2_")
  hist(d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage==un_mig_stage_ordered[i],]$jul,breaks=seq(0,370,by=10),
       xlab="Julian day",
       freq=F,
       main=paste0(temp_ss[[1]][1]," to ",temp_ss[[1]][2]))
}
dev.off()
par(mfrow=c(1,1)) #To negate unusual layout


#Both plots together
png("H:/JGD H drive backup 160320/Cuckoos/Plots/Results 230522/fig_Ax1_fixes.png",18,12,res=800,units="in")
# layout.matrix <- matrix(c(1,1,4,4,7,7,10,10,
#                           2,3,5,6,8,9,11,12,
#                           28,28,0,0,0,0,13,13,
#                           29,30,0,0,0,0,14,15,
#                           25,25,22,22,19,19,16,16,
#                           26,27,23,24,20,21,17,18),
#                         byrow=T, nrow = 6, ncol = 8)
# layout.matrix <- matrix(c(1,1,4,4,7,7,10,10,
#                           2,3,5,6,8,9,11,12,
#                           28,28,0,31,33,0,13,13,
#                           29,30,0,32,34,0,14,15,
#                           25,25,22,22,19,19,16,16,
#                           26,27,23,24,20,21,17,18),
#                         byrow=T, nrow = 6, ncol = 8)
layout.matrix <- matrix(c(25,25,28,28,1,1,4,4,
                          26,27,29,30,2,3,5,6,
                          22,22,0,31,33,0,7,7,
                          23,24,0,32,34,0,8,9,
                          19,19,16,16,13,13,10,10,
                          20,21,17,18,14,15,11,12),
                        byrow=T, nrow = 6, ncol = 8)
layout(layout.matrix,heights=c(1,3,1,3,1,3))
for(i in 1:length(un_mig_stage_ordered)){
  print(i)
  temp_ss<-strsplit(un_mig_stage_ordered[i],"_2_")
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="")
  text(0,0,paste0(temp_ss[[1]][1]," to ",temp_ss[[1]][2]),cex=2,font=2)
  par(mar=c(5.1, 3.1, 0.1, 1.1))
  plot(0,0,
       xlim=range(d1_spdf_2_with_mig_stage$location.long),
       ylim=range(d1_spdf_2_with_mig_stage$location.lat),
       col="white",
       xlab="Longitude",
       ylab="Latitude",
       main=paste0(temp_ss[[1]][1]," to ",temp_ss[[1]][2]))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col="lightblue1")
  plot(world_shp,add=T,col="light green",border="gray50")
  plot(sah,col="orange",border=NA,add=T)
  plot(ca_rf,col="dark green",border=NA,add=T)
  points(location.lat~location.long,data=d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage==un_mig_stage_ordered[i],],pch=20)
  hist(d1_spdf_2_with_mig_stage[d1_spdf_2_with_mig_stage$mig_stage==un_mig_stage_ordered[i],]$jul,breaks=seq(0,370,by=10),
       xlab="Julian day",
       freq=F,
       main="")
}

#Northbound migration map
par(mar=c(0.1, 0.1, 0.1, 0.1))
plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="")
text(0,0,paste0("Northbound"),cex=2,font=2)
par(mar=c(5.1, 3.1, 0.1, 1.1))

plot(0,0,
     xlim=range(d1_spdf_2_with_mig_stage$location.long),
     ylim=range(d1_spdf_2_with_mig_stage$location.lat),
     col="white",
     xlab="Longitude",
     ylab="Latitude",
     main="Northbound migration")
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
legend("bottomleft",c("SW","SE"),lwd=c(2,2),col=c("yellow","red"))

#Southbound migration map
par(mar=c(0.1, 0.1, 0.1, 0.1))
plot(0,0,xaxt="n",yaxt="n",col="white",bty="n",xlab="",ylab="")
text(0,0,paste0("Southbound"),cex=2,font=2)
par(mar=c(5.1, 3.1, 0.1, 1.1))

plot(0,0,
     xlim=range(d1_spdf_2_with_mig_stage$location.long),
     ylim=range(d1_spdf_2_with_mig_stage$location.lat),
     col="white",
     xlab="Longitude",
     ylab="Latitude",
     main="Southbound migration")
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
legend("bottomleft",c("SW","SE"),lwd=c(2,2),col=c("yellow","red"))

dev.off()
par(mfrow=c(1,1)) #To negate unusual layout

