###################################################################
### R script to infer feeding interactions from species traits   ##
### and obtain the results presented in Pecuchet et al.2020 GCB  ##
### 26/06/20 by Laurene Pecuchet, for questions: laurene.pecuchet@uit.no


## libraries
library(reshape)
library(reshape2)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(igraph)
library(NetIndices)
library(gbm)

# Function
CalcTSS=function(inter) # to calculate True Statistical Skill
{ aW=which(inter$link==1 & inter$pred1 ==1)
bW=which(inter$link==0 & inter$pred1 ==1) 
cW=which(inter$link==1 & inter$pred1 ==0) 
dW=which(inter$link==0 & inter$pred1 ==0) 
a=as.numeric(length(aW))
b=as.numeric(length(bW))
c=as.numeric(length(cW))
d=as.numeric(length(dW))
Py=a/(a+b)*100
Pn=d/(c+d)*100
TSS=(a*d-b*c)/((a+c)*(b+d))
TSSd=data.frame(a=a,b=b,c=c,d=d,TSS=TSS)
return(TSSd)}

#### Datasets ####
load("Datasets_gcb.15196.RData")

#sp_poly<-sp_poly[,colnames(sp_poly)%in%traits$Scientific_name_accepted]

#save(sp_poly, traits, fw, file="FW_Poly_Traits_061020.RData")

### Metaweb
head(fw)

# Pairwise list of all species interactions (minus pairwise species that do not co-occur in the Barents Sea)
head(links)

### MetaTraits
head(traits)

## Species names, taxonomic info, and convertion between foodweb sp. names <-> trait sp. names
head(f_names)

### List of species present in Arctic / Boreal / only (typical) Arctic / only (typical) Boreal
head(sp_list)
sp_boreal<-sp_list$sp_boreal # Species documented in boreal region
sp_arctic<-sp_list$sp_arctic[1:165] # Species documented in Arctic region
typarctic<-sp_list$typical_Arctic[1:13] # Arctic species not present in boreal region
typboreal<-sp_list$typical_boreal[1:70] # boreal species not present in Arctic region
new_sp<-sp_list$sp_new[1:11] # boreal species observed in Arctic in 2014-2017, while not observed there in 2004-2007


#### Boosted Regression Trees ####
## Merge interactions and predator traits / prey traits
# Prey traits
links<-merge(links, traits, by.x = "prey", by.y = "Scientific_name_accepted", all.x=T)
# Predator traits
links<-merge(links, traits, by.x = "predator", by.y = "Scientific_name_accepted", all.x=T)
head(links)#  y predator  /  x prey

## Input dataset for BRT analysis ##
links_all<-links 
## Traits predicting Arctic food web structure:  (get rid of taxonomic information and feeding traits of preys)
links_all<-links_all[,!(names(links_all) %in% c("AphiaID.x","Group.x","Species.x", "Class.x","Order.x","Family.x","Genus.x",
                                                "Body_defense.x","Body_defense.y", ## use continuous variable instead ("defence_cont")
                                                "AphiaID.y","Group.y","Species.y", "Class.y","Order.y","Family.y","Genus.y",
                                                "defence_cont.y","Kingdom.x","Phylum.x",  "Kingdom.y" ,"Phylum.y",
                                                "Living.x", "Non_living.x","X1mm.x","X0.1_1cm.x" ,"X1_10cm.x","X10cm.x" ,"Plant.x","Animal.x",
                                                "Suspension_feeder.x", "Deposit_feeder.x",   "Predator.x","occ"))]
links_all$Type_cat.x<-as.factor(links_all$Type_cat.x)
links_all$Type_cat.y<-as.factor(links_all$Type_cat.y)
summary(links_all)
table(links_all$link)

links_all$link<-as.numeric(as.character(links_all$link)) # For BRT: use numerical response, but bernouilli distribution

links_fw<-links_all[,-c(1,2)] # Eliminate species name
summary(links_fw)

### Rename variables
colnames(links_fw)<-gsub('.x', '.R', colnames(links_fw)) # Resource
colnames(links_fw)<-gsub('.y', '.C', colnames(links_fw)) # Consumer
names(links_fw) <- make.names(names(links_fw))

positive_link<-subset(links_fw, link==1) ## All documented links between species in the Barents Sea metaweb

## Dataset for Prediction ##
links_pred<-reshape2::melt(as.matrix(fw)) # Different dataset that input dataset, because  for prediction 
                             # purposes, we include species that do not co-occur in the Barents Sea (while excluded in modelling dataset)
colnames(links_pred)<-c("prey","predator","link")
summary(links_pred)

# Prey traits
links_pred<-merge(links_pred, traits, by.x = "prey", by.y = "Scientific_name_accepted", all.x=T)
# Predator traits
links_pred<-merge(links_pred, traits, by.x = "predator", by.y = "Scientific_name_accepted", all.x=T)
head(links_pred)#  y predator  /  x prey

links_name<-links_pred[,c(1:3)] # store the pairwise info on presence/absence link

## Traits predicting Arctic food web structure:
# Delete taxonomic info/traits that won't be use in the analysis (as done before)
links_pred<-links_pred[,!(names(links_pred) %in% c("AphiaID.x","Group.x","Species.x", "Class.x","Order.x","Family.x","Genus.x",
                                                   "Body_defense.x","Body_defense.y", ## use continuous variable instead ("defence_cont")
                                                   "AphiaID.y","Group.y","Species.y", "Class.y","Order.y","Family.y","Genus.y",
                                                   "defence_cont.y","Kingdom.x","Phylum.x",  "Kingdom.y" ,"Phylum.y",
                                                   "Living.x", "Non_living.x","X1mm.x","X0.1_1cm.x" ,"X1_10cm.x","X10cm.x" ,"Plant.x","Animal.x",
                                                   "Suspension_feeder.x", "Deposit_feeder.x",   "Predator.x","occ"))]
links_pred$Type_cat.x<-as.factor(links_pred$Type_cat.x)
links_pred$Type_cat.y<-as.factor(links_pred$Type_cat.y)
links_pred$link<-as.numeric(as.character(links_pred$link))
colnames(links_pred)<-gsub('.x', '.R', colnames(links_pred))
colnames(links_pred)<-gsub('.y', '.C', colnames(links_pred))
names(links_pred) <- make.names(names(links_pred))

niter<-vector()
TSS<-data.frame()
VI<-data.frame()
probabilities<-data.frame()

#for (j in 1:10){   ## for sensitivity test of the subsampling factor
#print(j)
j=2 # sub-sampling factor

for (i in 1:100){ # OBS! This loop can take very long
  print(i)
  
  set.seed(1+5*i) #for replicating results
  # subsample the non-documented links to obtain a better balance in the dataset
  true_zero<-sample_n(subset(links_fw, link==0),dim(positive_link)[1]*j) # 2461*2 zero-links
  
  links_all<-rbind(true_zero,positive_link)
  
  samp <- sample(nrow(links_all), 0.75 * nrow(links_all))
  traina <- links_all[samp, ]
  testa <- links_all[-samp, ]
  
  gbm1 <- gbm(link ~ ., data = traina,
              distribution = "bernoulli", n.trees = 1100, shrinkage = 0.05,             
              interaction.depth = 10, bag.fraction = 0.5, train.fraction = 0.75,  
              n.minobsinnode = 10, cv.folds = 5, keep.data = TRUE, 
              verbose = FALSE, n.cores = 1) 

  VI<-rbind(VI, summary(gbm1))
  
  # Check performance using 5-fold cross-validation
  niter<-c(niter,gbm.perf(gbm1, method = "cv"))
  
  ### Model performance
  Predictions<-predict(gbm1, newdata = testa, n.trees = gbm.perf(gbm1, method = "cv"), type="response")
  testa$pred1<-as.factor(ifelse(Predictions>0.5,1,0))
  TSS<-rbind(TSS,CalcTSS(testa))
  
  ## and add pred
  probabilities<-rbind(probabilities,data.frame(links_name,prob=predict(gbm1, newdata = links_pred, n.trees = gbm.perf(gbm1, method = "cv"), type="response")))
}
#}

# Save the results
#save(probabilities,TSS,niter, VI, file="GBM_results.Rdata")
# Load the results
#load("GBM_results.Rdata")

colMeans(TSS)
#Accuracy (portion of predicted correctly)
(mean(TSS$d)+mean(TSS$a))/(mean(TSS$d)+mean(TSS$a)+mean(TSS$b)+mean(TSS$c))
#Sensitivty (portion of links predicted correclty)
mean(TSS$a)/(mean(TSS$a)+mean(TSS$c))
#Specificity (portion of no-links predicted correclty)
mean(TSS$d)/(mean(TSS$b)+mean(TSS$d))


###### Figure 3: Relative trait importance ####
# Plot relative influence of each variable
VI<-ddply(VI, .(var), summarise,
           sdrel.inf=sd(rel.inf, na.rm=T),
           rel.inf=mean(rel.inf, na.rm=T))

VI <- VI[order(VI$var),]

VI$var<-as.factor(c("Metabolism Type (C)","Metabolic Type (R)","High mobility (C)","High mobility (R)","Carnivore (C)","Benthic (C)",
                     "Benthic (R)","Bentho-pelagic (C)","Bentho-pelagic (R)","Body Size (C)","Body Size (R)","Body toughness (R)",
                     "Deposit Feeder (C)","Ice Dependency (C)","Ice Dependency (R)","Living Resource (C)","Dead Resource (C)","Low Mobility (C)",
                     "Low Mobility (R)","Pelagic (C)","Pelagic (R)","Herbivore (C)","Predator (C)","No Mobility (C)","No Mobility (R)", 
                     "Suspension Feeder (C)","Preferred resource size 0.1-1cm (C)","Preferred resource size 1-10cm (C)",
                     "Preferred resource size >10cm (C)","Preferred resource size <1mm (C)"))

VI$var <- factor(VI$var, levels = VI$var[order(VI$rel.inf)])

  ggplot(data=VI, aes(x=rel.inf, y=var))+
  geom_point()+
  theme_bw()+
  theme(axis.text.y = element_text(face = "bold",size = 12), legend.title = element_blank())+
  ylab("Trait Consumer (C) and Resource (R)")+
  xlab("Variable Importance (BRT)")


############# Fig supp. Degree distribution documented - predicted ########
colnames(probabilities)<-c("Consumer","Resource","link","prob")
probabilities<-data.table(probabilities)  
  
probabilities<-probabilities[,list(prob=mean(prob),
                       sdprob=sd(prob),
                       link=mean(link)),
                 by=c("Consumer","Resource")]
  
probabilities$link<-as.numeric(as.character(probabilities$link))
# number of feeding links predicted to be positive with probability of link >0.5
# in function of whether this link was documented (1) or not (0)
table(subset(probabilities, prob>0.5)$link) 
table(subset(probabilities, prob<0.5)$link) # predicted to be negative (probability of link p<0.5)
  
table(subset(probabilities, prob>0.75)$link)
table(subset(probabilities, prob<0.25)$link)

# Histogram of predicted probability of documented link (link=1)  
ggplot(subset(probabilities, link==1), aes(x=prob)) +
  geom_density(aes(x=prob), fill="black", alpha=0.8, position="identity")+
  theme_bw()+
  xlab("Probability of link")

# Histogram of predicted probability of non-documented link (link=1)
ggplot(subset(probabilities, link==0), aes(x=prob)) +
  geom_density(aes(x=prob), fill="black", alpha=0.8, position="identity")+
  theme_bw()+
  xlab("Probability of link")

probabilities$link_obsprob<-probabilities$link
probabilities[probabilities$link==0&probabilities$prob>0.5,]$link_obsprob<-1
probabilities$link_obsprob75<-probabilities$link
probabilities[probabilities$link==0&probabilities$prob>0.75,]$link_obsprob75<-1

# Summary of observed vs predicted for for each species as consumers
consumer<-ddply(probabilities, .(Consumer), summarise,
                nR_obs=sum(link),
                nR_prob=sum(prob),
                nR_obsprob=sum(link_obsprob),
                nR_obsprob75=sum(link_obsprob75))

# Summary of observed vs predicted for each species as resources
resource<-ddply(probabilities, .(Resource), summarise,
                nC_obs=sum(link),
                nC_prob=sum(prob),
                nC_obsprob=sum(link_obsprob),
                nC_obsprob75=sum(link_obsprob75))

# Degree distribution for each speices (=number of links as resources and consumers)
degree_distrib<-cbind(consumer, resource)

#Number of documented and predicted (p>0.5 and p>0.75) resources per species
consumer$Consumer <- factor(consumer$Consumer, levels = consumer$Consumer[order(consumer$nR_obs, decreasing = T)])
consumer$diff<-consumer$nR_obsprob75-consumer$nR_obs
consumer[order(consumer$diff, decreasing = T),]

ggplot(data=consumer, aes(y=nR_obsprob, x=Consumer))+
  geom_bar(aes(y=nR_obsprob, x=Consumer),stat="identity", fill="grey60", alpha=0.6)+
  geom_bar(aes(y=nR_obsprob75, x=Consumer),stat="identity", fill="grey40")+
  geom_point(aes(y=nR_obs, x=Consumer))+
  theme_bw()+
  theme(axis.text.x = element_blank(), legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylab("Number of resources")+
  ggtitle("(a) ranked by: number of documented links")


#Number of documented and predicted (p>0.5 and p>0.75) consumers per species
resource$Resource <- factor(resource$Resource, levels = resource$Resource[order(resource$nC_obs, decreasing = T)])
resource$diff<-resource$nC_obsprob75-resource$nC_obs
resource[order(resource$diff, decreasing = T),]

ggplot(data=resource, aes(y=nC_obsprob, x=Resource))+
  geom_bar(aes(y=nC_obsprob, x=Resource),stat="identity", fill="grey60", alpha=0.6)+
  geom_bar(aes(y=nC_obsprob75, x=Resource),stat="identity", fill="grey40")+
  geom_point(aes(y=nC_obs, x=Resource))+
  theme_bw()+
  theme(axis.text.x = element_blank(), legend.title = element_blank(),
        axis.title.x = element_blank())+
  ylab("Number of consumers")+
  ggtitle("(b) ranked by: number of documented links")


########### Figure 4: trophic position of "typical boreal species ##########
########### in the Arctic food web 

## "typical"boreal species has consumers in the Arctic food web: 
## Number of documented links + not documented but predicted & from which predicted with "typical" Arctic species 

# documented
typboreal_consumers_documented<-subset(probabilities, Consumer%in%typboreal&Resource%in%sp_arctic&link==1)
typboreal_consumers_documented<-ddply(typboreal_consumers_documented, .(Consumer), summarise,nlinks=length(unique(Resource)))
typboreal[!typboreal%in%typboreal_consumers_documented$Consumer] # Species with no documented resource
typboreal_consumers_documented<-rbind(typboreal_consumers_documented, data.frame(Consumer="Arctozenus risso", nlinks=0))

# not documented but predicted p>0.75 and p>0.5
typboreal_consumers_notdocumented<-subset(probabilities, Consumer%in%typboreal&Resource%in%sp_arctic&link==0&prob>0.5)
typboreal_consumers_notdocumented1<-ddply(subset(typboreal_consumers_notdocumented, prob>0.75), .(Consumer), summarise,
                                          nlinks_nodoc75=length(unique(Resource)))
typboreal_consumers_notdocumented<-ddply(typboreal_consumers_notdocumented, .(Consumer), summarise, # p>0.5
                                         nlinks_nodoc=length(unique(Resource)))
typboreal_consumers_notdocumented<-merge(typboreal_consumers_notdocumented,typboreal_consumers_notdocumented1, all.x=T)
typboreal_consumers_notdocumented[is.na(typboreal_consumers_notdocumented$nlinks_nodoc75),]$nlinks_nodoc75<-0

# not documented but predicted with "typical" Arctic species p>0.75 and p>0.5
typboreal_consumers_notdocumented_typarctic<-subset(probabilities, Consumer%in%typboreal&Resource%in%typarctic&link==0&prob>0.5)
typboreal_consumers_notdocumented_typarctic1<-ddply(subset(typboreal_consumers_notdocumented_typarctic, prob>0.75), .(Consumer), summarise,
                                                  nlinks_nodoc_arc75=length(unique(Resource)))
typboreal_consumers_notdocumented_typarctic<-ddply(typboreal_consumers_notdocumented_typarctic, .(Consumer), summarise,
                                                 nlinks_nodoc_arc=length(unique(Resource)))
typboreal_consumers_notdocumented_typarctic<-merge(typboreal_consumers_notdocumented_typarctic,typboreal_consumers_notdocumented_typarctic1, all.x=T)
typboreal_consumers_notdocumented_typarctic[is.na(typboreal_consumers_notdocumented_typarctic$nlinks_nodoc_arc75),]$nlinks_nodoc_arc75<-0
typboreal[!typboreal%in%typboreal_consumers_notdocumented_typarctic$Consumer] # Species with no predicted links with "typical" Arctic
typboreal_consumers_notdocumented_typarctic<-rbind(typboreal_consumers_notdocumented_typarctic,
                                                   data.frame(Consumer=typboreal[!typboreal%in%typboreal_consumers_notdocumented_typarctic$Consumer], nlinks_nodoc_arc=0, nlinks_nodoc_arc75=0))

## Merging all of them
typboreal_consumers<-merge(typboreal_consumers_documented,typboreal_consumers_notdocumented, by="Consumer")
typboreal_consumers<-merge(typboreal_consumers,typboreal_consumers_notdocumented_typarctic, by="Consumer")

typboreal_consumers<-merge(typboreal_consumers,f_names[,c("GROUP","Trophospecies")], by.x="Consumer", by.y="Trophospecies")

## Add a row for Arctic food web average species role
fw_arctic<-fw[rownames(fw)%in%sp_arctic,colnames(fw)%in%sp_arctic]
typboreal_consumers<-rbind(typboreal_consumers, data.frame(Consumer="Arctic",nlinks=mean(mean(colSums(fw_arctic[,!colSums(fw_arctic)==0]))),
                                                         nlinks_nodoc=0,nlinks_nodoc75=0, nlinks_nodoc_arc=0,nlinks_nodoc_arc75=0,GROUP=NA))
## Reorder by number of links (from highest to lowest)
typboreal_consumers$Consumer <- factor(typboreal_consumers$Consumer, 
                                      levels = typboreal_consumers$Consumer[order(typboreal_consumers$nlinks, decreasing = T)])
typboreal_consumers<-typboreal_consumers[order(typboreal_consumers$nlinks, decreasing = T),]
typboreal_consumers$addpred<-typboreal_consumers$nlinks+typboreal_consumers$nlinks_nodoc_arc75

p1_75<- ggplot(data=typboreal_consumers)+
  geom_bar(aes(y=(nlinks+nlinks_nodoc), x=Consumer),stat="identity", fill="grey40")+
  geom_bar(aes(y=(nlinks+nlinks_nodoc75), x=Consumer),stat="identity", fill="black")+
  geom_bar(aes(y=(nlinks+nlinks_nodoc_arc75), x=Consumer),stat="identity", fill="lightblue")+
  geom_point(aes(y=max(colSums(fw_arctic)), x="Arctic"),shape=8, size=2, inherit.aes = F)+
  geom_bar(aes(y=nlinks, x=Consumer),stat="identity",fill="white", colour="black")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.4,size = 7, angle = 90), legend.title = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  ggtitle("(a) Number of resources (generality) in the Arctic food web")+
  scale_x_discrete(labels=c("Amblyraja radiata"=expression(bold("Amblyraja radiata")), "Anarhichas denticulatus"=expression(bold("Anarhichas denticulatus")),
                            "Anarhichas lupus"=expression(bold("Anarhichas lupus")),"Anarhichas minor"=expression(bold("Anarhichas minor")),
                            "Clupea harengus"=expression(bold("Clupea harengus")),"Gadus morhua"=expression(bold("Gadus morhua")),
                            "Meganyctiphanes norvegica"=expression(bold("Meganyctiphanes norvegica")),"Pollachius virens"=expression(bold("Pollachius virens")),
                            "Melanogrammus aeglefinus"=expression(bold("Melanogrammus aeglefinus")),
                            "Sebastes"=expression(bold("Sebastes sp.")),"Sebastes mentella"=expression(bold("Sebastes mentella")),
                            "Arctic"=expression(bolditalic("Arctic food web - mean and max")), parse=TRUE))

### "typical"boreal species has resources (preys) in the Arctic food web: 
### Number of documented links + not documented but predicted & from which predicted with "typical" Arctic species 

## Documented
typboreal_resourc_documented<-subset(probabilities, Resource%in%typboreal&Consumer%in%sp_arctic&link==1)
typboreal_resourc_documented<-ddply(typboreal_resourc_documented, .(Resource), summarise,
                                    nlinks=length(unique(Consumer)))
typboreal[!typboreal%in%typboreal_resourc_documented$Resource]
typboreal_resourc_documented<-rbind(typboreal_resourc_documented, data.frame(Resource=typboreal[!typboreal%in%typboreal_resourc_documented$Resource], nlinks=0))

## not documented but predicted p>0.75 and p>0.5
typboreal_resourc_notdocumented<-subset(probabilities, Resource%in%typboreal&Consumer%in%sp_arctic&link==0&prob>0.5)
typboreal_resourc_notdocumented1<-ddply(subset(typboreal_resourc_notdocumented, prob>0.75), .(Resource), summarise,
                                        nlinks_nodoc75=length(unique(Consumer)))
typboreal_resourc_notdocumented<-ddply(typboreal_resourc_notdocumented, .(Resource), summarise,
                                       nlinks_nodoc=length(unique(Consumer)))
typboreal_resourc_notdocumented<-merge(typboreal_resourc_notdocumented,typboreal_resourc_notdocumented1, all.x=T)
typboreal_resourc_notdocumented[is.na(typboreal_resourc_notdocumented$nlinks_nodoc75),]$nlinks_nodoc75<-0

## not documented but predicted only Arctic p>0.75 and p>0.5
typboreal_resourc_notdocumented_typarctic<-subset(probabilities, Resource%in%typboreal&Consumer%in%typarctic&link==0&prob>0.5)
typboreal_resourc_notdocumented_typarctic1<-ddply(subset(typboreal_resourc_notdocumented_typarctic, prob>0.75), .(Resource), summarise,
                                                  nlinks_nodoc_arc75=length(unique(Consumer)))
typboreal_resourc_notdocumented_typarctic<-ddply(typboreal_resourc_notdocumented_typarctic, .(Resource), summarise,
                                                 nlinks_nodoc_arc=length(unique(Consumer)))
typboreal_resourc_notdocumented_typarctic<-merge(typboreal_resourc_notdocumented_typarctic,typboreal_resourc_notdocumented_typarctic1, all.x=T)
typboreal_resourc_notdocumented_typarctic[is.na(typboreal_resourc_notdocumented_typarctic$nlinks_nodoc_arc75),]$nlinks_nodoc_arc75<-0

typboreal[!typboreal%in%typboreal_resourc_notdocumented_typarctic$Resource] # Species not preyed upon by "typical" Arctic
typboreal_resourc_notdocumented_typarctic<-rbind(typboreal_resourc_notdocumented_typarctic, data.frame(Resource=typboreal[!typboreal%in%typboreal_resourc_notdocumented_typarctic$Resource], nlinks_nodoc_arc=0, nlinks_nodoc_arc75=0))

## Merging all of them
typboreal_resourc<-merge(typboreal_resourc_documented,typboreal_resourc_notdocumented, by="Resource")
typboreal_resourc<-merge(typboreal_resourc,typboreal_resourc_notdocumented_typarctic, by="Resource")
typboreal_resourc<-merge(typboreal_resourc,f_names[,c(4,6)], by.x="Resource", by.y="Trophospecies")

## Add a row for Arctic food web average species role
fw_arctic<-fw[rownames(fw)%in%sp_arctic,colnames(fw)%in%sp_arctic]
typboreal_resourc<-rbind(typboreal_resourc, data.frame(Resource="Arctic",nlinks=mean(mean(colSums(fw_arctic[,!colSums(fw_arctic)==0]))),
                                                     nlinks_nodoc=0,nlinks_nodoc75=0, nlinks_nodoc_arc=0,nlinks_nodoc_arc75=0, GROUP=NA))

## Reorder by number of links (from highest to lowest)
typboreal_resourc$Resource <- factor(typboreal_resourc$Resource, 
                                    levels = typboreal_resourc$Resource[order(typboreal_resourc$nlinks, decreasing = T)])
typboreal_resourc<-typboreal_resourc[order(typboreal_resourc$nlinks, decreasing = T),]

p2_75<-ggplot(data=typboreal_resourc)+
  geom_bar(aes(y=(nlinks+nlinks_nodoc), x=Resource),stat="identity", fill="grey40")+
  geom_bar(aes(y=(nlinks+nlinks_nodoc75), x=Resource),stat="identity", fill="black")+
  geom_bar(aes(y=(nlinks+nlinks_nodoc_arc75), x=Resource),stat="identity", fill="lightblue")+
  geom_bar(aes(y=nlinks, x=Resource),stat="identity",fill="white", colour="black")+
  geom_point(aes(y=max(colSums(fw_arctic)), x="Arctic"),shape=8, size=2)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.4,size = 7, angle = 90), legend.title = element_blank(),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  ggtitle("(b) Number of consumers (vulnerability) in the Arctic food web")+
  scale_x_discrete(labels=c("Amblyraja radiata"=expression(bold("Amblyraja radiata")), "Anarhichas denticulatus"=expression(bold("Anarhichas denticulatus")),
                            "Anarhichas lupus"=expression(bold("Anarhichas lupus")),"Anarhichas minor"=expression(bold("Anarhichas minor")),
                            "Clupea harengus"=expression(bold("Clupea harengus")),"Gadus morhua"=expression(bold("Gadus morhua")),
                            "Meganyctiphanes norvegica"=expression(bold("Meganyctiphanes norvegica")),"Pollachius virens"=expression(bold("Pollachius virens")),
                            "Melanogrammus aeglefinus"=expression(bold("Melanogrammus aeglefinus")),
                            "Sebastes"=expression(bold("Sebastes sp.")),"Sebastes mentella"=expression(bold("Sebastes mentella")),
                            "Arctic"=expression(bolditalic("Arctic food web - mean and max")), parse=TRUE))

grid.arrange(p1_75,p2_75,nrow=2)


########### Figure 5: Example of shifting species COD ####

fw_oldarctic<-fw[rownames(fw)%in%c(sp_arctic,"Gadus morhua"),colnames(fw)%in%c(sp_arctic,"Gadus morhua")]
trophind<-TrophInd(fw_oldarctic)
trophind$species<-rownames(trophind)
set.seed(134)
trophind$xx<-runif(length(trophind$species))
trophind[trophind$TL<1.5,]$TL<-1.5
  
cod<-subset(probabilities, Consumer =="Gadus morhua"|Resource =="Gadus morhua")
cod<-subset(cod, Consumer %in%sp_arctic|Resource%in%sp_arctic)
cod<-subset(cod, link ==1|prob >0.5)
cod<-merge(cod, trophind[,-2], by.x="Consumer", by.y="species")
cod<-merge(cod, trophind[,-2], by.x="Resource", by.y="species")

# Put cod in the center of the x-axis
cod[cod$Resource=="Gadus morhua",]$xx.y<-0.5
cod[cod$Consumer=="Gadus morhua",]$xx.x<-0.5

## Cod as Prey in the Arctic food web
  ggplot()+
  geom_point(data=trophind, aes(xx,TL), alpha=0.5, size=1)+
  geom_segment(data=subset(cod, link==1&Resource=="Gadus morhua"), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=0.5, col="black", linetype="dotted") +
  geom_segment(data=subset(cod, link==0&prob>0.5&prob<0.75&Resource=="Gadus morhua"), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y),
               size=0.7, col="grey50") +
  geom_segment(data=subset(cod, link==0&prob>0.75&Resource=="Gadus morhua"), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y),
               size=1, col="black") +
  scale_size_identity()+
  geom_point(data=subset(cod, link==0&prob>0.75&Resource=="Gadus morhua"&Consumer%in%typarctic), mapping=aes(x=xx.x, y=TL.x),col="lightblue", size=3) +
  geom_label(data=subset(cod, link==0&prob>0.75&Resource=="Gadus morhua"&Consumer%in%typarctic), mapping=aes(x=xx.x, y=TL.x, label=Consumer),col="black", fill="lightblue") +
  scale_y_continuous(breaks=c(1.5,2,3,4,5),labels=c("1", "2", "3","4","5"))+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks = element_blank(), legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(),
        plot.title = element_text(size=12))+
  labs(y="Trophic level")+
  ggtitle("(b) Number of consumers (Vulnerability)")

## Cod as Predator in the Arctic food web
  ggplot()+
  geom_point(data=trophind, aes(xx,TL), alpha=0.5, size=1)+
  geom_segment(data=subset(cod, link==1&Consumer=="Gadus morhua"), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=0.5, col="black", linetype="dotted") +
  geom_segment(data=subset(cod, link==0&prob>0.5&Consumer=="Gadus morhua"), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y),
               size=0.7, col="grey50") +
  geom_segment(data=subset(cod, link==0&prob>0.75&Consumer=="Gadus morhua"), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y),
               size=1, col="black") +
  scale_size_identity()+
  geom_point(data=subset(cod, link==0&prob>0.75&Consumer=="Gadus morhua"&Resource%in%typarctic), mapping=aes(x=xx.y, y=TL.y),col="lightblue", size=3) +
  geom_label(data=subset(cod, link==0&prob>0.75&Consumer=="Gadus morhua"&Resource%in%typarctic), mapping=aes(x=xx.y, y=TL.y, label=Resource),col="black", fill="lightblue") +
  scale_y_continuous(breaks=c(1.5,2,3,4,5),labels=c("1", "2", "3","4","5"))+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks = element_blank(), legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(),
        plot.title = element_text(size=12))+
  labs(y="Trophic level")+
  ggtitle("(a) Number of resources (Generality)")


########### Figure 6: Rewiring Arctic Warm vs Very_warm ####
  
  ### Plot 1: Warm food web
  fw_oldarctic<-fw[rownames(fw)%in%sp_arctic,colnames(fw)%in%sp_arctic]
  fw_oldarctici<-graph.adjacency(as.matrix(fw_oldarctic))
  spingc <- spinglass.community(fw_oldarctici, spins = 3) ## Cut species into 3 clusters
  
  sping_mod <- spinglass.community(fw_oldarctici) ## Calculate modularity
  mod_old<-sping_mod$modularity 
  links_old<-sum(fw_oldarctic) ## Calculate number of links
  con_old<-sum(fw_oldarctic)/(dim(fw_oldarctic)[1]^2) ## Calculate connectance

  # If needed for nicer visualision, reassign grouping
  # spingc$membership<-dplyr::recode(spingc$membership,`1`=3,`2`=2,`3`=1)
  
  trophind_oldarctic<-TrophInd(fw_oldarctic)
  trophind_oldarctic$species<-rownames(trophind_oldarctic)
  trophind_oldarctic$xx<-runif(length(trophind_oldarctic$species))+(spingc$membership*1.2)
  trophind_oldarctic[trophind_oldarctic$TL<1.5,]$TL<-1.5
  
  probabilities_oldarctic<-subset(probabilities, link==1&Resource%in%sp_arctic&Consumer%in%sp_arctic)
  
  probabilities_oldarctic<-merge(probabilities_oldarctic, trophind_oldarctic, by.x="Consumer", by.y="species")
  probabilities_oldarctic<-merge(probabilities_oldarctic, trophind_oldarctic, by.x="Resource", by.y="species")
  
  trophind_oldarctic<-merge(trophind_oldarctic, f_names[,c("GROUP","Trophospecies","Phylum")], by.x="species", by.y="Trophospecies", all.x=T)
 
    ggplot()+
    geom_segment(data=probabilities_oldarctic, mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=0.2, col="grey80") +
     geom_point(data=trophind_oldarctic, aes(xx,TL, col=GROUP), size=1.5)+
    scale_colour_manual( values = c("Benthos"="orange","Zooplankton"="cyan", "Fish"="blue", "Basal"="darkgrey","Mammals"="lightpink","Seabirds"="magenta")) +
    annotate("text", label =c("(a) Arctic food web 2004-2007"), x = 2.5, y = 5.5, size = 4)+
    annotate("text", label =paste0( "nLinks = ",links_old), x = 2, y = 5.2, size = 3)+
    annotate("text", label =paste0( "Connectance = ", round(con_old*100,1),"%"), x = 2, y = 5, size = 3)+
    annotate("text", label =paste0( "Modularity = ", round(mod_old,3)), x = 2, y = 4.8, size = 3)+
    scale_y_continuous(breaks=c(1.5,2,3,4,5),labels=c("1", "2", "3","4","5"))+
    theme_bw()+
    theme(panel.border = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks = element_blank(), legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_text())+
    labs(y="Trophic level")
  
  ### Plot 2: Observed + new sp in Arctic FW
  fw_newarctic<-fw[rownames(fw)%in%c(sp_arctic,as.vector(new_sp)),colnames(fw)%in%c(sp_arctic,as.vector(new_sp))]
  fw_newarctici<-graph.adjacency(as.matrix(fw_newarctic))
  spingc <- spinglass.community(fw_newarctici, spins=3)
  
  sping_mod <- spinglass.community(fw_newarctici)
  mod_new<-sping_mod$modularity
  links_new<-sum(fw_newarctic)
  con_new<-sum(fw_newarctic)/(dim(fw_newarctic)[1]^2)

  trophind_newarctic<-TrophInd(fw_newarctic)
  trophind_newarctic$species<-rownames(trophind_newarctic)
  trophind_newarctic$xx<-runif(length(trophind_newarctic$species))+(spingc$membership*1.2)
  trophind_newarctic[trophind_newarctic$TL<1.5,]$TL<-1.5
  
  probabilities_newarctic<-subset(probabilities, link==1&Resource%in%c(sp_arctic,new_sp)&Consumer%in%c(sp_arctic,new_sp))
  probabilities_newarctic<-merge(probabilities_newarctic, trophind_newarctic, by.x="Consumer", by.y="species")
  probabilities_newarctic<-merge(probabilities_newarctic, trophind_newarctic, by.x="Resource", by.y="species")
  
  trophind_newarctic<-merge(trophind_newarctic, f_names[,c("GROUP","Trophospecies","Phylum")], by.x="species", by.y="Trophospecies", all.x=T)
  
  ggplot()+
    geom_segment(data=probabilities_newarctic, mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=0.2, col="grey80") +
    geom_segment(data=subset(probabilities_newarctic, Consumer%in%new_sp), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=0.4, col="grey40") +
    geom_point(data=trophind_newarctic, aes(xx,TL, col=GROUP), size=1)+
    geom_point(data=subset(trophind_newarctic, species%in%new_sp), aes(xx,TL, col=GROUP), size=3)+
    scale_colour_manual( values = c("Benthos"="orange","Zooplankton"="cyan", "Fish"="blue", "Basal"="darkgrey","Mammals"="lightpink","Seabirds"="magenta")) +
    geom_label(data=subset(trophind_newarctic, species%in%new_sp), aes(xx,TL, label=species), size=2)+
    annotate("text", label =c("(b) Arctic food web 2014-2017"), x = 2.5, y = 5.5, size = 4)+
    annotate("text", label =paste0( "nLinks = ",links_new), x = 2, y = 5.2, size = 3)+
    annotate("text", label =paste0( "Connectance = ", round(con_new*100,1),"%"), x = 2, y = 5, size = 3)+
    annotate("text", label =paste0( "Modularity = ", round(mod_new,3)), x = 2, y = 4.8, size = 3)+
    scale_y_continuous(breaks=c(1.5,2,3,4,5),labels=c("1", "2", "3","4","5"))+
    theme_bw()+
    theme(panel.border = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks = element_blank(), legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_text())+
    labs(y="Trophic level")
  
  ### Plot 3: Observed + new sp with predicted links in Arctic FW (all pred >0.75)
  
  probabilities_newarctic2<-subset(probabilities, link==1&Resource%in%c(sp_arctic,new_sp)&Consumer%in%c(sp_arctic,new_sp))
  probabilities_newarcticlinks<-subset(probabilities, prob>0.75&link==0&Resource%in%c(new_sp)&Consumer%in%c(typarctic))
  probabilities_newarcticlinks<-rbind(probabilities_newarcticlinks,subset(probabilities, prob>0.75&link==0&Consumer%in%c(new_sp)&Resource%in%c(typarctic)))
  ##Keep all new potential interaction with species present in Arctic and not Boreal
  probabilities_newarcticlinks[probabilities_newarcticlinks$link==0&probabilities_newarcticlinks$Consumer%in%c(typarctic,new_sp)&
                          probabilities_newarcticlinks$Resource%in%c(typarctic,new_sp),]$link<-1
  
  
  fw_newarctic<-fw[rownames(fw)%in%c(sp_arctic,as.vector(new_sp)),colnames(fw)%in%c(sp_arctic,as.vector(new_sp))]
  
  linksT<-melt(t(fw_newarctic))
  colnames(linksT)<-c("Consumer","Resource","link")
  summary(linksT)
  linksT<-rbind(linksT, probabilities_newarcticlinks[,c("Consumer","Resource","link")])
  
  fw_newlinksarctic<-reshape2::dcast(Consumer~Resource, value.var="link", sum, data=linksT, fill=0)
  rownames(fw_newlinksarctic)<-fw_newlinksarctic$Consumer
  fw_newlinksarctic<-fw_newlinksarctic[,-1]
  fw_newlinksarctic<-t(fw_newlinksarctic)
  
  fw_newlinksarctici<-graph.adjacency(as.matrix(fw_newlinksarctic))
  spingc <- spinglass.community(fw_newlinksarctici, spins=3)
  
  sping_mod <- spinglass.community(fw_newlinksarctici)
  mod_newlinks<-sping_mod$modularity
  links_newlinks<-sum(fw_newlinksarctic)
  con_newlinks<-sum(fw_newlinksarctic)/(dim(fw_newlinksarctic)[1]^2)
   
  trophind_newlinksarctic<-TrophInd(fw_newlinksarctic)
  trophind_newlinksarctic$species<-rownames(trophind_newlinksarctic)
  trophind_newlinksarctic$xx<-runif(length(trophind_newlinksarctic$species))+(spingc$membership*1.2)
  trophind_newlinksarctic[trophind_newlinksarctic$TL<1.5,]$TL<-1.5
  
  probabilities_newarctic3<-merge(probabilities_newarctic2, trophind_newlinksarctic, by.x="Consumer", by.y="species")
  probabilities_newarctic3<-merge(probabilities_newarctic3, trophind_newlinksarctic, by.x="Resource", by.y="species")
  
  probabilities_newarcticlinks2<-merge(probabilities_newarcticlinks, trophind_newlinksarctic, by.x="Consumer", by.y="species")
  probabilities_newarcticlinks2<-merge(probabilities_newarcticlinks2, trophind_newlinksarctic, by.x="Resource", by.y="species")
  
  trophind_newlinksarctic<-merge(trophind_newlinksarctic, f_names[,c("GROUP","Trophospecies","Phylum")], by.x="species", by.y="Trophospecies", all.x=T)
  
  ggplot()+
    geom_segment(data=probabilities_newarctic3, mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=0.2, col="grey80") +
    geom_segment(data=subset(probabilities_newarctic3, Consumer%in%new_sp), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=0.4, col="grey40") +
    geom_segment(data=subset(probabilities_newarcticlinks2, Consumer%in%new_sp), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=1, col="lightblue") +
    geom_segment(data=subset(probabilities_newarcticlinks2, Resource%in%new_sp), mapping=aes(x=xx.x, y=TL.x, xend=xx.y, yend=TL.y), size=1, col="lightblue") +
    geom_point(data=trophind_newlinksarctic, aes(xx,TL, col=GROUP), size=1)+
    geom_point(data=subset(trophind_newlinksarctic, species%in%new_sp), aes(xx,TL, col=GROUP), size=3)+
    scale_colour_manual( values = c("Benthos"="orange","Zooplankton"="cyan", "Fish"="blue", "Basal"="darkgrey","Mammals"="lightpink","Seabirds"="magenta")) +
    geom_label(data=subset(trophind_newlinksarctic, species%in%new_sp), aes(xx,TL, label=species), size=2)+
    annotate("text", label =c("(c) Arctic food web 2014-2017 - predicted links"), x = 3, y = 5.5, size = 4)+
    annotate("text", label =paste0( "nLinks = ",links_newlinks), x = 2, y = 5.2, size = 3)+
    annotate("text", label =paste0( "Connectance = ", round(con_newlinks*100,1),"%"), x = 2, y = 5, size = 3)+
    annotate("text", label =paste0( "Modularity = ", round(mod_newlinks,3)), x = 2, y = 4.8, size = 3)+
    scale_y_continuous(breaks=c(1.5,2,3,4,5),labels=c("1", "2", "3","4","5"))+
    theme_bw()+
    theme(panel.border = element_blank(),
          axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks = element_blank(), legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_text())+
    labs(y="Trophic level")

  
  