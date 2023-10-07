######################################
######### PREDICTING HF ##############
######################################
library(car)
library(mvtnorm)
library(MVN)
library(mvnormtest)
library(heplots)
library(MASS)
library(class)
library(e1071)
options(rgl.debug = TRUE)
library(rgl)
library(dbscan)
library(cluster)
library(fields)
library(glmnet)
library(leaps)
library(nlmeU)
library(nlme) 
library(corrplot)
library(lattice)
library(plot.matrix)
library(lme4)
library(insight)
library(sp)         
library(lattice)     
library(gstat) 
library(fda)
library(KernSmooth)
library(fields)  
library(dplyr)
install.packages("klaR")
library(klaR)
library(pROC)
library(ROCR)
library(cluster)
library(fields)
library(dplyr)
library(tidyverse)


##### DATA PRE PROCESSING #####
data_final <- read.csv('data_final.csv',header=T)
data_final <- subset(data_final, labelOUT != 'PERSO' ) # togliamo i pazienti persi
data_final <- dplyr::filter(data_final, tipo_prest == 41)   # teniamo solo pazienti ospedalizzati
data_final$labelOUT<-ifelse(data_final$labelOUT == "DECEDUTO", 1, 0)
data_final = data_final[,-1]

data_clean <- data[!duplicated(data$COD_REG),] # data senza duplicati (prima osservazione)

dati_maschi <- subset(data_clean, SESSO == "M")
dati_femmine <- subset(data_clean, SESSO == "F")

dev.new()
hist(dati_femmine$eta_Min, breaks = 40, col = rgb(1,0,0,0.5), xlab = "Age", main = "Age of first heart failure", ylab = "Occurrencies")
hist(dati_maschi$eta_Min, breaks = 40, col = rgb(1,1,0,0.5), add = T)
legend("topleft", legend=c("Females","Males"), col=c(rgb(1,0,0,0.5), rgb(1,1,0,0.5)), pt.cex=2, pch=15)
graphics.off()

# tengo solo le prime osservazioni per ogni paziente 
library(dplyr)
data2<- data%>%distinct(COD_REG, .keep_all = TRUE)  

##### PIE CHART #####
# vediamo numero deceduti, persi, troncati
attach(data2)
library(RColorBrewer)
labelOUT <- factor(labelOUT, levels=c('TRONCATO', 'DECEDUTO', 'PERSO'))
pie(table(data$labelOUT), col=brewer.pal(n = length(labelOUT), name = 'Pastel2'))
# deceduto 27.18%
# troncato 72.42%
# perso 0.4%

##### TRY PCA (INCONCLUSIVE) #####
# è possibile fare un PCA sulle variabili binarie? Sarebbe interessante vedere quali sono le malattie con più variabilità
data_binary1<-data[,13:32]
data_binary<-subset(data_binary1, tipo_prest=='41')
pc.binary<- princomp(data_binary,scores = TRUE)
pc.binary
summary(pc.binary)

layout(matrix(c(2,3,1,3),2,byrow=T))
barplot(pc.binary$sdev^2, las=2, main='Principal Components', ylab='Variances')
barplot(sapply(data_binary,sd)^2, las=2, main='Original variables', ylab='Variances')
plot(cumsum(pc.binary$sdev^2)/sum(pc.binary$sde^2), type='b', axes=F, xlab='number of components', ylab='contribution to the total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(data_binary),labels=1:ncol(data_binary),las=2)

scores.binary <- pc.binary$scores
scores.binary

detach(data2)

##### MORE DESCRIPTIVE ANALYSIS #####
attach(data)

codici <-factor(COD_REG)
table(codici)
length(table(codici)) # 187499 pazienti diversi

table(SESSO)

#data_reducted<- data[is.na(data)==FALSE]

ospedali <- factor(strutt_id)
table(ospedali)
length(table(ospedali))

table(labelOUT)

#colSums(is.na(HF.2006_2012))  11237234

stato <- factor(labelOUT)
table(stato)/length(COD_REG)

detach(data)

attach(data2)
num_pazienti = length(SESSO)  # numero pazienti: 187499
fr_sesso = table(SESSO)/num_pazienti
fr_sesso  # frequenze F: 0.5079867, M: 0.4920133 

table(SESSO) # numero F: 95247 M: 92252

fr_ASL = table(ASL_RESIDENZA)/num_pazienti  # divisione dei pazienti nelle varie ASL
fr_ASL 
table(ASL_RESIDENZA)  # numero pazienti in ogni ASL

# Riduzione dati
data_reducted <- subset(data2, labelOUT != 'PERSO' ) # togliamo i dati persi dal dataset delle prime osservazioni: rimangono 186055 pazienti
data_reductedTOT <- subset(data, labelOUT != 'PERSO' )  # togliamo i dati persi dal dataset generale: rimangono 11772721 osservazioni

codici <- data_reductedTOT$COD_REG
head(codici)

summary(codici)
frequenze_codici <- table(codici)
head(frequenze_codici)
summary(frequenze_codici)
cond <- frequenze_codici > 200
head(cond)
table(cond)

codici_unici <- data_reducted$COD_REG
tenere <- subset(codici_unici, cond == 'TRUE')
data_tenuti <- subset(data_reductedTOT, data_reductedTOT$COD_REG %in% tenere)

codici2 <- data_reducted$COD_REG

summary(codici)
frequenze_codici2 <- table(codici2)
head(frequenze_codici2)
summary(frequenze_codici2)
cond2 <- frequenze_codici2 > 50
head(cond2)
table(cond2)

detach(data2)

##### ANALYSIS BY ASL #####
data_reducted$death <- ifelse(data_reducted$labelOUT == 'DECEDUTO', 1,0)
attach(data_reducted)
pazientiASL <- table(data_reducted$ASL_RESIDENZA)
freq_death_ASL <-rep(1,15)

freq_death_ASL <-tapply(death, ASL_RESIDENZA, mean)
detach(data_reducted)

# data_final<- data_tenuti%>%distinct(COD_REG, .keep_all = TRUE)
data_final<- data_reductedTOT%>%distinct(COD_REG, .keep_all = TRUE)
data_final$death <- ifelse(data_final$labelOUT == 'DECEDUTO', 1,0)

attach(data_final)
pazientiASL <- table(data2$ASL_RESIDENZA)
pazientiASL
freq_death_ASL <-tapply(death, ASL_RESIDENZA, mean)
freq_death_ASL

library(RColorBrewer)
x11(width = 60, height = 30)
barplot(freq_death_ASL, col=brewer.pal(n = 12, name = 'Set3'), legend = TRUE, xlim = c(0,20))
# barplot(freq_death_ASL, col=cm.colors(n = length(freq_death_ASL)))

detach(data_final)

##### PIE CHARTS BINARY VARIABLES #####
library(tidyverse)
data_osp <- data_reductedTOT %>% dplyr::group_by(COD_REG) %>% slice_tail() # tengo l'ultima osservazione di ogni paziente
data_osp <- dplyr::filter(data_osp, tipo_prest == 41)  # tengo solo i pazienti OSPEDALIZZATI
data_osp$death <- NULL   # rimuovo la colonna death (non mi serve)
library(RColorBrewer)
names <- colnames(data_osp)
x11(width = 40, height = 30)
par(mfrow=c(4,6))
for(var in c(1:24)){
  uno <- sum(data_osp[,var+12])  # numero di 1 nelle colonne
  pie(c(uno, 145980-uno), col=brewer.pal(n = 3, name = 'Pastel2'), main = names[var+12], labels = c('yes','no'))
  print(names[var+12])
  print(uno/145980*100)
}

data_osp$msi <- rowSums( data_osp[,13:32] )/20  # colonna con la percentuale di malattie

# pie charts of variables to keep
library(tidyverse)
data_osp$wtloss <- NULL   # rimuovo le colonne che non mi servono
data_osp$hemiplegia <- NULL
data_osp$alcohol <- NULL
data_osp$coagulopathy <- NULL
data_osp$compdiabetes<- NULL
data_osp$liver <- NULL
data_osp$psychosis <- NULL
data_osp$pulmcirc <- NULL
data_osp$hivaids <- NULL
data_osp$ICD <- NULL
data_osp$CABG <- NULL
data_osp$PTCA <- NULL
library(RColorBrewer)
names <- colnames(data_osp)
x11(width = 40, height = 30)
par(mfrow=c(2,5))
for(var in c(1:10)){
  uno <- sum(data_osp[,var+12])  # numero di 1 nelle colonne
  pie(c(uno, 145980-uno), col=brewer.pal(n = 3, name = 'Pastel2'), main = names[var+12], labels = c(uno,145980-uno))
  print(names[var+12])
  print(uno/145980*100)
}

x11()
for(var in c(1:10)){
  uno[var] <- sum(data_osp[,var+12])  # numero di 1 nelle colonne
}
barplot(uno/145980, ylim = c(0,1))

##### CLUSTERING #####
# dataset ospedalizzati
ospedalizzati <- read.csv('ospedalizzati.csv', header = T,sep = ";")
# Aggiunta media patologie
ospedalizzati$msi <- rowSums( ospedalizzati[,13:32] )/20  # colonna con la percentuale di malattie

# Hierarchical clustering --> not feasible, computationally too heavy
patologie <- ospedalizzati[, c(15,16,20,21,22,25,26,28,32)] # teniamo le più presenti
pat.e <- dist(patologie,  method='euclidian')

x11()
image(1:150,1:150,as.matrix(pat.e), main='metrics: Euclidean', asp=1, xlab='i', ylab='j') # troppo grande

pat.es <- hclust(pat.e, method='single')
pat.ea <- hclust(pat.e, method='average')
pat.ec <- hclust(pat.e, method='complete')


## OTHER POSSIBLE APPROACH: KEEP ONLY PATIENTS WITH LONG HISTORY (>1 yrs)
# parto da data_reduced (dati non PERSI)
data_reducted <- subset(data,labelOUT!='PERSO')
# tengo solo le ultime osservazioni
library(dplyr)
data_final<- data_reducted%>%distinct(COD_REG, .keep_all = TRUE)

# vedere la differenza tra se è più di un anno
data_final <- subset(data_final, difftime(data_final$data_studio_out,data_final$data_rif_ev)>365)

#write.csv(data_final, "data_final.csv") # data set di partenza
data_final <- read.csv('data_final.csv', header = T)

##### CLUSTERING WITH KMODES #####
# Con kmodes, è quello da usare per il progetto

data_final <- read.csv('data_final.csv',header=T)
data_final <- subset(data_final, labelOUT != 'PERSO' ) # togliamo i pazienti persi
data_final <- dplyr::filter(data_final, tipo_prest == 41)   # teniamo solo pazienti ospedalizzati
data_final$labelOUT<-ifelse(data_final$labelOUT == "DECEDUTO", 1, 0)
data_final = data_final[,-1]

patologie <- data_final[, c(15,16,20,21,22,25,26,28,32)] # solo variabili più rilevanti
# K-means per variabili dicotomiche
set.seed(1)
ss2<-kmodes(patologie, modes=4, iter.max = 30, weighted = FALSE, fast = TRUE)
data_final$cluster <- factor(ss2$cluster) 

n1 <- sum(data_final$cluster == 1)
n1 # 18660
n2 <- sum(data_final$cluster == 2)
n2 # 35959
n3 <- sum(data_final$cluster == 3)
n3 # 64245
n4 <- sum(data_final$cluster == 4)
n4 # 11576


freq_pat1 <- colMeans(patologie[data_final$cluster == 1,])
freq_pat1 # più presenti: pumonarydz, hypertension
freq_pat2 <- colMeans(patologie[data_final$cluster == 2,])
freq_pat2 # più presenti: arrhythmia
freq_pat3 <- colMeans(patologie[data_final$cluster == 3,])
freq_pat3 # più presenti: nessuna
freq_pat4 <- colMeans(patologie[data_final$cluster == 4,])
freq_pat4 # più presenti: renal



# Frequenza patologie nei cluster:
x11()
par(mfrow = c(4,1))
barplot(freq_pat1, col = 'green', main = 'Cluster 1', ylim = c(0,1))
barplot(freq_pat2, col = 'red', main = 'Cluster 2', ylim = c(0,1))
barplot(freq_pat3, col = 'blue', main = 'Cluster 3', ylim = c(0,1))
barplot(freq_pat4, col = 'yellow', main = 'Cluster 4', ylim = c(0,1))

# Distribuzione dei cluster nelle ASL [ASL_residenza è la colonna 4]

my_palette <- c(#names had to be removed, data is protected by NDA)

my_palette_graph <- c(#names had to be removed, data is protected by NDA)

my_rainbow <- function(n, start = 0, end = 1, alpha = 1) {
  color_indices <- floor(seq(start, end, length.out = n + 1) * length(my_palette_graph))
  colors <- my_palette_graph[color_indices]
  if (alpha < 1) {
    colors <- adjustcolor(colors, alpha.f = alpha)
  }
  colors
}


x11()
par(mfrow = c(2,2))
plot(table(data_final$ASL_RESIDENZA[which(data_final$cluster == 1)]),ylab = '' , 
     xlab = '' , main = 'Cluster 1', col = '#NDA', ylim = c(0,10000), lwd=9)
plot(table(data_final$ASL_RESIDENZA[which(data_final$cluster == 2)]),ylab = '' , 
     xlab = '' , main = 'Cluster 2', col = '#NDA', ylim = c(0,10000), lwd=9)
plot(table(data_final$ASL_RESIDENZA[which(data_final$cluster == 3)]),ylab = '' , 
     xlab = '' , main = 'Cluster 3', col = '#NDA', ylim = c(0,10000), lwd=9)
plot(table(data_final$ASL_RESIDENZA[which(data_final$cluster == 4)]),ylab = '' , 
     xlab = '' , main = 'Cluster 4', col = '#NDA', ylim = c(0,10000), lwd=9)



##### Logistic Regression #####


# convertire le variabili sesso, cluster (già fatto prima) e ASL come factor
data_final$SESSO = factor(data_final$SESSO)
data_final$SESSO_bin = ifelse(data_final$SESSO == "M", 1, 0)
data_final$ASL_RESIDENZA = factor(data_final$ASL_RESIDENZA)

# riordino in maniera casuale il dataframe e lo divido in due parti, una per il traning dei modelli
# e l'altra per testarli

# Riordina in maniera casuale il dataframe e divido training data e test data

n = dim(data_final)[1]
data_final <- data_final %>% arrange(sample(n))
train = data_final[1:(n-40000),]
test = data_final[-c(1:(n-40000)),]


logreg <- glm(labelOUT ~ SESSO + cluster + ASL_RESIDENZA  +  eta_Min
              + eta_Min:cluster + eta_Min:ASL_RESIDENZA + SESSO:cluster,
              data = train, family='binomial')

# trattare l'ASL come factor è controverso perché l'ASL induce una dipendenza all'interno delle 
# osservazioni: non sarebbe corretto trattarla in questo modo

summary(logreg)
# Come leggere il summary (cos' come per tutti gli altri modelli):
# qui non abbiamo usato dummy variables per identificare l'appartenenza ad un cluster o ad una ASL
# quindi tutti i valori sono fatti rispetto al cluster 1 e ASL 301.
# E' come se il cluster 1 abbia un beta pari a 0, da quel che vedo io il cluster 2 ha un beta negativo
# quindi per il cluster due ci sarà una minore probabilità di morte (essendo che la variabile risposta
# è = 1 nel caso di paziente deceduto)


##Probabilità in funzione dell'età del primo arresto cardiaco

fittizio_c1_A305 = data.frame(SESSO = 'M', ASL_RESIDENZA = 305, cluster = 1, eta_Min = 20)
for( i in 21:100)
  fittizio_c1_A305 = rbind(fittizio_c1_A305, c(SESSO = 'M', ASL_RESIDENZA = 305, cluster = 1, eta_Min = i))
fittizio_c1_A305
fittizio_c1_A305$ASL_RESIDENZA = factor(fittizio_c1_A305$ASL_RESIDENZA)
fittizio_c1_A305$SESSO = factor(fittizio_c1_A305$SESSO)
fittizio_c1_A305$cluster= factor(fittizio_c1_A305$cluster)
fittizio_c1_A305$eta_Min = as.numeric(fittizio_c1_A305$eta_Min)
previsione_c1_A305 = predict(logreg, fittizio_c1_A305)
prob1 = exp(previsione_c1_A305)/(1+exp(previsione_c1_A305) )

fittizio_c2_A305 = data.frame(SESSO = 'M', ASL_RESIDENZA = 305, cluster = 2, eta_Min = 20)
for( i in 21:100)
  fittizio_c2_A305 = rbind(fittizio_c2_A305, c(SESSO = 'M', ASL_RESIDENZA = 305, cluster =2, eta_Min = i))
fittizio_c2_A305$ASL_RESIDENZA = factor(fittizio_c2_A305$ASL_RESIDENZA)
fittizio_c2_A305$SESSO = factor(fittizio_c2_A305$SESSO)
fittizio_c2_A305$cluster= factor(fittizio_c2_A305$cluster)
fittizio_c2_A305$eta_Min = as.numeric(fittizio_c2_A305$eta_Min)
previsione_c2_A305 = predict(logreg, fittizio_c2_A305)
prob2 = exp(previsione_c2_A305)/(1+exp(previsione_c2_A305) )

fittizio_c3_A305 = data.frame(SESSO = 'M', ASL_RESIDENZA = 305, cluster = 3, eta_Min = 20)
for( i in 21:100)
  fittizio_c3_A305 = rbind(fittizio_c3_A305, c(SESSO = 'M', ASL_RESIDENZA = 305, cluster = 3, eta_Min = i))
fittizio_c3_A305$ASL_RESIDENZA = factor(fittizio_c3_A305$ASL_RESIDENZA)
fittizio_c3_A305$SESSO = factor(fittizio_c3_A305$SESSO)
fittizio_c3_A305$cluster= factor(fittizio_c3_A305$cluster)
fittizio_c3_A305$eta_Min = as.numeric(fittizio_c3_A305$eta_Min)
previsione_c3_A305 = predict(logreg, fittizio_c3_A305)
prob3 = exp(previsione_c3_A305)/(1+exp(previsione_c3_A305) )

fittizio_c4_A305 = data.frame(SESSO = 'M', ASL_RESIDENZA = 305, cluster = 4, eta_Min = 20)
for( i in 21:100)
  fittizio_c4_A305 = rbind(fittizio_c4_A305, c(SESSO = 'M', ASL_RESIDENZA = 305, cluster = 4, eta_Min = i))
fittizio_c4_A305$ASL_RESIDENZA = factor(fittizio_c4_A305$ASL_RESIDENZA)
fittizio_c4_A305$SESSO = factor(fittizio_c4_A305$SESSO)
fittizio_c4_A305$cluster= factor(fittizio_c4_A305$cluster)
fittizio_c4_A305$eta_Min = as.numeric(fittizio_c4_A305$eta_Min)
previsione_c4_A305 = predict(logreg, fittizio_c4_A305)
prob4 = exp(previsione_c4_A305)/(1+exp(previsione_c4_A305) )

# look for example at ASL 311
fittizio_c1_A311 = data.frame(SESSO = 'M', ASL_RESIDENZA = 311, cluster = 1, eta_Min = 20)
for( i in 21:100)
  fittizio_c1_A311 = rbind(fittizio_c1_A311, c(SESSO = 'M', ASL_RESIDENZA = 311, cluster = 1, eta_Min = i))
fittizio_c1_A311
fittizio_c1_A311$ASL_RESIDENZA = factor(fittizio_c1_A311$ASL_RESIDENZA)
fittizio_c1_A311$SESSO = factor(fittizio_c1_A311$SESSO)
fittizio_c1_A311$cluster= factor(fittizio_c1_A311$cluster)
fittizio_c1_A311$eta_Min = as.numeric(fittizio_c1_A311$eta_Min)
previsione_c1_A311 = predict(logreg, fittizio_c1_A311)
prob12 = exp(previsione_c1_A311)/(1+exp(previsione_c1_A311) )

fittizio_c2_A311 = data.frame(SESSO = 'M', ASL_RESIDENZA = 311, cluster = 2, eta_Min = 20)
for( i in 21:100)
  fittizio_c2_A311 = rbind(fittizio_c2_A311, c(SESSO = 'M', ASL_RESIDENZA = 311, cluster =2, eta_Min = i))
fittizio_c2_A311$ASL_RESIDENZA = factor(fittizio_c2_A311$ASL_RESIDENZA)
fittizio_c2_A311$SESSO = factor(fittizio_c2_A311$SESSO)
fittizio_c2_A311$cluster= factor(fittizio_c2_A311$cluster)
fittizio_c2_A311$eta_Min = as.numeric(fittizio_c2_A311$eta_Min)
previsione_c2_A311 = predict(logreg, fittizio_c2_A311)
prob22 = exp(previsione_c2_A311)/(1+exp(previsione_c2_A311) )

fittizio_c3_A311 = data.frame(SESSO = 'M', ASL_RESIDENZA = 311, cluster = 3, eta_Min = 20)
for( i in 21:100)
  fittizio_c3_A311 = rbind(fittizio_c3_A311, c(SESSO = 'M', ASL_RESIDENZA = 311, cluster = 3, eta_Min = i))
fittizio_c3_A311$ASL_RESIDENZA = factor(fittizio_c3_A311$ASL_RESIDENZA)
fittizio_c3_A311$SESSO = factor(fittizio_c3_A311$SESSO)
fittizio_c3_A311$cluster= factor(fittizio_c3_A311$cluster)
fittizio_c3_A311$eta_Min = as.numeric(fittizio_c3_A311$eta_Min)
previsione_c3_A311 = predict(logreg, fittizio_c3_A311)
prob32 = exp(previsione_c3_A311)/(1+exp(previsione_c3_A311) )

fittizio_c4_A311 = data.frame(SESSO = 'M', ASL_RESIDENZA = 311, cluster = 4, eta_Min = 20)
for( i in 21:100)
  fittizio_c4_A311 = rbind(fittizio_c4_A311, c(SESSO = 'M', ASL_RESIDENZA = 311, cluster = 4, eta_Min = i))
fittizio_c4_A311$ASL_RESIDENZA = factor(fittizio_c4_A311$ASL_RESIDENZA)
fittizio_c4_A311$SESSO = factor(fittizio_c4_A311$SESSO)
fittizio_c4_A311$cluster= factor(fittizio_c4_A311$cluster)
fittizio_c4_A311$eta_Min = as.numeric(fittizio_c4_A311$eta_Min)
previsione_c4_A311 = predict(logreg, fittizio_c4_A311)
prob42 = exp(previsione_c4_A311)/(1+exp(previsione_c4_A311) )


x11()
par(mfrow = c(1,2))
plot(20:100,prob1, type = 'l', col = '#NDA', xlab = 'Age'
     , ylab = 'Probability of death', ylim = c(0,1), main = '###', lwd=3)
points(20:100,prob2, type = 'l', col = '#NDA', lwd=3)
points(20:100,prob3, type = 'l', col = '#NDA', lwd=3)
points(20:100,prob4, type = 'l', col = '#NDA', lwd=3)
abline(h = 0, col = 'grey')
abline(h = 1, col = 'grey')
abline(v = 95)
abline(h=prob1[76],col='#NDA',lty=2)
abline(h=prob2[76],col='#NDA',lty=2)
abline(h=prob3[76],col='#NDA',lty=2)
abline(h=prob4[76],col='#NDA',lty=2)
legend("bottomright", legend=c('Cluster 1',
                               'Cluster 2',
                               'Cluster 3',
                               'Cluster 4'), fill=c('#...'), cex=.7)
plot(20:100,prob12,type = 'l', col = '', xlab = 'Age'
     , ylab = 'Probability of death', ylim = c(0,1), main = '###', lwd=3)
points(20:100,prob22, type = 'l', col = '#NDA', lwd=3)
points(20:100,prob32, type = 'l', col = '#NDA', lwd=3)
points(20:100,prob42, type = 'l', col = '#NDA', lwd=3)
abline(h = 0, col = 'grey')
abline(h = 1, col = 'grey')
abline(v = 95)
abline(h=prob12[76],col='#NDA',lty=2)
abline(h=prob22[76],col='#NDA',lty=2)
abline(h=prob32[76],col='#NDA',lty=2)
abline(h=prob42[76],col='#NDA',lty=2)
legend("bottomright", legend=c('Cluster 1',
                               'Cluster 2',
                               'Cluster 3',
                               'Cluster 4'), fill=c(#...), cex=.7)

#####Tabella con le probabilità di morte rispetto a ogni ASL####
####30 ANNI
#CLUSTER 1
fittizio_c1_M_30 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 1, eta_Min = 30)
for( i in 1:14)
  fittizio_c1_M_30= rbind(fittizio_c1_M_30, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 1, eta_Min = 30))
fittizio_c1_M_30
fittizio_c1_M_30$ASL_RESIDENZA = factor(fittizio_c1_M_30$ASL_RESIDENZA)
fittizio_c1_M_30$SESSO = factor(fittizio_c1_M_30$SESSO)
fittizio_c1_M_30$cluster= factor(fittizio_c1_M_30$cluster)
fittizio_c1_M_30$eta_Min = as.numeric(fittizio_c1_M_30$eta_Min)
previsione_c1_M_30 = predict(logreg, fittizio_c1_M_30)
prob1_M_30 = exp(previsione_c1_M_30)/(1+exp(previsione_c1_M_30) )

#CLUSTER 2
fittizio_c2_M_30 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 2, eta_Min = 30)
for( i in 1:14)
  fittizio_c2_M_30= rbind(fittizio_c2_M_30, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 2, eta_Min = 30))
fittizio_c2_M_30
fittizio_c2_M_30$ASL_RESIDENZA = factor(fittizio_c2_M_30$ASL_RESIDENZA)
fittizio_c2_M_30$SESSO = factor(fittizio_c2_M_30$SESSO)
fittizio_c2_M_30$cluster= factor(fittizio_c2_M_30$cluster)
fittizio_c2_M_30$eta_Min = as.numeric(fittizio_c2_M_30$eta_Min)
previsione_c2_M_30 = predict(logreg, fittizio_c2_M_30)
prob2_M_30 = exp(previsione_c2_M_30)/(1+exp(previsione_c2_M_30) )

#CLUSTER 3
fittizio_c3_M_30 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 3, eta_Min = 30)
for( i in 1:14)
  fittizio_c3_M_30= rbind(fittizio_c3_M_30, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 3, eta_Min = 30))
fittizio_c3_M_30
fittizio_c3_M_30$ASL_RESIDENZA = factor(fittizio_c3_M_30$ASL_RESIDENZA)
fittizio_c3_M_30$SESSO = factor(fittizio_c3_M_30$SESSO)
fittizio_c3_M_30$cluster= factor(fittizio_c3_M_30$cluster)
fittizio_c3_M_30$eta_Min = as.numeric(fittizio_c3_M_30$eta_Min)
previsione_c3_M_30 = predict(logreg, fittizio_c3_M_30)
prob3_M_30 = exp(previsione_c3_M_30)/(1+exp(previsione_c3_M_30) )

#CLUSTER 4
fittizio_c4_M_30 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 4, eta_Min = 30)
for( i in 1:14)
  fittizio_c4_M_30= rbind(fittizio_c4_M_30, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 4, eta_Min = 30))
fittizio_c4_M_30
fittizio_c4_M_30$ASL_RESIDENZA = factor(fittizio_c4_M_30$ASL_RESIDENZA)
fittizio_c4_M_30$SESSO = factor(fittizio_c4_M_30$SESSO)
fittizio_c4_M_30$cluster= factor(fittizio_c4_M_30$cluster)
fittizio_c4_M_30$eta_Min = as.numeric(fittizio_c4_M_30$eta_Min)
previsione_c4_M_30 = predict(logreg, fittizio_c4_M_30)
prob4_M_30 = exp(previsione_c4_M_30)/(1+exp(previsione_c4_M_30) )


####60 ANNI
#CLUSTER 1
fittizio_c1_M_60 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 1, eta_Min = 60)
for( i in 1:14)
  fittizio_c1_M_60= rbind(fittizio_c1_M_60, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 1, eta_Min = 60))
fittizio_c1_M_60
fittizio_c1_M_60$ASL_RESIDENZA = factor(fittizio_c1_M_60$ASL_RESIDENZA)
fittizio_c1_M_60$SESSO = factor(fittizio_c1_M_60$SESSO)
fittizio_c1_M_60$cluster= factor(fittizio_c1_M_60$cluster)
fittizio_c1_M_60$eta_Min = as.numeric(fittizio_c1_M_60$eta_Min)
previsione_c1_M_60 = predict(logreg, fittizio_c1_M_60)
prob1_M_60 = exp(previsione_c1_M_60)/(1+exp(previsione_c1_M_60) )

#CLUSTER 2
fittizio_c2_M_60 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 2, eta_Min = 60)
for( i in 1:14)
  fittizio_c2_M_60= rbind(fittizio_c2_M_60, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 2, eta_Min = 60))
fittizio_c2_M_60
fittizio_c2_M_60$ASL_RESIDENZA = factor(fittizio_c2_M_60$ASL_RESIDENZA)
fittizio_c2_M_60$SESSO = factor(fittizio_c2_M_60$SESSO)
fittizio_c2_M_60$cluster= factor(fittizio_c2_M_60$cluster)
fittizio_c2_M_60$eta_Min = as.numeric(fittizio_c2_M_60$eta_Min)
previsione_c2_M_60 = predict(logreg, fittizio_c2_M_60)
prob2_M_60 = exp(previsione_c2_M_60)/(1+exp(previsione_c2_M_60) )

#CLUSTER 3
fittizio_c3_M_60 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 3, eta_Min = 60)
for( i in 1:14)
  fittizio_c3_M_60= rbind(fittizio_c3_M_60, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 3, eta_Min = 60))
fittizio_c3_M_60
fittizio_c3_M_60$ASL_RESIDENZA = factor(fittizio_c3_M_60$ASL_RESIDENZA)
fittizio_c3_M_60$SESSO = factor(fittizio_c3_M_60$SESSO)
fittizio_c3_M_60$cluster= factor(fittizio_c3_M_60$cluster)
fittizio_c3_M_60$eta_Min = as.numeric(fittizio_c3_M_60$eta_Min)
previsione_c3_M_60 = predict(logreg, fittizio_c3_M_60)
prob3_M_60 = exp(previsione_c3_M_60)/(1+exp(previsione_c3_M_60) )

#CLUSTER 4
fittizio_c4_M_60 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 4, eta_Min = 60)
for( i in 1:14)
  fittizio_c4_M_60= rbind(fittizio_c4_M_60, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 4, eta_Min = 60))
fittizio_c4_M_60
fittizio_c4_M_60$ASL_RESIDENZA = factor(fittizio_c4_M_60$ASL_RESIDENZA)
fittizio_c4_M_60$SESSO = factor(fittizio_c4_M_60$SESSO)
fittizio_c4_M_60$cluster= factor(fittizio_c4_M_60$cluster)
fittizio_c4_M_60$eta_Min = as.numeric(fittizio_c4_M_60$eta_Min)
previsione_c4_M_60 = predict(logreg, fittizio_c4_M_60)
prob4_M_60 = exp(previsione_c4_M_60)/(1+exp(previsione_c4_M_60) )

####80 ANNI
#CLUSTER 1
fittizio_c1_M_80 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 1, eta_Min = 80)
for( i in 1:14)
  fittizio_c1_M_80= rbind(fittizio_c1_M_80, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 1, eta_Min = 80))
fittizio_c1_M_80
fittizio_c1_M_80$ASL_RESIDENZA = factor(fittizio_c1_M_80$ASL_RESIDENZA)
fittizio_c1_M_80$SESSO = factor(fittizio_c1_M_80$SESSO)
fittizio_c1_M_80$cluster= factor(fittizio_c1_M_80$cluster)
fittizio_c1_M_80$eta_Min = as.numeric(fittizio_c1_M_80$eta_Min)
previsione_c1_M_80 = predict(logreg, fittizio_c1_M_80)
prob1_M_80 = exp(previsione_c1_M_80)/(1+exp(previsione_c1_M_80) )

#CLUSTER 2
fittizio_c2_M_80 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 2, eta_Min = 80)
for( i in 1:14)
  fittizio_c2_M_80= rbind(fittizio_c2_M_80, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 2, eta_Min = 80))
fittizio_c2_M_80
fittizio_c2_M_80$ASL_RESIDENZA = factor(fittizio_c2_M_80$ASL_RESIDENZA)
fittizio_c2_M_80$SESSO = factor(fittizio_c2_M_80$SESSO)
fittizio_c2_M_80$cluster= factor(fittizio_c2_M_80$cluster)
fittizio_c2_M_80$eta_Min = as.numeric(fittizio_c2_M_80$eta_Min)
previsione_c2_M_80 = predict(logreg, fittizio_c2_M_80)
prob2_M_80 = exp(previsione_c2_M_80)/(1+exp(previsione_c2_M_80) )

#CLUSTER 3
fittizio_c3_M_80 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 3, eta_Min = 80)
for( i in 1:14)
  fittizio_c3_M_80= rbind(fittizio_c3_M_80, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 3, eta_Min = 80))
fittizio_c3_M_80
fittizio_c3_M_80$ASL_RESIDENZA = factor(fittizio_c3_M_80$ASL_RESIDENZA)
fittizio_c3_M_80$SESSO = factor(fittizio_c3_M_80$SESSO)
fittizio_c3_M_80$cluster= factor(fittizio_c3_M_80$cluster)
fittizio_c3_M_80$eta_Min = as.numeric(fittizio_c3_M_80$eta_Min)
previsione_c3_M_80 = predict(logreg, fittizio_c3_M_80)
prob3_M_80 = exp(previsione_c3_M_80)/(1+exp(previsione_c3_M_80) )

#CLUSTER 4
fittizio_c4_M_80 = data.frame(SESSO = 'M', ASL_RESIDENZA = 301, cluster = 4, eta_Min = 80)
for( i in 1:14)
  fittizio_c4_M_80= rbind(fittizio_c4_M_80, c(SESSO = 'M', ASL_RESIDENZA = 301+i, cluster = 4, eta_Min = 80))
fittizio_c4_M_80
fittizio_c4_M_80$ASL_RESIDENZA = factor(fittizio_c4_M_80$ASL_RESIDENZA)
fittizio_c4_M_80$SESSO = factor(fittizio_c4_M_80$SESSO)
fittizio_c4_M_80$cluster= factor(fittizio_c4_M_80$cluster)
fittizio_c4_M_80$eta_Min = as.numeric(fittizio_c4_M_80$eta_Min)
previsione_c4_M_80 = predict(logreg, fittizio_c4_M_80)
prob4_M_80 = exp(previsione_c4_M_80)/(1+exp(previsione_c4_M_80) )
# 4° cluster Maschi 80 anni


#Error on test data

label_stimati_reg = predict(logreg, test)
prob_stimata_reg = exp(label_stimati_reg)/(1+exp(label_stimati_reg) )
test$label_stimato_reg = ifelse(prob_stimata_reg > 0.5, 1, 0)
misc = table('label reali' = test$labelOUT, 'label stimati' = test$label_stimato_reg)
misc
ERR <- 0
for(g in 1:2){
  ERR <- ERR + sum(misc[g,-g])
}
ERR/dim(test)[1] # sembra aggirarsi intorno al 30/35%
# 0.34365

##### ROC curve
# grafico usato per valutare i classificatori binari (la nostra regressione logistica)
# come valutarlo è spiegato in fondo

# Ottenere le previsioni della probabilità di morte dal modello di regressione logistica
predictions <- predict(logreg, newdata = test, type = "response")

# Creare la curva ROC
roc_obj <- roc(test$labelOUT, predictions, curve = T) 

# Visualizzare la curva ROC
x11()
plot(roc_obj, main = "ROC curve", print.auc = T)
abline(h = 0, v = 1, col = 'grey')
abline(h = 1, v = 0, col = 'grey')

# x11()
# plot(roc_obj, main = "ROC  curve", xlab = "False positive rate",
#      ylab = "True positive rate", print.thres = "best", print.auc = TRUE)
# abline(h = 0, v = 1, col = 'grey')
# abline(h = 1, v = 0, col = 'grey')
# Calcolare l'AUC (area sotto la curva)
auc <- auc(roc_obj)
print(auc)

# Valutazione del modello
# - l'area sotto la curva ROC è data dall'AUC, se è sopra 0.5 si può ritenere un buon modello


##### Mixed effects Linear Models
PVRE <- function(lmmodel) {
  sigma2_eps <- as.numeric(get_variance_residual(lmmodel))
  sigma2_b <- as.numeric(get_variance_random(lmmodel))
  PVRE <- sigma2_b/(sigma2_b+sigma2_eps)
  PVRE
}

# Provo a fare 3 tipi di modelli:
# - Modello 1: labelOUT ~ eta + sesso + cluster + (1|strutt_id)
# - Modello 2: labelOUT ~ eta + sesso + cluster + (1|ASL_RESIDENZA)
# - Modello 3: labelOUT ~ eta + sesso + cluster + (1|strutt_id/ASL_RESIDENZA)
# Per vedere quale è il migliore in termini di PVRE (quale mi spiega meglio la variabilità)

# Modello 1:
fm16.1mer1 <- lmer(labelOUT ~ eta_Min + SESSO + cluster + (1|strutt_id),
                   data = data_final)
summary(fm16.1mer1)

# PVRE modello 1
PVRE(fm16.1mer1) # 0.02594834

# Modello 2:
fm16.1mer2 <- lmer(labelOUT ~ eta_Min + SESSO + cluster + (1|ASL_RESIDENZA),
                   data = data_final)
summary(fm16.1mer2)

# PVRE modello 2
PVRE(fm16.1mer2) # 0.0005681652 -> la ASL ci spiega molto meno quindi scarterei questo modello

# Modello 3:
fm16.1mer3 <- lmer(labelOUT ~ eta_Min + SESSO + cluster + (1|strutt_id/ASL_RESIDENZA),
                   data = data_final)
summary(fm16.1mer3)

# PVRE modello 3
PVRE(fm16.1mer3) # 0.02819235 -> simile a quello che prende in considerazione solo gli ospedali

# L'idea è quella di scegliere il modello più semplice perché anche se è meno performante 
# ha una robustezza maggiore (giustificarlo in vari modi) 

# Quindi alla fine scelgo il modello 1

# Una volta scelto il modello aggiungo la slope per cluster così posso vedere come 
# l’ospedale agisce sui diversi clusters di pazienti senza dover fare 4 modelli differenti.

# Provo a fare LMM annidati: trattare cluster come covariata random per vedere qual
# è l’effetto di essere in quel cluster all’interno di quell’ospedale rispetto all’effetto 
# medio che ha essere in quel cluster (e questo è dato dal cluster anche fuori dalla random intercept)
# sulla probabilità di morte -> può essere che la probabilità di morte aumenta/diminuisce rispetto
# all’andamento medio.

# Se nella distribuzione delle slope a livello ospedale vediamo che si distanziano dallo 
# zero identifichiamo gli ospedali che hanno un effetto nettamente 
# positivo/negativo per quel cluster di pazienti 
# (prendere magari quelli con IC nettamente diverso da zero).


n = dim(data_final)[1]
data_final <- data_final %>% arrange(sample(n))
train = data_final[1:(n-40000),]
test = data_final[-c(1:(n-40000)),]



fm16.2mer <- lmer(labelOUT ~ eta_Min + cluster + SESSO + (1+cluster|strutt_id),
                  data = data_final, control=lmerControl(optimizer="bobyqa",
                                                         optCtrl=list(maxfun=2e5)))
summary(fm16.2mer)

PVRE(fm16.2mer) # 0.0266838



x11()
dotplot(ranef(fm16.2mer, condVar=T))
# confidence interval for each of b0i -> faraway from zero means that is significant

# Validation sul modello:
label_stimati <- predict(fm16.2mer, newdata=test, allow.new.levels = T)
length(label_stimati)
sum(label_stimati < 0) # probabilità < 0? Ignorare
test$label_stimato = ifelse(label_stimati > 0.5, 1, 0)

misc = table('label reali' = test$labelOUT, 'label stimati' = test$label_stimato)
misc
ERR <- 0
for(g in 1:2){
  ERR <- ERR + sum(misc[g,-g])
}
ERR/dim(test)[1] # 0.3381


# The dotplot shows the point and interval estimates for the random effects, 
# ordering them and highlighting which are significantly different from the mean (0)

# IDEA: per ogni parametro (quindi cluster1, cluster2, cluster3, cluster4) seleziono gli intervalli
# di confidenza che non contengono lo zero.
# Questo mi permette di vedere quali ospedali hanno effetto positivo/negativo rispetto ai vari cluster
# quindi in pratica sto vedendo quale ospedale cura meglio pazienti di un determinato cluster

# computazionalmente ci mette una vita e mezza
# confintervals <- confint(fm16.2mer,oldNames=TRUE)

# Codice per selezionare magari il primo e l'ultimo ospedale per ogni cluster

intercepts <- data.frame(ranef(fm16.2mer, condVar=T)$strutt_id)

# CLUSTER 1
intercepts[which(intercepts$X.Intercept.==min(intercepts$X.Intercept.)),]
intercepts[which(intercepts$X.Intercept.==max(intercepts$X.Intercept.)),]
# in this case the mininum is NDA with an intercept of -0.1700568 
# and the maximum is NDA  with an intercept of 0.2552814
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA


# CLUSTER 2
intercepts[which(intercepts$cluster2==min(intercepts$cluster2)),]
intercepts[which(intercepts$cluster2==max(intercepts$cluster2)),]
# in this case the mininum is NDA  with an intercept of -0.03033882  
# and the maximum is NDA  with an intercept of 0.01871143  
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA


# CLUSTER 3
intercepts[which(intercepts$cluster3==min(intercepts$cluster3)),]
intercepts[which(intercepts$cluster3==max(intercepts$cluster3)),]
# in this case the mininum is NDA  with an intercept of -0.05051491   
# and the maximum is NDA  with an intercept of 0.0376465   
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA


# CLUSTER 4
intercepts[which(intercepts$cluster4==min(intercepts$cluster4)),]
intercepts[which(intercepts$cluster4==max(intercepts$cluster4)),]
# in this case the mininum is NDA with an intercept of -0.02753709  
# and the maximum is NDA with an intercept of 0.03115628  
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA
data_final[which(data_final$strutt_id == 'NDA'),4] # ASL NDA

# Idea: 4 cartine della lombardia con i 4 cluster e asl migliore in verde e peggiore in rosso
# scrivere che gli ospedali dei dati sono anonimi e si può risalire solo all'asl
# la cosa importante che emerge da questa analisi è che abbiamo trovato che esistono ospedali
# in grado di curare meglio gli hearth failure e che potrebbero essere presi come modello dagli altri
# per migliorare la qualità del servizio complessivo


##### A FEW TABLE FOR THE PRESENTATION
data_final <- read.csv('data_final.csv',header=T)
data_final <- subset(data_final, labelOUT != 'PERSO' ) # togliamo i pazienti persi
data_final <- dplyr::filter(data_final, tipo_prest == 41)   # teniamo solo pazienti ospedalizzati
data_final$labelOUT<-ifelse(data_final$labelOUT == "DECEDUTO", 1, 0)
data_final = data_final[,-1]

data_final$SESSO = factor(data_final$SESSO)
data_final$SESSO_bin = ifelse(data_final$SESSO == "M", 1, 0)
data_final$ASL_RESIDENZA = factor(data_final$ASL_RESIDENZA)

patologie <- data_final[, c(15,16,20,21,22,25,26,28,32)] # solo variabili più rilevanti
# K-means per variabili dicotomiche
set.seed(1)
ss2<-kmodes(patologie, modes=4, iter.max = 30, weighted = FALSE, fast = TRUE)
data_final$cluster <- factor(ss2$cluster) 

n = dim(data_final)[1]
data_final <- data_final %>% arrange(sample(n))
train = data_final[1:(n-40000),]
test = data_final[-c(1:(n-40000)),]

logreg <- glm(labelOUT ~ SESSO + cluster + ASL_RESIDENZA  +  eta_Min
              + eta_Min:cluster + eta_Min:ASL_RESIDENZA + SESSO:cluster,
              data = train, family='binomial')

prob_stimata <- function(sex, clust, eta) {
  fittizio_c1_M_30 = data.frame(SESSO = sex, ASL_RESIDENZA = 301, cluster = clust, eta_Min = eta)
  for( i in 1:14)
    fittizio_c1_M_30= rbind(fittizio_c1_M_30, c(SESSO = sex, ASL_RESIDENZA = 301+i, cluster = clust, eta_Min = eta))
  fittizio_c1_M_30
  fittizio_c1_M_30$ASL_RESIDENZA = factor(fittizio_c1_M_30$ASL_RESIDENZA)
  fittizio_c1_M_30$SESSO = factor(fittizio_c1_M_30$SESSO)
  fittizio_c1_M_30$cluster= factor(fittizio_c1_M_30$cluster)
  fittizio_c1_M_30$eta_Min = as.numeric(fittizio_c1_M_30$eta_Min)
  previsione_c1_M_30 = predict(logreg, fittizio_c1_M_30)
  prob1_M_30 = exp(previsione_c1_M_30)/(1+exp(previsione_c1_M_30) )
  prob1F30 <- 0
  for(i in 1:15){
    prob1F30 <- prob1F30 + prob1_M_30[i]*table(data_final$ASL_RESIDENZA)[i]
  }
  prob1F30 <- prob1F30/sum(table(data_final$ASL_RESIDENZA))
  prob1F30
}

## FEMMINE
# 30 ANNI
prob1F30 <- prob_stimata('F', 1, 30)
prob2F30 <- prob_stimata('F', 2, 30)
prob3F30 <- prob_stimata('F', 3, 30)
prob4F30 <- prob_stimata('F', 4, 30)

# 60 ANNI
prob1F60 <- prob_stimata('F', 1, 60)
prob2F60 <- prob_stimata('F', 2, 60)
prob3F60 <- prob_stimata('F', 3, 60)
prob4F60 <- prob_stimata('F', 4, 60)

# 80 ANNI
prob1F80 <- prob_stimata('F', 1, 80)
prob2F80 <- prob_stimata('F', 2, 80)
prob3F80 <- prob_stimata('F', 3, 80)
prob4F80 <- prob_stimata('F', 4, 80)

## MASCHI
# 30 ANNI
prob1M30 <- prob_stimata('M', 1, 30)
prob2M30 <- prob_stimata('M', 2, 30)
prob3M30 <- prob_stimata('M', 3, 30)
prob4M30 <- prob_stimata('M', 4, 30)

# 60 ANNI
prob1M60 <- prob_stimata('M', 1, 60)
prob2M60 <- prob_stimata('M', 2, 60)
prob3M60 <- prob_stimata('M', 3, 60)
prob4M60 <- prob_stimata('M', 4, 60)

# 80 ANNI
prob1M80 <- prob_stimata('M', 1, 80)
prob2M80 <- prob_stimata('M', 2, 80)
prob3M80 <- prob_stimata('M', 3, 80)
prob4M80 <- prob_stimata('M', 4, 80)


# femmine: '#aa5b4d'
# maschi: '#0f64a5'



x11()
par(mfrow=c(2,2))
plot(c(30,60,80),rbind(1-prob1F30,1-prob1F60,1-prob1F80), type='l',xlim=c(20,90), ylim=c(0,1),
     col='', xlab = 'Age', ylab = 'Survival probability', lwd = 3)
lines(c(30, 60, 80), rbind(1 - prob1M30, 1 - prob1M60, 1 - prob1M80),
      col = '', lwd = 3)
abline(v=30, lty=2)
abline(v=60, lty=2)
abline(v=80, lty=2)
abline(h=0.5, lty=2, col='')
legend('bottomleft', legend=c('Female', 'Male'), lwd = 3, col=c('',''))
title('Cluster 1')
#title(main="How the probability of survival evolves with respect to the patient's age")
#legend(x = 'left', legend = factor(levels(data_final$ASL_RESIDENZA)), col = my_rainbow(15), lwd = 2)

plot(c(30,60,80),rbind(1-prob2F30,1-prob2F60,1-prob2F80), type='l',xlim=c(20,90), ylim=c(0,1),
     col='', xlab = 'Age', ylab = 'Survival probability', lwd = 3)
lines(c(30, 60, 80), rbind(1 - prob2M30, 1 - prob2M60, 1 - prob2M80),
      col = '', lwd = 3)
abline(v=30, lty=2)
abline(v=60, lty=2)
abline(v=80, lty=2)
abline(h=0.5, lty=2, col='')
legend('bottomleft', legend=c('Female', 'Male'), lwd = 3, col=c('',''))
title('Cluster 2')

plot(c(30,60,80),rbind(1-prob3F30,1-prob3F60,1-prob3F80), type='l',xlim=c(20,90), ylim=c(0,1),
     col='', xlab = 'Age', ylab = 'Survival probability', lwd = 3)
lines(c(30, 60, 80), rbind(1 - prob3M30, 1 - prob3M60, 1 - prob3M80),
      col = '', lwd = 3)
abline(v=30, lty=2)
abline(v=60, lty=2)
abline(v=80, lty=2)
abline(h=0.5, lty=2, col='')
legend('bottomleft', legend=c('Female', 'Male'), lwd = 3, col=c('',''))
title('Cluster 3')

plot(c(30,60,80),rbind(1-prob4F30,1-prob4F60,1-prob4F80), type='l',xlim=c(20,90), ylim=c(0,1),
     col='', xlab = 'Age', ylab = 'Survival probability', lwd = 3)
lines(c(30, 60, 80), rbind(1 - prob4M30, 1 - prob4M60, 1 - prob4M80),
      col = '', lwd = 3)
abline(v=30, lty=2)
abline(v=60, lty=2)
abline(v=80, lty=2)
abline(h=0.5, lty=2, col='#7d3253')
legend('bottomleft', legend=c('Female', 'Male'), lwd = 3, col=c('',''))
title('Cluster 4')

# percentuale maschi e femmine in ogni cluster
cl1 <- data_final[which(data_final$cluster==1),]
cl2 <- data_final[which(data_final$cluster==2),]
cl3 <- data_final[which(data_final$cluster==3),]
cl4 <- data_final[which(data_final$cluster==4),]

table(cl1$SESSO) # -> 0.4823151 di femmine, 0.5176849 di maschi
# F    M 
# 9000 9660 
table(cl2$SESSO) # -> 0.5437026 di femmine, 0.4562974 di maschi
# F     M 
# 19551 16408 
table(cl3$SESSO) # -> 0.4912289 di femmine, 0.5087711 di maschi
# F     M 
# 31559 32686 
table(cl4$SESSO) # -> 0.4406531 di femmine, 0.5593469 di maschi
# F    M 
# 5101 6475 






