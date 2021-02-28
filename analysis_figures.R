#title: Analysis and figures for Knight et al. 2017 Recommendations for acoustic recognizer performance assessment with application to five common automated signal recognition programs. Avian Conservation and Ecology 12(2):14
#author: Elly C. Knight
#date: Feb 24, 2018

library(tidyverse)
library(stringi)
library(ggplot2)
library(sciplot)
library(reshape2)
library(splitstackshape)
library(lme4)
library(gridExtra)
library(grid)
library(pROC)
library(PRROC)
library(ROCR)
library(survival)
library(AICcmodavg)
library(ResourceSelection)
library(MuMIn)
library(unmarked)
library(readxl)
library(RColorBrewer)

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

score <- rep(seq(from=0.00, to=0.99, by=0.01), 6)
recognizer <- append(rep("CNN", 100), rep("Kaleidoscope", 100)) %>% 
  append(rep("MonitoR", 100)) %>% 
  append(rep("RavenPro", 100)) %>% 
  append(rep("SongScope", 100)) %>% 
  append(rep("Human", 100))

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="center"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

gg_color_hue <- function(n) {
  hues = seq(0, 360, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 6
cols = gg_color_hue(n)
cols

#"#F8766D" "#B79F00" "#00BA38" "#00BFC4" "#619CFF" "#F564E3"

col.h <- "#FF6C91"
col.c <- "#CD9600"
col.k <- "#49B500"
col.m <- "#00C1A9"
col.r <- "#00A9FF" 
col.s <- "#E36EF6"

lines <- c(1,2,5,4,3,6)

#DATA WRANGLING----

#1. Get list of files----
setwd("/Users/ellyknight/Documents/UoA/Projects/RecognizerVariation/Processing/ONData/")

files <- as.data.frame(list.files())
colnames(files) <- c("ID")
files$ID <- as.character((files$ID))
files$site <- stri_sub(files$ID, from=1, length=6)
files$month <- stri_sub(files$ID, from=12, length=2)
files$day <- stri_sub(files$ID, from=14, length=2)
#files$hour <- stri_sub(files$ID, from=-10, length=2)
#files$min <- stri_sub(files$ID, from=-8, length=2)
files$site <- as.factor(files$site)


#2. SongScope ----
setwd("/Users/ellyknight/Documents/UoA/Projects/RecognizerVariation/Analysis/")

s <- read.table("ON-CONIPeent-SS_20_0_results_validated.txt", header=FALSE, sep="\t")
colnames(s) <- c("filepath", "start", "duration", "level", "quality", "score", "template", "CONI")
s$recognizer <- "SongScope"
s$min <- as.numeric(stri_sub(s$start, from=2, length=1))
s$sec <- as.numeric(stri_sub(s$start, from=4, length=6))
s$score.r <- (s$score-min(s$score))/(max(s$score)-min(s$score))

#3. CNN----
c <- read.table("ON_CONIPeent-CNN_0.001_0.1_results_validated.txt", header=FALSE, sep="\t")
colnames(c) <- c("filepath", "start", "duration", "level", "quality", "score", "template", "CONI")
c$recognizer <- "CNN"
c$min <- as.numeric(stri_sub(c$start, from=5, length=1))
c$sec <- as.numeric(stri_sub(c$start, from=7, length=6))
c$score.r <- (c$score-min(c$score))/(max(c$score)-min(c$score))

#4. MonitoR----
m <- read.table("ON_CONIMonitoR0.1_100_0.1_results_validated.txt", header=FALSE, sep="\t")
colnames(m) <- c("filepath", "start", "duration", "level", "quality", "score", "template", "CONI")
m$recognizer <- "MonitoR"
m$min <- as.numeric(stri_sub(m$start, from=5, length=1))
m$sec <- as.numeric(stri_sub(m$start, from=7, length=6))
m$score.r <- (m$score-min(m$score))/(max(m$score)-min(m$score))

#5. Kaleidoscope----
k <- read.csv("ON_CONIPeent_BC100_cluster_results_validated.csv", header=TRUE)
k <- subset(k, TOP1MATCH.=="y")
k <- subset(k, CHANNEL==1)
k$start.min <- k$OFFSET/60 
k$sec <- round(as.numeric(stri_sub(k$start.min, from=2, length=8))*60,2)
k$sec[is.na(k$sec)] <- 0
k$min <- round((k$OFFSET-k$sec)/60)
k$min[is.na(k$min)] <- 0
k$CONI <- k$MANUAL.ID
k$CONI <- gsub("Noise", 0, k$CONI)
k$CONI <- gsub("y", 1, k$CONI)
k$CONI[is.na(k$CONI)] <- 0
k$score <- 2-k$TOP1DIST
k$score.r <- (k$score-min(k$score))/(max(k$score)-min(k$score))
k$recognizer <- "Kaleidoscope"
k$ID <- k$IN.FILE
k.1 <- dplyr::select(k, ID, recognizer, score, score.r, min, sec, CONI)

#6. RavenPro----
r <- read.csv("ON_CONIPeent-Raven_results_validated.csv")
r <- subset(r, use=="y")
r$recognizer <- "RavenPro"
r$start.min <- r$timediff_ms/60000
r$sec <- round(as.numeric(stri_sub(r$start.min, from=2, length=8))*60)
r$sec[is.na(r$sec)] <- 0
r$min <- round((r$timediff_ms/1000-r$sec)/60)
r$min[is.na(r$min)] <- 0
r$score.r <- (r$score-min(r$score))/(max(r$score)-min(r$score))
r$CONI <- r$targetspp
r$CONI <- gsub("y", 1, r$CONI)
r$CONI <- gsub("n", 0, r$CONI)
r$CONI <- as.numeric(r$CONI)
r.1 <- dplyr::select(r, ID, recognizer, score, score.r, min, sec, CONI)

#7. Human----
h.1 <- read.csv("ON_CONIPeent-Human_results_validated.csv", header=TRUE)
h.1$recognizer <- "Human"
h.1$start.min <- h.1$timediff_ms/60000
h.1$sec <- round(as.numeric(stri_sub(h.1$start.min, from=2, length=8))*60)
h.1$sec[is.na(h.1$sec)] <- 0
h.1$min <- round((h.1$timediff_ms/1000-h.1$sec)/60)
h.1$min[is.na(h.1$min)] <- 0
h.1$CONI <- 1
h.1$pres <- 1
h.1$score <- 1
h.1$score.r <- 1
h.1$ID <- h.1$filename
h <- dplyr::select(h.1, ID, recognizer, score, score.r, min, sec, CONI)

#8. Collate----

#Bind 3 SS validated dataframes together
all.1 <- rbind(s,c,m)
all.1$ID <- stri_sub(all.1$filepath, from=-34, length=34)
all.2 <- dplyr::select(all.1, ID, recognizer, score, score.r, min, sec, CONI)
all.2$CONI <- gsub("y", 1, all.2$CONI)
all.2$CONI <- gsub("n", 0, all.2$CONI)
all.2$CONI <- as.numeric(all.2$CONI)

#Bind all recognizer data together
all.3 <- rbind(all.2, k.1, r.1, h)
all.3$pres <- as.numeric(ifelse(all.3$CONI > 0, 1, 0))
all.3$recognizer <- as.factor(all.3$recognizer)
all.3$site <- stri_sub(all.3$ID, from=1, length=6)

#remove files that were not listener processed
all <- subset(all.3, ID != "SM4706_20120613_205500S_215500.wav")
all <- subset(all, ID != "SM4706_20120613_205500S_225500.wav")

#9. Histogram----
#all$pres <- as.factor(all$pres)
#all$pres <- factor(all$pres, c(0,1))

all.s <- all %>%
  filter(recognizer=="SongScope")
all.c <- all %>%
  filter(recognizer=="CNN")
all.r <- all %>%
  filter(recognizer=="RavenPro")
all.k <- all %>%
  filter(recognizer=="Kaleidoscope")
all.m <- all %>%
  filter(recognizer=="MonitoR")

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.position="none",
        axis.title=element_blank())

legend.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

hist.c <- ggplot(all.c) +
  geom_histogram(aes(x=score, fill=pres), bins=50) +
  scale_fill_manual(values=c("grey30", "grey70")) +
  annotate("text", x=0.5, y=1600, hjust=0.5, vjust=1, label="CNN", size=5) +
  my.theme

hist.k <- ggplot(all.k) +
  geom_histogram(aes(x=score.r, fill=pres), bins=50) +
  scale_fill_manual(values=c("grey30", "grey70")) +
  annotate("text", x=0.55, y=600, hjust=0.55, vjust=1, label="Kaleidoscope", size=5) +
  my.theme

hist.r <- ggplot(all.r) +
  geom_histogram(aes(x=score, fill=pres), bins=50) +
  scale_fill_manual(values=c("grey30", "grey70")) +
  annotate("text", x=0.6, y=900, hjust=0.5, vjust=1, label="RavenPro", size=5) +
  my.theme

hist.s <- ggplot(all.s) +
  geom_histogram(aes(x=score, fill=pres), bins=50) +
  scale_fill_manual(values=c("grey30", "grey70")) +
  annotate("text", x=50, y=650, hjust=0.5, vjust=1, label="Song Scope", size=5) +
  my.theme

hist.m <- ggplot(all.m) +
  geom_histogram(aes(x=score, group=pres, fill=pres), bins=50) +
  scale_fill_manual(values=c("grey30", "grey70"), labels=c("False positive", "True positive")) +
  annotate("text", x=11, y=1100, hjust=0.5, vjust=1, label="MonitoR", size=5) +
  labs(fill="") +
  legend.theme

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


legend <- g_legend(hist.m)

grid.arrange(arrangeGrob(hist.k, hist.r, hist.s, hist.c, hist.m + my.theme, legend,
             nrow=2, ncol=3),
             bottom=textGrob("Recognizer score", gp=gpar(fontsize=18)),
             left=textGrob("Recognizer hits", gp=gpar(fontsize=18), rot=90))

#10. Summary stats----
all$pres <- as.numeric(all$pres)

all.site <- all %>%
  dplyr::group_by(ID) %>% 
  dplyr::summarize(dets=sum(pres)) %>% 
  mutate(site=substr(ID, 1, 6))

all.rec.pres <- all.site %>% 
  mutate(pres = ifelse(dets> 0, 1, 0)) %>% 
  dplyr::summarize(pres = sum(pres))

all.site.pres <- all.site %>% 
  group_by(site) %>% 
  summarize(dets=sum(dets)) %>% 
  mutate(pres=ifelse(dets>0,1,0)) %>% 
  dplyr::summarize(pres=sum(pres))
all.site.pres

mean(all.site$dets)
sd(all.site$dets)

#OVERALL PERFORMANCE----
all$pres <- as.numeric(all$pres)

#1. Gold standard----
gold <- all %>% 
  dplyr::group_by(ID, recognizer) %>% 
  dplyr::summarize(hit=n(), det=sum(pres)) %>% 
  right_join(files, by="ID") %>% 
  dplyr::select(ID, recognizer, det) %>% 
  spread(key=recognizer, value=det) %>% 
  mutate(CNN = ifelse(is.na(CNN),0,CNN)) %>% 
  mutate(Human = ifelse(is.na(Human),0,Human)) %>% 
  mutate(Kaleidoscope = ifelse(is.na(Kaleidoscope),0,Kaleidoscope)) %>% 
  mutate(MonitoR = ifelse(is.na(MonitoR),0,MonitoR)) %>% 
  mutate(RavenPro = ifelse(is.na(RavenPro),0,RavenPro)) %>% 
  mutate(SongScope = ifelse(is.na(SongScope),0,SongScope)) %>% 
  mutate(gold = pmax(CNN, Human, Kaleidoscope, MonitoR, RavenPro, SongScope)) %>% 
  dplyr::select(ID, gold)

gold.all <- sum(gold$gold)
gold.all

#2. Recognizer data frame all----

score <- rep(seq(from=0.00, to=0.99, by=0.01), 6)
recognizer <- append(rep("CNN", 100), rep("Kaleidoscope", 100)) %>% 
  append(rep("MonitoR", 100)) %>% 
  append(rep("RavenPro", 100)) %>% 
  append(rep("SongScope", 100)) %>% 
  append(rep("Human", 100))

all.all = data.frame()
for(i in 1:length(score)){
  score.i <- paste0(score[i])
  recognizer.i <- paste0(recognizer[i])
  data.1 <- all %>%
    dplyr::filter(score.r > as.numeric(score.i), recognizer==recognizer.i) %>% 
    dplyr::summarize(det=sum(pres), hit=n()) %>% 
    mutate(r=(det)/(gold.all), p=det/hit, fp=(hit-det)/hit, tp=det/hit, fn=(gold.all-det)/gold.all, acc=(det-(hit-det))/hit) %>% 
    mutate(f=(2*p*r)/(p+r)) %>% 
    ungroup()
  data.1$recognizer <- as.factor(recognizer.i)
  data.1$threshold <- as.numeric(score.i)
  all.all <- rbind(all.all, data.1)
}

#3. Precision, Recall, F-score plot----

all.all$recognizer <- factor(all.all$recognizer, c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"))

pr.theme <- theme_classic() +
  theme(text=element_text(size=14, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0), size=18),
        axis.title.y=element_text(margin=margin(0,10,0,0), size=18),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.position="none")

f.theme <- theme_classic() +
  theme(text=element_text(size=14, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0), size=18),
        axis.title.y=element_text(margin=margin(0,10,0,0), size=18),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.justification=c(1,0), legend.position=c(1,0),
        legend.background = element_rect(fill=alpha('white', 0.01)),
        legend.key.width=unit(2, "cm"))

gg.p <- ggplot(all.all) +
  geom_line(aes(x=threshold, y=p, linetype=recognizer, colour=recognizer), lwd=1.2) +
  f.theme +
  labs(x="Score Threshold", y="Precision") +
  scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                      values=cols,
                      labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                        values=c(1,2,5,4,3,6),
                        labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
  ylim(0,1) + xlim(0,1)


gg.r <- ggplot(all.all, aes(group=recognizer)) +
  geom_line(aes(x=threshold, y=r, linetype=recognizer, colour=recognizer), lwd=1.2) +
  pr.theme +
  labs(x="Score Threshold", y="Recall") +
    scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                        values=cols,
                        labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                        values=c(1,2,5,4,3,6),
                        labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
  ylim(0,1) + xlim(0,1)


gg.f <- ggplot(all.all, aes(group=recognizer)) +
  geom_line(aes(x=threshold, y=f, linetype=recognizer, colour=recognizer), lwd=1.2) +
  pr.theme +
  labs(x="Score Threshold", y=expression(paste("F-score (",beta,"=1)"))) +
  scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                      values=cols,
                      labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                        values=c(1,2,5,4,3,6),
                        labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
  ylim(0,1) + xlim(0,1)


grid.arrange(gg.p, gg.r, gg.f, ncol=3)

#4. ROC & Precision recall----
all.s <- all %>%
  filter(recognizer=="SongScope")
fg.s <- subset(all.s, pres==1)
bg.s <- subset(all.s, pres==0)

all.c <- all %>%
  filter(recognizer=="CNN")
fg.c <- subset(all.c, pres==1)
bg.c <- subset(all.c, pres==0)

all.r <- all %>%
  filter(recognizer=="RavenPro")
fg.r <- subset(all.r, pres==1)
bg.r <- subset(all.r, pres==0)

all.k <- all %>%
  filter(recognizer=="Kaleidoscope")
fg.k <- subset(all.k, pres==1)
bg.k <- subset(all.k, pres==0)

all.m <- all %>%
  filter(recognizer=="MonitoR")
fg.m <- subset(all.m, pres==1)
bg.m <- subset(all.m, pres==0)

#PR predictions
pr.s <- pr.curve(scores.class0=fg.s$score.r, scores.class1=bg.s$score.r, sorted=FALSE, dg.compute=TRUE, curve=TRUE)
pr.s
pr.c <- pr.curve(scores.class0=fg.c$score.r, scores.class1=bg.c$score.r, sorted=FALSE, dg.compute=TRUE, curve=TRUE)
pr.c
pr.r <- pr.curve(scores.class0=fg.r$score.r, scores.class1=bg.r$score.r, sorted=FALSE, dg.compute=TRUE, curve=TRUE)
pr.r
pr.k <- pr.curve(scores.class0=fg.k$score.r, scores.class1=bg.k$score.r, sorted=FALSE, dg.compute=TRUE, curve=TRUE)
pr.k
pr.m <- pr.curve(scores.class0=fg.m$score.r, scores.class1=bg.m$score.r, sorted=FALSE, dg.compute=TRUE, curve=TRUE)
pr.m

#ROC predictions
roc.s <- roc.curve(scores.class0=fg.s$score.r, scores.class1=bg.s$score.r, sorted=FALSE, curve=TRUE)
roc.s
roc.c <- roc.curve(scores.class0=fg.c$score.r, scores.class1=bg.c$score.r, sorted=FALSE, curve=TRUE)
roc.c
roc.k <- roc.curve(scores.class0=fg.k$score.r, scores.class1=bg.k$score.r, sorted=FALSE, curve=TRUE)
roc.k
roc.r <- roc.curve(scores.class0=fg.r$score.r, scores.class1=bg.r$score.r, sorted=FALSE, curve=TRUE)
roc.r
roc.m <- roc.curve(scores.class0=fg.m$score.r, scores.class1=bg.m$score.r, sorted=FALSE, curve=TRUE)
roc.m

#PR curves
par(oma=c(0,0,0,0), mar=c(4.5,4.5,0.5,0.5), mfrow=c(1,2))
plot(pr.s, col=col.s, lwd=3, auc.main=FALSE, legend=FALSE, main="", cex.lab=1.5, axes=FALSE, lty=2)
axis(1)
axis(2)
plot(pr.c, col=col.c, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=3)
plot(pr.r, col=col.r, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=4)
plot(pr.k, col=col.k, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=5)
plot(pr.m, col=col.m, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=6)
legend.text <- c("CNN (0.94)", "Kaleidoscope (0.77)", "MonitoR (0.88)" , "RavenPro (0.82)", "SongScope (0.87)")
legend.lty <- c(2,3,4,5,6)
legend.colour <- c(col.c, col.k, col.m, col.r, col.s)
legend(0.4, 0.32, col = legend.colour, legend = legend.text, lty=legend.lty, lwd=3, title="Processing approach (AUC)", bty="o", cex=1.2, box.col="white")

#ROC curves
plot(roc.s, col=col.s, lwd=3, auc.main=FALSE, legend=FALSE, main="", xlab="False positive rate", ylab="True positive rate", cex.lab=1.5, axes=FALSE, lty=2)
axis(1)
axis(2)
plot(roc.c, col=col.c, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=3)
plot(roc.r, col=col.r, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=4)
plot(roc.k, col=col.k, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=5)
plot(roc.m, col=col.m, lwd=3, auc.main=FALSE, legend=FALSE, main="", add=TRUE, lty=6)
legend.text <- c("CNN (0.88)", "Kaleidoscope (0.53)", "MonitoR (0.78)" , "RavenPro (0.72)", "SongScope (0.90)")
legend.colour <- c(col.c, col.k, col.m, col.r, col.s)
legend(0.4, 0.32, col = legend.colour, legend = legend.text, lty=legend.lty, lwd=3, title="Processing approach (AUC)", bty="o", cex=1.2, box.col="white")

#PRESENCE ABSENCE----

#1. File list----
xtabs(~site+day, files)

id.pres <- c("SM4620_20120613_205500S_215500.wav", "SM4706_20120625_205800S_215800.wav", "SM4793_20120613_205500S_215500.wav", "SM4820_20120613_205500S_225500.wav", "SM4821_20120613_205500S_215500.wav", "SM4822_20120613_205500S_215500.wav", "SM4823_20120613_205500S_215500.wav", "SM4824_20120613_205500S_225500.wav", "SM4825_20120613_205500S_215500.wav", "SM4826_20120613_205500S_215500.wav", "SM4827_20120613_205500S_215500.wav", "SM4828_20120613_205500S_215500.wav", "SM4830_20120613_205500S_225500.wav", "SM4831_20120613_205500S_215500.wav", "SM4832_20120613_205500S_225500.wav", "SM4833_20120625_205800S_215800.wav", "SM4837_20120613_205500S_215500.wav", "SM4840_20120613_205500S_215500.wav", "SM4841_20120613_205500S_215500.wav", "SM4843_20120613_205500S_215500.wav", "SM4844_20120613_205500S_215500.wav", "SM5069_20120613_205500S_215500.wav", "SM5074_20120613_205500S_215500.wav", "SM5090_20120613_205500S_215500.wav", "SM5091_20120613_205500S_225500.wav", "SM7517_20120613_205500S_215500.wav", "SM7518_20120613_205500S_215500.wav", "SM7522_20120625_205800S_225800.wav", "SM7524_20120613_205500S_215500.wav", "SM7594_20120613_205500S_215500.wav", "SM7615_20120613_205500S_215500.wav", "SM7616_20120613_205500S_215500.wav", "SM7617_20120613_205500S_215500.wav", "SM7621_20120613_205500S_215500.wav", "SM7673_20120613_205500S_215500.wav", "SM7679_20120613_205500S_215500.wav", "SM7682_20120613_205500S_215500.wav", "SM7687_20120613_205500S_215500.wav", "SM7689_20120613_205500S_225500.wav", "SM7697_20120613_205500S_215500.wav", "SM7698_20120613_205500S_225500.wav", "SM7699_20120613_205500S_215500.wav", "SM7701_20120613_205500S_225500.wav", "SM7708_20120625_205800S_215800.wav", "SM7714_20120625_205800S_215800.wav")

files.pres <- files[files$ID %in% id.pres,]

#2. Gold standard----
gold.pres <- all %>% 
  group_by(ID, recognizer) %>% 
  dplyr::summarize(hit=n(), det=sum(pres)) %>% 
  right_join(files.pres, by="ID") %>% 
  dplyr::select(ID, recognizer, det) %>% 
  spread(key=recognizer, value=det) %>% 
  mutate(CNN = ifelse(is.na(CNN),0,CNN)) %>% 
  mutate(Human = ifelse(is.na(Human),0,Human)) %>% 
  mutate(Kaleidoscope = ifelse(is.na(Kaleidoscope),0,Kaleidoscope)) %>% 
  mutate(MonitoR = ifelse(is.na(MonitoR),0,MonitoR)) %>% 
  mutate(RavenPro = ifelse(is.na(RavenPro),0,RavenPro)) %>% 
  mutate(SongScope = ifelse(is.na(SongScope),0,SongScope)) %>% 
  mutate(gold = pmax(CNN, Human, Kaleidoscope, MonitoR, RavenPro, SongScope)) %>% 
  mutate(pres = ifelse(gold > 0, 1, 0)) %>% 
  dplyr::select(ID, gold, pres)

#3. Recognizer data frame by recording----

all.pres <- data.frame()
for(i in 1:length(score)){
  score.i <- paste0(score[i])
  recognizer.i <- paste0(recognizer[i])
  data.1 <- all %>%
    dplyr::filter(score.r > as.numeric(score.i), recognizer==recognizer.i) %>%
    group_by(ID) %>% 
    dplyr::summarize(det=sum(pres), hit=n()) %>% 
    right_join(gold.pres, by = "ID") %>% 
    mutate(det = ifelse(is.na(det), 0, det)) %>% 
    mutate(hit = ifelse(is.na(hit), 0, hit)) %>% 
    mutate(pres.rec = ifelse(det > 0, 1, 0)) %>% 
    ungroup()
  data.1$recognizer <- recognizer.i
  data.1$threshold <- as.numeric(score.i)
  all.pres <- rbind(all.pres, data.1)
}

all.pres$correct <- ifelse(all.pres$pres==all.pres$pres.rec, 1, 0)

all.pres.s <- all.pres %>%
  filter(recognizer=="SongScope")
all.pres.c <- all.pres %>%
  filter(recognizer=="CNN")
all.pres.r <- all.pres %>%
  filter(recognizer=="RavenPro")
all.pres.k <- all.pres %>%
  filter(recognizer=="Kaleidoscope")
all.pres.m <- all.pres %>%
  filter(recognizer=="MonitoR")
all.pres.h <- all.pres %>%
  filter(recognizer=="Human")

#4. GLM----

newdata.pres.1 <- with(all.pres.s, expand.grid(threshold=unique(threshold)))

#SongScope
fn.s.0 <- glm(correct ~ 1, data=all.pres.s, family="binomial")
fn.s.1 <- glm(correct ~ threshold, data=all.pres.s, family="binomial")
fn.s.2 <- glm(correct ~ threshold + I(threshold^2), data=all.pres.s, family="binomial")
fn.s.3 <- glm(correct ~ threshold + I(threshold^2) + I(threshold^3), data=all.pres.s, family="binomial")
anova(fn.s.3, test="Chisq")
fn.s.list <- c(fn.s.1, fn.s.2, fn.s.3)
model.sel(fn.s.0, fn.s.1, fn.s.2, fn.s.3, rank=AIC)
hoslem.test(all.pres.s$correct, fitted(fn.s.2))

newdata.pres <- as.data.frame(predict(fn.s.3, newdata.pres.1, se.fit=TRUE, type="response")) %>% 
  mutate(lwr=fit-1.96*se.fit, upr=fit+1.96*se.fit, recognizer="SongScope") %>% 
  dplyr::select(fit, lwr, upr, recognizer) %>% 
  cbind(newdata.pres.1)

#CNN
fn.c.0 <- glm(correct ~ 1, data=all.pres.c, family="binomial")
fn.c.1 <- glm(correct ~ threshold, data=all.pres.c, family="binomial")
fn.c.2 <- glm(correct ~ threshold + I(threshold^2), data=all.pres.c, family="binomial")
fn.c.3 <- glm(correct ~ threshold + I(threshold^2) + I(threshold^3), data=all.pres.c, family="binomial")
anova(fn.c.3, test="Chisq")
model.sel(fn.c.0, fn.c.1, fn.c.2, fn.c.3, rank=AIC)
hoslem.test(all.pres.s$correct, fitted(fn.c.3))

newdata.pres <- as.data.frame(predict(fn.c.3, newdata.pres.1, se.fit=TRUE, type="response")) %>% 
  mutate(lwr=fit-1.96*se.fit, upr=fit+1.96*se.fit, recognizer="CNN") %>% 
  dplyr::select(fit, lwr, upr, recognizer) %>% 
  cbind(newdata.pres.1) %>% 
  rbind(newdata.pres)

#RavenPro
fn.r.0 <- glm(correct ~ 1, data=all.pres.r, family="binomial")
fn.r.1 <- glm(correct ~ threshold, data=all.pres.r, family="binomial")
fn.r.2 <- glm(correct ~ threshold + I(threshold^2), data=all.pres.r, family="binomial")
fn.r.3 <- glm(correct ~ threshold + I(threshold^2) + I(threshold^3), data=all.pres.r, family="binomial")
model.sel(fn.r.0, fn.r.1, fn.r.2, fn.r.3, rank=AIC)
anova(fn.r.3, test="Chisq")

newdata.pres <- as.data.frame(predict(fn.r.0, newdata.pres.1, se.fit=TRUE, type="response")) %>% 
  mutate(lwr=fit-1.96*se.fit, upr=fit+1.96*se.fit, recognizer="RavenPro") %>% 
  dplyr::select(fit, lwr, upr, recognizer) %>% 
  cbind(newdata.pres.1) %>% 
  rbind(newdata.pres)

#Kaleidoscope
fn.k.0 <- glm(correct ~ 1, data=all.pres.k, family="binomial")
fn.k.1 <- glm(correct ~ threshold, data=all.pres.k, family="binomial")
fn.k.2 <- glm(correct ~ threshold + I(threshold^2), data=all.pres.k, family="binomial")
fn.k.3 <- glm(correct ~ threshold + I(threshold^2) + I(threshold^3), data=all.pres.k, family="binomial")
model.sel(fn.k.0, fn.k.1, fn.k.2, fn.k.3, rank=AIC)
anova(fn.k.3, test="Chisq")

newdata.pres <- as.data.frame(predict(fn.k.3, newdata.pres.1, se.fit=TRUE, type="response")) %>% 
  mutate(lwr=fit-1.96*se.fit, upr=fit+1.96*se.fit, recognizer="Kaleidoscope") %>% 
  dplyr::select(fit, lwr, upr, recognizer) %>% 
  cbind(newdata.pres.1) %>% 
  rbind(newdata.pres)

#MonitoR
fn.m.0 <- glm(correct ~ 1, data=all.pres.m, family="binomial")
fn.m.1 <- glm(correct ~ threshold, data=all.pres.m, family="binomial")
fn.m.2 <- glm(correct ~ threshold + I(threshold^2), data=all.pres.m, family="binomial")
fn.m.3 <- glm(correct ~ threshold + I(threshold^2) + I(threshold^3), data=all.pres.m, family="binomial")
model.sel(fn.m.0, fn.m.1, fn.m.2, fn.m.3, rank=AIC)
anova(fn.m.3, test="Chisq")
hoslem.test(all.pres.s$correct, fitted(fn.m.2))

newdata.pres <- as.data.frame(predict(fn.m.2, newdata.pres.1, se.fit=TRUE, type="response")) %>% 
  mutate(lwr=fit-1.96*se.fit, upr=fit+1.96*se.fit, recognizer="MonitoR") %>% 
  dplyr::select(fit, lwr, upr, recognizer) %>% 
  cbind(newdata.pres.1) %>% 
  rbind(newdata.pres)

#Human
fn.h.0 <- glm(correct ~ 1, data=all.pres.h, family="binomial")
fn.h.1 <- glm(correct ~ threshold, data=all.pres.h, family="binomial")
fn.h.2 <- glm(correct ~ threshold + I(threshold^2), data=all.pres.h, family="binomial")
fn.h.3 <- glm(correct ~ threshold + I(threshold^2) + I(threshold^3), data=all.pres.h, family="binomial")
summary(fn.h.3)
anova(fn.h.3, test="Chisq")

newdata.pres <- as.data.frame(predict(fn.h.0, newdata.pres.1, se.fit=TRUE, type="response")) %>% 
  mutate(lwr=fit-1.96*se.fit, upr=fit+1.96*se.fit, recognizer="Human") %>% 
  dplyr::select(fit, lwr, upr, recognizer) %>% 
  cbind(newdata.pres.1) %>% 
  rbind(newdata.pres)

#5. Plot-----
newdata.pres$recognizer <- as.factor(newdata.pres$recognizer)
newdata.pres$recognizer <- factor(newdata.pres$recognizer, c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"))

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        axis.title.x=element_text(margin=margin(10,0,0,0), size=18),
        axis.title.y=element_text(margin=margin(0,10,0,0), size=18),
        legend.justification=c(0,0), legend.position=c(0,0),
        legend.background = element_rect(fill=alpha('white', 0.01)),
        legend.key.width=unit(2, "cm"))

ggplot(data=newdata.pres, aes(group=recognizer)) +
  geom_ribbon(aes(x=threshold, ymin=lwr, ymax=upr), alpha=0.3, col="grey") +
  geom_line(aes(x=threshold, y=fit, colour=recognizer, linetype=recognizer), lwd=1.2, linejoin="round") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  my.theme+
  labs(x="Score Threshold", y="Presence-absence recall") +
  scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                      values=cols,
                      labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"),
                      name="Processing\napproach") +
  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                        values=lines,
                        labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"),
                        name="Processing\napproach")


#OCCUPANCY----

#1. File list----
xtabs(~site+day, files)

files.occ  <- files[!(files$site %in% c("SM4620", "SM7517", "SM7522", "SM7616", "SM7699", "SM7708", "SM7714")), ]

xtabs(~site, files.occ)

files.occ$visit <- c(1,2,1,2,3,1,2,3,1,2,1,2,1,2,3,4,1,2,3,1,2,3,4,1,2,1,2,3,1,2,3,4,1,2,1,2,3,1,2,1,2,1,2,3,4,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,3,1,2,3,1,2,1,2,1,2,1,2,3,4,1,2,1,2,3,1,2,3,4,1,2,3,4,1,2,1,2,1,2,3,1,2,3,1,2,3,4,1,2,1,2,3)

#2. Gold standard----
gold.occ <- all %>% 
  group_by(ID, recognizer) %>% 
  dplyr::summarize(hit=n(), det=sum(pres)) %>% 
  right_join(files.occ, by="ID") %>% 
  dplyr::select(ID, site, recognizer, det, visit) %>% 
  spread(key=recognizer, value=det) %>% 
  mutate(CNN = ifelse(is.na(CNN),0,CNN)) %>% 
  mutate(Human = ifelse(is.na(Human),0,Human)) %>% 
  mutate(Kaleidoscope = ifelse(is.na(Kaleidoscope),0,Kaleidoscope)) %>% 
  mutate(MonitoR = ifelse(is.na(MonitoR),0,MonitoR)) %>% 
  mutate(RavenPro = ifelse(is.na(RavenPro),0,RavenPro)) %>% 
  mutate(SongScope = ifelse(is.na(SongScope),0,SongScope)) %>% 
  mutate(gold = pmax(CNN, Human, Kaleidoscope, MonitoR, RavenPro, SongScope)) %>% 
  mutate(pres = ifelse(gold > 0, 1, 0)) %>% 
  dplyr::select(ID, gold, pres, visit, site)

gold.site <- gold.occ %>% 
  group_by(site) %>% 
  summarize(sum = sum(pres)) %>% 
  mutate(pres = ifelse(sum > 0, 1, 0))

sum(gold.site$pres)
count(gold.site$pres)

#3. Occupancy data frame----

#lists of options to subset data by
score <- rep(seq(from=0.00, to=0.99, by=0.01), 6)
recognizer <- append(rep("CNN", 100), rep("Kaleidoscope", 100)) %>% 
  append(rep("MonitoR", 100)) %>% 
  append(rep("RavenPro", 100)) %>% 
  append(rep("SongScope", 100)) %>% 
  append(rep("Human", 100))

all.occ <- data.frame()

for(i in 1:length(score)){
  
  score.i <- paste0(score[i])
  recognizer.i <- paste0(recognizer[i])
  data.1 <- all %>%
    dplyr::filter(score.r > score.i, recognizer==recognizer.i) %>% 
    group_by(ID) %>% 
    dplyr::summarize(det=sum(pres), hit=n()) %>% 
    right_join(gold.occ, by = "ID") %>% 
    mutate(det = ifelse(is.na(det), 0, det)) %>% 
    mutate(hit = ifelse(is.na(hit), 0, hit)) %>% 
    mutate(pres.rec = ifelse(det > 0, 1, 0)) %>% 
    dplyr::select(site, visit, pres.rec) %>% 
    ungroup()
  data.1$recognizer <- recognizer.i
  data.1$threshold <- as.numeric(score.i)
  data.2 <- data.1 %>% 
    dplyr::group_by(recognizer, threshold) %>% 
    spread(key=visit, value=pres.rec) %>% 
    ungroup()
  data.2 <- as.data.frame(data.2)
  
  umf <- unmarkedFrameOccu(y = data.2[,4:7], siteCovs = data.2[,1:3])
  thresh1 <- occu(~1 ~1, umf)
  state <- as(backTransform(thresh1, 'state'), "data.frame") %>% 
    select(Estimate, SE) %>% 
    mutate(lower=Estimate-1.96*SE, upper=Estimate+1.96*SE)
  colnames(state) <- c("Occu", "OccuSE", "OccuLower", "OccuUpper")
  det <- as(backTransform(thresh1, type="det"), "data.frame") %>% 
    select(Estimate, SE) %>% 
    mutate(lower=Estimate-1.96*SE, upper=Estimate+1.96*SE)
  colnames(det) <- c("Det", "DetSE", "DetLower", "DetUpper")
  occ <- cbind(state, det)
  occ$recognizer <- recognizer.i
  occ$threshold <- score.i
  all.occ <- rbind(all.occ, occ)

}

all.occ <- all.occ %>% 
  mutate(OccuLower = ifelse(OccuLower < 0, 0, OccuLower),
         OccuLower = ifelse(OccuLower > 1, 1, OccuLower),
         OccuUpper = ifelse(OccuUpper < 0, 0, OccuUpper),
         OccuUpper = ifelse(OccuUpper > 1, 1, OccuUpper),
         DetLower = ifelse(DetLower < 0, 0, DetLower),
         DetLower = ifelse(DetLower > 1, 1, DetLower),
         DetUpper = ifelse(DetUpper < 0, 0, DetUpper),
         DetUpper = ifelse(DetUpper > 1, 1, DetUpper))

#4. Plot----

all.occ$recognizer <- as.factor(all.occ$recognizer)
all.occ$recognizer <- factor(all.occ$recognizer, c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"))
all.occ$threshold <- as.numeric(all.occ$threshold)

Occu <- ggplot(data=all.occ) +
  geom_ribbon(aes(x=threshold, ymin=OccuLower, ymax=OccuUpper, group=recognizer), alpha=0.3, col="grey") +
  geom_line(aes(x=threshold, y=Occu, colour=recognizer, linetype=recognizer), lwd=1.2)+
  my.theme +
  labs(x="Score Threshold", y="Occupancy", colour="Recognizer")+
  theme(legend.position="none") +
  scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                      values=cols,
                      labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"),
                      name="Processing\napproach") +
  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                        values=lines,
                        labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"),
                        name="Processing\napproach") +
  scale_y_continuous(lim=c(0,1))

Occu.wrap <- Occu + facet_wrap(~recognizer, nrow=3, ncol=2)

Det <- ggplot(data=all.occ) +
  geom_ribbon(aes(x=threshold, ymin=DetLower, ymax=DetUpper, group=recognizer), alpha=0.3, col="grey") +
  geom_line(aes(x=threshold, y=Det, colour=recognizer, linetype=recognizer), lwd=1.2)+
  my.theme +
  labs(x="Score Threshold", y="Detection")+
  theme(legend.position="none")+
  scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                      values=cols,
                      labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"),
                      name="Processing\napproach") +
  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"),
                        values=lines,
                        labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"),
                        name="Processing\napproach") +
  scale_y_continuous(lim=c(0,1))

Det.wrap <- Det + facet_wrap(~recognizer, nrow=3, ncol=2)

grid.arrange(Det.wrap, Occu.wrap, ncol=2)

#CALL RATE----

#1. File list----
files.rate <- rbind(files, files, files, files, files)
min <- append(rep(0, 117), rep(1, 117)) %>% 
  append(rep(2, 117)) %>% 
  append(rep(3, 117)) %>% 
  append(rep(4, 117))
files.rate <- cbind(files.rate, min)

#2. Gold standard----
gold.rate <- all %>% 
  group_by(ID, recognizer, min) %>% 
  dplyr::summarize(hit=n(), det=sum(pres)) %>% 
  right_join(files.rate, by=c("ID", "min")) %>% 
  dplyr::select(ID, recognizer, det, min) %>% 
  group_by(ID, min) %>% 
  tidyr::spread(key=recognizer, value=det) %>% 
  mutate(CNN = ifelse(is.na(CNN),0,CNN)) %>% 
  mutate(Human = ifelse(is.na(Human),0,Human)) %>% 
  mutate(Kaleidoscope = ifelse(is.na(Kaleidoscope),0,Kaleidoscope)) %>% 
  mutate(MonitoR = ifelse(is.na(MonitoR),0,MonitoR)) %>% 
  mutate(RavenPro = ifelse(is.na(RavenPro),0,RavenPro)) %>% 
  mutate(SongScope = ifelse(is.na(SongScope),0,SongScope)) %>% 
  mutate(gold = pmax(CNN, Human, Kaleidoscope, MonitoR, RavenPro, SongScope)) %>% 
  mutate(pres = ifelse(gold > 0, 1, 0)) %>% 
  dplyr::select(ID, min, gold, pres) %>% 
  ungroup()

#3. Recognizer data frame----

#lists of options to subset data by
score <- rep(seq(from=0.00, to=0.99, by=0.01), 6)
recognizer <- append(rep("CNN", 100), rep("Kaleidoscope", 100)) %>% 
  append(rep("MonitoR", 100)) %>% 
  append(rep("RavenPro", 100)) %>% 
  append(rep("SongScope", 100)) %>% 
  append(rep("Human", 100))

#summarize data by options

all.rate <- data.frame()

for(i in 1:length(score)){
  score.i <- paste0(score[i])
  recognizer.i <- paste0(recognizer[i])
  data.1 <- all %>%
    dplyr::filter(score.r > as.numeric(score.i), recognizer==recognizer.i) %>% 
    group_by(ID, min) %>% 
    dplyr::summarize(det=sum(pres), hit=n()) %>% 
    right_join(gold.rate, by = c("ID", "min")) %>% 
    mutate(det = ifelse(is.na(det), 0, det)) %>% 
    ungroup()
  cor.1 <- as.data.frame(cor(data.1$det, data.1$gold))
  colnames(cor.1) <- "cor"
  cor.1$recognizer <- recognizer.i
  cor.1$threshold <- as.numeric(score.i)
  all.rate <- rbind(all.rate, cor.1)
}

#4. Plot----

all.rate$recognizer <- factor(all.rate$recognizer, c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"))

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        axis.title.x=element_text(margin=margin(10,0,0,0), size=18),
        axis.title.y=element_text(margin=margin(0,10,0,0), size=18),
        legend.justification=c(0,0), legend.position=c(0,0),
        legend.background = element_rect(fill=alpha('white', 0.01)),
        legend.key.width=unit(2, "cm"))

ggplot(data=all.rate, aes(group=recognizer)) +
  geom_line(aes(x=threshold, y=cor, colour=recognizer, linetype=recognizer), lwd=1.2) +
  scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"), values=cols,
                      labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach")+
  labs(x="Score Threshold", y="Call rate Spearman correlation", colour="Processing\napproach") +
  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"), values=lines,
                      labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach")+
  labs(x="Score Threshold", y="Call rate Spearman correlation", colour="Processing\napproach") +
  scale_y_continuous(lim=c(0,1)) +
  my.theme

#EFFICIENCY----

#1. Determine cutoffs----

opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = x+y
    ind = which(d == max(d, na.rm=TRUE))
    c(precision = y[[ind]], recall = x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

#predictions
pred.s <- prediction(all.s$score.r, all.s$pres)
pred.c <- prediction(all.c$score.r, all.c$pres)
pred.r <- prediction(all.r$score.r, all.r$pres)
pred.k <- prediction(all.k$score.r, all.k$pres)
pred.m <- prediction(all.m$score.r, all.m$pres)

perf.c <- performance(pred.c, measure="prec", x.measure="rec")
plot(perf.c, col="red", xlim = c(0,1), ylim=c(0,1))
perf.s <- performance(pred.s, measure="prec", x.measure="rec")
plot(pr.s, col="blue",  add=TRUE)
perf.m <- performance(pred.m, measure="prec", x.measure="rec")
plot(perf.m, col="green", add=TRUE)
perf.k <- performance(pred.k, measure="prec", x.measure="rec")
plot(perf.k, col="purple", add=TRUE)
perf.r <- performance(pred.r, measure="prec", x.measure="rec")
plot(perf.r, col="orange", add=TRUE)
legend.text <- c("CNN", "SongScope", "MonitoR", "Kaleidoscope", "RavenPro")
legend.colour <- c("red", "blue", "green", "purple", "orange")
legend("topright", col = legend.colour, legend = legend.text, lty=1)

print(opt.cut(perf.c, pred.c))
print(opt.cut(perf.k, pred.k))
print(opt.cut(perf.m, pred.m))
print(opt.cut(perf.r, pred.r))
print(opt.cut(perf.s, pred.s))

eff.c <- subset(all.c, score.r>0.214)
count(eff.c)
eff.k <- subset(all.k, score.r>0)
count(eff.k)
eff.m <- subset(all.m, score.r>0.0005)
count(eff.m)
eff.r <- subset(all.r, score.r>0.067)
count(eff.r)
eff.s <- subset(all.s, score.r>0.220)
count(eff.s)



#2. Import data----

e <- read.csv("Efficiency2.csv")

hours <- as.data.frame(c(1:1000))
recognizer <- as.character(e$Recognizer)
build <- e$Build
run.0 <- e$Run.0
run.acc <- e$Run.Accuracy

e.long <- data.frame()

for(i in 1:length(recognizer)){
  
  time.0 <- (build[i] + hours*run.0[i])/hours
  time.acc <- (build[i] + hours*run.acc[i])/hours
  data <- cbind(hours, time.0, time.acc)
  colnames(data) <- c("Hours", "Time.0", "Time.Acc")
  data$Recognizer <- recognizer[i]
  e.long <- rbind(e.long, data)
  
}


e.long$Recognizer <- as.factor(e.long$Recognizer)
e.long$Hours <- as.numeric(e.long$Hours)

e.long$Recognizer <- factor(e.long$Recognizer, c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"))

my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        axis.title.x=element_text(margin=margin(10,0,0,0), size=18),
        axis.title.y=element_text(margin=margin(0,10,0,0), size=18),
        legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill=alpha('white', 0.01)),
        legend.key.width=unit(2, "cm"))

ggplot(data=e.long, aes(group=Recognizer)) +
  geom_line(aes(x=Hours, y=Time.Acc, colour=Recognizer, linetype=Recognizer), lwd=1.2) +
  ylim(0,5) +
#  scale_colour_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"), values=cols, labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
#  labs(x="Database size (GB)", y="Processing time (hr/GB)", colour="Processing\napproach") +
#  scale_linetype_manual(breaks=c("Human", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "SongScope"), values=lines,labels=c("Human listening", "CNN", "Kaleidoscope", "MonitoR", "RavenPro", "Song Scope"), name="Processing\napproach") +
#  labs(x="Database size (GB)", y="Processing time (hr/GB)", colour="Processing\napproach") +
  my.theme

#SPECTOGRAM FIGURE----

library(tuneR)
library(signal)
library(seewave)


rec <- readWave("RD-207-03-0_0+1_20150629_230000.wav", from=130, to=137.5, units="seconds")
rec.down <- downsample(rec, 11025)
#spectro(rec.down)
spectro(rec, wl=510, wn="blackman", flim=c(0,9))

