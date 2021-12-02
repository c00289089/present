#morphology assay

#load in mRNA transport excel doc.
library(readxl)
morph <- read_xlsx(path = "NMorph_R.xlsx", 
                   sheet = 1, col_names = T,skip = 1) 
head(morph)

#subset the nuclear and cytoplasmic signals as well as the ratio columns
morph.sub<-morph[,c('Percent Abnormal Y',
                    'Percent Abnormal O')] 
head(morph.sub)


co.names <- c('perc.ab.y',
              'perc.ab.o')

colnames(morph.sub) <- co.names 
head(morph.sub)

##################################################################
library(plyr)

summarize(.data=morph.sub,mean.ab.young=mean(perc.ab.y), sd.per.ab.y=sd(perc.ab.y)) 

summarize(.data=morph.sub,mean.ab.young=mean(perc.ab.o), sd.per.ab.y=sd(perc.ab.o)) 


morph.sub.desc<-arrange(df = morph.sub,desc(perc.ab.o))
head(morph.sub.desc)

library(tidyverse)

p <- morph.sub.desc


outlier.perc <- if_else(condition = p > 70, true = "outlier", 
                        false = "normal")
outlier.perc



no.outliers <- subset(x = morph.sub.desc, 
                      perc.ab.o < 70 & perc.ab.y < 70) 

no.outliers
#removes the previously identified outliers

library(wesanderson)

str(no.outliers)
hist1<-hist(no.outliers$perc.ab.o)
hist2<-hist(no.outliers$perc.ab.y)
lines(density(no.outliers$perc.ab.y,adjust = .1))

t.test(x=no.outliers$perc.ab.o,y=no.outliers$perc.ab.y)

#slight evidence of a significant difference.  


####################################################################

#load in mRNA transport excel doc. --> checking cellprofiler pipeline
library(readxl)
mrna <- read_xlsx(path = "cell.prof.fish.xlsx", 
                  sheet = 1, col_names = T,skip = 2) 
head(mrna)



nuc2cyt.1 <- function(x){ 
  
  y <- x/mrna$cyto.1 
  return(y) }

mrna$r.ratio.1<-nuc2cyt.1(mrna$nuc.1)



nuc2cyt.3 <- function(a){
  
  b <- a/mrna$cyto.3 
  return(b) 
}


mrna$r.ratio.3<-nuc2cyt.3(mrna$nuc.3)
head(mrna)
#calculating ratios using R to check against other calculations



#making histogram or density plot

hist.1<-hist(mrna$r.ratio.1)
lines(density(mrna$r.ratio.1))

hist.2<-hist(x=mrna$r.ratio.1)
lines(density(mrna$r.ratio.3))



#visually comparring data to see if any prominent features occur.
##Image 3 appears to shifted to the right, indicating inconsistency 
##in the Cell Profiler pipeline made.  

#making Whisker plots

mrna.sub<-mrna[c('r.ratio.1', 'r.ratio.3')]

head(mrna.sub)

library(ggplot2)

#making box-whicker plots

ratio.1<-ggplot(data = mrna.sub, aes(x = , y=r.ratio.1, fill=r.ratio.1)) +
  geom_boxplot() +
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        panel.background = element_rect(fill="pink",colour="black"))

ratio.1

ratio.3<-ggplot(data = mrna.sub, aes(x = , y=r.ratio.3, fill=r.ratio.3)) +
  geom_boxplot() +
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        panel.background = element_rect(fill="tan",colour="black"))

ratio.3

bpmrna<-boxplot(mrna.sub) #additional whisker boxplot with both

summary(mrna.sub)

#t-test - slide 1 v slide 3
t.test(x=mrna.sub$r.ratio.1,y=mrna.sub$r.ratio.3)


#p value <<<0.05 --> suggests values significantly differ-->pipline isn't 
#functioning properly

boxplot(mrna.sub)

w.1a<-ggplot(data = mrna.sub, aes(x = r.ratio.1, y = 'y')) +
  geom_boxplot()
w.1a
png('w.1a.png')

dev.off()

z<-ggplot(data = mrna.sub, aes(x = r.ratio.3, y = )) +
  geom_boxplot(outlier.colour = "red")
z

png('z.png')

dev.off()

################################################

#qPCR

library(readxl)
qpcr <- read_xlsx(path = "qpcr.r.copy.xlsx", 
                  sheet = 1, col_names = T,skip = 0,na = "0") 
head(qpcr)




#subset the nuclear and cytoplasmic signals as well as the ratio columns
sub.qpcr<-qpcr[,c('num','protein',
                  'old','young')] 
head(sub.qpcr)

pd<-arrange(.data = sub.qpcr,desc(sub.qpcr,c(num)))
head(pd)



p <- pd
outlier.gene <- if_else(condition = p > 5, true = "outlier", 
                        false = "normal")


no.out.qpcr <- subset(x = sub.qpcr, 
                      old < 5 & young < 5) 
library(ggplot2)
p<-ggplot(data = no.out.qpcr, aes(x = young, y = num)) +
  geom_point()

p

png('p.png')
dev.off()

plot.1<-ggplot(data = no.out.qpcr, aes(x = young,y=num)) +
  geom_line()+ geom_point(color='blue',size=2)+
  scale_fill_manual(values=wes_palette("GrandBudapest1",n=2)+
                      stat_smooth(method = 'lm'))

plot.1

##########################################################
#ddply to execute function over df
no.out.qpcr

old.10 <- ddply(.data = no.out.qpcr, .variables = c("old","young"), function(r){
  
  no.out.qpcr$old.10 <- no.out.qpcr$old*10
  
}, .progress = "text", .inform = T)

head(old.10)

##################################################################

#exporting data sets

write.table(no.out.qpcr, file = "export.one.txt", sep = "\t",
            row.names = TRUE, col.names = NA)

#################################################################

#exporting plots

ggsave('plot.1.png')
dev.off()
##############################################################################
#Making line plot

library(readxl)
qpcr.b <- read_xlsx(path = "qpcr.r.copyb.xlsx", 
                    sheet = 1, col_names = T,skip = 0,na = "0") 
head(qpcr.b)

qpcr.bsub<-qpcr.b[,c("age","avg")]

head(qpcr.bsub)

library("RColorBrewer")

p <- ggplot(data = qpcr.bsub, aes(x = age, y = avg, color='red')) +
  geom_point() +
  geom_smooth(method = "lm")+
  scale_color_brewer(palette = "Set2")
p

png("p.png")
dev.off()

#Testing linear model

y <- qpcr.bsub$age
x <- qpcr.bsub$avg


#linear model
fit <- lm(formula = y~x)
plot(fit)
summary(fit)

#linear model is not likely a good fit for the data


#Check polynomial model
fit2 <- lm(formula = age ~ poly(x = avg, degree = 2), data = qpcr.bsub)
plot(fit2)
summary(fit2)
#not a good model but some indcation of outliers skewing results; R squared=.4%

fit3 <- lm(formula =  age~ poly(x = avg, degree = 3), data = qpcr.bsub)
summary(fit3)



#not a good model but some indication of outliers skewing results; R squared=.34%


