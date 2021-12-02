#morphology assay

#load in mRNA transport excel doc.
library(readxl)
morph <- read_xlsx(path = "NMorph_R.xlsx", 
                   sheet = 1, col_names = T,skip = 1) #satisfies number 3
head(morph)

#subset the nuclear and cytoplasmic signals as well as the ratio columns
morph.sub<-morph[,c('Percent Abnormal Y',
                    'Percent Abnormal O')] #satisfies numbers 5 and 6
head(morph.sub)

co.names <- c('perc.ab.y',
              'perc.ab.o')

colnames(morph.sub) <- co.names 
head(morph.sub)

##################################################################
library(plyr)

summarize(.data=morph.sub,mean.ab.young=mean(perc.ab.y), sd.per.ab.y=sd(perc.ab.y)) #satisfies number 8

summarize(.data=morph.sub,mean.ab.young=mean(perc.ab.o), sd.per.ab.y=sd(perc.ab.o)) #satisfies number 8


morph.sub.desc<-arrange(df = morph.sub,desc(perc.ab.o)) #satisfies numbers 7
head(morph.sub.desc)

summary(morph.sub.desc)

library(tidyverse)

p <- morph.sub.desc


outlier.perc <- if_else(condition = p > 70, true = "outlier", 
                        false = "normal")
#checks to see if any outliers are present in the data set

outlier.perc

no.outliers <- subset(x = morph.sub.desc, 
                      perc.ab.o < 70 & perc.ab.y < 70) 
#removes the previously identified outliers


ggplot(data = no.outliers, aes(x = perc.ab.y,y= num)) +
  geom_col()+ geom_point(color='blue',size=2)+
  scale_fill_manual(values=wes_palette("GrandBudapest1",n=2)+
                      stat_smooth(method = 'loess'))








####################################################################

#load in mRNA transport excel doc.
library(readxl)
mrna <- read_xlsx(path = "cell.prof.fish.xlsx", 
                  sheet = 1, col_names = T,skip = 2) #satisfies number 3
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
#calculating ratios using R instead of using pre-existing ratios

#######################################################################

#making histogram or density plot

hist.1<-hist(mrna$r.ratio.1) 

hist.2<-hist(mrna$r.ratio.3)
#visually comparring data to see if any prominent features occur.
##Image 3 appears to shifted to the right, indicating inconsistency 
##in the Cell Profiler pipeline made.  

#making Whisker plots

mrna.sub<-mrna[c('r.ratio.1', 'r.ratio.3')]

head(mrna.sub)
old
library(ggplot2)

#making box-whicker plots

ratio.1<-ggplot(data = mrna.sub, aes(x = , y=r.ratio.1, fill=r.ratio.1)) +
  geom_boxplot() +
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        panel.background = element_rect(fill="pink",colour="black"))


ratio.3<-ggplot(data = mrna.sub, aes(x = , y=r.ratio.3, fill=r.ratio.3)) +
  geom_boxplot() +
  theme(axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        panel.background = element_rect(fill="tan",colour="black"))

boxplot(mrna.sub)

#############################################################################

#t-test - old v young
summary(mrna.sub)

library(stats)
#onesided
t.test(x = mrna.sub,conf.level = 0.95)
#suggests that r.ratio.1 is larger than r.ratio.3 --> pipline isn't functioning
##properly



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












ggplot(data = mrna, aes(x = r.ratio.1, y = )) +
  geom_histogram(bins = 50)

ggplot(data = mrna, aes(x = r.ratio.3, y = )) +
  geom_histogram()












#subset the nuclear and cytoplasmic signals as well as the ratio columns
mrna.sub<-mrna[,c('ONUC:CYT','YNUC:CYT')] #satisfies numbers 5 and 6
head(mrna.sub)


library(plyr)

mrna.desc<-arrange(.data = mrna.sub,desc(mrna.sub,c(YNUC:CYT))) #satisfies numbers 7
head(mrna.desc)

rm.na<-rm(mrna.asc)

summarize(mrna.desc$c(1:2))

################################################

#qPCR

library(readxl)
qpcr <- read_xlsx(path = "qpcr.r.copy.xlsx", 
                  sheet = 1, col_names = T,skip = 0,na = "0") #satisfies number 3
head(qpcr)




#subset the nuclear and cytoplasmic signals as well as the ratio columns
sub.qpcr<-qpcr[,c('num','protein',
                  'old','young')] 
head(sub.qpcr)


p <- sub.qpcr
outlier.gene <- if_else(condition = p > 5, true = "outlier", 
                        false = "normal")


no.out.qpcr <- subset(x = sub.qpcr, 
                      old < 5 & young < 5) 
library(ggplot2)
p<-ggplot(data = no.out.qpcr, aes(x = young, y = num)) +
  geom_point()

png('p.png')
dev.off()

plot.1<-ggplot(data = no.out.qpcr, aes(x = young,y=num)) +
  geom_col()+ geom_point(color='blue',size=2)+
  scale_fill_manual(values=wes_palette("GrandBudapest1",n=2)+
                      stat_smooth(method = 'loess'))

##########################################################
#ddply to execute function over df

old.10 <- ddply(.data = no.out.qpcr, .variables = c("old","young"), function(r){
  
  no.out.qpcr$old.10x <- no.out.qpcr$old*10 #says to add parcel.length.km and 
  #fill it with values of the parcel.length.m divided by 1000
  return(r) #return the new field, 'x'
  
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


# Write data to csv files:  
# decimal point = "." and value separators = comma (",")
write.csv(mtcars, file = "mtcars.csv")
# Write data to csv files: 
# decimal point = comma (",") and value separators = semicolon (";")
write.csv2(mtcars, file = "mtcars.csv")





#Making line plot

library(readxl)
qpcr.b <- read_xlsx(path = "qpcr.r.copyb.xlsx", 
                    sheet = 1, col_names = T,skip = 0,na = "0") #satisfies number 3
head(qpcr.b)

qpcr.bsub<-qpcr.b[,c("age","avg")]

head(qpcr.bsub)

library("RColorBrewer")

p <- ggplot(data = qpcr.bsub, aes(x = age, y = avg, color='red')) +
  #use the 'color=" argument for elements like points and lines
  #use the 'fill =" argument for elements like bars or box-whiskers or areas
  geom_point() +
  geom_smooth(method = "lm")+
  scale_color_brewer(palette = "Set2")
p

png("p.png")
dev.off()


#Testing linear model
