library(ggplot2)

## Scatter plots
# Actual data
stat <- read.table('PATH\\Cd11b_MK_EM.txt',header = T)
stat <- stat/1250*400 # 1250 pix = 400 um
stat$grp <- stat$CXCR4_high_MK-stat$CXCR4_low_MK
stat$grp[stat$grp>0] <- "low"
stat$grp[stat$grp<0] <- "high"
p1 <- ggplot(stat,aes(x=CXCR4_low_MK,y=CXCR4_high_MK,color = grp))+
  geom_point(alpha = 0.6)+
  xlim(0,100)+
  ylim(0,100)+
  xlab("To CXCR4 low MKs (um)")+
  ylab("To CXCR4 high MKs (um)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank())+
  scale_color_manual(values = c("#0000FF","#FF0000")[2:1])
p1
ggsave(p1,filename = "PATH\\Cd11b_MK_actual_sactter.pdf",width = 4,height = 3)

# Random data
stat <- read.table('PATH\\example_low_high_dist.txt',header = T)
stat <- stat/1250*400 # 1250 pix = 400 um

stat$grp <- stat$CXCR4_high_MK-stat$CXCR4_low_MK
stat$grp[stat$grp>0] <- "low"
stat$grp[stat$grp<0] <- "high"
p1 <- ggplot(stat,aes(x=CXCR4_low_MK,y=CXCR4_high_MK,color = grp))+
  geom_point(alpha = 0.6)+
  xlim(0,300)+
  ylim(0,300)+
  xlab("To CXCR4 low MKs (um)")+
  ylab("To CXCR4 high MKs (um)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank())+
  scale_color_manual(values = c("#0000FF","#FF0000")[2:1])
p1
ggsave(p1,filename = "PATH\\Cd11b_MK_random_sactter.pdf",width = 4,height = 3)

## Density plots
# Density plot for actual data
stat <- read.table('PATH\\Cd11b_MK_EM.txt',header = T)
stat2 <- data.frame(grp = c(rep("low",85),rep("high",85)),
                    data = c(stat$CXCR4_low_MK,stat$CXCR4_high_MK),
                    fill = c(rep("blue",85),rep("red",85)))
mean(stat$CXCR4_low_MK)
summary(stat$CXCR4_low_MK)
mean(stat$CXCR4_high_MK)
summary(stat$CXCR4_high_MK)

ks.test(stat$CXCR4_low_MK,stat$CXCR4_high_MK) # KS test

p2 <- ggplot(stat2,aes(x=data,fill=fill))+
  geom_density(highition = 'identity', stat = 'density',alpha = 0.6)+
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank())+
  geom_vline(aes(xintercept=25.62))+
  geom_vline(aes(xintercept=15.35))+
  xlim(0,80)
p2
ggsave(p2,filename = "PATH\\Cd11b_MK_actual_density.pdf",width = 9,height = 4)

# Random distribution mean density
hi_mean <- read.table("PATH\\closer_high_mean_dist.txt")
lo_mean <- read.table("PATH\\closer_low_mean_dist.txt")
hi_mean <- hi_mean/1250*400 # mean distance 1250 px = 400 um
lo_mean <- lo_mean/1250*400 # mean distance 1250 px = 400 um
stat3 <- data.frame(grp = c(rep("low",500),rep("high",500)),
                    data = c(lo_mean$V1,hi_mean$V1),
                    fill = c(rep("blue",500),rep("red",500)))
mean(hi_mean$V1) # 35.3686
mean(lo_mean$V1) # 27.76236
ks.test(hi_mean$V1,lo_mean$V1) # p = 2.2e-16

pnorm(15.36,mean = mean(hi_mean$V1),sd = sd(hi_mean$V1),lower.tail = T) # p = 1.77e-10
pnorm(25.62,mean = mean(lo_mean$V1),sd = sd(lo_mean$V1),lower.tail = T) # p = 0.1399

p3 <- ggplot(stat3,aes(x=data,fill=fill))+
  geom_density(highition = 'identity', stat = 'density',alpha = 0.2)+
  scale_fill_manual(values = c("blue","red"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank())+
  geom_histogram(aes(y=..density..), alpha=0.5, color = "grey",
                 position="identity")+
  geom_vline(aes(xintercept=25.62))+
  geom_vline(aes(xintercept=15.35))+
  geom_vline(aes(xintercept=35.37))+
  geom_vline(aes(xintercept=27.76))+
  xlim(0,80)
p3
ggsave(p3,filename = "PATH\\Cd11b_MK_random_density.pdf",width = 9,height = 4)
