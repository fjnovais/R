########################################
#           POOL-ANALYSIS             #
########################################
#
#By Francisco Novais
#
#PLINK_assoc.file + Z_test_pool.file = compair.file
#
#1-Concatenate two plink and z_test_pool files
order_qui_z <- gpsob_qui_Z [order (gpsob_qui_Z[,1]), ]
order_plink <- gpsob [order(gpsob[,2]), ]
order_plink$CHR <-NULL
compair <- merge(order_plink,order_qui_z, by=1)
#
#2-Correlations between allelotyping and genotyping
MAF_A <- subset.data.frame(compair, subset= compair$A1 == "A")
MAF_B <- subset.data.frame(compair, subset= compair$A1 == "B")
MAF_B$freq.B.case <- MAF_B$F_A
MAF_B$freq.B.control <- MAF_B$F_U
MAF_A$freq.B.case <- 1- MAF_A$F_A
MAF_A$freq.B.control <- 1- MAF_A$F_U
compair_freq <- rbind.data.frame(MAF_A,MAF_B)
corfreq.case <- cor (compair_freq$p_top, compair_freq$freq.B.case, use= "pairwise.complete.obs")
corfreq.control <- cor (compair_freq$p_low, compair_freq$freq.B.control, use= "pairwise.complete.obs")
cor.rept1 <- cor (replicatas$GPT1.B.Allele.Freq,replicatas$GPT1.B.Allele.Freq.1,use= "pairwise.complete.obs" )
#
#3-Plots between pool frequency and indivudial frequency
plot (compair_freq$p_top~compair_freq$freq.B.case,xlab = "individual frequency", ylab = "pool frequency")
abline(lm(compair_freq$p_top~compair_freq$freq.B.case), col="red")
summary(lm(compair_freq$p_top~compair_freq$freq.B.case))
#
plot (compair_freq$p_low~compair_freq$freq.B.case,xlab = "individual frequency", ylab = "pool frequency")
abline(lm(compair_freq$p_low~compair_freq$freq.B.case), col="red")
summary(lm(compair_freq$p_low~compair_freq$freq.B.case))
#
pop_strat <- mean(compair_freq$p.value)/0.456
#
assoc.ind1 <- merge(rank.ind1,compair_freq,by=1)
cor(assoc.ind1$p_low, assoc.ind1$freq.B.control, use= "pairwise.complete.obs")
cor (assoc.ind1$p_top, assoc.ind1$freq.B.case, use= "pairwise.complete.obs")
assoc.pool2 <- merge(rank.pool2, compair_freq, by=1)
cor(assoc.pool2$p_low, assoc.pool2$freq.B.control, use= "pairwise.complete.obs")
cor (assoc.pool2$p_top, assoc.pool2$freq.B.case, use= "pairwise.complete.obs")
assoc.general3 <- merge(rank.general3, compair_freq, by=1)
cor(assoc.general3$p_low, assoc.general3$freq.B.control, use= "pairwise.complete.obs")
cor (assoc.general3$p_top, assoc.general3$freq.B.case, use= "pairwise.complete.obs")
#
#4-GWAS:files montage
pool.manh <- data.frame (compair_freq$SNP,compair_freq$Chr,compair_freq$BP,compair_freq$p.value)
names (pool.manh)[1] <- "SNP"
names (pool.manh)[2] <- "CHR"
names (pool.manh)[3] <- "BP"
names (pool.manh)[4] <- "P"
pool.manh <- pool.manh [order(pool.manh[,2]), ]
#
ind.manh <- data.frame (compair_freq$SNP,compair_freq$Chr,compair_freq$BP,compair_freq$P)
names (ind.manh)[1] <- "SNP"
names (ind.manh)[2] <- "CHR"
names (ind.manh)[3] <- "BP"
names (ind.manh)[4] <- "P"
ind.manh <- ind.manh [order(ind.manh[,2]), ]
#
pool.manh$P[pool.manh$P == 0] <- NA
manh.pool <- na.exclude(pool.manh)
ind.manh$P[ind.manh$P == 0] <- NA
manh.ind <- na.exclude (ind.manh)
#
#5-Manhattan Plot (General and per chromossome)
manhattan(manh.pool, col=c("black", "chocolate2"),suggestiveline = F,genomewideline = F)
par(new=T)
manhattan (manh.ind,col=c("black", "chocolate2"), suggestiveline = F,genomewideline = F)
#Same SNPs with significance
same.signif.pool<-data.frame(rank.p.pool$SNP , rank.p.pool$CHR.x , rank.p.pool$BP.x ,rank.p.pool$P.y )
names (same.signif.pool)[1] <- "SNP"
names (same.signif.pool)[2] <- "CHR"
names (same.signif.pool)[3] <- "BP"
names (same.signif.pool)[4] <- "P"
same.signif.pool <- same.signif.pool [order(same.signif.pool[,2]), ]
same.signif.ind<-data.frame(rank.p.ind$SNP , rank.p.ind$CHR.x , rank.p.ind$BP.x ,rank.p.ind$P.x )
names (same.signif.ind)[1] <- "SNP"
names (same.signif.ind)[2] <- "CHR"
names (same.signif.ind)[3] <- "BP"
names (same.signif.ind)[4] <- "P"
same.signif.ind <- same.signif.ind [order(same.signif.ind[,2]), ]
listsnps <- same.signif.ind$SNP
#chr 01
manhattan(subset(rank.pool, CHR==1), suggestiveline = F,genomewideline = F, cex=0.93, highlight =listsnps)
manhattan(subset(rank.ind, CHR==1), suggestiveline = F,genomewideline = F, cex=0.93,  ylim=range(1e-11:10), highlight =listsnps)
manhattan(subset(same.signif.pool, CHR==1), suggestiveline = F,genomewideline = F, cex=0.93)
manhattan(subset(same.signif.ind, CHR==1), suggestiveline = F,genomewideline = F, cex=0.93,  ylim=range(1e-6:5))
#
#6-Ranking p-values
rank.ind1<- subset.data.frame(manh.ind, subset= manh.ind$P < 0.01)
rank.pool2 <- subset.data.frame(manh.pool,subset= manh.pool$P < 0.05)
rank.general3<-merge.data.frame(rank.ind1,rank.pool, by= "SNP")
rank.p.ind <- rank.general [order(rank.general[,4]),]
rank.p.pool <- rank.general [order(rank.general[,7]),]
#
#
#OBS:Can't apply bonferroni or fdr corretion
vector.p <- manh.pool$P
teste<- p.adjust(vector.p,method = "fdr", n=21060)
sum(teste < 0.05)
vector.i <- manh.ind$P
teste2 <-p.adjust(vector.i, method = "fdr", n=21060)
sum(teste2 < 0.05)
#
#7-Genotyped population GWAS
geno1 <- data.frame (geno$SNP,geno$CHR,geno$BP,geno$P)
names (geno1)[1] <- "SNP"
names (geno1)[2] <- "CHR"
names (geno1)[3] <- "BP"
names (geno1)[4] <- "P"
geno1 <- geno1 [order(geno1[,2]), ]
#
geno1$P[geno1$P == 0] <- NA
geno1 <- na.exclude(geno1)
manhattan(geno1, col=c("black", "chocolate2"),suggestiveline = F,genomewideline = F)
rank.geno<- subset.data.frame(geno1, subset= geno1$P < 0.01)
teste<-merge.data.frame(rank.geno,rank.pool, by= "BP")
#
# Basic plot of x and y :
plot(x,y,col=rgb(0.4,0.4,0.8,0.6), pch=16 , cex=1.3 , xlab="" , ylab="") 
#
# Can we find a polynome that fit this function ?
model=lm(y ~ x + I(x^2) + I(x^3))
#
# I can get the features of this model :
summary(model)
model$coefficients
summary(model)$adj.r.squared
#
#For each value of x, I can get the value of y estimated by the model, and the confidence interval around this value.
myPredict <- predict( model , interval="predict" )

#Finally, I can add it to the plot using the line and the polygon function with transparency.
ix <- sort(x,index.return=T)$ix
lines(x[ix], myPredict[ix , 1], col=2, lwd=2 )  
polygon(c(rev(x[ix]), x[ix]), c(rev(myPredict[ ix,3]), myPredict[ ix,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)

