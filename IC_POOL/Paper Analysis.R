##############################
#      PAPER ANALYSIS        #
##############################
#
#by Francisco Novais
#
#
#1-Pool Quality Control 
getwd()
#1-PP14
setwd("C:/Users/Francisco Novais/Documents/Posgraduacao/IC/Dados allelotipagem/PP14_pools")                     
getwd()
PP14_GS <- read.delim("~/Posgraduacao/IC/Dados allelotipagem/PP14_pools/PP14_pool_GS.txt")
PP14_GS1 <-subset(PP14_GS,subset = PP14_GS$Cluster_Sep >0.40)#Cluster separation > .40
PP14_GS2 <-subset(PP14_GS1,subset = PP14_GS1$GENERAL_MEAN <0.98)#Exclude fixed SNP
PP14_GS3 <-subset(PP14_GS2,subset = PP14_GS2$GENERAL_MEAN >0.02)#Exclude fixed SNP
levels(PP14_GS$Chr)#see Chr names
PP14_GS4 <-subset(PP14_GS3,subset = PP14_GS3$Chr !="X")#Exclude non autossomal SNPs
PP14_GS5 <-subset(PP14_GS4,subset = PP14_GS4$Chr !="Y")#Exclude non autossomal SNPs
PP14_GS6 <-subset(PP14_GS5,subset = PP14_GS5$Chr != 0)#Exclude non autossomal SNPS
#Dos 54.609 ficaram 28.581 (cluster separation > 40 = 928 , snps fixados= 24.120 (falar da eficacia pra Nelore), non autossomicos =980)=total 26.028. Lembrando que apenas 25,492 sao polimorficos em Nelore
unique(PP14_GS5$Chr)#see Chr names
rm(PP14_GS1,PP14_GS2,PP14_GS3,PP14_GS4,PP14_GS5)#GS6 final first QC
#
#TOP and LOW frequency pool mean 
PP14_FREQ <- edit(PP14_GS6)
PP14_FREQ$FREQ_TOP <- ((PP14_FREQ$PP14TOP_1M + PP14_FREQ$PP14TOP_2M)/2)
PP14_FREQ$FREQ_LOW <- ((PP14_FREQ$PP14LOW_1M + PP14_FREQ$PP14LOW_2M)/2)
#
#
#2-Correlation between sub-pool and individual genotyping
#2.2.1-PP14-building files
PP14_LOWA <- read.table("C:/Users/Francisco Novais/Desktop/TRABALHO GMABT/PROJETO IC/Dados experimento/GWAS/PP14_LOWA", header=TRUE, quote="\"")
PP14_LOWB <- read.table("C:/Users/Francisco Novais/Desktop/TRABALHO GMABT/PROJETO IC/Dados experimento/GWAS/PP14_LOWB", header=TRUE, quote="\"")
PP14_TOPA <- read.table("C:/Users/Francisco Novais/Desktop/TRABALHO GMABT/PROJETO IC/Dados experimento/GWAS/PP14_TOPA", header=TRUE, quote="\"")
PP14_TOPB <- read.table("C:/Users/Francisco Novais/Desktop/TRABALHO GMABT/PROJETO IC/Dados experimento/GWAS/PP14_TOPB", header=TRUE, quote="\"")
PP14_TOPB$QUARTZO <- NULL #Exlude due low call rate
#TOP A individual Frequency
t10 <- apply(PP14_TOPA==0,1,sum,na.rm=T)
t11 <- apply(PP14_TOPA==1,1,sum,na.rm=T)
t12 <- apply(PP14_TOPA==2,1,sum,na.rm=T)
t <- t10 + t11 + t12
pta <- ((2*t12)+t11)/(2*t2)#calculate allele frequencies
PP14_TOPA$TOPA_FREQ <- pta
#TOP B individual Frequency
t20 <- apply(PP14_TOPB==0,1,sum,na.rm=T)
t21 <- apply(PP14_TOPB==1,1,sum,na.rm=T)
t22 <- apply(PP14_TOPB==2,1,sum,na.rm=T)
t2 <- t20 + t21 + t22
ptb <- ((2*t22)+t21)/(2*t2)#calculate allele frequencies
PP14_TOPB$TOPB_FREQ <- ptb
#LOW A individual Frequency
l10 <- apply(PP14_LOWA==0,1,sum,na.rm=T)
l11 <- apply(PP14_LOWA==1,1,sum,na.rm=T)
l12 <- apply(PP14_LOWA==2,1,sum,na.rm=T)
l1 <- l10 + l11 + l12
pla <- ((2*l12)+l11)/(2*l1)#calculate allele frequencies
PP14_LOWA$LOWA_FREQ <- pla
#LOW B individual Frequency
l20 <- apply(PP14_LOWB==0,1,sum,na.rm=T)
l21 <- apply(PP14_LOWB==1,1,sum,na.rm=T)
l22 <- apply(PP14_LOWB==2,1,sum,na.rm=T)
l2 <- l20 + l21 + l22
plb <- ((2*l22)+l21)/(2*l2)#calculate allele frequencies LOW B
PP14_LOWB$LOWB_FREQ <- plb
rm(l1,l10,l11,l12,l2,l20,l21,l22,t,t10,t11,t12,t2,t20,t21,t22,pla,plb,pta,ptb)
#
#2.2.2-Correlations sub-pool
pset_topa <- merge(x = PP14_FREQ,PP14_TOPA, by=1)#subset_topa
pset_topb <- merge(x = PP14_FREQ,PP14_TOPB, by=1)#subset_topb
pset_topb <- pset_topb[complete.cases (pset_topb[,28]),]
pset_topb <- edit(pset_topb)
write.table(pset_topb,file = "C:/Users/Francisco Novais/Documents/Posgraduacao/IC/Dados allelotipagem/PP14_pools/pset_topb.txt")
pset_topb <- pset_topb[-21164,]
pset_lowa <- merge(x = PP14_FREQ,PP14_LOWA, by=1)#subset_lowa
pset_lowb <- merge(x = PP14_FREQ,PP14_LOWB, by=1)#subset_lowb
#
corr.ta <- cor (pset_topa$PP14TOP_1M,pset_topa$TOPA_FREQ,use= "pairwise.complete.obs" )#correlation TOP A
corr.tb <- cor (pset_topb$PP14TOP_2M,pset_topb$TOPB_FREQ, use= "pairwise.complete.obs")#correlation TOP B
corr.la <- cor (pset_lowa$PP14LOW_1M,pset_lowa$LOWA_FREQ,use= "pairwise.complete.obs" )#correlation LOW A
corr.lb <- cor (pset_lowb$PP14LOW_2M,pset_lowb$LOWB_FREQ,use= "pairwise.complete.obs" )#correlation LOW B
#
#2.2.3-Correlations pools
#Building files with individual allele frequency
PP14_TOP <- merge(PP14_TOPA,PP14_TOPB, by=1)
PP14_TOP$TOPA_FREQ <- PP14_TOP$TOPB_FREQ <- NULL
t10 <- apply(PP14_TOP==0,1,sum,na.rm=T)
t11 <- apply(PP14_TOP==1,1,sum,na.rm=T)
t12 <- apply(PP14_TOP==2,1,sum,na.rm=T)
t1 <- t10 + t11 + t12
pta <- ((2*t12)+t11)/(2*t1)#calculate allele frequencies
PP14_TOP$TOP_FREQ <- pta
PP14_TOP <- data.frame (PP14_TOP$SNP,PP14_TOP$TOP_FREQ)
#
PP14_LOW <- merge(PP14_LOWA,PP14_LOWB, by=1)
PP14_LOW$LOWA_FREQ <- PP14_LOW$LOWB_FREQ <- NULL
l10 <- apply(PP14_LOW==0,1,sum,na.rm=T)
l11 <- apply(PP14_LOW==1,1,sum,na.rm=T)
l12 <- apply(PP14_LOW==2,1,sum,na.rm=T)
l1 <- l10 + l11 + l12
pla <- ((2*l12)+l11)/(2*l1)#calculate allele frequencies
PP14_LOW$LOW_FREQ <- pla
rm(l1,l10,l11,l12,l2,l20,l21,l22,t1,t10,t11,t12,t2,t20,t21,t22,pla,plb,pta,ptb)
PP14_LOW <- data.frame (PP14_LOW$SNP,PP14_LOW$LOW_FREQ)
#creating subset with freq pool and individual
pset_top <- data.frame (pset_topa$Name,pset_topa$FREQ_TOP)
pset_top <- merge (pset_top,PP14_TOP, by=1)
pset_low <- data.frame (pset_lowa$Name,pset_lowa$FREQ_LOW)
pset_low <- merge (pset_low,PP14_LOW, by=1)
#replicate correlation
repta <- cor (PP14_GS6$PP14TOP_1,PP14_GS6$PP14TOP_1.1, "pairwise.complete.obs")
> reptb <- cor (PP14_GS6$PP14TOP_2,PP14_GS6$PP14TOP_2.1, "pairwise.complete.obs")
> repla <- cor (PP14_GS6$PP14LOW_1,PP14_GS6$PP14LOW_1.1, "pairwise.complete.obs")
> replb <- cor (PP14_GS6$PP14LOW_2,PP14_GS6$PP14LOW_2.1, "pairwise.complete.obs")
#correlation 0.95
corr.t <- cor (pset_top$pset_topa.FREQ_TOP,pset_top$PP14_TOP.TOP_FREQ,use= "pairwise.complete.obs" )#correlation TOP
corr.l <- cor (pset_low$pset_lowa.FREQ_LOW,pset_low$PP14_LOW.LOW_FREQ,use= "pairwise.complete.obs" )#correlation LOW 
#
#3-Plots between pool frequency and indivudial frequency
#3.2-PP14
#3.2.1-TOP_A
plot (pset_topa$PP14TOP_1~pset_topa$TOPA_FREQ,xlab ="individual allele frequency" , ylab = "pooled allele frequency", main= "TOP_A", cex.main = 0.9)
abline(lm(pset_topa$PP14TOP_1~pset_topa$TOPA_FREQ), col="red")
par(mfrow = c(1,2))
mtext("Low group", side=3, outer=TRUE, line=-1)
summary(lm(pset_topa$PP14TOP_1~pset_topa$TOPA_FREQ))
#3.2.1-TOP_B
plot (pset_top$pset_topa.FREQ_TOP~pset_top$PP14_TOP.TOP_FREQ,xlab ="individual allele frequency" , ylab = "pooled allele frequency", main= "TOP_B",cex.main = 0.9)
abline(lm(pset_top$pset_topa.FREQ_TOP~pset_top$PP14_TOP.TOP_FREQ), col="red")
summary(lm(pset_topb$PP14TOP_2M~pset_topb$TOPB_FREQ))
#3.2.2-LOW_A
plot (pset_lowa$PP14LOW_1M~pset_lowa$LOWA_FREQ,xlab ="individual allele frequency" , ylab = "pooled allele frequency", main= "LOW_A",cex.main = 0.9)
abline(lm(pset_lowa$PP14LOW_1M~pset_lowa$LOWA_FREQ), col="red")
summary(lm(pset_lowa$PP14LOW_1M~pset_lowa$LOWA_FREQ))
#3.2.2-LOW_B
plot (pset_lowb$PP14LOW_2M~pset_lowb$LOWB_FREQ,xlab ="individual allele frequency" , ylab = "pooled allele frequency", main= "LOW_B",cex.main = 0.9)
abline(lm(pset_lowb$PP14LOW_2M~pset_lowb$LOWB_FREQ), col="red")
summary(lm(pset_lowb$PP14LOW_2M~pset_lowb$LOWB_FREQ))
#
##4.2-PP14-GWAS POOL
pp14_assoc<- pp14_assoc[-c(19261),]
pp14_gwas <- data.frame(pp14_assoc$Name,pp14_assoc$Chr,pp14_assoc$Position, pp14_assoc$p.value)
pp14_gwas <- pp14_gwas [order(pp14_gwas[,2]), ]
pp14_gwas <-edit(pp14_gwas)
pp14_gwas<- na.exclude(pp14_gwas)
pp14_gwas$P[pp14_gwas$P == 0] <- NA
pool.significativo <- pp14_p001$SNP
manhattan(pp14_gwas,suggestiveline = F,genomewideline = -log10(1e-05), ylim=range(1e-9:8))
pp14_p001 <- subset (pp14_gwas,subset = pp14_gwas$P <0.01)
qq(pp14_gwas$P)#Q-Q plot PP14
pp14_sign <- subset (pp14_gwas,subset = pp14_gwas$P < 10^-5)
pp14_sug <- subset (pp14_gwas,subset = pp14_gwas$P < 10^-3)
#
pp14_gwas$fdr <- p.adjust(pp14_gwas$P, method = "fdr")
POOL.bonf <- subset.data.frame(pp14_gwas, pp14_gwas$fdr < 0.05)
pp14_gwas$bonf <- pp14_gwas$fdr <- NULL
#
manhattan(pp14_gwas,suggestiveline = F,genomewideline = F,highlight = c("Hapmap25090-BTA-157559","BTB-00282168", "ARS-BFGL-NGS-21116"))
##4.3-PP14-GWAS INDIVIDUAL same animals- CHI-SQUARE
PP14_PLINK.txt <- read.table("~/Posgraduacao/IC/Dados allelotipagem/PP14_pools/PP14_PLINK.txt.assoc", header=TRUE, quote="\"")
PP14_PLINK.txt$CHR <- NULL
pp14_plink <- merge.data.frame (PP14_GS6,PP14_PLINK.txt, by=1)
pp14_plink$Address <- pp14_plink$PP14TOP_1 <-pp14_plink$PP14TOP_1.1 <-pp14_plink$PP14TOP_1M <- pp14_plink$PP14TOP_2 <- pp14_plink$PP14TOP_2.1 <- pp14_plink$PP14TOP_2M <- pp14_plink$PP14LOW_1 <- pp14_plink$PP14LOW_1.1 <- pp14_plink$PP14LOW_1M <- pp14_plink$PP14LOW_2 <- pp14_plink$PP14LOW_2.1 <-pp14_plink$PP14LOW_2M <- pp14_plink$BP <-pp14_plink$A1 <-pp14_plink$F_A <-pp14_plink$F_U <-pp14_plink$A2 <-pp14_plink$CHISQ <-pp14_plink$OR <- NULL
pp14_plink$P[pp14_plink$P ==0] <-NA
pp14_plink <- na.exclude(pp14_plink)#24.999 SNP
pp14_plink <- edit(pp14_plink)
manhattan (pp14_plink,suggestiveline = F, genomewideline = F,ylim=range(1e-8:7), highlight = qui.plink.significativo, main= "same animals GWAS-no adjust (snp sign <0.01= red point)")
qq(pp14_plink$P)
plink_p001 <- subset (pp14_plink,subset = pp14_plink$P <0.01)
qui.plink.significativo <- plink_p001$SNP
pp14_plink$bonf <- p.adjust(pp14_plink$P, method = "bonferroni")
TESTE <- val.FDR
val.FDR <- subset(TESTE, TESTE < 0.01)
#
pp14_p001 <- pp14_p001 [order(pp14_p001[,2]), ]
plink_p001 <- plink_p001 [order(plink_p001[,2]), ]
teste <- merge.data.frame(plink_p001,pp14_p001,by=1)
#
##4.4-PP14-GWAS validation 268 animal- CHI-SQUARE
chi.268.gwas <- data.frame (validacao_freq4$SNP,validacao_freq4$CHR,validacao_freq4$BP, validacao_freq4$P)#509.501SNP
chi.268.gwas <-edit(chi.268.gwas)
chi.268.gwas <- chi.268.gwas [order(chi.268.gwas[,2]), ]
chi.268.gwas$P[chi.268.gwas$P ==0] <-NA
chi.268.gwas <- na.exclude(chi.268.gwas)#509.501SNP
val.sign <- subset (validacao_freq4, subset = validacao_freq4$P < 0.01)
chi.268.significativo <- val.sign$SNP
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = chi.268.significativo, main= "286 chi-square GWAS-no adjust (snp sign <0.01= red point)")
qq(pp14_plink$P)
#pool signif snps
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = pool.significativo, main= "286 chi-square GWAS-no adjust (pool snp sign <0.01= red point)")
#individual signif snps
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = qui.plink.significativo, main= "286 chi-square GWAS-no adjust (ind snp sign <0.01= red point)")
#bonferroni correction
chi.268.bonferroni <- p.adjust(validacao_freq4$P, method = "bonferroni")
validacao_freq4$bonferroni <- chi.268.bonferroni
chi.268.bonferroni <- subset(validacao_freq4, validacao_freq4$bonferroni < 0.01)
chi.268.bonferroni <- chi.268.bonferroni$SNP
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = chi.268.bonferroni, main= "286 chi-square GWAS-bonferroni (snp sign <0.01= red point)")
#fdr correction
chi.268.fdr <- p.adjust(validacao_freq4$P, method = "fdr")
validacao_freq4$fdr <- chi.268.fdr
chi.268.fdr <- subset(validacao_freq4, validacao_freq4$fdr < 0.01)
chi.268.fdr <- chi.268.fdr$SNP
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = chi.268.fdr, main= "286 chi-square GWAS-fdr (snp sign <0.01= red point)")
#
##4.4-PP14-GWAS validation 268 animal- wald teste
#data montage
fdr268 <- p.adjust(qassoc.268animal$P,method = "fdr")
qassoc.268animal$fdr <- fdr268
bonferroni268<- p.adjust(qassoc.268animal$P, method ="bonferroni")
qassoc.268animal$bonf <- bonferroni268
q268.p <- subset(qassoc.268animal, subset= qassoc.268animal$P <0.01)
q268.fdr <- subset(qassoc.268animal, subset= qassoc.268animal$fdr <0.01)
q268.bonferroni <- subset(qassoc.268animal, subset= qassoc.268animal$bonf <0.01)
q268.p <- q268.p [order(q268.p[,2]), ]
snps.qplink.p <- snps.qplink.p [order(snps.qplink.p [,2]),]
q268.fdr <- q268.fdr [order(q268.fdr[,2]), ]
snps.qplink.fdr <- snps.qplink.fdr [order(snps.qplink.fdr [,2]),]
comum.268.1301.animal.fdr <- merge.data.frame (snps.qplink.fdr,q268.fdr, by=2)
q268.fdr$CHR <- NULL
pool.validacao.268 <- merge.data.frame (q268.fdr,pp14_p001, by=1)
#gwas
walt.268.gwas <- data.frame (qassoc.268animal$SNP,qassoc.268animal$CHR,qassoc.268animal$BP, qassoc.268animal$P)#735.291SNP
walt.268.gwas <-edit(walt.268.gwas)
walt.268.gwas <- walt.268.gwas [order(walt.268.gwas[,2]), ]
walt.268.gwas$P[walt.268.gwas$P ==0] <-NA
walt.268.gwas <- na.exclude(walt.268.gwas)#654.872SNP
watl.268significativo <- subset (qassoc.268animal, subset = qassoc.268animal$P < 0.01)
manhattan (walt.268.gwas, suggestiveline = F, genomewideline = F, highlight = pool.significativo, main= "286 walt-test GWAS-no adjust (pool snp sign <0.01= red point)")
#bonferroni correction
walt.268.bonferroni <- p.adjust(qassoc.268animal$P, method = "bonferroni")
qassoc.268animal$bonferroni <- walt.268.bonferroni
walt.268.bonferroni <- subset(qassoc.268animal, qassoc.268animal$bonferroni < 0.01)
walt.268.bonferroni <- walt.268.bonferroni$SNP
manhattan (walt.268.gwas, suggestiveline = F, genomewideline = F, highlight = walt.268.bonferroni, main= "286 walt-test GWAS-bonferroni (snp sign <0.01= red point)")
#fdr correction
walt.268.fdr <- p.adjust(qassoc.268animal$P, method = "fdr")
qassoc.268animal$fdr <- walt.268.fdr
walt.268.fdr <- subset(qassoc.268animal, qassoc.268animal$fdr < 0.01)
walt.268.fdr <- walt.268.fdr$SNP
manhattan (walt.268.gwas, suggestiveline = F, genomewideline = F, highlight = walt.268.fdr, main= "286 walt-square GWAS-fdr (snp sign <0.01= red point)")
#
##4.4-PP14-GWAS validation 1301 animal- wald teste
#data montage
fdr1301 <- p.adjust(qassoc.1301animal$P,method = "fdr")
qassoc.1301animal$fdr <- fdr1301
bonferroni1301<- p.adjust(qassoc.1301animal$P, method ="bonferroni")
qassoc.1301animal$bonf <- bonferroni1301
q1301.p <- subset(qassoc.1301animal, subset= qassoc.1301animal$P <0.01)
q1301.fdr <- subset(qassoc.1301animal, subset= qassoc.1301animal$fdr <0.01)
q1301.bonferroni <- subset(qassoc.1301animal, subset= qassoc.1301animal$bonf <0.01)
q1301.p <- q1301.p [order(q1301.p[,2]), ]
snps.qplink.p <- snps.qplink.p [order(snps.qplink.p [,2]),]
q1301.fdr <- q1301.fdr [order(q1301.fdr[,2]), ]
snps.qplink.fdr <- snps.qplink.fdr [order(snps.qplink.fdr [,2]),]
comum.1301.1301.animal.fdr <- merge.data.frame (snps.qplink.fdr,q1301.fdr, by=2)
q1301.fdr$CHR <- NULL
pool.validacao.1301 <- merge.data.frame (q1301.fdr,pp14_p001, by=1)
#gwas
walt.1301.gwas <- data.frame (qassoc.1301animal$SNP,qassoc.1301animal$CHR,qassoc.1301animal$BP, qassoc.1301animal$P)#735.291SNP
walt.1301.gwas <-edit(walt.1301.gwas)
walt.1301.gwas <- walt.1301.gwas [order(walt.1301.gwas[,2]), ]
walt.1301.gwas$P[walt.1301.gwas$P ==0] <-NA
walt.1301.gwas <- na.exclude(walt.1301.gwas)#693.911SNP
watl.1301significativo <- subset (qassoc.1301animal, subset = qassoc.1301animal$P < 0.01)
manhattan (walt.1301.gwas, suggestiveline = F, genomewideline = F, highlight = c("Hapmap25090-BTA-157559","BTB-00282168", "ARS-BFGL-NGS-21116"), main= "1301 walt-test GWAS-no adjust (pool snp sign <0.01= red point)")
#bonferroni correction
walt.1301.bonferroni <- p.adjust(qassoc.1301animal$P, method = "bonferroni")
qassoc.1301animal$bonferroni <- walt.1301.bonferroni
walt.1301.bonferroni <- subset(qassoc.1301animal, qassoc.1301animal$bonferroni < 0.01)
walt.1301.bonferroni <- walt.1301.bonferroni$SNP
manhattan (walt.1301.gwas, suggestiveline = F, genomewideline = F, highlight = walt.1301.bonferroni, main= "1301 walt-test GWAS-bonferroni (snp sign <0.01= red point)")
#fdr correction
walt.1301.fdr <- p.adjust(qassoc.1301animal$P, method = "fdr")
qassoc.1301animal$fdr <- walt.1301.fdr
walt.1301.fdr <- subset(qassoc.1301animal, qassoc.1301animal$fdr < 0.01)
walt.1301.fdr <- walt.1301.fdr$SNP
manhattan (walt.1301.gwas, suggestiveline = F, genomewideline = F, highlight = walt.1301.fdr, main= "1301 walt-square GWAS-fdr (snp sign <0.01= red point)")
#
###4.5-PP14-GWAS subset 1301 animal- wald teste
#data montage
subset.gwas <- data.frame (plink.subset$SNP,plink.subset$CHR,plink.subset$BP, plink.subset$P)#47782SNP
subset.gwas <-edit(subset.gwas)
subset.gwas <- subset.gwas [order(subset.gwas[,2]), ]
subset.gwas$P[subset.gwas$P ==0] <-NA
subset.gwas <- na.exclude(subset.gwas)#40.573SNP
subset.sign <- subset (subset.gwas, subset = subset.gwas$P < 0.01)
subset.significativo <- subset.sign$SNP
manhattan (subset.gwas, suggestiveline = F, genomewideline = F, highlight = subset.significativo, main= "1309 subset GWAS-no adjust (snp sign <0.01= red point)")
qq(pp14_plink$P)
#pool signif snps
manhattan (subset.gwas, suggestiveline = F, genomewideline = F, highlight = c ("BTB-00282168", "ARS-BFGL-NGS-21116"), main= "1309 subset GWAS-no adjust (pool snp sign <0.01= red point)")
subset.gwas$bonf <- p.adjust(subset.gwas$P, method = "bonferroni")
val.bonf <- subset.data.frame(subset.gwas, subset.gwas$bonf < 0.01)
subset.gwas$bonf <- NULL
#
manhattan(subset.gwas,suggestiveline = F,genomewideline = F,highlight = val.FDR, main="POOL GWAS-fdr (snp sign <0.01= red point)")

#individual signif snps
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = qui.plink.significativo, main= "286 chi-square GWAS-no adjust (ind snp sign <0.01= red point)")
#bonferroni correction
chi.268.bonferroni <- p.adjust(validacao_freq4$P, method = "bonferroni")
validacao_freq4$bonferroni <- chi.268.bonferroni
chi.268.bonferroni <- subset(validacao_freq4, validacao_freq4$bonferroni < 0.01)
chi.268.bonferroni <- chi.268.bonferroni$SNP
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = chi.268.bonferroni, main= "286 chi-square GWAS-bonferroni (snp sign <0.01= red point)")
#fdr correction
chi.268.fdr <- p.adjust(validacao_freq4$P, method = "fdr")
validacao_freq4$fdr <- chi.268.fdr
chi.268.fdr <- subset(validacao_freq4, validacao_freq4$fdr < 0.01)
chi.268.fdr <- chi.268.fdr$SNP
manhattan (chi.268.gwas, suggestiveline = F, genomewideline = F, highlight = chi.268.fdr, main= "286 chi-square GWAS-fdr (snp sign <0.01= red point)")
#
#
#5-VALIDATION (file validacao 1301-chi-square)
MAF_A <- subset.data.frame(validacao, subset= validacao$A1 == "A")#frequency of case and control
MAF_B <- subset.data.frame(validacao, subset= validacao$A1 == "B")
MAF_B$freq.B.case <- MAF_B$F_A
MAF_B$freq.B.control <- MAF_B$F_U
MAF_A$freq.B.case <- 1- MAF_A$F_A
MAF_A$freq.B.control <- 1- MAF_A$F_U
validacao_freq <- rbind.data.frame(MAF_A,MAF_B)
validacao_freq <- validacao_freq [order(validacao_freq[,2]), ]
validacao_freq1 <-subset(validacao_freq,subset = validacao_freq$freq.B.case <0.98)#Exclude fixed SNP
validacao_freq2 <-subset(validacao_freq1,subset = validacao_freq1$freq.B.control <0.98)#Exclude fixed SNP
validacao_freq3 <-subset(validacao_freq2,subset = validacao_freq2$freq.B.case >0.02)#Exclude fixed SNP
validacao_freq4 <-subset(validacao_freq3,subset = validacao_freq3$freq.B.control >0.02)#Exclude fixed SNP
rm(validacao_freq1,validacao_freq2, validacao_freq3)
val.sign <- subset (validacao_freq4, subset = validacao_freq4$P < 0.01)
val.sign$CHR <- NULL
val <- merge.data.frame (PP14_GS6,val.sign, by=1)#subset between array 50k and HD SNP
val.GWAS <- merge.data.frame(PP14_GS6,validacao_freq4, by=1)
val.sign <- val.sign [order(val.sign[,8]), ]
pp14_p001 <- pp14_p001 [order(pp14_p001[,1]), ]
pub <- merge.data.frame (val.sign,pp14_p001, by=1)
pub <- pub [order(pub[,3]), ]
val.gwas <- data.frame (val.GWAS$Name,val.GWAS$Chr,val.GWAS$BP,val.GWAS$P)
val.gwas <- edit (val.gwas)
val.gwas$P[val.gwas$P == 0] <- NA
val.gwas<- na.exclude(val.gwas)
manhattan(val.gwas,suggestiveline = F,genomewideline = F)
qq(val.gwas$P)
#
##5-VALIDATION (file validacao 1301-chi-square)
MAF_A <- subset.data.frame(validacao, subset= validacao$A1 == "A")#frequency of case and control
MAF_B <- subset.data.frame(validacao, subset= validacao$A1 == "B")
MAF_B$freq.B.case <- MAF_B$F_A
MAF_B$freq.B.control <- MAF_B$F_U
MAF_A$freq.B.case <- 1- MAF_A$F_A
MAF_A$freq.B.control <- 1- MAF_A$F_U
validacao_freq <- rbind.data.frame(MAF_A,MAF_B)
validacao_freq <- validacao_freq [order(validacao_freq[,2]), ]
validacao_freq1 <-subset(validacao_freq,subset = validacao_freq$freq.B.case <0.98)#Exclude fixed SNP
validacao_freq2 <-subset(validacao_freq1,subset = validacao_freq1$freq.B.control <0.98)#Exclude fixed SNP
validacao_freq3 <-subset(validacao_freq2,subset = validacao_freq2$freq.B.case >0.02)#Exclude fixed SNP
validacao_freq4 <-subset(validacao_freq3,subset = validacao_freq3$freq.B.control >0.02)#Exclude fixed SNP
rm(validacao_freq1,validacao_freq2, validacao_freq3)
val.sign <- subset (validacao_freq4, subset = validacao_freq4$P < 0.01)
val.sign$CHR <- NULL
val <- merge.data.frame (PP14_GS6,val.sign, by=1)#subset between array 50k and HD SNP
val.GWAS <- merge.data.frame(PP14_GS6,validacao_freq4, by=1)
val.sign <- val.sign [order(val.sign[,8]), ]
pp14_p001 <- pp14_p001 [order(pp14_p001[,1]), ]
pub <- merge.data.frame (val.sign,pp14_p001, by=1)
pub <- pub [order(pub[,3]), ]
val.gwas <- data.frame (val.GWAS$Name,val.GWAS$Chr,val.GWAS$BP,val.GWAS$P)
val.gwas <- edit (val.gwas)
val.gwas$P[val.gwas$P == 0] <- NA
val.gwas<- na.exclude(val.gwas)
manhattan(val.gwas,suggestiveline = F,genomewideline = F)
qq(val.gwas$P)
##5-VALIDATION USING PP14 EBV's(file plink.qassoc 1300 animals)
fdr <- p.adjust(plink.qassoc$P,method = "fdr")
plink.qassoc$fdr <- fdr
bonferroni <- p.adjust(plink.qassoc$P, method ="bonferroni")
plink.qassoc$bonf <- bonferroni
snps.qplink.p <- subset(plink.qassoc, subset= plink.qassoc$P <0.01)
snps.qplink.fdr <- subset(plink.qassoc, subset= plink.qassoc$fdr <0.01)
snps.qplink.bonferroni <- subset(plink.qassoc, subset= plink.qassoc$bonf <0.01)
snps.qplink.bonferroni <- snps.qplink.bonferroni [order(snps.qplink.bonferroni[,2]), ]
snps.qplink.bonferroni$CHR <- NULL
pool.validacao.todos <- merge.data.frame (snps.qplink.bonferroni,pp14_p001, by=1)

#
plink_p001 <- plink_p001 [order(plink_p001 [,1]),]
q268.p$CHR <- NULL
merge.quiq.t <- merge(q268.p,plink_p001, by=1)
hhh <- read.table(file = "../Desktop/hhh.txt",header = T, skip = 3)
trace ("manhattan", edit=T)
