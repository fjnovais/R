##############################
#  PAPER ANALYSIS POOL       #
##############################
#
#by Francisco Novais
#
pp14_qui <- edit(pp14_qui)
pp14_qui$P[pp14_qui$P == 0] <- NA
pp14_qui<- na.exclude(pp14_qui)
pp14_p001 <- subset(pp14_qui,subset= pp14_qui$P < 0.01)
pool.sign <- pp14_p001$SNP
manhattan(pp14_qui,suggestiveline = F,genomewideline = F, main="POOL GWAS-no adjust")
#correction
POOL.FDR <- p.adjust(pp14_qui$P, method = "fdr")
TESTE <- POOL.FDR
POOL.FDR <- subset(TESTE, subset= TESTE < 0.01 )
qq(pp14_qui$P, main="Q-Q plot pool")#Q-Q plot PP14
#
#1309 animais-subset (40k snps)
subset_plink <- edit(subset_plink)
subset_plink$P[subset_plink$P == 0] <- NA
subset_plink<- na.exclude(subset_plink)
manhattan(subset_plink,suggestiveline = F,genomewideline = F,highlight = pool.sign, main="subset 1309/GWAS-no adjust (pool sign <0.01= red point)")
#correction
subset.bonf <- subset_plink
subset.bonf$bonferroni <- p.adjust(subset_plink$P, method = "bonferroni")
subset.bonf <- subset(subset.bonf, subset= subset.bonf$bonferroni < 0.01 )
subset.bonf <- subset.bonf$SNP
manhattan(subset_plink,suggestiveline = F,genomewideline = F,highlight = subset.bonf, main="subset 1309/GWAS-bonferroni (snp sign <0.01= red point)")
qq(subset_plink$P, main="Q-Q plot pool")#Q-Q plot PP14
#
#1309 animais-HD (700k snps)
hd_plink <- edit(hd_plink)
hd_plink$P[hd_plink$P == 0] <- NA
hd_plink<- na.exclude(hd_plink)
manhattan(hd_plink,suggestiveline = F,genomewideline = F,highlight = pool.sign, main="hd 1309/GWAS-no adjust (pool sign <0.01= red point)")
#correction
hd.bonf <- hd_plink
hd.bonf$bonferroni <- p.adjust(hd_plink$P, method = "bonferroni")
hd.bonf <- subset(hd.bonf, subset= hd.bonf$bonferroni < 0.01 )
hd.bonf <- hd.bonf$SNP
manhattan(hd_plink,suggestiveline = F,genomewideline = F,highlight = hd.bonf, main="hd 1309/GWAS-bonferroni (snp sign <0.01= red point)")
qq(hd_plink$P, main="Q-Q plot pool")#Q-Q plot PP14
#
#268 animais-HD (700k snps)
chi268_plink <- edit(chi268_plink)
chi268_plink$P[chi268_plink$P == 0] <- NA
chi268_plink<- na.exclude(chi268_plink)
manhattan(chi268_plink,suggestiveline = F,genomewideline = F,highlight = pool.sign, main="chi268 1309/GWAS-no adjust (pool sign <0.01= red point)")
#correction
chi268.bonf <- chi268_plink
chi268.bonf$bonferroni <- p.adjust(chi268_plink$P, method = "bonferroni")
chi268.bonf <- subset(chi268.bonf, subset= chi268.bonf$bonferroni < 0.01 )
chi268.bonf <- chi268.bonf$SNP
manhattan(chi268_plink,suggestiveline = F,genomewideline = F,highlight = chi268.bonf, main="chi268 1309/GWAS-bonferroni (snp sign <0.01= red point)")
qq(chi268_plink$P, main="Q-Q plot pool")#Q-Q plot PP14
#
#268 animais-usando wald-test-HD (700k snps)
walt268_plink <- edit(walt268_plink)
walt268_plink$P[walt268_plink$P == 0] <- NA
walt268_plink<- na.exclude(walt268_plink)
manhattan(walt268_plink,suggestiveline = F,genomewideline = F,highlight = pool.sign, main="walt268/GWAS-no adjust (pool sign <0.01= red point)")
#correction
walt268.bonf <- walt268_plink
walt268.bonf$bonferroni <- p.adjust(walt268_plink$P, method = "bonferroni")
walt268.bonf <- subset(walt268.bonf, subset= walt268.bonf$bonferroni < 0.01 )
walt268.bonf <- walt268.bonf$SNP
manhattan(walt268_plink,suggestiveline = F,genomewideline = F,highlight = walt268.bonf, main="walt268/GWAS-bonferroni (snp sign <0.01= red point)")
qq(walt268_plink$P, main="Q-Q plot pool")#Q-Q plot PP14
#
#FOR ANNOTATION
#POOL
write.table(pp14_p001,"../annotation.pool.txt")
#1309 subset
write.table(subset.bonf,"../annotation.subset1309.txt")
#1309 hd
write.table(hd.bonf,"../annotation.hd1309.txt")