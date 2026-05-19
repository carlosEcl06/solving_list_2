#IBI5086 - 2026
#Lista 2

#Simulando dados com efeito de SNPs

#Modelo Uniloco: Gerar Y com efeito de um único loco de SNP

n=400
x1<-rbinom(n,2,0.32)
x2<-rbinom(n,2,0.25)
table(x1,x2)

gera_phen <- function (loci,j){
 y <- rep(NA,nrow(loci))
 for(i in 1:nrow(loci)){
 if (loci[i,j] == 0) y[i]<-rnorm(1,mi-a,s2)
 if (loci[i,j] == 1) y[i]<-rnorm(1,mi+d,s2)
 if (loci[i,j] == 2) y[i]<-rnorm(1,mi+a,s2)}
 cbind(loci[,j],y)
}

mi=100; a=10; d=2; s2=4
phen1<-gera_phen(cbind(x1,x2),1)
mi=100; a=2; d=8; s2=6
phen2<-gera_phen(cbind(x1,x2),2)

boxplot(phen1[,2],phen2[,2])

tapply(phen1[,2],phen1[,1],mean)
tapply(phen1[,2],phen1[,1],sd)
boxplot(phen1[,2]~ phen1[,1]) #ef. aditivo do SNP


tapply(phen2[,2],phen2[,1],mean)
tapply(phen2[,2],phen2[,1],sd)
boxplot(phen2[,2]~ phen2[,1]) #há ef. de dominância do SNP

#Construa o gráfico de perfis de médias
#Ajustar modelos lineares para os efeitos individuais desses SNPs
#Interpretar os resultados do ajuste do modelo linear
#Os parâmetros foram recuperados sob este ajuste?


##Modelo biloco: Gerar Y com efeito dos dois SNPs
library(car)
xa1<-Recode(x1,"0=(-1);1=0;2=1")
xd1<-Recode(x1,"0=0;1=1;2=0")
xa2<-Recode(x2,"0=(-1);1=0;2=1")
xd2<-Recode(x2,"0=0;1=1;2=0")

mi=110
a1=20; d1=0
a2=5; d2=4
a1a2=20; a1d2=22
d1a2=5; d1d2=0
s2=25

phen <- function(xa1,xd1,xa2,xd2){
 y<-matrix(NA,length(xa1),1)
 for(i in 1:length(xa1)){
 y[i]<-
rnorm(1,mi+a1*xa1[i]+d1*xd1[i]+a2*xa2[i]+d2*xd2[i]
       +a1a2*(xa1[i]*xa2[i])+a1d2*(xa1[i]*xd2[i])
       +d1a2*(xd1[i]*xa2[i])+d1d2*(xd1[i]*xd2[i]))
}
return(y)
}

ymod=phen(xa1,xd1,xa2,xd2)

fit1<- lm(ymod~xa1+xd1+xa2+xd2+xa1*xa2+xa1*xd2+xd1*xa2+xd1*xd2)
summary(fit1)

#Interprete
#Construa o gráfico com os perfis de médias
#Compare os par6ametros e as estimativas: o ajuste permitiu a recuperação dos parâmetros?


##Lista 2 
#Gerando Dados POPULACAO 1
r<-25
sigma<-6

#Pop1
mi01<-80
mi02<-90
mi11<-80
mi12<-90
mi21<-80
mi22<-90

#Gerando Dados
set.seed=9519
resp01<-rnorm(r,mi01,sigma)
set.seed=9521
resp02<-rnorm(r,mi02,sigma)
set.seed=9523
resp11<-rnorm(r,mi11,sigma)
set.seed=9525
resp12<-rnorm(r,mi12,sigma)
set.seed=9527
resp21<-rnorm(r,mi21,sigma)
set.seed=9529
resp22<-rnorm(r,mi22,sigma)

resp<-stack(data.frame(cbind(resp01,resp02,resp11,resp12,resp21,resp22)))
resp
fg<-rep(0:2,each=50)
fg
ft<-rep(rep(1:2,each=25),3)
ft
dat1<-data.frame(cbind(resp,fg,ft))
dat1
str(dat1)
names(dat1)<-c("y","grup","fg","ft")
ft<-factor(ft)
fg<-factor(fg)

dat<-dat1
attach(dat)

boxplot(y~fg,main="Efeito Principal Genético - P1")
boxplot(y~ft,main="Efeito Principal de Trat - P1")
interaction.plot(fg,ft,y,main="Perfis de Médias: Interação - P1")
#Comentar (descritivamente) o efeito de interação
#Lembrar que, se há efeito de interação, os efeitos principais ficam comprometidos

fit1<-aov(y~fg*ft)
summary(fit1)
fit1$coefficients
#Interpretar a significância dos efeitos testados

#Equivalentemente
fit1<-lm(y~fg*ft)
anova(fit1)
summary(fit1)

#Analise de resíduos do modelo
plot(fit1$fit, fit1$res, xlab="Valores Ajustados", ylab="Resíduos")
abline(h=0)
title("Resíduos vs Preditos")

qqnorm(fit1$residuals,ylab="Residuos", main=NULL)
qqline(fit1$residuals)
title("Gráfico Normal de Probabilidade dos Resíduos")

# Teste de normalidade de Shapiro-Wilk
shapiro.test(fit1$residuals)

#Teste da homogeneidade das variâncias
bartlett.test(y,grup) 

#Ajuste do modelo reduzido
fit2<-aov(y~ft)
summary(fit2) #existe efeito principal de tratamento: miT1 < miT2
fit2$coefficients

fit2<-lm(y~ft)
summary(fit2)

###Gerando Dados POPULACAO 2
r<-25
sigma<-6

#Pop2
mi01<-75
mi02<-80
mi11<-80
mi12<-86
mi21<-90
mi22<-97

set.seed=9519
resp01<-rnorm(r,mi01,sigma)
set.seed=9521
resp02<-rnorm(r,mi02,sigma)
set.seed=9523
resp11<-rnorm(r,mi11,sigma)
set.seed=9525
resp12<-rnorm(r,mi12,sigma)
set.seed=9527
resp21<-rnorm(r,mi21,sigma)
set.seed=9529
resp22<-rnorm(r,mi22,sigma)

resp<-stack(data.frame(cbind(resp01,resp02,resp11,resp12,resp21,resp22)))
resp
fg<-rep(0:2,each=50)
fg
ft<-rep(rep(1:2,each=25),3)
ft
dat2<-data.frame(cbind(resp,fg,ft))
dat2
str(dat2)
names(dat2)<-c("y","grup","fg","ft")

ft<-factor(ft)
fg<-factor(fg)
dat2

dat<-dat2
attach(dat)

boxplot(y~fg,main="Efeito Genético - P2")
boxplot(y~ft,main="Efeito de Trat - P2")
interaction.plot(fg,ft,y,main="Perfis de Médias - P2")

fit1<-aov(y~fg*ft)
summary(fit1)
fit1$coefficients

#Equivalentemente
fit1<-lm(y~fg*ft)
anova(fit1)
summary(fit1)

#Analise de resíduos do modelo
plot(fit1$fit, fit1$res, xlab="Valores Ajustados", ylab="Resíduos")
abline(h=0)
title("Resíduos vs Preditos")

qqnorm(fit1$residuals,ylab="Residuos", main=NULL)
qqline(fit1$residuals)
title("Gráfico Normal de Probabilidade dos Resíduos")

# Teste de normalidade de Shapiro-Wilk
shapiro.test(fit1$residuals)

#Teste da homogeneidade das variâncias
bartlett.test(y,grup) 

#Ajuste do modelo reduzido
fit2<-aov(y~fg+ft)
summary(fit2) 
#existe efeito principal genético e de tratamento
fit2$coefficients

fit2<-lm(y~fg+ft)
summary(fit2)

#Estudo dos efeitos principais do fator Trat
#Tukey somente para o comando "aov"
fit2.t <- TukeyHSD(fit2, "fg")
fit2.t
plot(fit2.t)

fit2.t <- TukeyHSD(fit2, "ft")
fit2.t
plot(fit2.t)
#Interprete as diferenças e os IC

##Gerando dados de P3
r<-25
sigma<-6

#Pop3
mi01<-75
mi02<-75
mi11<-80
mi12<-85
mi21<-85
mi22<-100

set.seed=9519
resp01<-rnorm(r,mi01,sigma)
set.seed=9521
resp02<-rnorm(r,mi02,sigma)
set.seed=9523
resp11<-rnorm(r,mi11,sigma)
set.seed=9525
resp12<-rnorm(r,mi12,sigma)
set.seed=9527
resp21<-rnorm(r,mi21,sigma)
set.seed=9529
resp22<-rnorm(r,mi22,sigma)

resp<-stack(data.frame(cbind(resp01,resp02,resp11,resp12,resp21,resp22)))
resp
fg<-rep(0:2,each=50)
fg
ft<-rep(rep(1:2,each=25),3)
ft
dat3<-data.frame(cbind(resp,fg,ft))
dat3
str(dat3)
names(dat3)<-c("y","grup","fg","ft")
ft<-factor(ft)
fg<-factor(fg)
dat3

dat<-dat3
attach(dat)

boxplot(y~fg,main="Efeito Genético - P3")
boxplot(y~ft,main="Efeito de Trat - P3")
interaction.plot(fg,ft,y,main="Perfis de Médias - P3")

fit1<-aov(y~fg*ft)
summary(fit1)
fit1$coefficients

#Equivalentemente
fit1<-lm(y~fg*ft)
anova(fit1)
summary(fit1)

#Analise de resíduos do modelo
plot(fit1$fit, fit1$res, xlab="Valores Ajustados", ylab="Resíduos")
abline(h=0)
title("Resíduos vs Preditos")

qqnorm(fit1$residuals,ylab="Residuos", main=NULL)
qqline(fit1$residuals)
title("Gráfico Normal de Probabilidade dos Resíduos")

# Teste de normalidade de Shapiro-Wilk
shapiro.test(fit1$residuals)

#Teste da homogeneidade das variâncias
bartlett.test(y,grup) 

#Estudo do efeito de interação entre os fatores

#Método de Tukey (somente para aov)
fit1.t <- TukeyHSD(fit1, "fg:ft")
fit1.t
plot(fit1.t)

#Neste caso, Tukey pode estar testando pares de grupos que não são de interesse
#Ex:Vamos então testar somente diferença entre Tratamentos para cada genótipo

#Função para calcular IC e valor-p para contrastes entre médias
ci=function(dif,qmres,gl,n1,n2,alpha) {   
   se=sqrt(qmres*(1/n1 + 1/n2))   
   tval=dif/se
   dt <- qt(1-alpha/2, gl) * se
   pval <- 1-pt(abs(tval),gl)
   low=dif-dt
   up=dif+dt
   m=cbind(dif,se,low,up,pval)
   dimnames(m)[[2]]=c("estimate","se",paste(100*(1-alpha),"% Conf.",sep=""),"limits","p-val")
   m 
   } 

mig<-tapply(y, list(fg,ft), mean)
mig
qmres<-anova(fit1)$"Mean Sq"[4]
qmres
gl<-anova(fit1)$"Df"[4]
gl

#fg=0: miT1-miT2
dif <- mig[1,1]-mig[1,2]
n1 <- length(dat[fg==0 & ft==1,]$y)
n2 <- length(dat[fg==0 & ft==2,]$y)
c1<-ci(dif,qmres,gl,n1,n2,0.05)
c1
p1<-c1[5]

#fg=1: miT1-miT2
dm <- mig[2,1]-mig[2,2]
n1 <- length(dat[fg==1 & ft==1,]$y)
n2 <- length(dat[fg==1 & ft==2,]$y)
tt <- dm/sqrt(qmres*(1/n1 + 1/n2))
c2<-ci(dm,qmres,gl,n1,n2,alpha=0.05)
c2
p2<-c2[5]

#fg=2: miT1-miT2
dm <- mig[3,1]-mig[3,2]
n1 <- length(dat[fg==2 & ft==1,]$y)
n2 <- length(dat[fg==2 & ft==2,]$y)
c3<-ci(dm,qmres,gl,n1,n2,0.05)
c3
p3<-c3[5]

results <- c(p1,p2,p3)
results

#install.packages("gtools")
library(gtools)

results.o <- mixedsort(results, scientific=TRUE, decreasing=TRUE)
results.o  #útil para ordenar sob notação científica dos valores

# Obtendo valores-p ajustados 
p.adjust.methods
adjustb <- p.adjust(results.o,method="bonferroni")
cbind(results.o,adjustb) 

adjustfdr <- p.adjust(results.o,method="fdr")
cbind(results.o,adjustfdr)

adjusth <- p.adjust(results.o,method="holm")
cbind(results.o,adjusth)

cbind(results.o,adjustb,adjusth,adjustfdr)

##Qual é a conclusão para alfa global=5%?
#Ajuste por Bonferroni: mi01=mi02 mi11<mi12 mi21<mi22
#Ajuste por Holm:       mi01=mi02 mi11<mi12 mi21<mi22
#Ajuste por FDR:        mi01=mi02 mi11<mi12 mi21<mi22



#Combinando os dados das diferentes Populações
datp<-data.frame(rbind(dat1,dat2,dat3))
dim(datp)
pop<-c(rep(1,150),rep(2,150),rep(3,150))
pop
datp<-data.frame(cbind(datp,pop))
head(datp)
str(datp)
attach(datp)
datp$fg<-factor(datp$fg)
datp$ft<-factor(datp$ft)
datp$pop<-factor(datp$pop)
str(datp)


##Combinando os dados das diferentes Populações: Alternativa 1
fitp1<-aov(y~pop+fg+ft+fg*ft,data=datp)
anova(fitp1)

#Note que se os fatores fg,ft e pop não foram definidos como "factor"
#o R considerará como variáveos quantitativas e na anova eles entrarão com
#1 grau de liberdade (efeito linear), o que não é possível para ft e pop
#Para fg, podemos considerar como uma variável ordinal (com 2 graus de liberdade), 
#ou podemos também considerar quantitativa discreta, e neste caso,
#o efeito linear do SNP pode ser uma alternativa viável (1 grau de liberdade)

#Analise de resíduos do modelo
plot(fitp1$fit, fitp1$res, xlab="Valores Ajustados", ylab="Resíduos")
abline(h=0)
title("Resíduos vs Preditos")

qqnorm(fitp1$residuals,ylab="Residuos", main=NULL)
qqline(fitp1$residuals)
title("Gráfico Normal de Probabilidade dos Resíduos")

# Teste de normalidade de Shapiro-Wilk
shapiro.test(fitp1$residuals)

#Teste da homogeneidade das variâncias
bartlett.test(fitp1$residuals ~ factor(round(fitp1$fit,3)))
#alfa=1%: não rejeitar

#Neste caso, o efeito principal de pop não é significante (a 5%) MAS
#mesmo não sendo significante é importante ficar no modelo para realizarmos
#o ajuste por população ao combinar os dados

#Mas, nos dados combinados podemos eliminar o efeito de interação que 
#não é significante e ajustar o modelo mais reduzido

fitp2<-aov(y~pop+fg+ft,data=datp)
anova(fitp2)

#Estudo do efeito principal do fator fg para os dados combinados e ajustados por pop
#Tukey para o efeito do SNP
fit2.t <- TukeyHSD(fitp2, "fg")
fit2.t
plot(fit2.t) #Conclusão: mi(aa) < mi(Aa) < mi(AA)

fit2.t <- TukeyHSD(fitp2, "ft")
fit2.t
plot(fit2.t) #Conclusão: mi(T1) < mi(T2)

#O efeito de ft nos dados combinados nem precisa ser testado por Tukey já que 
#o teste pode ser feito na própria tabela de anova
#a não ser para conseguirmos o IC


boxplot(datp$y~datp$fg,main="Efeito Genético: P1+P2+P3")
boxplot(datp$y~datp$ft,main="Efeito de Trat: P1+P2+P3")
interaction.plot(datp$fg,datp$ft,datp$y,main="Perfis de Médias: P1+P2+P3")


#Combinando os dados das diferentes Populações: Alternativa 2
fitp3<-aov(y~pop+fg+ft+fg*ft*pop,data=datp)
anova(fitp3) 
#O efeito de interação de ordem 2 (entre os 3 fatores) é significante
#Neste caso, não há como reduzir o modelo (denominado saturado)

#Analise de resíduos do modelo
plot(fitp3$fit, fitp2$res, xlab="Valores Ajustados", ylab="Resíduos")
abline(h=0)
title("Resíduos vs Preditos")

qqnorm(fitp3$residuals,ylab="Residuos", main=NULL)
qqline(fitp3$residuals)
title("Gráfico Normal de Probabilidade dos Resíduos")

#Agora, para estudar o efeito de interação pop:fg:ft
#Vamos construir constrastes entre grupos de interesse
#Parametrização adotada pelo R no modelo anova ajustado

#Função para calcular IC e valor-p para contrastes entre médias
#a partir de uma ANOVA
ci=function(dif,qmres,gl,n1,n2,alpha) {   
   se=sqrt(qmres*(1/n1 + 1/n2))   
   tval=dif/se
   dt <- qt(1-alpha/2, gl) * se
   pval <- 1-pt(abs(tval),gl)
   low=dif-dt
   up=dif+dt
   m=cbind(dif,se,low,up,pval)
   dimnames(m)[[2]]=c("estimate","se",paste(100*(1-alpha),"% Conf.",sep=""),"limits","p-val")
   m 
   } 


fitp3$coefficients
mi123<-tapply(datp$y, list(datp$fg,datp$ft,datp$pop), mean)
mi123
qmres<-anova(fitp3)$"Mean Sq"[8]
qmres
gl<-anova(fitp3)$"Df"[8]
gl

#P1 fg=0,1,2: miT1-miT2
dif <- mi123[1,1,1]-mi123[1,2,1]
n1 <- length(datp[pop==1 & fg==0 & ft==1,]$y)
n2 <- length(datp[pop==1 & fg==0 & ft==2,]$y)
c1<-ci(dif,qmres,gl,n1,n2,0.05)
c1
p1<-c1[5]

dif <- mi123[2,1,1]-mi123[2,2,1]
n1 <- length(datp[pop==1 & fg==1 & ft==1,]$y)
n2 <- length(datp[pop==1 & fg==1 & ft==2,]$y)
c2<-ci(dif,qmres,gl,n1,n2,0.05)
c2
p2<-c2[5]

dif <- mi123[3,1,1]-mi123[3,2,1]
n1 <- length(datp[pop==1 & fg==2 & ft==1,]$y)
n2 <- length(datp[pop==1 & fg==2 & ft==2,]$y)
c3<-ci(dif,qmres,gl,n1,n2,0.05)
c3
p3<-c3[5]

#P2 fg=0,1,2: miT1-miT2
dif <- mi123[1,1,2]-mi123[1,2,2]
n1 <- length(datp[pop==2 & fg==0 & ft==1,]$y)
n2 <- length(datp[pop==2 & fg==0 & ft==2,]$y)
c4<-ci(dif,qmres,gl,n1,n2,0.05)
c4
p4<-c4[5]

dif <- mi123[2,1,2]-mi123[2,2,2]
n1 <- length(datp[pop==2 & fg==1 & ft==1,]$y)
n2 <- length(datp[pop==2 & fg==1 & ft==2,]$y)
c5<-ci(dif,qmres,gl,n1,n2,0.05)
c5
p5<-c5[5]

dif <- mi123[3,1,2]-mi123[3,2,2]
n1 <- length(datp[pop==2 & fg==2 & ft==1,]$y)
n2 <- length(datp[pop==2 & fg==2 & ft==2,]$y)
c6<-ci(dif,qmres,gl,n1,n2,0.05)
c6
p6<-c6[5]

#P3 fg=0,1,2: miT1-miT2
dif <- mi123[1,1,3]-mi123[1,2,3]
n1 <- length(datp[pop==3 & fg==0 & ft==1,]$y)
n2 <- length(datp[pop==3 & fg==0 & ft==2,]$y)
c7<-ci(dif,qmres,gl,n1,n2,0.05)
c7
p7<-c7[5]

dif <- mi123[2,1,3]-mi123[2,2,3]
n1 <- length(datp[pop==3 & fg==1 & ft==1,]$y)
n2 <- length(datp[pop==3 & fg==1 & ft==2,]$y)
c8<-ci(dif,qmres,gl,n1,n2,0.05)
c8
p8<-c8[5]

dif <- mi123[3,1,3]-mi123[3,2,3]
n1 <- length(datp[pop==3 & fg==2 & ft==1,]$y)
n2 <- length(datp[pop==3 & fg==2 & ft==2,]$y)
c9<-ci(dif,qmres,gl,n1,n2,0.05)
c9
p9<-c9[5]

results <- c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
results

#install.packages("gtools")
library(gtools)

results.o <- mixedsort(results, scientific=TRUE, decreasing=TRUE)
results.o  #útil para ordenar sob notação científica dos valores

# Obtendo valores-p ajustados 
p.adjust.methods
adjustb <- p.adjust(results.o,method="bonferroni")
cbind(results.o,adjustb) 

adjustfdr <- p.adjust(results.o,method="fdr")
cbind(results.o,adjustfdr)

adjusth <- p.adjust(results.o,method="holm")
cbind(results.o,adjusth)

cbind(results.o,adjustb,adjusth,adjustfdr)


##Combinando os dados das diferentes Populações: Alternativa 3
##Considerando População como um Fator Aleatório (não vimos no curso): 
## Neste caso, P1, P2 e P3 são consideradas uma amostra aleatória das possíveis populações que poderiam fazer parte do estudo 
##Ajuste do Modelo Linear Misto: fatores fixos (fg e ft) e aleatórios (pop)
head(datp)
dim(datp)
grupn<-rep(rep(seq(1,6),each=25),3)
grupn #efeito combinado: 3x2=6 grupos no total
interac<- rep(seq(1,18),each=25)
interac
datp<-data.frame(cbind(datp,grupn,interac))
head(datp)
str(datp)
attach(datp)
with(datp, interaction.plot(factor(grupn),factor(pop),y,main="Gráfico de interação"))

install.packages("lme4")
library(lme4)

fit.mlm<- lmer(y ~ factor(grupn) + (1|pop), data=datp)
fit.mlm
#modela efeito aleaório somente para o fator pop (uma amostra aleatório de 3 níveis das possíveis pop do estudo)
confint(fit.mlm, method="Wald") #IC só para efeito fixos
#nenhum IC dos ef fixos inclui o zero: há efeito significante
confint(fit.mlm) ##solução perfilada para cálculo do IC
confint(fit.mlm, method="boot", nsim=1000) #IC Bootstrap
#este efeito aleatório pode ser considerado nulo: a dependência entre respostas dentro de pop é desprezível

#Vamos investigar um modelo mais sofisticado: além de um componente aleatório de pop há um ef adicional devido a pop*grupn
fit.mlm<- lmer(y ~ factor(grupn) + (1|pop) + (1|interac), data=datp)
fit.mlm
summary(fit.mlm)
anova(fit.mlm)

#Intervalos de Confiança aproximados
confint(fit.mlm, method="Wald") #IC só para efeito fixos
confint(fit.mlm) ##solução perfilada
confint(fit.mlm, method="boot", nsim=1000) #IC Bootstrap

fixef(fit.mlm) #estimativa de mi e tau
#Variância dos ef. fixos (mihat)
RX <- getME(fit.mlm, "RX") 
sigma2 <- sigma(fit.mlm)^2  #sigma_e
vbeta<-sigma2*chol2inv(RX)
vbeta  ##na diagonal: var de mihat e tauhat

ranef(fit.mlm) #preditores dos efeitos aleatórios
VarCorr(fit.mlm)  #sd dos ef. aleatórios
print(vc <- VarCorr(fit.mlm), comp = c("Variance","Std.Dev."))

#Modelo estrutural (matricial): Y = Xbeta + Zu + e
#X e Z são amtrizes de planejamento dos efeitos fixos e aleatórios
#beta é o vetor de parâmetros (efeitos fixos)
#u é vetor dos efeitos aleatórios

model.matrix(fit1,data=dat)
Y<-as.vector(getME(fit1,"y"))  #vetor de resposta
Y
X<-as.matrix(getME(fit1,"X"))  #matriz de planejamento dos ef fixos
X
Z<-(as.matrix(getME(fit1,"Z"))) #matriz de planejamento dos ef aleatórios
Z

#Análise de Resíduos
plot(fit.mlm) ##residuo condicional ehat
plot(fit.mlm,type=c("p","smooth")) ## resíduo x ajustado
qqnorm(residuals(fit.mlm))
qqline(residuals(fit.mlm))
#veja outros gráficos no material de Aula

##Comparações múltiplas 
install.packages("emmeans")
library(emmeans)
emm1<-emmeans(fit.mlm,pairwise ~ factor(grupn),data=datp)
emm1

#comparação dos modelos via o ajuste de máxima verossimilanças
#Note que
fit0.mlm<- lmer(y ~ factor(grupn) + (1|pop), data=datp)
fit1.mlm<- lmer(y ~ factor(grupn) + (1|pop) + (1|interac), data=datp)
anova(fit0.mlm,fit1.mlm)
#Rej H0: o modelo fit1 se ajusta melhor que o fit0

#Ajuste também:
fit.mlm<- lmer(y ~ factor(grupn) +  (1|interac), data=datp)
fit.mlm
summary(fit.mlm)