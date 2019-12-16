# *******************************************************************
# Procedimento:  PCR = Regressão via Componentes principais
# Entrada: C:\\dados\\ScriptsR\\saidas\trabalhoF$$$.dat
#          C:\\dados\\ScriptsR\\saidas\trabalhoG$$$.dat
#  
# Saída: C:\\dados\\ScriptsR\\saidas\saida$$$.dat
#        C:\dados\ScriptsR\Saidas\grafico.pdf
# Autor: C.D.Cruz
# Data: 06/2013
# *******************************************************************


# *******************************************************************
# 1. Definição da trilha de dados
# *******************************************************************
library(rrBLUP)
library (HapEstXXR)
library(pls)
library(stringr)
setwd("C:\\dados\\ScriptsR\\Saidas")


# *******************************************************************
# 2. Leitura do arquivo auxiliar e espectros
# *******************************************************************
Aux<-read.table("auxiliar$$$.dat",header=F)
valork= Aux[1,1]
ndim= Aux[2,1]


NomeAqdados<-read.table("trabalho$$$.dat",header=F)
aqT =NomeAqdados[1,1]

nir<-read.table(levels(aqT)[1],header=T, check.names = TRUE)
dados<-read.table(levels(aqT)[2],header=T, check.names = TRUE)
colnames(nir) <- str_sub(colnames(nir), start = 2)




#transformar os dados do espectro  em uma matriz
nir=as.matrix(nir)


#head(snp)
#tail(snp)





# *******************************************************************
# 4. Análise das variáveis
# *******************************************************************
sink("saida$$$.dat")
pdf(file='grafico.pdf')



cat("**********************************************************", "\n")
cat("**                                                      **  ", "\n")
cat("**            Programa Genes - Programa R               **", "\n")
cat("**                                                      **  ", "\n")
cat("**          NIR = Regressão via Componentes principais  **" ,"\n")
cat("**********************************************************", "\n")



AuxBD<-read.table("basedado$$$.dat",header=T)
cat("Data da Analise:", date(),"\n")
cat("Base de dados : ", colnames(AuxBD)[1] , "\n")

cat("**********************************************************", "\n")


valorR2t <- matrix(rep(NA),ncol(dados)-1,valork)
valorR2v <- matrix(rep(NA),ncol(dados)-1,valork)
valorCP <- matrix(rep(NA),ncol(dados)-1,valork)
valorbeta1 <- matrix(rep(NA),ncol(dados)-1,valork)
valorh2gen <- matrix(rep(NA), ncol(dados)-1,valork)
valorrggl <- matrix(rep(NA), ncol(dados)-1,valork)
valorrggr <- matrix(rep(NA), ncol(dados)-1,valork)
valorerrot <- matrix(rep(NA),ncol(dados)-1,valork)
valorerrov <- matrix(rep(NA),ncol(dados)-1,valork)

valorR2ajt <- matrix(rep(NA), ncol(dados)-1,valork)
valorR2ajv <- matrix(rep(NA), ncol(dados)-1,valork)

for(i in 2:ncol(dados)){


nindTotal = nrow(dados)
intervalo = round(nindTotal / valork)
inix = 1
#********************************************
for(kkuu in 1:valork){
iniy = inix
fimy = iniy + intervalo - 1
inix = fimy + 1
if (fimy > nindTotal) {
fimy =  nindTotal
}

sorteio <- NULL
ll=0
for(j1 in 1:nindTotal){
	if (j1 < iniy  || j1 > fimy){
	ll=ll+1
	sorteio[ll]= j1
	}
}
falta=nindTotal-ll
for(j2 in 1:falta){
	sorteio[ll+j2]= iniy +j2-1
}
#********************************************

#print(sorteio)

M<-nir[t(sorteio),]

fenotipo<-dados[t(sorteio),]

nindVal = fimy-iniy +1
nindTeste = nindTotal -  nindVal


varsel <- i
variavel<-(fenotipo[1: nindTeste,varsel])
ini = nindTeste+1
variavelval<-(fenotipo[ini:nindTotal,varsel])


M1 <- as.matrix(M[1:nindTeste,])
ini = nindTeste+1
M2 <- as.matrix(M[ini:nindTotal,])



# *******************************************************************
# 4.1 Estatísticas preliminares
# *******************************************************************
cat("  ", "\n")
cat("  ", "\n")

cat("======================================================","\n")
cat("Variável = ", colnames(fenotipo)[varsel], "\n")
cat("Validação = ", kkuu, "\n")


cat("======================================================","\n")
cat("Média =", mean(variavel), "\n")
cat("Desvio-padrão", sd(variavel), "\n")
cat("variância = ", var(variavel), "\n")
cat("Num. obs = ", length(variavel), "\n")
cat("DP média =" , sd(variavel)/sqrt(length(variavel)), "\n")


cat("======================================================","\n")
cat("  ", "\n")
cat("(*) Identificar o número de componentes que explica 80% da variação de X   ", "\n")
cat("Número de componentes estabelecidos pelo usuário  = ", ndim, "\n")
cat("  ", "\n")


pcr=pcr(variavel~M1)
names(pcr)
summary(pcr)

# *******************************************************************
# 4.2 Vetor de efeito de marcadores 
# *******************************************************************
# número de componentes - min(n_obs,n_par)-1    Estabelecer o número que explica a var. na variavel em 80%
nc_pcr=ndim

# coeficientes associados as marcas head(m_pcr)
m_pcr=as.matrix(pcr$coefficients[,,nc_pcr])


cat("Efeito das ondas = PCR","\n")
print(m_pcr)

par(mfrow=c(3,1))  # divide a tela gráfica em 3
seq <-(1:1:length(m_pcr))

plot(seq,m_pcr,xlab="Marcas",ylab="Efeito",main = colnames(fenotipo)[i], type = "l")
plot(seq,abs(m_pcr),xlab="Ondas",ylab="Abs(Efeito)",type = "l")
# visualização dos mais informativos 
plot(density(m_pcr))

# Ondas odenadas pela importância
cat("  ", "\n")
cat("Comprimento de ondas ordenados pela importância ","\n")
print(order(abs(m_pcr), decreasing=TRUE))

# *******************************************************************
# 4.3 Vetor de valores genômicos
# *******************************************************************

gbv_pcr=M1%*%m_pcr # valor genético genômico (GEBVs) dos indivíduos



cat("  ", "\n")
cat("Valores preditos ","\n")
print(gbv_pcr)


# R2 vfen e vgenômico
seq <-(1:1:length(variavel))
r2= cor(gbv_pcr,variavel)^2
cat("  ", "\n")
cat("R²(VPred,Média) =", r2, "\n")
valorR2t[i-1,kkuu] <- r2


  erro<-gbv_pcr-variavel
  erro2<-erro^2
  soma<-sum(erro2)
  nn<-length(erro)
  reqmt<-sqrt((1/nn)*soma)   
  valorerrot[i-1,kkuu] <- reqmt
cat("  ", "\n")
cat("Raiz(EQM) =", reqmt, "\n")


r2aj = 1- ((nn-1)/(nn-ncol(M1)-1))*(1-r2)
valorR2ajt[i-1,kkuu] <- r2aj
cat("  ", "\n")
cat("R2aj =", r2aj, "\n")


 
# identificação dos fenotipos no arquivo fenotipico
rownames(gbv_pcr)<-c(dados$id[1:nindTeste])

#ordenamento dos valores genômicos
gbv1=sort(as.matrix(gbv_pcr)[,],decreasing=TRUE) #ordenamento dos valores
plot(hist(gbv1))
gbv_20=quantile(gbv1, probs =c(80)/100) #identificação dos 20% melhores
top_20=gbv1[gbv1>=gbv_20]
cat("  ", "\n")
cat("Valores preditos ordenados - 20% ","\n")
print(top_20)



par(mfrow=c(3,1))  # divide a tela gráfica em 3
plot(seq,gbv_pcr,xlab="Ind - Teste",ylab="VPred", main = colnames(fenotipo)[i], type = "l")
plot(seq,variavel,xlab="Ind - Teste",ylab="Média fen",  type = "l")
plot(variavel,gbv_pcr,xlab="Média fen",ylab="VPred")




# *******************************************************************
# 4.4 Vetor de valores preditos - dados de validação
# *******************************************************************
gbv_pcr_val=M2%*%m_pcr
cat("  ", "\n")
cat("Valores preditos - Indivíduos de validação ","\n")
print(gbv_pcr_val)


# R2 vfen e vgenômico
seq <-(1:1:length(variavelval))
r2= cor(gbv_pcr_val,variavelval)^2
#nlin = nrow(variavelval)
#r2 = 1 - ((nlin -1)/(nlin - nvar -1))*(1-r2)
cat("  ", "\n")
cat("R²(EGBV,Média) =", r2, "\n")
valorR2v[i-1,kkuu] <- r2



  erro<-gbv_pcr_val-variavelval
  erro2<-erro^2
  soma<-sum(erro2)
  nn<-length(erro)
  reqmv<-sqrt((1/nn)*soma) 
valorerrov[i-1,kkuu] <- reqmv
cat("  ", "\n")
cat("Raiz(EQM) =", reqmv, "\n")


r2aj = 1- ((nn-1)/(nn-ncol(M1)-1))*(1-r2)
valorR2ajv[i-1,kkuu] <- r2aj
cat("  ", "\n")
cat("R2aj =", r2aj, "\n")

par(mfrow=c(3,1))  # divide a tela gráfica em 3
plot(seq,gbv_pcr_val,xlab="Ind - Validação",ylab="VPred", main = colnames(fenotipo)[i], type = "l")
plot(seq,variavelval,xlab="Ind - Validação",ylab="Média fen",  type = "l")
plot(variavelval,gbv_pcr_val,xlab="Média fen",ylab="VPred")

# *******************************************************************
# 4.5 Estimativas do Modelo 
# *******************************************************************

# Capacidade Preditiva
seq <-(1:1:length(variavelval))
CP= cor(gbv_pcr_val,variavelval)
#cat("  ", "\n")
cat("CP(VPred,Média) =", CP, "\n")
valorCP[i-1,kkuu] <- CP

# Viés
seq <-(1:1:length(variavel))
beta1= ((cov(gbv_pcr,variavel))/(var(variavel)))
#cat("  ", "\n")
cat("Beta1 - Cov(VPred,Média)/var(variavel) =", beta1, "\n")
valorbeta1[i-1,kkuu] <- beta1



}
}


res.par <- cbind(valorR2t,valorR2v, valorCP, valorbeta1, valorR2ajt, valorR2ajv)

nvar = ncol(fenotipo)-1

valorR2t1 <- matrix(rep(NA),nvar,1)
valorR2t2 <- matrix(rep(NA),nvar,1)
valorR2v1 <- matrix(rep(NA),nvar,1)
valorR2v2 <- matrix(rep(NA),nvar,1)
valorCP1 <- matrix(rep(NA),nvar,1)
valorCP2 <- matrix(rep(NA),nvar,1)
valorbeta11 <- matrix(rep(NA),nvar,1)
valorbeta12 <- matrix(rep(NA),nvar,1)


valorerrot1 <- matrix(rep(NA),nvar,1)
valorerrot2 <- matrix(rep(NA),nvar,1)
valorerrov1 <- matrix(rep(NA),nvar,1)
valorerrov2 <- matrix(rep(NA),nvar,1)



valorR2ajt1 <- matrix(rep(NA),nvar,1)
valorR2ajt2 <- matrix(rep(NA),nvar,1)
valorR2ajv1 <- matrix(rep(NA),nvar,1)
valorR2ajv2 <- matrix(rep(NA),nvar,1)



for(i1 in 1:nvar)
{
temp <- valorR2t[i1,]
valorR2t1[i1,1] <- mean(temp)
valorR2t2[i1,1] <- sd(temp)

temp <- valorR2v[i1,]
valorR2v1[i1,1] <- mean(temp)
valorR2v2[i1,1] <- sd(temp)


temp <- valorCP[i1,]
valorCP1[i1,1] <- mean(temp)
valorCP2[i1,1] <- sd(temp)


temp <- valorbeta1[i1,]
valorbeta11[i1,1] <- mean(temp)
valorbeta12[i1,1] <- sd(temp)


temp <-valorerrot[i1,]
valorerrot1[i1,1] <- mean(temp)
valorerrot2[i1,1] <- sd(temp)

temp <-valorerrov[i1,]
valorerrov1[i1,1] <- mean(temp)
valorerrov2[i1,1] <- sd(temp)




temp <-valorR2ajt[i1,]
valorR2ajt1[i1,1] <- mean(temp)
valorR2ajt2[i1,1] <- sd(temp)


temp <-valorR2ajv[i1,]
valorR2ajv1[i1,1] <- mean(temp)
valorR2ajv2[i1,1] <- sd(temp)

}


res.par1 <- cbind(valorR2t1,valorR2v1, valorCP1, valorbeta11,  valorerrot1,  valorerrov1, valorR2ajt1,valorR2ajv1 )
colnames(res.par1) <- c("R²Trein.","R²Val.", "Cap.Pred.", "Viés", "REQMt", "REQMv",  "R2ajt", "R2ajv")

res.par2 <- cbind(valorR2t2,valorR2v2, valorCP2, valorbeta12,  valorerrot2, valorerrov2,  valorR2ajt2,valorR2ajv2)
colnames(res.par2) <- c("R²Trein.","R²Val.", "Cap.Pred.", "Viés",  "REQMt", "REQMv", "R2ajt", "R2ajv")


cat("**********************************************************", "\n")
cat("    ","\n")
cat("---------------------------------------------------- ","\n")
cat("     Resumo  das Estimativas  ","\n")
cat("---------------------------------------------------- ","\n")
cat("  ","\n")
cat("Colunas: Estatística repetida k vezes   Linhas: variáveis", "\n")
cat("R²Trein.  R²Val.  Cap.Pred.  Viés  REQMt REQMv R2ajt R2ajv", "\n")
print(res.par)
cat("    ","\n")


cat("Médias das estatísticas  ","\n")
print(res.par1)
cat("    ","\n")

cat("Desvio padrão das estatísticas  ","\n")
print(res.par2)
cat("    ","\n")


cat("---------------------------------------------------- ","\n")


sink()


# =======================================
# Saida de arquivo resumo
# =======================================


sink("saida1$$$.dat")

cat("**********************************************************", "\n")
cat("**                                                      **  ", "\n")
cat("**            Programa Genes - Programa R               **", "\n")
cat("**                                                      **  ", "\n")
cat("**         NIR = Regressão via Componentes principais   **" ,"\n")
cat("**********************************************************", "\n")



AuxBD<-read.table("basedado$$$.dat",header=T)
cat("Data da Analise:", date(),"\n")
cat("Base de dados : ", colnames(AuxBD)[1] , "\n")

cat("**********************************************************", "\n")
cat("    ","\n")
cat("---------------------------------------------------- ","\n")
cat("     Resumo  das Estimativas  ","\n")
cat("---------------------------------------------------- ","\n")
cat("  ","\n")
cat("Colunas: Estatística repetida k vezes   Linhas: variáveis", "\n")
cat("R²Trein.  R²Val.  Cap.Pred.  Viés  REQMt REQMv R2ajt R2ajv", "\n")
print(res.par)
cat("    ","\n")


cat("Médias das estatísticas  ","\n")
print(res.par1)
cat("    ","\n")

cat("Desvio padrão das estatísticas  ","\n")
print(res.par2)
cat("    ","\n")


cat("---------------------------------------------------- ","\n")

sink()

dev.off()


