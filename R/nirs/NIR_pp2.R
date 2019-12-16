# *******************************************************************
# Procedimento: NIR Pré-processamento
# Entrada: C:\\dados\\ScriptsR\\saidas\trabalho$$$.dat
#        : C:\dados\ScriptsR\Saidas\auxiliar$$$.dat
# Saída:   C:\\dados\\ScriptsR\\saidas\saida$$$.dat
# Data: 18-05-2015
# *******************************************************************



# *******************************************************************
# 1. Definição da trilha de dados
# *******************************************************************
setwd("C:\\dados\\ScriptsR\\saidas")
library(prospectr)
library(stringr)


# *******************************************************************
# 2. Leitura dos dados
# *******************************************************************
NomeAqdados<-read.table("trabalho$$$.dat",header=F)
aqT =NomeAqdados[1,1]
dados<-read.table(levels(aqT)[1],header=T, check.names = TRUE)
colnames(dados) <- str_sub(colnames(dados), start = 2)


Aux<-read.table("auxiliar$$$.dat",header=F)

acessoi =Aux[1,1]
acessof =Aux[2,1]
pacote =Aux[3,1]
valorm =Aux[4,1]
graupol =Aux[5,1]
lagx =Aux[6,1]



# *******************************************************************
# 3. Cabeçalho
# *******************************************************************
sink("saida$$$.dat")
pdf(file='grafico.pdf')
cat("**********************************************************", "\n")
cat("**                                                      **  ", "\n")
cat("**            Programa Genes - Programa R               **", "\n")
cat("**                                                      **  ", "\n")
cat("**      Procedimento: NIR Pré-processamento             **" ,"\n")
cat("**                                                      **  ", "\n")
cat("**********************************************************", "\n")
AuxBD<-read.table("basedado$$$.dat",header=T)
cat("Data da Analise:", date(),"\n")
cat("Base de dados : ", colnames(AuxBD)[1] , "\n")
Sys.time()
cat("**********************************************************", "\n")
cat("  ", "\n")
#cat("Informações Preliminares  ", "\n")
#summary(dados)
cat("  ", "\n")



# *******************************************************************
# 1. Plotagem de todos os spectruns originais
# *******************************************************************
par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados)),t(dados[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância original",main =  "Todos acessos")



dados.out <- dados

# *******************************************************************
# 2. Inicio das análises
# *******************************************************************


for(j in 7:18){
Modsel= Aux[j,1]
if (Modsel != 0) {
cat(" ","\n")
cat("****************************************************************", "\n")
cat("Pré-processamento realizado :", Modsel,"\n")
dados.in <- dados.out
cat("Número de linhas ",nrow(dados.in),"\n")
cat("Número de coluna ",ncol(dados.in),"\n")







if (Modsel == 1) {
# *******************************************************************
# 2. Binning
# After averaging, the spectrum can be further resampled (binning)
# We keep here one 1 out every 10 data points
# *******************************************************************
dados.out <- binning(dados.in,bin.size=pacote)
cat("Binning - Média ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")

}




if (Modsel == 2) {
# *******************************************************************
# 3. Média móvel
# window size of 2m+1 bands
# Note that the m(valorm) first and last bands are lost in the process
# *******************************************************************
dados.out <- movav(dados.in, w =  2*valorm+1)
cat("Média móvel ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")
}

if (Modsel == 3) {
# *******************************************************************
# 4. Savitzky-Golay filtering
# p = polynomial order w = window size (must be odd) m = m-th derivative (0
# = smoothing) The function accepts vectors, data.frames or matrices. For a
# matrix input, observations should be arranged row-wise
# *******************************************************************
dados.out <- savitzkyGolay(dados.in, p = graupol, w =  2*valorm+1, m = 0)
cat("Savitzky-Golay filtering  ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")
}




if (Modsel == 4) {
# *******************************************************************
# 5. 1a Derivada
# X = wavelength Y = spectral matrix n = order
# *******************************************************************
dados.out <- t(diff(t(dados.in), differences = 1)) # first derivative
cat("1a. derivada  ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")
}


if (Modsel == 5) {
# *******************************************************************
# 6. 2a. Derivada
# X = wavelength Y = spectral matrix n = order
# *******************************************************************
dados.out <- t(diff(t(dados.in), differences = 2)) # second derivative
cat("2a. derivada  ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")

}




if (Modsel == 6) {
# *******************************************************************
# 7. 1a. derivative with a gap 
# m = order of the derivative w = window size ( = {2 * gap size} + 1) s =
# segment size first derivative with a gap of 10 bands
# *******************************************************************
dados.out <- gapDer(X = dados.in, m = 1, w = 2*lagx+1, s = lagx)
cat("1a. derivada com gap  ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")

}

if (Modsel == 7) {
# *******************************************************************
# 8. 2a. derivative with a gap of 10 bands
# m = order of the derivative w = window size ( = {2 * gap size} + 1) s =
# segment size first derivative with a gap of 10 bands
# *******************************************************************
dados.out <- gapDer(X = dados.in, m = 2, w = 2*lagx+1, s = lagx)
cat("2a. derivada com gap  ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")


}


if (Modsel == 8) {
# *******************************************************************
# 9. Normalização
# *******************************************************************
dados.out <- standardNormalVariate(X = dados.in)
cat("Normalização ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")


}





if (Modsel == 9) {
# *******************************************************************
# 10. detrend signal
# X = input spectral matrix wav = band centers
# *******************************************************************
dados.out <- detrend(X = dados.in, wav = as.numeric(colnames(dados)))
cat("Detrend signal ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")


}


if (Modsel == 10) {
# *******************************************************************
# 11. Centering and scaling - blockScale
# X = spectral matrix
# type = "soft" or "hard"
# The ouptut is a list with the scaled matrix (Xscaled) and the divisor (f)
# *******************************************************************
dados.out <- blockScale(X = dados.in,type="hard")$Xscaled
cat("Centering and scaling - blockScale ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")


}



if (Modsel == 11) {
# *******************************************************************
# 12. Centering and scaling - blockNorm
# X = spectral matrix
# type = "soft" or "hard"
# The ouptut is a list with the scaled matrix (Xscaled) and the divisor (f)
# *******************************************************************
dados.out <- blockNorm(X = dados.in, targetnorm = 1)$Xscaled
cat("Centering and scaling -blockNorm ","\n")
cat("****************************************************************", "\n")
cat(" ","\n")
}




if (Modsel == 12) {
# *******************************************************************
# 13 Continuum removal
# type of data: 'R' for reflectance (default), 'A' for absorbance
# *******************************************************************
dados.out <- continuumRemoval(X = dados.in, type = "A")
cat("Continuum removal","\n")
cat("****************************************************************", "\n")
cat(" ","\n")
}

}
}

write.table(dados.out, file="dados_out.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Resultado final  ","\n")
cat("****************************************************************", "\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dados_out.dat ","\n")
cat("Número de linhas ",nrow(dados.out),"\n")
cat("Número de coluna ",ncol(dados.out),"\n")
cat("****************************************************************", "\n")
par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.out)),t(dados.out[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância transformadal",main =  "Todos acessos")




sink()
dev.off()