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



# *******************************************************************
# 2. Leitura dos dados
# *******************************************************************
NomeAqdados<-read.table("trabalho$$$.dat",header=F)
aqT =NomeAqdados[1,1]
dadosOrig<-read.table(levels(aqT)[1],header=F)

eixox <- dadosOrig[1:2,]
dados=dadosOrig[2:nrow(dadosOrig),]


Aux<-read.table("auxiliar$$$.dat",header=F)

acessoi =Aux[1,1]
acessof =Aux[2,1]
Modsel =Aux[3,1]
pacote =Aux[4,1]
valorm =Aux[5,1]
graupol =Aux[6,1]
lagx =Aux[7,1]

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



if (Modsel == 0) {
# *******************************************************************
# 1. Plotagem dos spectruns 
# *******************************************************************
cat("****************************************************************", "\n")
cat("Spetrum ","\n")
cat("****************************************************************", "\n")
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3
for(j in acessoi:acessof){
#par(new=TRUE)
plot(as.numeric(eixox[1,]),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)
}
}

if (Modsel == 1) {
# *******************************************************************
# 2. Binning
# After averaging, the spectrum can be further resampled (binning)
# We keep here one 1 out every 10 data points
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3
dados.bin <- binning(dados,bin.size=pacote)
eixox.bin <- binning(eixox,bin.size=pacote)



for(j in acessoi:acessof){
#par(new=TRUE)
plot(as.numeric(eixox[1,]),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)

# new data points
points(as.numeric(eixox.bin[1,]),dados.bin[j,],pch=2,col = 2)
#legend("topleft",legend=c("bin.size = "),pch = 2:1, col = 2:1)

}
write.table(dados.bin, file="dadosmed.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Binning - redução por janela = ",pacote,"\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosmed.dat ","\n")
cat("Número de linhas ",nrow(dados.bin),"\n")
cat("Número de coluna ",ncol(dados.bin),"\n")

cat("****************************************************************", "\n")

}




if (Modsel == 2) {
# *******************************************************************
# 3. Média móvel
# window size of 2m+1 bands
# Note that the m(valorm) first and last bands are lost in the process
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3

dados.mm <- movav(dados, w =  2*valorm+1)
i1=valorm+1
i2=ncol(dados)-valorm
eixow <- eixox[,i1:i2]

for(j in acessoi:acessof){

plot(as.numeric(eixox[1,]),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)
lines(as.numeric(eixow[1,]), dados.mm[j, ], col = "red")
legend("topright", legend = c("Original", "média móvel"), lty = c(1, 1), col = 1:2)

}


write.table(dados.mm, file="dadosmedmov.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Média móvel:  Valor m = ",valorm,"\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosmedmov.dat ","\n")
cat("Número de linhas ",nrow(dados.mm),"\n")
cat("Número de coluna ",ncol(dados.mm),"\n")

cat("****************************************************************", "\n")

}

if (Modsel == 3) {
# *******************************************************************
# 4. Savitzky-Golay filtering
# p = polynomial order w = window size (must be odd) m = m-th derivative (0
# = smoothing) The function accepts vectors, data.frames or matrices. For a
# matrix input, observations should be arranged row-wise
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3


dados.sg <- savitzkyGolay(dados, p = graupol, w =  2*valorm+1, m = 0)


i1=valorm+1
i2=ncol(dados)-valorm
eixow <- eixox[,i1:i2]

for(j in acessoi:acessof){

plot(as.numeric(eixox[1,]),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)
lines(as.numeric(eixow[1,]), dados.sg[j, ], col = "red")
legend("topright", legend = c("Original", "Savitzky-Golay"), lty = c(1, 1), col = 1:2)

}


write.table(dados.mm, file="dadossg.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Média móvel:  Valor m = ",valorm,"\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadossg.dat ","\n")
cat("Número de linhas ",nrow(dados.sg),"\n")
cat("Número de coluna ",ncol(dados.sg),"\n")

cat("****************************************************************", "\n")

}




if (Modsel == 4) {
# *******************************************************************
# 5. Derivadas
# X = wavelength Y = spectral matrix n = order
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3

dados.d1 <- t(diff(t(dados), differences = 1)) # first derivative
dados.d2 <- t(diff(t(dados), differences = 2)) # second derivative


i1=1
i2=ncol(dados)-1
i3=ncol(dados)-2
eixow1 <- eixox[,i1:i2]
eixow2 <- eixox[,i1:i3]

for(j in acessoi:acessof){
plot(as.numeric(eixow1[1,]), dados.d1[j, ], type = "l", xlab = "Comprimento de onda", ylab = "")
lines(as.numeric(eixow2[1,]), dados.d2[j, ], col = "red")
legend("topright", legend = c("1a der", "2a der"), lty = c(1, 1), col = 1:2)
}

write.table(dados.d1, file="dadosd1.dat", row.names = F,col.names = F)
write.table(dados.d2, file="dadosd2.dat", row.names = F,col.names = F)

cat("****************************************************************", "\n")
cat("Primeira derivada","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosd1.dat ","\n")
cat("Número de linhas ",nrow(dados.d1),"\n")
cat("Número de coluna ",ncol(dados.d1),"\n")

cat("Segunda derivada","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosd2.dat ","\n")
cat("Número de linhas ",nrow(dados.d2),"\n")
cat("Número de coluna ",ncol(dados.d2),"\n")
cat("****************************************************************", "\n")
}




if (Modsel == 5) {
# *******************************************************************
# 6. derivative with a gap of 10 bands
# m = order of the derivative w = window size ( = {2 * gap size} + 1) s =
# segment size first derivative with a gap of 10 bands
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3

dados.gd1 <- t(diff(t(dados), differences = 1, lag = lagx)) # first derivative
dados.gd2 <- t(diff(t(dados), differences = 2, lag = lagx)) # second derivative



# m = order of the derivative w = window size ( = {2 * gap size} + 1) s =
# segment size first derivative with a gap of 10 bands
gsd1 <- gapDer(X = NIRsoil$spc, m = 1, w = 11, s = 10)
gsd1 <- gapDer(X = NIRsoil$spc, m = 1, w = 11, s = 10)




i1=1
i2=ncol(dados)-1
i3=ncol(dados)-2
eixow1 <- eixox[,i1:i2]
eixow2 <- eixox[,i1:i3]

for(j in acessoi:acessof){
plot(as.numeric(eixow1[1,]), dados.d1[j, ], type = "l", xlab = "Comprimento de onda", ylab = "")
lines(as.numeric(eixow2[1,]), dados.d2[j, ], col = "red")
legend("topright", legend = c("1a der", "2a der"), lty = c(1, 1), col = 1:2)
}

write.table(dados.d1, file="dadosd1.dat", row.names = F,col.names = F)
write.table(dados.d2, file="dadosd2.dat", row.names = F,col.names = F)

cat("****************************************************************", "\n")
cat("Primeira derivada","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosd1.dat ","\n")
cat("Número de linhas ",nrow(dados.d1),"\n")
cat("Número de coluna ",ncol(dados.d1),"\n")

cat("Segunda derivada","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosd2.dat ","\n")
cat("Número de linhas ",nrow(dados.d2),"\n")
cat("Número de coluna ",ncol(dados.d2),"\n")
cat("****************************************************************", "\n")
}





sink()
dev.off()