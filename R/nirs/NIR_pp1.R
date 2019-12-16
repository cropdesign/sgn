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
par(mfrow=c(ng+1,1))  # divide a tela gráfica 
for(j in acessoi:acessof){
#par(new=TRUE)
plot(as.numeric(colnames(dados)),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)
}
par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados)),t(dados[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  "Todos acessos")


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



for(j in acessoi:acessof){
#par(new=TRUE)
plot(as.numeric(colnames(dados)),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)
# new data points
points(as.numeric(colnames(dados.bin)) ,dados.bin[j,],pch=2,col = 2)
#legend("topleft",legend=c("bin.size = "),pch = 2:1, col = 2:1)
}
par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.bin)),t(dados.bin[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - média",main =  "Todos acessos")

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


for(j in acessoi:acessof){
plot(as.numeric(colnames(dados)),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)
lines(as.numeric(colnames(dados.mm)), dados.mm[j, ], col = "red")
legend("topright", legend = c("Original", "média móvel"), lty = c(1, 1), col = 1:2)
}

par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.mm)),t(dados.mm[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - média móvel",main =  "Todos acessos")


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

for(j in acessoi:acessof){
plot(as.numeric(colnames(dados)),dados[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorvância",main =  j)
lines(as.numeric(colnames(dados.sg)), dados.sg[j, ], col = "red")
legend("topright", legend = c("Original", "Savitzky-Golay"), lty = c(1, 1), col = 1:2)
}


par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.sg)),t(dados.sg[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - Savitzky-Golay",main =  "Todos acessos")


write.table(dados.sg, file="dadossg.dat", row.names = F,col.names = F)
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

for(j in acessoi:acessof){
plot(as.numeric(colnames(dados.d1)), dados.d1[j, ], type = "l", xlab = "Comprimento de onda", ylab = "")
lines(as.numeric(colnames(dados.d2)), dados.d2[j, ], col = "red")
legend("topright", legend = c("1a der", "2a der"), lty = c(1, 1), col = 1:2)
}


par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.d1)),t(dados.d1[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - 1a. derivada",main =  "Todos acessos")



par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.d2)),t(dados.d2[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - 2a. derivada",main =  "Todos acessos")



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

dados.d1 <- t(diff(t(dados), differences = 1)) # first derivative
dados.d2 <- t(diff(t(dados), differences = 2)) # second derivative


dados.gd1 <- t(diff(t(dados), differences = 1, lag = lagx)) # first derivative
dados.gd2 <- t(diff(t(dados), differences = 2, lag = lagx)) # second derivative

dados.gsd1 <- gapDer(X = dados, m = 1, w = 2*lagx+1, s = lagx)
dados.gsd2 <- gapDer(X = dados, m = 2, w = 2*lagx+1, s = lagx)




for(j in acessoi:acessof){
plot(as.numeric(colnames(dados.d1)), dados.d1[j, ], type = "l", xlab = "Comprimento de onda", ylab = "")
lines(as.numeric(colnames(dados.gsd1)), dados.gsd1[j, ], col = "red")
legend("topright", legend = c("1a der", "1a der com gap"), lty = c(1, 1),col = 1:2)
}

write.table(dados.gsd1, file="dadosgsd1.dat", row.names = F,col.names = F)
write.table(dados.gsd2, file="dadosgsd2.dat", row.names = F,col.names = F)


par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.gsd1)),t(dados.gsd1[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - 1a. derivada com gap",main =  "Todos acessos")



par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.gsd2)),t(dados.gsd2[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - 2a. derivada com gap",main =  "Todos acessos")




cat("****************************************************************", "\n")
cat("Primeira derivada com gap","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosgd1.dat ","\n")
cat("Número de linhas ",nrow(dados.gsd1),"\n")
cat("Número de coluna ",ncol(dados.gsd1),"\n")

cat("Segunda derivada com gap","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosgd2.dat ","\n")
cat("Número de linhas ",nrow(dados.gsd2),"\n")
cat("Número de coluna ",ncol(dados.gsd2),"\n")
cat("****************************************************************", "\n")
}



if (Modsel == 6) {
# *******************************************************************
# 7. Normalização
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3
dados.snv <- standardNormalVariate(X = dados)


for(j in acessoi:acessof){
#par(new=TRUE)
plot(as.numeric(colnames(dados.snv)),dados.snv[j, ], type = "l", xlab = "Comprimento de onda",ylab = "Absorv. Normatizada",main =  j)
}



par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.snv)),t(dados.snv[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - Normalizado",main =  "Todos acessos")



write.table(dados.snv, file="dadossnv.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Normalização = ","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadossnv.dat ","\n")
cat("Número de linhas ",nrow(dados.snv),"\n")
cat("Número de coluna ",ncol(dados.snv),"\n")
cat("****************************************************************", "\n")

}





if (Modsel == 7) {
# *******************************************************************
# 8. detrend signal
# X = input spectral matrix wav = band centers
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3

dados.dt <- detrend(X = dados, wav = as.numeric(colnames(dados)))

for(j in acessoi:acessof){
#par(new=FALSE)
plot(as.numeric(colnames(dados)),dados[j,], type = "l", xlab = "Band number", ylab = "",main =  j)
par(new = T)
plot(as.numeric(colnames(dados)),dados.dt[j,], xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "red", type = "l")
axis(4, col = "red")
legend("topright", legend = c("raw", "detrend signal"), lty = c(1, 1), col = 1:2)

}

par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.dt)),t(dados.dt[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - Detrend signal",main =  "Todos acessos")



write.table(dados.dt, file="dadosdt.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Detrend signal ","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosdt.dat ","\n")
cat("Número de linhas ",nrow(dados.dt),"\n")
cat("Número de coluna ",ncol(dados.dt),"\n")
cat("****************************************************************", "\n")

}






if (Modsel == 8) {
# *******************************************************************
# 9. Centering and scaling - blockScale
# X = spectral matrix
# type = "soft" or "hard"
# The ouptut is a list with the scaled matrix (Xscaled) and the divisor (f)
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3

dados.bs <- blockScale(X = dados,type="hard")$Xscaled
sum(apply(dados.bs,2,var)) # this works!


for(j in acessoi:acessof){
plot(as.numeric(colnames(dados.bs)), dados.bs[j, ], type = "l", xlab = "Comprimento de onda", ylab = "Absorv. Normatizada",main =  j)
}


par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.bs)),t(dados.bs[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - Centralizado e normatizado 0 a 1",main =  "Todos acessos")


write.table(dados.bs, file="dadosbs.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Centering and scaling - blockScale ","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosbs.dat ","\n")
cat("Número de linhas ",nrow(dados.bs),"\n")
cat("Número de coluna ",ncol(dados.bs),"\n")
cat("****************************************************************", "\n")

}





if (Modsel == 9) {
# *******************************************************************
# 10. Centering and scaling - blockNorm
# X = spectral matrix
# type = "soft" or "hard"
# The ouptut is a list with the scaled matrix (Xscaled) and the divisor (f)
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3

dados.bn <- blockNorm(X = dados, targetnorm = 1)$Xscaled
sum(dados.bn^2)



for(j in acessoi:acessof){
plot(as.numeric(colnames(dados.bn)), dados.bn[j, ], type = "l", xlab = "Comprimento de onda", ylab = "Absorv. Normatizada",main =  j)
}

par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.bn)),t(dados.bn[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - Centralizado e normatizado",main =  "Todos acessos")


write.table(dados.bs, file="dadosbs.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Centering and scaling -  blockNorm ","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadosbn.dat ","\n")
cat("Número de linhas ",nrow(dados.bn),"\n")
cat("Número de coluna ",ncol(dados.bn),"\n")
cat("****************************************************************", "\n")

}




if (Modsel == 10) {
# *******************************************************************
# 11 Continuum removal
# type of data: 'R' for reflectance (default), 'A' for absorbance
# *******************************************************************
ng= acessof - acessoi +1
par(mfrow=c(ng,1))  # divide a tela gráfica em 3

dados.cr <- continuumRemoval(X = dados, type = "A")

# plot of the 10 first abs spectra
matplot(as.numeric(colnames(dados)), t(dados[acessoi:acessof, ]), type = "l",
ylim = c(0, 0.6), xlab = "Comprimento de onda/nm", ylab = "Absorvância")
matlines(as.numeric(colnames(dados)), t(dados.cr[acessoi:acessof, ]))

par(new=FALSE)
par(mfrow=c(1,1))
matplot(as.numeric(colnames(dados.cr)),t(dados.cr[, ]), type = "l", xlab = "Comprimento de onda",ylab = "Absorvância - Continuum removal",main =  "Todos acessos")



write.table(dados.cr, file="dadoscr.dat", row.names = F,col.names = F)
cat("****************************************************************", "\n")
cat("Continuum removal ","\n")
cat("Arquivo salvo em =  c: dados scripstR saidas dadoscr.dat ","\n")
cat("Número de linhas ",nrow(dados.cr),"\n")
cat("Número de coluna ",ncol(dados.cr),"\n")
cat("****************************************************************", "\n")

}




sink()
dev.off()