# eegkit sample R files, revised from the document written by Nathaniel E. Helwig <helwig@umn.edu>

library(eegkit)  # in Mac OS X, eegkit needs to install XQuartz for rpl package http://xquartz.org .

############## eegcap function ##############
eegcap()

data(eegcoord)
mycols <- rep("white",87)
enames <- rownames(eegcoord) 
mycols[enames=="A1"] <- "green" 
mycols[enames=="A2"] <- "light blue" 
mycols[enames=="NZ"] <- "pink" 

eegcap(col.point = mycols)

eegcap("10-20")
eegcap("10-20", plotlabels = FALSE)

myelectrodes <- c("FP1","FP2","FPZ","F7","F3","FZ", "F4","F8","T7","C3","CZ","C4","T8", "P7","P3","PZ","P4","P8","O1","O2")
eegcap(myelectrodes)

############## eegcap2d function ##############
eegcap2d()
eegcap2d(col.point = mycols)

eegcap2d("10-20", cex.label = -1)
myelectrodes <- c("FP1","FP2","FPZ","F7","F3","FZ", "F4","F8","T7","C3","CZ","C4","T8", "P7","P3","PZ","P4","P8","O1","O2")
eegcap2d(myelectrodes)

############## eegcapdense function ##############

eegcapdense()
eegcapdense("10-20", plotlabels = FALSE)
myelectrodes <- c("FP1","FP2","FPZ","F7","F3","FZ", "F4","F8","T7","C3","CZ","C4","T8", "P7","P3","PZ","P4","P8","O1","O2")
eegcapdense(myelectrodes)

############## eegcoord data ##############

plot3d(eegcoord$x,eegcoord$y,eegcoord$z,size=10,col="green")
text3d(eegcoord$x,eegcoord$y,eegcoord$z,texts=enames,col="blue") 
plot(eegcoord$xproj,eegcoord$yproj,cex=2,col="green",pch=19) 
text(eegcoord$xproj,eegcoord$yproj,labels=enames,col="blue")

############## eegdense ##############
data(eegdense)
plot3d(eegdense$x,eegdense$y,eegdense$z,size=2,col="green")
plot(eegdense$xproj,eegdense$yproj,cex=1,col="green",pch=19) 

########## eegfilter ##############
# parameters for signal 
Fs <- 1000                          # 1000 Hz signal
s <- 3                              # 3 seconds of data 
t <- seq(1, s * Fs) / Fs            # time sequence
n <- length(t)                      # number of data points 
freqs <- c(1, 5, 10, 20)            # frequencies
amp <- c(2, 1.5, 3, 1.75)           # strengths (amplitudes)
phs <- c(0, pi/6, pi/4, pi/2)       # phase shifts

# create data generating signals 
mu <- rep(0, n)
for(j in 1:length(freqs)){ 
  mu <- mu + amp[j] * sin(2*pi*t*freqs[j] + phs[j])
} 
set.seed(1)                         # set random seed
e <- rnorm(n)                       # Gaussian error
y <- mu + e                         # data = mean + error 

par(mfrow = c(2,1))
plot(t, mu, type = 'l',
     xlab = "Time (s)", ylab = expression("Strength (" * mu * "V)"), 
     main = "Time series of Noise-Free Data")
plot(t, y, type = 'l',
     xlab = "Time (s)", ylab = expression("Strength (" * mu * "V)"), 
     main = "Time series of Noise-Added Data")

# fft of noise-free data
ef <- eegfft(mu, Fs = Fs, upper = 40)
head(ef)
ef[ef$strength > 0.25, ]

par(mfrow = c(1,2)) 
plot(x = ef$frequency, y = ef$strength, t = "l", 
     xlab = "Frequency (Hz)", ylab = expression("Strength (" * mu * "V)"), 
     main = "FFT of Noise-Free Data")

cbind(amp, ef$strength[ef$strength > 0.25]) 
cbind(phs - pi/2, ef$phase[ef$strength > 0.25])

# fft of noise-added data
ef <- eegfft(y, Fs = Fs, upper = 40) 
head(ef)
ef[ef$strength > 0.25,]

plot(x = ef$frequency, y = ef$strength, t = "l", 
     xlab = "Frequency (Hz)", ylab = expression("Strength (" * mu * "V)"), 
     main = "FFT of Noise-Added Data")

cbind(amp, ef$strength[ef$strength > 0.25]) 
cbind(phs - pi/2, ef$phase[ef$strength > 0.25])

# Filtering

yf.but <- eegfilter(y, Fs = Fs, lower = 2, upper = 15, method = "butter", order = 4) # butterworth
yf.fir <- eegfilter(y, Fs = Fs, lower = 2, upper = 15, method = "fir1", order = 350) # 350-th order FIR filter

par(mfrow=c(2,1), mar = c(5, 4.5, 4, 2) + 0.1)
plot(t, yf.but, type = 'l',
     xlab = "Time (s)", ylab = expression("Strength (" * mu * "V)"), 
     main = "Time series data filtered by butterworth")
plot(t, yf.fir, type = 'l',
     xlab = "Time (s)", ylab = expression("Strength (" * mu * "V)"), 
     main = "Time series data filtered by FIR")

s5 <- 1.5 * sin(2*pi*t*5 + pi/6)
s10 <- 3 * sin(2*pi*t*10 + pi/4)
yftrue <- s5 + s10
mean((yf.but - yftrue)^2)
mean((yf.fir - yftrue)^2)

par(mfrow = c(1,1))
plot(t, yftrue, type = "l", lty = 1, lwd = 2, ylim = c(-3, 3)) 
lines(t, yf.but, col = "blue", lty = 2, lwd = 2) 
lines(t, yf.fir, col = "red", lty = 3, lwd = 2)
legend("topright", legend = c("Truth", "Butterworth", "FIR"), lty = 1:3, lwd = 2, col = c("black", "blue", "red"), bty = "n")

# power spectrum in db
par(mfrow=c(1,3), mar = c(5, 4.5, 4, 2) + 0.1) 
eegpsd(y, Fs = Fs, upper = 50, t = "b", main = "Before Filtering", lwd = 2)
rect(2, -63, 15, 1, col = rgb(0.5,0.5,0.5,1/4)) 
legend("topright", legend = "2-15 Hz Filter", fill = rgb(0.5,0.5,0.5,1/4), bty = "n")
eegpsd(yf.but, Fs = Fs, upper = 50, t = "b", main = "After Butterworth Filter", lwd = 2)
eegpsd(yf.fir, Fs = Fs, upper = 50, t = "b", main = "After FIR Filter", lwd = 2)

# power spectrum in mv^2
par(mfrow=c(1,3), mar = c(5, 4.5, 4, 2) + 0.1) 
eegpsd(y, Fs = Fs, upper = 50, unit = "mV^2", t = "b", main = "Before Filtering", lwd = 2)
rect(2, 0, 15, 1.05, col = rgb(0.5,0.5,0.5,1/4)) 
legend("topright", legend = "2-15 Hz Filter", fill = rgb(0.5,0.5,0.5,1/4), bty = "n")
eegpsd(yf.but, Fs = Fs, upper = 50, unit = "mV^2", t = "b", main = "After Butterworth Filter", lwd = 2)
eegpsd(yf.fir, Fs = Fs, upper = 50, unit = "mV^2", t = "b", main = "After FIR Filter", lwd = 2)

############## eeghead ##############
data(eeghead) 
shade3d(eeghead) 
eeghead$material$color <- rep("black",length(eeghead$material$color)) 
wire3d(eeghead)

############## eegica - ICA - independent component analysis ##############
# get "c" subjects of "eegdata" data
data(eegdata) 
idx <- which(eegdata$group=="c") 
eegdata <- eegdata[idx,]

# get average data (across subjects) 
eegmean <- tapply(eegdata$voltage,list(eegdata$channel,eegdata$time),mean)

# remove ears and nose 
acnames <- rownames(eegmean) 
idx <- c(which(acnames=="X"),which(acnames=="Y"),which(acnames=="nd")) 
eegmean <- eegmean[-idx,]

# get spatial coordinates (for plotting) 
data(eegcoord) 
cidx <- match(rownames(eegmean),rownames(eegcoord))

# temporal ICA with 4 components
icatime <- eegica(eegmean,4) 
icatime$vafs
quartz() 
par(mfrow=c(4,2)) 
tseq <- (0:255)*1000/255 
for(j in 1:4){ 
  par(mar=c(5.1,4.6,4.1,2.1)) 
  sptitle <- bquote("VAF: "*.(round(icatime$vafs[j],4))) 
  eegtime(tseq,icatime$S[,j],main=bquote("Component "*.(j)),cex.main=1.5) 
  eegspace(eegcoord[cidx,4:5],icatime$M[,j],main=sptitle) 
  }

# spatial ICA with 4 components 
icaspace <- eegica(eegmean,4,type="space") 
icaspace$vafs 
quartz() 
par(mfrow=c(4,2)) 
tseq <- (0:255)*1000/255 
for(j in 1:4){ 
  par(mar=c(5.1,4.6,4.1,2.1)) 
  sptitle <- bquote("VAF: "*.(round(icaspace$vafs[j],4))) 
  eegtime(tseq,icaspace$M[,j],main=bquote("Component "*.(j)),cex.main=1.5) 
  eegspace(eegcoord[cidx,4:5],icaspace$S[,j],main=sptitle) 
}

############## eegmesh ##############
data(eegmesh) 
wire3d(eegmesh) 
eegmesh$material$color <- rep("red",length(eegmesh$material$color)) 
shade3d(eegmesh)

############## eegpsd ##############

# create data generating signals 
n <- 1000                           # 1000 Hz signal 
s <- 2                              # 2 seconds of data
t <- seq(0, s, length.out = s * n)  # time vector 
s1 <- sin(2*pi*t)                   # 1 Hz sinusoid
s5 <- sin(2*pi*t*5)                 # 5 Hz sinusoid
s10 <- sin(2*pi*t*10)               # 10 Hz sinusoid
s20 <- sin(2*pi*t*20)               # 20 Hz sinusoid
# create data 
set.seed(1)                         # set random seed
e <- rnorm(s * n, sd = 0.25)        # Gaussian error
mu <- s1 + s5 + s10 + s20           # 1 + 5 + 10 + 20 Hz mean 
y <- mu + e                         # data = mean + error
# plot psd (single channel) 
eegpsd(y, Fs = n, upper = 30, t = 'l')
# plot psd (multi-channel) 
ym <- cbind(s1, s5, s10, s20) 
eegpsd(ym, Fs = n, upper = 30, units = "mV")

############## eegresample ##############
# create vector with N = 200 time points 
N <- 200 
x <- sin(4 * pi * seq(0, 1, length.out = N))
# down-sample (i.e., decrease sampling rate) to n = 100 
y <- eegresample(x, n = 100) 
mean((y - sin(4 * pi * seq(0, 1, length.out = 100)))^2)
# up-sample (i.e., increase sampling rate) to n = 500 
z <- eegresample(x, n = 500) 
mean((z - sin(4 * pi * seq(0, 1, length.out = 500)))^2)
# plot results 
par(mfrow = c(1,3)) 
plot(x, main = "Original (N = 200)") 
plot(y, main = "Down-sampled (n = 100)") 
plot(z, main = "Up-sampled (n = 500)")

# create matrix with N = 500 time points and 2 columns 
N <- 500 
x <- cbind(sin(2 * pi * seq(0, 1, length.out = N)), sin(4 * pi * seq(0, 1, length.out = N)))
# down-sample (i.e., decrease sampling rate) to n = 250 
y <- eegresample(x, n = 250) 
ytrue <- cbind(sin(2 * pi * seq(0, 1, length.out = 250)), sin(4 * pi * seq(0, 1, length.out = 250)))
mean((y - ytrue)^2)
# up-sample (i.e., increase sampling rate) to n = 1000 
z <- eegresample(x, n = 1000) 
ztrue <- cbind(sin(2 * pi * seq(0, 1, length.out = 1000)), sin(4 * pi * seq(0, 1, length.out = 1000)))
mean((z - ztrue)^2)
# plot results 
par(mfrow = c(1,3)) 
plot(x[,1], main = "Original (N = 500)", cex = 0.5) 
points(x[,2], pch = 2, col = "blue", cex = 0.5) 
plot(y[,1], main = "Down-sampled (n = 250)", cex = 0.5) 
points(y[,2], pch = 2, col = "blue", cex = 0.5) 
plot(z[,1], main = "Up-sampled (n = 1000)", cex = 0.5) 
points(z[,2], pch = 2, col = "blue", cex = 0.5)

################ eegsim ##############
data(eegcoord) 
chnames <- rownames(eegcoord) 
tseq <- seq(0,1,length.out=200)
quartz(width=18,height=6) 
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,11), 2, 6, byrow = TRUE))
eegspace(eegcoord[,4:5],p1s(chnames),cex.point=1,main=expression(psi[p1]),cex.main=2,vlim=c(-3,9)) 
eegtime(tseq,p1t(tseq),ylim=c(-1,1),asp=1/2,main=expression(tau[p1]),cex.main=2,
        xlab="Time After Stimulus (sec)")
eegspace(eegcoord[,4:5],p2s(chnames),cex.point=1,main=expression(psi[p2]),cex.main=2,vlim=c(-3,9)) 
eegtime(tseq,p2t(tseq),ylim=c(-1,1),asp=1/2,main=expression(tau[p2]),cex.main=2,
        xlab="Time After Stimulus (sec)")
eegspace(eegcoord[,4:5],p3s(chnames),cex.point=1,main=expression(psi[p3]),cex.main=2,vlim=c(-3,9)) 
eegtime(tseq,p3t(tseq),ylim=c(-1,1),asp=1/2,main=expression(tau[p3]),cex.main=2,
        xlab="Time After Stimulus (sec)")
eegspace(eegcoord[,4:5],n1s(chnames),cex.point=1,main=expression(psi[n1]),cex.main=2,vlim=c(-3,9)) 
eegtime(tseq,n1t(tseq),ylim=c(-1,1),asp=1/2,main=expression(tau[n1]),cex.main=2,
        xlab="Time After Stimulus (sec)")
eegspace(eegcoord[,4:5],n2s(chnames),cex.point=1,main=expression(psi[n2]),cex.main=2,vlim=c(-3,9)) 
eegtime(tseq,n2t(tseq),ylim=c(-1,1),asp=1/2,main=expression(tau[n2]),cex.main=2,
        xlab="Time After Stimulus (sec)")
plot(seq(-10,10),seq(-10,10),type="n",axes=FALSE,xlab="",ylab="") 
text(0,8,labels=expression(omega[p1]*" = "*psi[p1]*tau[p1]),cex=2) 
text(0,4,labels=expression(omega[n1]*" = "*psi[n1]*tau[n1]),cex=2) 
text(0,0,labels=expression(omega[p2]*" = "*psi[p2]*tau[p2]),cex=2) 
text(0,-4,labels=expression(omega[n2]*" = "*psi[n2]*tau[n2]),cex=2) 
text(0,-8,labels=expression(omega[p3]*" = "*psi[p3]*tau[p3]),cex=2)

quartz(width=15,height=3) 
tseq <- c(50,150,250,350,450)/1000 
par(mfrow=c(1,5)) 
for(j in 1:5){ 
  eegspace(eegcoord[,4:5],eegsim(chnames,rep(tseq[j],87)),vlim=c(-6.8,5.5),
  main=paste(tseq[j]*1000," ms"),cex.main=2) 
  }

################### eegsmooth ##############
########## EXAMPLE 1: Temporal ##########
# get "PZ" electrode of "c" subjects in "eegdata" data 
data(eegdata) 
idx <- which(eegdata$channel=="PZ" & eegdata$group=="c") 
eegdata <- eegdata[idx,]
# temporal smoothing 
eegmod <- eegsmooth(eegdata$voltage,time=eegdata$time)
# define data for prediction 
time <- seq(min(eegdata$time),max(eegdata$time),length.out=100) 
yhat <- predict(eegmod,newdata=time,se.fit=TRUE)
# plot results using eegtime 
eegtime(time*1000/255,yhat$fit,voltageSE=yhat$se.fit,ylim=c(-4,4),main="Pz")

########## EXAMPLE 2: Spatial ##########
data(eegdata)
idx <- which(eegdata$time==65L & eegdata$group=="c") 
eegdata <- eegdata[idx,]
# remove ears, nose, and reference (Cz) 
idx <- c(which(eegdata$channel=="X"),which(eegdata$channel=="Y"), which(eegdata$channel=="nd"),which(eegdata$channel=="Cz"))
eegdata <- eegdata[-idx,]
# match to eeg coordinates 
data(eegcoord) 
cidx <- match(eegdata$channel,rownames(eegcoord))
# spatial smoothing 
eegmod <- eegsmooth(eegdata$voltage,space=eegcoord[cidx,1:3])
# use dense cap for prediction 
mycap <- levels(factor(eegdata$channel)) 
ix <- eegcapdense(mycap,type="2d",index=TRUE) 
data(eegdense) 
space <- eegdense[ix,1:3] 
yhat <- predict(eegmod,newdata=space)
eegspace(eegdense[ix,4:5],yhat)
