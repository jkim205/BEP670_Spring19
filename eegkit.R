library(eegkit)  # in Mac OS X, eegkit needs to install XQuartz for rpl package http://xquartz.org .

# eegcap function
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

# eegcap2d function
eegcap2d()
eegcap2d(col.point = mycols)

eegcap2d("10-20", cex.label = -1)
myelectrodes <- c("FP1","FP2","FPZ","F7","F3","FZ", "F4","F8","T7","C3","CZ","C4","T8", "P7","P3","PZ","P4","P8","O1","O2")
eegcap(myelectrodes)
