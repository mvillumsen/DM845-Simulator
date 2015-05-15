library(ggplot2)
setwd("~/Dropbox/SDU/8. semester/DM845/Project/DM845-Simulator/src")

readLength <- c('500', '1000', '1500')
errorProb <- c('0.05', '0.10', '0.15')

unphasable100k <- rbind(
    c(22915/190722, 24452/184634, 26343/178859),
    c(227/344830, 568/334865, 954/325212),
    c(1/469897, 9/457648, 24/445837))
row.names(unphasable100k) <- readLength
colnames(unphasable100k) <- errorProb
unphasable500k <- rbind(
    c(NA, 65726/640715, NA),
    c(207/883588, 578/874403, 993/863965),
    c(1/963613, NA, 12/953189))
row.names(unphasable500k) <- readLength
colnames(unphasable500k) <- errorProb

png("../output/plots/unphasable100k.png", width = 1024, height = 1024)
par(mar = c(5,5,5,5))
plot(colnames(unphasable100k), unphasable100k[1,], type='b', xlab = 'Errors in %',
     ylab = 'Unphasable positions in %', ylim=c(2.128126e-06, 0.15), col="green",
     pch=16, cex=1.5, cex.lab=2, cex.axis=1.5)
lines(colnames(unphasable100k), unphasable100k[2,], type='b', col="red", pch=16)
lines(colnames(unphasable100k), unphasable100k[3,], type='b', col="blue", pch=16)
title("Performance of Phasing with 100k Reads", cex.main=2)
legend("topleft", rownames(unphasable100k), pch=16, col=c('green', 'red', 'blue'),
       bty='o', cex=2, title = "Number of reads")
dev.off()

png("../output/plots/unphasable500k.png", width = 1024, height = 1024)
par(mar = c(5,5,5,5))
plot(colnames(unphasable500k), unphasable500k[1,], type='b', xlab = 'Errors in %',
     ylab = 'Unphasable positions in %', ylim=c(1.037761e-06,0.1025823), col="green",
     pch=16, cex=1.5, cex.lab=2, cex.axis=1.5)
lines(colnames(unphasable500k), unphasable500k[2,], type='b', col="red", pch=16)
lines(colnames(unphasable500k), unphasable500k[3,], type='b', col="blue", pch=16)
title("Performance of Phasing with 500k Reads", cex.main=2)
legend("topleft", rownames(unphasable500k), pch=16, col=c('green', 'red', 'blue'),
       bty='o', cex=2, title = "Number of reads")
dev.off()

file_list <- list.files("../output/results/")

for (file in file_list) {
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")) {
        dataset <- read.table(paste("../output/results/", file, sep = ""),
                              header = T, sep="\t")
    }
    
    # if the merged dataset does exist, append to it
    else if (exists("dataset")) {
        temp_dataset <- read.table(paste("../output/results/", file, sep = ""),
                                   header = T, sep="\t")
        dataset <- rbind(dataset, temp_dataset)
        rm(temp_dataset)
    }
}

matchPercent <- vector()
for (i in 1:nrow(dataset)) {
    matchPercent <- c(matchPercent, ( dataset[i,1][1] / sum(dataset[i,] ) ) * 100 )
}

## Created manually
reads.100k <- rbind(
    c(83.372161, 84.153290, 84.938491),
    c(65.953661, 67.118105, 68.226144),
    c(53.678436, 55.195709, 56.615366) )
reads.500k <- rbind(
    c(NA, 43.609658, NA),
    c(13.076512, 14.768164, 16.572427),
    c(5.016365, NA, 7.673753) )

row.names(reads.100k) <- readLength
colnames(reads.100k) <- errorProb
row.names(reads.500k) <- readLength
colnames(reads.500k) <- errorProb

## Create and save plot for 100k reads
png("../output/plots/plot100k.png", width = 1024, height = 1024)
par(mar = c(5,5,5,5))
plot(rownames(reads.100k), reads.100k[,1], type = "b", col = "red",
     xlab = "Read Length", ylab = "Match in %", ylim=c(50, 85), pch=16, cex.lab=2,
     cex.axis=1.5)
lines(rownames(reads.100k), reads.100k[,2], type = "b", col = "blue", pch=16)
lines(rownames(reads.100k), reads.100k[,3], type = "b", col = "green", pch=16)
title("Match with 100k Reads", cex.main=2, cex.lab=2)
legend("topright", colnames(reads.100k), pch=16, col=c('red', 'blue', 'green'),
       bty='o', cex=2, title = "Error Probability")
dev.off()

## Create and save plot for 500k reads
png("../output/plots/plot500k.png", width = 1024, height = 1024)
par(mar = c(5,5,5,5))
plot(rownames(reads.500k), reads.500k[,2], type="b", col="red",
     xlab = "Read Length", ylab = "Match in %", ylim=c(0,50), pch=16, cex.lab=2,
     cex.axis=1.5)
lines(rownames(reads.500k), reads.500k[,1], type="b", col="blue", pch=16)
lines(rownames(reads.500k), reads.500k[,3], type="b", col="green", pch=16)
title("Match with 500k Reads", cex.main=2)
legend("topright", colnames(reads.500k), pch=16, col=c('red', 'blue', 'green'),
       bty='o', cex=2, title = "Error Probability")
dev.off()


## Plot memory/time performance
## Since values are very close for all three error probabilites we just take the average
perf.100k <- rbind(
    c( (78.21+72.76+77.58)/3, ((373268+369580+363296)/3)/1024),
    c( (87.47+86.90+87.52)/3, ((502804+496472+483564)/3)/1024),
    c( (101.27+90.92+95.59)/3, ((613572+603444+597664)/3)/1024))

## NotSome values are not available for 500k
perf.500k <- rbind(
    c( 106.38, 816752/1024),
    c( (152.98+145.79+145.15)/3, ((1144728+1134464+1123820)/3)/1024),
    c( (176.03+172.91)/2, ((1290584+1273172)/2)/1024))

colnames(perf.100k) <- c('Time in sec.', 'Memory usage in MB')
colnames(perf.500k) <- c('Time in sec.', 'Memory usage in MB')
row.names(perf.100k) <- c('500', '1000', '1500')
row.names(perf.500k) <- c('500', '1000', '1500')

png("../output/plots/plotPerformance.png", width = 1024, height = 1024)
par(mar = c(5,5,5,5))
plot(perf.100k[,2], perf.100k[,1], type="p", col="blue", xlab = "Memory in MB",
     ylab = "Time in Seconds", xlim=c(300,1300), ylim=c(70, 180), pch=16,
     cex=1.5, cex.lab=2, cex.axis=1.5)
points(perf.500k[,2], perf.500k[,1], type="p", col="red", pch=16, cex=1.5)
text(perf.100k[,2], perf.100k[,1], row.names(perf.100k), pos=3, cex=1, col="blue")
text(perf.500k[,2], perf.500k[,1], row.names(perf.500k), pos=3, cex=1, col="red")
title("Performance (Memory vs Time in sec.)", cex.main=2)
legend("topleft", c('100.000', '500.000'), col=c('red', 'blue'), pch=16,
       bty='o', cex=2, title = "Number of Reads")
dev.off()
