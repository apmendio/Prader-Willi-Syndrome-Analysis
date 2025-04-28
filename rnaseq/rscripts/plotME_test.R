alldata_f <- alldata
dat <- rbind(alldata_f, alldata_m)
dat$grp <- rep(factor(1:2), times = c(nrow(alldata_f), nrow(alldata_m)))

dat2 <- rbind(alldata_f, alldata_m)
x <- c("female", "male")
x
dat2$group <- rep(x, times = c(nrow(alldata_f), nrow(alldata_m)))
q <- ggplot(dat, aes(x = Timepoint, y = MEorange, color = grp)) +
  geom_point() +
  geom_smooth(method = "loess") +
  scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21))
q
##plot data from different dataframes: geom_smooth does not work for this method
p <- ggplot() + geom_point(data=alldata, aes(x = Timepoint, y = MEmidnightblue), color = "red") +
  geom_point(data=alldata_m, aes(x = Timepoint, y = MEmidnightblue), color = "blue") +
  scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21))
p + geom_smooth(method = "loess")
print(p)
p

q <- ggplot(alldata, aes(x = Timepoint, y = MEdarkred)) +
  geom_point() +
  geom_smooth(method = "loess") +
  scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21))
q

dat <- rbind(alldata_cm, alldata_cf)

dat$Genotype.Entrainment.Sex <- rep(Experiment, times = c(nrow(alldata_cm), nrow(alldata_cf)))

dat2 <- rbind(alldata_f, alldata_m)
Experiment <- c("CLAMS_Male", "CLAMS_Female")
x
dat2$group <- rep(x, times = c(nrow(alldata_f), nrow(alldata_m)))