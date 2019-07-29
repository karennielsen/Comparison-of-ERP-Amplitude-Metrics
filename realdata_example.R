# Real Data example for 
# "Comparison of Common Amplitude Metrics in Event-Related Potential Analysis"
# Karen E. Nielsen & Richard Gonzalez

# Data source: https://archive.ics.uci.edu/ml/datasets/eeg+database
# Acknowledgments to Henri Begleiter, 
# Neurodynamics Laboratory at the State University of New York Health Center at Brooklyn

library(dplyr)
library(tidyr)
library(ggplot2)


## Read in file(s) containing UCI data
# this can be done by reading in all files (e.g., with rbind) and keeping only control subjects
# or reading in control subjects only

## assume Step 2 of ERP analysis has been done already:
# data assumed to be band-pass filtered and artifact-rejected. 
# data is collected with Cz reference. 
# data is epoched
# no baseline data available, so assume data has been baselined

# note: time here is in samples at 256Hz - so time/256*1000 = time in milliseconds

## check that the example data has been filtered:
i=0
x.spec <- spectrum(fulldata$voltage[(1+i*256):((i+1)*256)],log="no",span=5,plot=FALSE)
spx <- x.spec$freq * 256
spy <- 2*x.spec$spec
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l", ylim=c(0,1200))
for(i in 1:20) {
x.spec <- spectrum(fulldata$voltage[(1+i*256):((i+1)*256)],log="no",span=5,plot=FALSE)
spx <- x.spec$freq * 256
spy <- 2*x.spec$spec
lines(spy~spx, col=i+1)
}
#Based on spectral density plot and Zhang et al. 1997 paper, seems to be 
# bandpass filtered between .02 and 50 Hz

#try out all three conditions:
salldata<-fulldata

## create subject-averaged ERP waveforms for each condition
salldata_grouped<-group_by(salldata, subject, channel, condition, time)
salldata_aerp<-summarize(salldata_grouped, mvoltage=mean(voltage), ntrials=n())


##get trial counts per condition
salldata_aerp_trials<-subset(salldata_aerp, channel=="CPZ" & time==0)
tapply(salldata_aerp_trials$ntrials, salldata_aerp_trials$condition, summary)

#use a prespecified window to isolate component of interest. 
#justification:
## the c247 from 220-260, surface energy 200-260... I'll use 200-300 for now. 
## luck says 50ms for avg, 150ms for peak, so split the difference at 100ms?
#We are expecting a positive deflection, but also compute min in case of negative deflection
sall_aerp_windowed<-subset(salldata_aerp, time>200/1000*256 & time <300/1000*256)
sall_aerp_grouped<-group_by(sall_aerp_windowed, subject, channel, condition)
sall_maxavg<-summarize(sall_aerp_grouped, max_v=max(mvoltage), avg_v=mean(mvoltage), min_v=min(mvoltage))

#pull the channel of interest (just max/avg summaries)
sall_maxavg_1ch<-subset(sall_maxavg, channel=="CPZ" & condition!="S2 match")
#reshape to setup for paired t-test
sall_maxavg_1ch_long<-gather(sall_maxavg_1ch, sumtype, value, c(max_v, min_v, avg_v))
sall_maxavg_1ch_long$cond_sum<-paste(sall_maxavg_1ch_long$condition, sall_maxavg_1ch_long$sumtype)
sall_maxavg_1ch_long<-sall_maxavg_1ch_long[,-which(colnames(sall_maxavg_1ch_long) %in% c("sumtype", "condition"))]
sall_maxavg_1ch_wide<-spread(sall_maxavg_1ch_long, cond_sum, value)

#paired t-tests
test_max<-t.test(sall_maxavg_1ch_wide$`S1 obj max_v`,
                 sall_maxavg_1ch_wide$`S2 nomatch max_v`,
                 paired=TRUE)
test_max

test_avg<-t.test(sall_maxavg_1ch_wide$`S1 obj avg_v`,
              sall_maxavg_1ch_wide$`S2 nomatch avg_v`,
              paired=TRUE)
test_avg



#pull the channel of interest  (all timepoints)
salldata_aerp_chs<-subset(salldata_aerp, channel %in% c("CPZ"))

salldata_aerp_chs$subjcond<-paste0(salldata_aerp_chs$subject,salldata_aerp_chs$condition)
#ggplot(subset(salldata_aerp_chs, condition!='S2 match'), aes(x=time/256*1000, y=mvoltage, color=condition)) +
#  geom_line(aes(group = subjcond)) +
#  ylim(-15,10) +
#  geom_vline(xintercept = 200) +
#  geom_vline(xintercept = 300) +
#  geom_hline(yintercept = 0) +
#  facet_grid(rows=vars(channel))

#get grand averages waveforms for this channel
salldata_aerp_chs<-group_by(salldata_aerp_chs, channel, condition, time)
salldata_grandaerp_chs<-summarize(salldata_aerp_chs, gm_voltage=mean(mvoltage))

#for publication, y axis and legend were edited to include units and simplify condition names:
ggplot(subset(salldata_grandaerp_chs, condition!='S2 match'), aes(x=time/256*1000, y=gm_voltage, group=condition)) +
  geom_line(aes(color=condition, lty=condition)) +
  ylim(-3,3) +
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 300) +
  geom_hline(yintercept = 0) + 
  theme_bw(base_size = 12)+
  xlab("Time (milliseconds)") +
  ylab("Voltage") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#find the max and average of the grand-averaged waveforms
salldata_grandaerp_chs_windowed<-subset(salldata_grandaerp_chs, condition!='S2 match' & time>200/1000*256 & time <300/1000*256)
#max s2 nomatch
max(salldata_grandaerp_chs_windowed$gm_voltage[which(salldata_grandaerp_chs_windowed$condition=="S2 nomatch")])
#max s1
max(salldata_grandaerp_chs_windowed$gm_voltage[which(salldata_grandaerp_chs_windowed$condition=="S1 obj")])
#avg s2 nomatch
mean(salldata_grandaerp_chs_windowed$gm_voltage[which(salldata_grandaerp_chs_windowed$condition=="S2 nomatch")])
#avg s1
mean(salldata_grandaerp_chs_windowed$gm_voltage[which(salldata_grandaerp_chs_windowed$condition=="S1 obj")])
