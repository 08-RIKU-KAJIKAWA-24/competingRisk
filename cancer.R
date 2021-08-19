library(survival)
library(cmprsk)

#データの読み込み
cancer_data <- read.csv("cancer_data.csv", header=T)

#エンドポイントと競合リスクの設定
cancer <- cancer_data
others <- cancer_data
cancer$status <- ifelse(cancer$status==1, 1, 0)
others$status <- ifelse(others$status==2, 1, 0)

#Kaplan-Meier推定
km_cancer <- survfit(Surv(time, status)~1, data=cancer)
km_others <- survfit(Surv(time, status)~1, data=others)

#Kaplan-Meier曲線の描画
dev.new(width=5, height=5)
plot(km_cancer ,xlab="day", ylab="probability", xlim=c(0, 2065), ylim=c(0,1), conf.int=F, mark.time=T)
par(new=T)
plot(km_others, col=2, conf.int=F, mark.time=T)
legend("bottomleft", paste(c("death.cancer", "death.others")), lty=1, col=1:2)

#累積発生割合(K-M法)
km_cancer <- cuminc(ftime=cancer$time, fstatus=cancer$status, cencode=0)
km_others<- cuminc(ftime=others$time, fstatus=others$status, cencode=0)
result_cancer <- timepoints(km_cancer, cancer$time)
result_others <- timepoints(km_others, others$time)
est_cancer <- result_cancer$"est"
est_others <- result_others$"est"

#Kaplan-Meier推定の表
km_table <- data.frame(time=c(0, cancer_data$time), at_risk=c(200, seq(200,1,-1)), death.cancer=c("-", cancer$status), death.others=c("-", others$status), eventFree_1=c(1, 1-est_cancer), eventFree_2=c(1, 1-est_others))

#累積発生割合(CIF)
cif <- cuminc(ftime=cancer_data$time, cancer_data$status)
result <- timepoints(cif, cancer$time)
est <- result$"est"
est_cancer <- est[1,]
est_others <- est[2, ]

#累積発生割合(CIF)の表
cif_table <- data.frame(time=c(0, cancer_data$time), at_risk=c(200, seq(200,1,-1)), death.cancer=c("-", cancer$status), death.others=c("-", others$status), CI_cancer=c(0, est_cancer), CI_others=c(0, est_others), total=c(0, est_cancer+est_others))

#========群間の検定==========
#Log-rank検定
logrank_cancer <- survdiff(Surv(time, status==1)~invasion, data=cancer_data)
logrank_others <- survdiff(Surv(time, status==2)~invasion, data=cancer_data)

#群間における累積発生割合の描画(Log-rank)
logrank_cancer_km <- survfit(Surv(time, status==1)~invasion, data=cancer_data)
logrank_others_km <- survfit(Surv(time, status==2)~invasion, dat=cancer_data)
dev.new(width=5, height=5)
plot(logrank_cancer_km, xlab="", ylab="", fun="event", lty=1:2)
par(new=T)
plot(logrank_others_km, xlab="day", ylab="probability", fun="event", lty=3:4)
legend("topleft", paste(c("10:death.cancer", "30:death.cancer", "10:death.others", "30:death.others")), lty=1:4, col=1)

#Gray検定
gray_test <- cuminc(ftime=cancer_data$time, fstatus=cancer_data$status, group=cancer_data$invasion)

#群間における累積発生割合の描画(Gray)
gray_graph <- cuminc(ftime=cancer_data$time, fstatus=cancer_data$status, group=cancer_data$invasion)
dev.new(width=5, height=5)
plot(gray_graph, xlab="day", ylab="probability")