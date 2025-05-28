### Numerical Studies RMSE ###
model <- c(rep("DNAmf",100),rep("RNAmf",100),rep("CoKriging",100),rep("BM",100),rep("FBM",100))
model <- factor(model, levels=c("DNAmf", "RNAmf", "CoKriging", "BM", "FBM"))

df.tuolinear.rmse <- c(result.tuolinear.rmse[,1], result.tuolinear.rmse[,2], result.tuolinear.rmse[,3], result.tuolinear.rmse[,5], result.tuolinear.rmse[,4])

padditive <- ggplot(data.frame(df.tuolinear.rmse, model), aes(x=model, y=df.tuolinear.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Additive") + labs(x="", y = "RMSE") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 5),
        panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

df.tuononlinear.rmse <- c(result.tuononlinear.rmse[,1], result.tuononlinear.rmse[,2], result.tuononlinear.rmse[,3], result.tuononlinear.rmse[,5], result.tuononlinear.rmse[,4])

pnonadditive <- ggplot(data.frame(df.tuononlinear.rmse, model), aes(x=model, y=df.tuononlinear.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Non-additive") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

df.currin.rmse <- c(result.currin.rmse[,1], result.currin.rmse[,2], result.currin.rmse[,3], result.currin.rmse[,5], result.currin.rmse[,4])

pcurrin <- ggplot(data.frame(df.currin.rmse, model), aes(x=model, y=df.currin.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Currin") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

### Numerical Studies CRPS ###

df.tuolinear.crps <- c(result.tuolinear.crps[,1], result.tuolinear.crps[,2], result.tuolinear.crps[,3], result.tuolinear.crps[,5], result.tuolinear.crps[,4])

padditive2 <- ggplot(data.frame(df.tuolinear.crps, model), aes(x=model, y=df.tuolinear.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + labs(x="", y = "CRPS") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 5),
        panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

df.tuononlinear.crps <- c(result.tuononlinear.crps[,1], result.tuononlinear.crps[,2], result.tuononlinear.crps[,3], result.tuononlinear.crps[,5], result.tuononlinear.crps[,4])

pnonadditive2 <- ggplot(data.frame(df.tuononlinear.crps, model), aes(x=model, y=df.tuononlinear.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

df.currin.crps <- c(result.currin.crps[,1], result.currin.crps[,2], result.currin.crps[,3], result.currin.crps[,5], result.currin.crps[,4])

pcurrin2 <- ggplot(data.frame(df.currin.crps, model), aes(x=model, y=df.currin.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

figure4 <- ggarrange(padditive, pnonadditive, pcurrin, padditive2, pnonadditive2, pcurrin2, 
                     ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

