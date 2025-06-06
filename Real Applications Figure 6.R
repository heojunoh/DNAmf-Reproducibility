### Real Applications RMSE ###
model <- c(rep("DNAmf",100),rep("BM",100),rep("FBM",100))
model <- factor(model, levels=c("DNAmf", "BM", "FBM"))

df.poisson.rmse <- c(result.poisson.rmse[,1], result.poisson.rmse[,3], result.poisson.rmse[,2])

ppoisson <- ggplot(data.frame(df.poisson.rmse, model), aes(x=model, y=df.poisson.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Poisson") + labs(x="", y = "RMSE") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 5),
        panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))
df.plate.rmse <- c(result.plate.rmse[,1], result.plate.rmse[,5], result.plate.rmse[,4])

pplate <- ggplot(data.frame(df.plate.rmse, model), aes(x=model, y=df.plate.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Plate") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

df.heat.rmse <- c(result.heat.rmse[,1], result.heat.rmse[,5], result.heat.rmse[,4])

pheat <- ggplot(data.frame(df.heat.rmse, model), aes(x=model, y=df.heat.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Heat") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

### Real Applications CRPS ###

df.poisson.crps <- c(result.poisson.crps[,1], result.poisson.crps[,3], result.poisson.crps[,2])

ppoisson2 <- ggplot(data.frame(df.poisson.crps, model), aes(x=model, y=df.poisson.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + labs(x="", y = "CRPS") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 5),
        panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))
df.plate.crps <- c(result.plate.crps[,1], result.plate.crps[,5], result.plate.crps[,4])

pplate2 <- ggplot(data.frame(df.plate.crps, model), aes(x=model, y=df.plate.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

df.heat.crps <- c(result.heat.crps[,1], result.heat.crps[,5], result.heat.crps[,4])

pheat2 <- ggplot(data.frame(df.heat.crps, model), aes(x=model, y=df.heat.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


figure6 <- ggarrange(ppoisson, pplate, pheat, ppoisson2, pplate2, pheat2, 
                     ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

