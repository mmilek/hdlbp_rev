#figS4

dat<-read.delim("data/s4d.txt", header=T)
mel<-melt(dat[,1:3])
ggbarplot(mel, x = "protein", y = "value",
          ylab= "% total TC transitions", xlab = "protein", color = "protein", fill = "protein",
          add = c("mean_sd"), palette = "jco",
          position = position_dodge(.8)) +
  geom_jitter(aes(protein, value, fill = protein), shape = 21, color = "black",  position = position_jitterdodge(jitter.height = -2, jitter.width = 1.5))
