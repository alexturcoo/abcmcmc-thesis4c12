###IMPORTING THE TAB DELIM FILES FOR 10 DISTANCES
df_10_im0_ii1 <- read.delim("10_distances_plots/10_distances_march14_N0.0-1.0_1000sims_1mutation_im0_ii1.txt", header = FALSE)
df_10_im1_ii0 <- read.delim("10_distances_plots/10_distances_march15_N0.0-1.0_1000sims_1mutation_im1_ii0.txt", header = FALSE)
df_10_im1_ii1 <- read.delim("10_distances_plots/10_distances_march15_N0.0-1.0_1000sims_1mutation_im1_ii1.txt", header = FALSE)

###IMPORTING THE TAB DELIM FILES FOR 100 DISTANCES
df_100_im0_ii1 <- read.delim("100_distances_plots/100_distances_march15_N0.0-1.0_1000sims_1mutation_im0_ii1.txt", header = FALSE)
df_100_im1_ii0 <- read.delim("100_distances_plots/100_distances_march15_N0.0-1.0_1000sims_1mutation_im1_ii0.txt", header = FALSE)
df_100_im1_ii1 <- read.delim("100_distances_plots/100_distances_march15_N0.0-1.0_1000sims_1mutation_im1_ii1.txt", header = FALSE)

###IMPORTING THE TAB DELIM FILES FOR 1000 DISTANCES
df_1000_im0_ii1 <- read.delim("1000_distances_plots/1000_distances_im0_ii1.txt", header = FALSE)
df_1000_im1_ii1 <- read.delim("1000_distances_plots/1000_distances_im1_ii1.txt", header = FALSE)

### ADDING A COLUMN OF EITHER 1:10 OR 1:100 FOR PLOTTING
col_10 <- c(1:10)
col_100 <- c(1:100)
col_1000 <- c(1:1000)

df_10_im0_ii1 <- cbind(df_10_im0_ii1, col_10)
df_10_im1_ii0 <- cbind(df_10_im1_ii0, col_10)
df_10_im1_ii1 <- cbind(df_10_im1_ii1, col_10)

df_100_im0_ii1 <- cbind(df_100_im0_ii1, col_100)
df_100_im1_ii0 <- cbind(df_100_im1_ii0, col_100)
df_100_im1_ii1 <- cbind(df_100_im1_ii1, col_100)

df_1000_im0_ii1 <- cbind(df_1000_im0_ii1, col_1000)
df_1000_im1_ii1 <- cbind(df_1000_im1_ii1, col_1000)

### GRAPHING THE SHIT
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())
### 10 Distances, IM0 II1
ggplot(df_10_im0_ii1, aes(x=col_10, y=V5)) + 
  geom_point() +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 0 Indel Rate 1")

### 10 Distances, IM1 II0
ggplot(df_10_im1_ii0, aes(x=col_10, y=V5)) + 
  geom_point() +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 1 Indel Rate 0")

### 10 Distances, IM1 II1
ggplot(df_10_im1_ii1, aes(x=col_10, y=V5)) + 
  geom_point() +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 1 Indel Rate 1")

### 100 Distances, IM0 II1
ggplot(df_100_im0_ii1, aes(x=col_100, y=V5)) + 
  geom_point() +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 0 Indel Rate 1")

### 100 Distances, IM1 II0
ggplot(df_100_im1_ii0, aes(x=col_100, y=V5)) + 
  geom_point() +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 1 Indel Rate 0")

### 100 Distances, IM1 II1
ggplot(df_100_im1_ii1, aes(x=col_100, y=V5)) + 
  geom_point() +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 1 Indel Rate 1")

### 1000 Distances, IM0 II1
x <- ggplot(df_1000_im0_ii1, aes(x=col_1000, y=V5)) + 
  geom_point(color="red") +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 0 Indel Rate 1")

### 1000 Distances, IM1 II1
y <- ggplot(df_1000_im1_ii1, aes(x=col_1000, y=V5)) + 
  geom_point(color="red") +
  xlab("Iteration") +
  ylab("Distance") +
  ggtitle("Mutation Rate 1 Indel Rate 1")

jpeg("distances.jpeg", units="in", width=12, height=6, res=300)
ggarrange(x,y)
dev.off()



