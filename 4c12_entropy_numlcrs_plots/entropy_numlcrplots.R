###IMPORTING THE TAB DELIM FILES FOR 1000 DISTANCES
df_1000_im1_ii0.01 <- read.delim("10_distances_im1_ii0.01.txt", header = FALSE)
df_1000_im1_ii0.1 <- read.delim("10_distances_im1_ii0.1.txt", header = FALSE)
df_1000_im1_ii0.5 <- read.delim("10_distances_im1_ii0.5.txt", header = FALSE)
df_1000_im1_ii1 <- read.delim("10_distances_im1_ii1.txt", header = FALSE)
df_1000_im1_ii2 <- read.delim("10_distances_im1_ii2.txt", header = FALSE)
df_1000_im1_ii10 <- read.delim("10_distances_im1_ii10.txt", header = FALSE)

### ADDING A COLUMN OF EITHER 1:10 OR 1:100 FOR PLOTTING
col_1000 <- c(1:1000)

df_1000_im1_ii0.01 <- cbind(df_1000_im1_ii0.01, col_1000)
df_1000_im1_ii0.1 <- cbind(df_1000_im1_ii0.1, col_1000)
df_1000_im1_ii0.5 <- cbind(df_1000_im1_ii0.5, col_1000)
df_1000_im1_ii1 <- cbind(df_1000_im1_ii1, col_1000)
df_1000_im1_ii2 <- cbind(df_1000_im1_ii2, col_1000)
df_1000_im1_ii10 <- cbind(df_1000_im1_ii10, col_1000)


### GRAPHING THE SHIT
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
theme_set(theme_pubr())

###V7 is avg num of LCRs
###V8 is avg entropy of LCRs
###V9 is Sequence Entropy
### col_1000 is the indexes
### 1000 Distances, IM0 II0.01
###V7 is avg num of LCRs
###V8 is avg entropy
### col_1000 is the indexes
### 1000 Distances, IM1 II0.01
a1 <- ggplot(df_1000_im1_ii0.01, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  xlab("Iteration") +
  ylab("# of LCRs")

b1 <- ggplot(df_1000_im1_ii0.01, aes(x=col_1000)) + 
  geom_line(aes(y = V8)) +
  xlab("Iteration") +
  ylab("Avg Entropy of LCRs")

c1 <- ggplot(df_1000_im1_ii0.01, aes(x=col_1000)) + 
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("Avg Entropy of Sequence")

figure1 <- ggarrange(a1, b1, c1,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)
figure1

#####PUTTING LINES ON SAME PLOT
ggplot(df_1000_im1_ii0.01, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  geom_line(aes(y = V8)) +
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("# of LCRs")

###V7 is avg num of LCRs
###V8 is avg entropy
### col_1000 is the indexes
### 1000 Distances, IM0 II0.01
a2 <- ggplot(df_1000_im1_ii0.1, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  xlab("Iteration") +
  ylab("# of LCRs")

b2 <- ggplot(df_1000_im1_ii0.1, aes(x=col_1000)) + 
  geom_line(aes(y = V8)) +
  xlab("Iteration") +
  ylab("Avg Entropy of LCRs")

c2 <- ggplot(df_1000_im1_ii0.1, aes(x=col_1000)) + 
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("Avg Entropy of Sequence")

figure2 <- ggarrange(a2, b2, c2,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)
figure2

#####PUTTING LINES ON SAME PLOT
ggplot(df_1000_im1_ii0.1, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  geom_line(aes(y = V8)) +
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("# of LCRs")

###V7 is avg num of LCRs
###V8 is avg entropy
### col_1000 is the indexes
### 1000 Distances, IM0 II0.01
a3 <- ggplot(df_1000_im1_ii0.5, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  xlab("Iteration") +
  ylab("# of LCRs")

b3 <- ggplot(df_1000_im1_ii0.5, aes(x=col_1000)) + 
  geom_line(aes(y = V8)) +
  xlab("Iteration") +
  ylab("Avg Entropy of LCRs")

c3 <- ggplot(df_1000_im1_ii0.5, aes(x=col_1000)) + 
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("Avg Entropy of Sequence")

figure3 <- ggarrange(a3, b3, c3,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)
figure3

#####PUTTING LINES ON SAME PLOT
ggplot(df_1000_im1_ii0.5, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  geom_line(aes(y = V8)) +
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("# of LCRs")

###V7 is avg num of LCRs
###V8 is avg entropy
### col_1000 is the indexes
### 1000 Distances, IM0 II0.01
a4 <- ggplot(df_1000_im1_ii1, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  xlab("Iteration") +
  ylab("# of LCRs")

b4 <- ggplot(df_1000_im1_ii1, aes(x=col_1000)) + 
  geom_line(aes(y = V8)) +
  xlab("Iteration") +
  ylab("Avg Entropy of LCRs")

c4 <- ggplot(df_1000_im1_ii1, aes(x=col_1000)) + 
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("Avg Entropy of Sequence")

figure4 <- ggarrange(a4, b4, c4,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)
figure4

#####PUTTING LINES ON SAME PLOT
ggplot(df_1000_im1_ii1, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  geom_line(aes(y = V8)) +
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("# of LCRs")

###V7 is avg num of LCRs
###V8 is avg entropy
### col_1000 is the indexes
### 1000 Distances, IM0 II0.01
a5 <- ggplot(df_1000_im1_ii2, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  xlab("Iteration") +
  ylab("# of LCRs")

b5 <- ggplot(df_1000_im1_ii2, aes(x=col_1000)) + 
  geom_line(aes(y = V8)) +
  xlab("Iteration") +
  ylab("Avg Entropy of LCRs")

c5 <- ggplot(df_1000_im1_ii2, aes(x=col_1000)) + 
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("Avg Entropy of Sequence")

figure5 <- ggarrange(a5, b5, c5,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)
figure5

#####PUTTING LINES ON SAME PLOT
ggplot(df_1000_im1_ii2, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  geom_line(aes(y = V8)) +
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("# of LCRs")

###V7 is avg num of LCRs
###V8 is avg entropy
### col_1000 is the indexes
### 1000 Distances, IM1 II10
a6 <- ggplot(df_1000_im1_ii10, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  xlab("Iteration") +
  ylab("# of LCRs")

b6 <- ggplot(df_1000_im1_ii10, aes(x=col_1000)) + 
  geom_line(aes(y = V8)) +
  xlab("Iteration") +
  ylab("Avg Entropy of LCRs")

c6 <- ggplot(df_1000_im1_ii10, aes(x=col_1000)) + 
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("Avg Entropy of Sequence")

figure6 <- ggarrange(a6, b6, c6,
                      labels = c("A", "B", "C"),
                      ncol = 2, nrow = 2)
figure6

#####PUTTING LINES ON SAME PLOT
ggplot(df_1000_im1_ii10, aes(x=col_1000)) + 
  geom_line(aes(y = V7)) + 
  geom_line(aes(y = V8)) +
  geom_line(aes(y = V9)) +
  xlab("Iteration") +
  ylab("# of LCRs")

pBIG <- df_1000_im1_ii10 %>% ggplot(aes(x = col_1000, y = V9)) +
  geom_line(colour = "blue")

pSmall <- df_1000_im1_ii10 %>% select(col_1000, V7, V8) %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = col_1000, y = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_line() +
  scale_colour_manual(values = c("red","green"))

ggarrange(pBIG, pSmall, ncol = 1)
