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
###MAKING NEW PLOTS WITH THE NOT NORMALIZED VALUES

pSmall1 <- df_1000_im1_ii0.01 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Number of LCRs", "Avg Entropy of LCRs", "Sequence Entropy") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = col_1000, y = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_line() +
  #geom_smooth() +
  xlab("Iteration") +
  scale_colour_manual(values = c("black","red", "blue"))

pSmall1

##IM1 II 0.1
pSmall2 <- df_1000_im1_ii0.1 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Number of LCRs", "Avg Entropy of LCRs", "Sequence Entropy") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = col_1000, y = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_line() +
  #geom_smooth() +
  xlab("Iteration") +
  scale_colour_manual(values = c("black","red", "blue"))

pSmall2

##IM1 II 0.5
pSmall3 <- df_1000_im1_ii0.5 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Number of LCRs", "Avg Entropy of LCRs", "Sequence Entropy") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = col_1000, y = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_line() +
  #geom_smooth() +
  xlab("Iteration") +
  scale_colour_manual(values = c("black","red", "blue"))

pSmall3

##IM1 II 1
pSmall4 <- df_1000_im1_ii1 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Number of LCRs", "Avg Entropy of LCRs", "Sequence Entropy") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = col_1000, y = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_line() +
  #geom_smooth() +
  xlab("Iteration") +
  scale_colour_manual(values = c("black","red", "blue"))

pSmall4

##IM1 II 2
pSmall5 <- df_1000_im1_ii2 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Number of LCRs", "Avg Entropy of LCRs", "Sequence Entropy") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = col_1000, y = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_line() +
  #geom_smooth() +
  xlab("Iteration") +
  scale_colour_manual(values = c("black","red", "blue"))

pSmall5

##IM1 II 10
pSmall6 <- df_1000_im1_ii10 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Number of LCRs", "Avg Entropy of LCRs", "Sequence Entropy") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = col_1000, y = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_line() +
  #geom_smooth() +
  xlab("Iteration") +
  #geom_smooth() +
  scale_colour_manual(values = c("black","red", "blue"))

pSmall6

jpeg("im0.01-0.1-0.5.jpeg", units="in", width=12, height=6, res=300)
ggarrange( pSmall2, pSmall3,
          labels = c("ii0.1", "ii0.5"))
dev.off()

jpeg("im1-2-10.jpeg", units="in", width=12, height=6, res=300)
ggarrange(ggarrange(pSmall4, pSmall5, ncol = 2, labels = c("ii1", "ii2")), pSmall6, 
          labels = c("","ii10"), nrow = 2)
dev.off()
### Testing out histograms
pSmall7a <- df_1000_im1_ii0.01 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Avg Entropy of LCRs") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_histogram() +
  xlab("Value") +
  #geom_smooth() +
  scale_colour_manual(values = c("Orange","orange", "purple"))

pSmall7a

pSmall7b <- df_1000_im1_ii0.01 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Sequence Entropy") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_histogram() +
  xlab("Value") +
  #geom_smooth() +
  scale_colour_manual(values = c("Orange","orange", "purple"))

pSmall7b

pSmall7c <- df_1000_im1_ii0.01 %>% rename("Number of LCRs" = "V7", "Avg Entropy of LCRs" = "V8", "Sequence Entropy" = "V9") %>% select(col_1000, "Number of LCRs") %>%
  pivot_longer(-col_1000, names_to = "LCRInfo", values_to = "Values") %>%
  ggplot(aes(x = Values, group = LCRInfo, colour = LCRInfo)) +
  geom_histogram() +
  xlab("Value") +
  #geom_smooth() +
  scale_colour_manual(values = c("Orange","orange", "purple"))

pSmall7c

ggarrange(pSmall7a, pSmall7b, pSmall7c)




