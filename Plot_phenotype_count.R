# Import library
library(ggplot2)
library(tidyverse)

# Import data (csv format)
df = read.csv("Leoni_27h.csv")
View(df)   
# NOTE:  adjust the column name if necessary, column names should be Dead, Interphase, Prometaphase, and Well
colnames(df)
rownames(df) = df$Well
df$Well <-  NULL
df_count <- df %>% mutate(cell_count = rowSums(df))
df_count$cell_count
# Normalization: divide each cell type with the total number of cells detected
df_norm = data.frame(Dead = df_count$Dead/df_count$cell_count,
                     Interphase = df_count$Interphase/df_count$cell_count,
                     Prometaphase = df_count$Prometaphase/df_count$cell_count)
df_norm$Well = rownames(df_norm)

# NOTE: adjust this if there was some pipetting mistakes when pipetting the siRNA!
df_norm = df_norm %>% mutate(treatment = case_when(
  Well == "G10" | Well == "H13" | Well == "I15" ~ "siPLK1",
  Well == "G11" | Well == "H15" | Well == "I12" ~ "siKIF11",
  Well == "G12" | Well == "H14" | Well == "I10" ~ "siNC",
  Well == "G13" | Well == "H10" | Well == "I14" ~ "Monastrol",
  Well == "G14" | Well == "H12" | Well == "I11" ~ "DMSO",
  Well == "G15" | Well == "H11" | Well == "I13" ~ "Mock"
))
df_norm = na.omit(df_norm)
df_norm$Well <- NULL

df_norm_mean = df_norm %>%
  group_by(treatment) %>%
  summarise(Dead =  mean(Dead), Interphase = mean(Interphase), Prometaphase = mean(Prometaphase)) %>%
  pivot_longer(!treatment, names_to = "phenotype", values_to = "norm_count_mean")
df_norm_sd = df_norm %>%
  group_by(treatment) %>%
  summarise(Dead = sd(Dead), Interphase = sd(Interphase), Prometaphase = sd(Prometaphase)) %>%
  pivot_longer(!treatment, names_to = "phenotype", values_to = "norm_count_sd")
# Get individual data for plotting
df_rep1 = df_norm[1:6,] %>%
  pivot_longer(!treatment, names_to = "phenotype", values_to = "rep1")
df_rep2 = df_norm[7:12,] %>% pivot_longer(!treatment, names_to = "phenotype", values_to = "rep2")
df_rep3 = df_norm[13:18,] %>% pivot_longer(!treatment, names_to = "phenotype", values_to = "rep3")
df_plot = df_norm_mean %>%
  left_join(df_rep1,by=c("treatment","phenotype")) %>%
  left_join(df_rep2,by=c("treatment","phenotype")) %>%
  left_join(df_rep3,by=c("treatment","phenotype"))
df_plot$treatment = factor(df_plot$treatment, levels = c("siPLK1", "siKIF11", "siNC","Monastrol","DMSO","Mock"))
# Plot
ggplot(df_plot, aes(x=treatment,y = norm_count_mean, fill = phenotype))+
  theme_classic() + 
  geom_bar(position = "dodge", stat = "identity")+
  geom_point(aes(y=rep1),position = position_dodge(width = .9))+
  geom_point(aes(y=rep2),position = position_dodge(width = .9))+
  geom_point(aes(y=rep3),position = position_dodge(width = .9))+
  scale_y_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), expand = c(0,0.01))+
  labs(x = "Treatment", y = "Cell count / total cell count per well")+
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size =14))
