setwd("D:/Feng/DIA_PASEF/20250409_PentaPASEF_Classic_DIA/By_Gradient/45min")

library(arrow)
library(diann)
library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(ggpubr)
library(stringr)
library(ggbreak)
library(ggplot2)
df <- read_parquet("report.parquet")
PentaPASEF_ClassicDIA_Files <- read_excel("45min.xlsx")
## change the run name to real sample name based on the list, to avoid the mistakes in the sample queue 
df$Run <- dplyr::recode(df$Run, !!!setNames(PentaPASEF_ClassicDIA_Files$NewSampleName, PentaPASEF_ClassicDIA_Files$FileName))
unique(df$Run) ## check if the rename is sucessful
df$File.Name <- df$Run

## filter out false assignments:contaminants
clean_data <- df %>%
  filter(!str_detect(Protein.Group, "cRAP-"))
##remove decoy
clean_data <- clean_data %>%
  filter(!str_detect(Decoy, "1"))

clean_data <- clean_data %>%
  mutate(species = case_when(
    str_detect(Protein.Names, "_9LACO") ~ "L.murinus",
    str_detect(Protein.Names, "_SALRD") ~ "S.ruber",
    TRUE ~ "background"
  ))

data1 <- clean_data %>%
  mutate(
    Gradient = case_when(
      grepl("^5min_FAST", Run) ~ "5min_FAST",
      grepl("^5min", Run) ~ "5min",
      grepl("^22min", Run) ~ "22min",
      grepl("^45min", Run) ~ "45min",
      TRUE ~ NA_character_
    )
  )


## remove the false assignments of labeled L.m, S.r, labeled background precursors
precursors_to_exclude <- data1 %>%
  filter(grepl("NSCtrl", Run)) %>%
  group_by(Gradient) %>%
  filter(
    (species == "L.murinus"  & grepl("UniMod:259|UniMod:267", Precursor.Id, ignore.case = TRUE)) |
      (species == "background" & grepl("UniMod:259|UniMod:267", Precursor.Id, ignore.case = TRUE)) |
      (species == "S.ruber")
  ) %>%
  select(Gradient, Precursor.Id, species) %>%
  distinct()

write.csv(precursors_to_exclude,"precursors_to_exclude_basedOn_NGCtrl.csv")

data2 <- data1 %>%
  left_join(precursors_to_exclude, by = c("Gradient", "Precursor.Id")) %>%
  filter(is.na(species.y)) %>%     # Keep only non-excluded rows
  select(-species.y) %>%
  rename(species = species.x)

## extract the precursors either based on the protein group qvalue or gene group qvalue
precursors <- diann_matrix(data2, id.header="Precursor.Id",q = 0.01, quantity.header = "Precursor.Normalised")
precursors <- as.data.frame(precursors)
precursors$Precursor <- rownames(precursors)


precursors$ProteinGroup <- data2$Protein.Group[match(row.names(precursors),data2$Precursor.Id)]
precursors$ProteinName <- data2$Protein.Names[match(row.names(precursors),data2$Precursor.Id)]
precursors$GeneNames <- data2$Genes[match(row.names(precursors),data2$Precursor.Id)]
precursors$StrippedSequence <- data2$Stripped.Sequence[match(row.names(precursors),data2$Precursor.Id)]
precursors$ModifiedSequence <- data2$Modified.Sequence[match(row.names(precursors),data2$Precursor.Id)]

write.csv(precursors[,c(37:42,1:36)], "Filtered_Precursors_all_FDR0.01.csv", na="NA",eol = "\n")

count_table_pr <- precursors[,-c(38:42)] %>%
  pivot_longer(
    cols = -Precursor,               
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Sample) %>%
  summarise(Precursor_Count = n_distinct(Precursor)) %>%
  ungroup()
library(stringr)
# Extract dilution replicate for coloring
count_table_pr <- count_table_pr %>%
  mutate(Dilution_Replicate = str_extract(Sample, "Dil\\d+"),
         Gradient = str_extract(Sample, ".*(?=_Dil)")
  )

count_table_pr_clean <- count_table_pr[!grepl("NSCtrl",count_table_pr$Sample),]

write.csv(count_table_pr, "Precursor_IDs.csv")

f1 <- ggplot() +
  geom_bar(
    data = count_table_pr_clean[grepl("45min",count_table_pr_clean$Gradient),], 
    mapping = aes(x = Sample, y = Precursor_Count, fill = Dilution_Replicate),
    alpha = 0.6,
    colour = "black",
    position = "dodge",
    stat = "identity"
  ) +
  theme_bw() +
  ggtitle("PrecursorIDs_DIA-PASEF_45min") +
  theme(
    plot.title = element_text(size = rel(1.2), lineheight = 0.9, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.title = element_text(size = 12, face = "plain"),
    panel.grid.minor = element_blank()
  ) +
  ylab("Number of precursors") +
  scale_fill_brewer(palette = "Dark2")+scale_y_continuous(
    limits = c(0, 55000),
    breaks = seq(0, 55000, by = 10000)
  )


# extract the peptide MaxLFQ quantity either based on stripped sequence or modified sequence. This is recommended by the author

peptides.maxlfq <- diann_maxlfq(data2[data2$Q.Value <= 0.01 & data2$PG.Q.Value <= 0.01,], sample.header = "Run",group.header="Modified.Sequence", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")


peptides.maxlfq <- as.data.frame(peptides.maxlfq)
peptides.maxlfq$ModifiedPeptide <- rownames(peptides.maxlfq)

peptides.maxlfq$ProteinGroup <- data2$Protein.Group[match(row.names(peptides.maxlfq),data2$Modified.Sequence)]
peptides.maxlfq$ProteinName <- data2$Protein.Names[match(row.names(peptides.maxlfq),data2$Modified.Sequence)]
peptides.maxlfq$GeneNames <- data2$Genes[match(row.names(peptides.maxlfq),data2$Modified.Sequence)]
peptides.maxlfq$StrippedSequence <- data2$Stripped.Sequence[match(row.names(peptides.maxlfq),data2$Modified.Sequence)]

write.csv(peptides.maxlfq[,c(37:41,1:36)], "Filtered_peptides_Quant_all_FDR0.01_45min.csv", na="NA",eol = "\n")

count_table_Pep <- peptides.maxlfq[,-c(38:41)] %>%
  pivot_longer(
    cols = -ModifiedPeptide,               
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Sample) %>%
  summarise(Peptide_Count = n_distinct(ModifiedPeptide)) %>%
  ungroup()
library(stringr)
# Extract dilution replicate for coloring
count_table_Pep <- count_table_Pep %>%
  mutate(Dilution_Replicate = str_extract(Sample, "Dil\\d+"),
         Gradient = str_extract(Sample, ".*(?=_Dil)")
  )

count_table_Pep_clean <- count_table_Pep[!grepl("NSCtrl",count_table_Pep$Sample),]
write.csv(count_table_Pep, "ModifiedPeptides_IDs_45min.csv")


##remove one file: 22min_Dil2_2_2 due to low IDs
##count_table_Pep_clean <- count_table_Pep_clean[!grepl("22min_Dil1_2_2", count_table_Pep_clean$Sample),]

f5 <- ggplot() +
  geom_bar(
    data = count_table_Pep_clean[grepl("45min",count_table_Pep_clean$Gradient),],
    mapping = aes(x = Sample, y = Peptide_Count, fill = Dilution_Replicate),
    alpha = 0.6,
    colour = "black",
    position = "dodge",
    stat = "identity"
  ) +
  theme_bw() +
  ggtitle("PeptideIDs_DIA-PASEF_45min-FAST") +
  theme(
    plot.title = element_text(size = rel(1.2), lineheight = 0.9, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.title = element_text(size = 12, face = "plain"),
    panel.grid.minor = element_blank()
  ) +
  ylab("Number of peptides") +
  scale_fill_brewer(palette = "Dark2")+scale_y_continuous(
    limits = c(0, 55000),
    breaks = seq(0, 55000, by = 10000)
  )


## peptide subset for L.m and S.r
peptides.maxlfq <- peptides.maxlfq %>%
  mutate(species = case_when(
    str_detect(ProteinName, "_9LACO") ~ "L.murinus",
    str_detect(ProteinName, "_SALRD") ~ "S.ruber",
    TRUE ~ "background"
  ))


Peptide_9LACO_filtered <- peptides.maxlfq %>%
  filter(
    (
      (grepl("L.murinus", species) & grepl("UniMod:259|UniMod:267", ModifiedPeptide))
      
    )
  )
# Subset for S.r
Peptide_SALRD_filtered <- peptides.maxlfq %>%
  filter(
    (
      (grepl("S.ruber", species) & !grepl("UniMod:259|UniMod:267", ModifiedPeptide))
      
    )
  )

count_table_Pep_LM <- Peptide_9LACO_filtered[,-c(38:42)] %>%
  pivot_longer(
    cols = -ModifiedPeptide,               
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Sample) %>%
  summarise(Peptide_Count = n_distinct(ModifiedPeptide)) %>%
  ungroup()
library(stringr)
# Extract dilution replicate for coloring
count_table_Pep_LM <- count_table_Pep_LM %>%
  mutate(Dilution_Replicate = str_extract(Sample, "Dil\\d+"),
         Gradient = str_extract(Sample, ".*(?=_Dil)")
  )

count_table_Pep_LM_clean <- count_table_Pep_LM[!grepl("NSCtrl",count_table_Pep_LM$Sample),]

write.csv(count_table_Pep_LM, "Filtered_ModifiedPeptides_IDs_LM_45min-FAST.csv")
##remove one file: 45min_Dil2_2_2 due to low IDs
##count_table_Pep_LM_clean <- count_table_Pep_LM_clean[!grepl("45min_Dil1_2_2", count_table_Pep_LM_clean$Sample),]


f9 <- ggplot() +
  geom_bar(
    data = count_table_Pep_LM_clean[grepl("45min",count_table_Pep_LM_clean$Gradient),],
    mapping = aes(x = Sample, y = Peptide_Count, fill = Dilution_Replicate),
    alpha = 0.6,
    colour = "black",
    position = "dodge",
    stat = "identity"
  ) +
  theme_bw() +
  ggtitle("L.murinus_PeptideIDs_DIA-PASEF_45min-FAST") +
  theme(
    plot.title = element_text(size = rel(1.2), lineheight = 0.9, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.title = element_text(size = 12, face = "plain"),
    panel.grid.minor = element_blank()
  ) +
  ylab("Number of peptides") +
  scale_fill_brewer(palette = "Dark2")+ylim(0,800)+ scale_y_cut(breaks=c(20,200), which = c(1,2,3),scales=c(0.5,0.5,0.5), expand = F)


count_table_Pep_SR <- Peptide_SALRD_filtered[,-c(38:42)] %>%
  pivot_longer(
    cols = -ModifiedPeptide,               
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Sample) %>%
  summarise(Peptide_Count = n_distinct(ModifiedPeptide)) %>%
  ungroup()
library(stringr)
# Extract dilution replicate for coloring
count_table_Pep_SR <- count_table_Pep_SR %>%
  mutate(Dilution_Replicate = str_extract(Sample, "Dil\\d+"),
         Gradient = str_extract(Sample, ".*(?=_Dil)")
  )
## remove NSCtrl samples
count_table_Pep_SR_clean <- count_table_Pep_SR[!grepl("NSCtrl",count_table_Pep_SR$Sample),]

write.csv(count_table_Pep_SR, "Filtered_ModifiedPeptides_IDs_SR_45min.csv")

##remove one file: 45min_Dil2_2_2 due to low IDs
#count_table_Pep_SR_clean <- count_table_Pep_SR_clean[!grepl("45min_Dil1_2_2", count_table_Pep_SR_clean$Sample),]

f13 <- ggplot() +
  geom_bar(
    data = count_table_Pep_SR_clean[grepl("45min",count_table_Pep_SR_clean$Gradient),],
    mapping = aes(x = Sample, y = Peptide_Count, fill = Dilution_Replicate),
    alpha = 0.6,
    colour = "black",
    position = "dodge",
    stat = "identity"
  ) +
  theme_bw() +
  ggtitle("S.ruber_PeptideIDs_DIA-PASEF_45min-FAST") +
  theme(
    plot.title = element_text(size = rel(1.2), lineheight = 0.9, face = "plain", colour = "black", vjust = 0.5, hjust = 0.5),
    axis.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12, face = "plain", colour = "black"),
    axis.title = element_text(size = 12, face = "plain"),
    panel.grid.minor = element_blank()
  ) +
  ylab("Number of peptides") +
  scale_fill_brewer(palette = "Dark2")+ylim(0,1800)+scale_y_cut(breaks=c(20, 200), which = c(1,2,3),scales=c(0.5,0.5,0.5), expand = F)


count_table_Pep_LM_clean <- count_table_Pep_LM_clean %>%
  mutate(Group = "L.murinus")

count_table_Pep_SR_clean <- count_table_Pep_SR_clean %>%
  mutate(Group = "S.ruber")
LM_SR_combined_count <- bind_rows(count_table_Pep_LM_clean, count_table_Pep_SR_clean)

write_tsv(LM_SR_combined_count, "Filtered_combined_peptide_count_LM_SR_45min.tsv")




##  injection replicate CV

library(dplyr)
library(stringr)
library(ggpubr)

data3 <- peptides.maxlfq

##remove one file: 22min_Dil2_2_2 due to low IDs
##data4 <- data3[,!grepl("22min_Dil1_2_2",colnames(data3))]

## remove NGCtrl samples

data4 <- data3[,!grepl("NSCtrl",colnames(data3))]


data5 <- data4[,-c(29:32)] %>%
  pivot_longer(
    cols = -c(ModifiedPeptide,species),               
    names_to = "Sample",
    values_to = "Intensity"
  )

data5 <- data5 %>%
  mutate(Dilution_Replicate = str_extract(Sample, ".*Dil\\d+"),
         Dilution_ID = str_extract(Sample, ".*Dil\\d+_\\d+"),
         Gradient = str_extract(Sample, ".*(?=_Dil)")
  )

## MS injection replicate
data6 <- data5 %>%
  group_by(Dilution_ID, Dilution_Replicate, ModifiedPeptide, species, Gradient) %>%
  summarise(
    Avg_Intensity = mean(Intensity, na.rm = TRUE),
    SD_Intensity = sd(Intensity, na.rm = TRUE),
    CV = (SD_Intensity / Avg_Intensity) * 100,
    .groups = "drop"
  )

## dilution replicated CV
data7 <- data5 %>%
  group_by(Dilution_Replicate, ModifiedPeptide, species, Gradient) %>%
  summarise(
    Avg_Intensity = mean(Intensity, na.rm = TRUE),
    SD_Intensity = sd(Intensity, na.rm = TRUE),
    CV = (SD_Intensity / Avg_Intensity) * 100,
    .groups = "drop"
  )

f17 <- ggplot(data6[data6$Gradient == "45min",], aes(x=Dilution_ID, ## to have the x-axis in a centrain order you prefer
                                                     y=CV, fill=Dilution_Replicate)) + geom_boxplot(alpha=0.5, width=0.8,outliers = F, outlier.shape = 21, outlier.size = 1.5, na.rm = T)+ ylim(0,max(data6$CV)) + ylab("CV")+
  scale_fill_brewer(palette = "Dark2")+geom_abline(intercept = 20, slope = 0, linetype=2, size=0.8,colour="red")## to change the box size ## specify different colors for differnt conditions

f17 <- f17 + theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 12, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), 
        axis.title.y = element_text(size=12, face = "plain"),axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks.x = element_blank())+ ggtitle("MS-replicate_CV_DIA-PASEF")+ylim(0,100)


f18 <- ggplot(data7[data7$Gradient == "45min",], aes(x=Dilution_Replicate, ## to have the x-axis in a centrain order you prefer
                                                     y=CV, fill=Dilution_Replicate)) + geom_boxplot(alpha=0.5, width=0.8,outliers = F, outlier.shape = 21, outlier.size = 1.5,na.rm = T)+ ylim(0,max(data7$CV)) + ylab("CV")+
  scale_fill_brewer(palette = "Dark2")+geom_abline(intercept = 20, slope = 0, linetype=2, size=0.8,colour="red")## to change the box size ## specify different colors for differnt conditions

f18 <- f18 + theme_bw() + theme(plot.title=element_text(size=rel(1.5), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(legend.position = "top", legend.text = element_text(size = 12, face = "plain"),legend.title = element_text(size = 12, face = "plain")) + 
  theme(axis.text.x = element_blank(),axis.text = element_text(size = 12, face = "plain", color = "black"), 
        axis.title.y = element_text(size=12, face = "plain"),axis.title.x = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks.x = element_blank())+ ggtitle("Dil-replicate_CV_DIA-PASEF")+ylim(0,100)
ggarrange(f17, f18, widths = c(1,0.6), common.legend = T)
write.csv(data6, "45min-FAST_peptide_averaged_intensity_by_injection_CV.csv")
write.csv(data7, "45min-FAST_peptide_averaged_intensity_by_Dilution_CV.csv")


## extracting dilution 1 as example to plot the dynamic range

data8 <- data7[data7$Dilution_Replicate == "45min_Dil1",]
data8<- data8 %>% arrange(desc(Avg_Intensity))
data8$LogIntensity <- log2(data8$Avg_Intensity)
data9 <- data8[data8$LogIntensity != "-Inf",]
data9 <- data8[data8$LogIntensity != "NaN",]
data9$Rank <- 1:nrow(data9)

data9_SR_LM <- data9 %>%
  filter(str_detect(species, "L.murinus|S.ruber"))

data9_SR_LM <- data9_SR_LM %>%
  mutate(
    Rank_nudged = case_when(
      species == "L.murinus" ~ Rank - 200,  # shift left
      species == "S.ruber" ~ Rank + 200,    # shift right
      TRUE ~ as.numeric(Rank)
    )
  )

data9_SR_LM <- data9_SR_LM %>%
  filter(
    !(
      (grepl("L.murinus", species) & !grepl("UniMod:259|UniMod:267",ModifiedPeptide))
    )
  )


ggplot()+ geom_line(data=data9, mapping=aes(x=Rank, y=LogIntensity), color="wheat3", size=2) + 
  geom_point(data=data9_SR_LM, mapping=aes(x=Rank_nudged, y=LogIntensity, fill = species, colour = species, shape = species), size=4, alpha=0.5) + 
  theme_bw()+
  theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), legend.position = "top",
        legend.text = element_text(size = 12, face = "plain", color = "black"),axis.title = element_text(size = 12,face = "plain", color="black")) + 
  xlab("Peptide rank")+ ylab("Log2(peptide intensity)")+ scale_shape_manual(values = c(21,24))+ ggtitle("DynamicRange_DIA-PASEF-45min-FAST-Dil1")+
  theme(panel.grid.minor = element_blank()) + scale_fill_manual(values = c("L.murinus"="lightblue","S.ruber"="salmon")) + scale_colour_manual(values = c("lightblue","salmon"))

## Bench mark quantitative accuracy
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
## remove the false assignments of L.murinus (non-labeled)
data10 <- data7 %>%
  filter(
    !(
      (grepl("L.murinus", species) & !grepl("UniMod:259|UniMod:267",ModifiedPeptide)) |
        (grepl("S.ruber", species) & grepl("UniMod:259|UniMod:267",ModifiedPeptide)) |
        (grepl("background", species) & grepl("UniMod:259|UniMod:267",ModifiedPeptide))
    )
  )


## Dil1 and Dil2
data11 <- data10 %>%
  select(Dilution_Replicate, ModifiedPeptide, species,Gradient,Avg_Intensity) %>%
  pivot_wider(names_from = Dilution_Replicate, values_from = Avg_Intensity) %>%
  filter(!is.na(`45min_Dil1`) & !is.na(`45min_Dil2`)) %>%
  mutate(
    Log2_A = log2(`45min_Dil1`),
    Log2_B = log2(`45min_Dil2`),
    log2ratio = Log2_A - Log2_B,
    log2B = Log2_B
  )
write.csv(data11, "DIA_45min_QuantitativeAccuracy_Dil1-Dil2.csv", row.names = F)
hline_df <- data.frame(
  species = c("background", "L.murinus", "S.ruber"),
  hline_y = c(0, log2(10), log2(10))
)

ggplot(data11, aes(x = log2B, y = log2ratio, color = species)) +
  geom_point(alpha = 0.5, size = 2.5) +
  geom_smooth(se = T, method = "gam", formula = y ~ s(x), linetype = "dashed", color = "gray6") +
  
  geom_hline(data = hline_df, aes(yintercept = hline_y, color = species),size=1,
             linetype = "dashed") +
  scale_color_manual(values = c("L.murinus" = "lightblue", "S.ruber" = "salmon", "background" = "tan1")) +
  labs(
    x = expression(Log[2](Dil2)),
    y = expression(Log[2](Dil1:Dil2)),
    color = "Species"
  ) + theme_bw() +facet_wrap(~species, nrow = 3, scales = "free_y")+theme(plot.title=element_text(size=rel(1.2), lineheight=.9,face="plain", colour="black", vjust = 0.5, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 12, face = "plain", color = "black"), legend.position = "top",
        legend.text = element_text(size = 12, face = "plain", color = "black"),axis.title = element_text(size = 12,face = "plain", color="black"), strip.text = element_text(size = 12,face = "italic", color="black")) + 
  scale_shape_manual(values = c(21,24))+theme(panel.grid.minor = element_blank())+ ggtitle("QuantAccuracy_DIA-PASEF-45min Dil1 Vs. Dil2")+ylim(-6,6)

save.image("D:/Feng/DIA_PASEF/20250409_PentaPASEF_Classic_DIA/By_Gradient/45min/Rdata_45min.RData")





## Prepare data for taxa annotation
## average the peptide intensity from injection replicates 
data10 <- data4 %>%
  pivot_longer(
    cols = -c(ModifiedPeptide,species, ProteinGroup, ProteinName, GeneNames, StrippedSequence),               
    names_to = "Sample",
    values_to = "Intensity"
  )

data10 <- data10 %>%
  mutate(Dilution_Replicate = str_extract(Sample, ".*Dil\\d+"),
         Dilution_ID = str_extract(Sample, ".*Dil\\d+_\\d+"),
         Gradient = str_extract(Sample, ".*(?=_Dil)")
  )

## MS injection replicate
data11 <- data10 %>%
  group_by(Dilution_ID,ModifiedPeptide, ProteinGroup,ProteinName,GeneNames,StrippedSequence,species, Gradient) %>%
  summarise(
    Avg_Intensity = mean(Intensity, na.rm = TRUE),
    SD_Intensity = sd(Intensity, na.rm = TRUE),
    CV = (SD_Intensity / Avg_Intensity) * 100,
    .groups = "drop"
  )

write.csv(data11, "45min_Injection_Averaged_Intensity.csv", row.names = F)

## remove human peptides
data12 <- data11[!grepl("HUMAN", data11$ProteinName),]

## spread the data into wide format
library(tidyr)
library(dplyr)

data13 <- data12 %>%
  pivot_wider(
    id_cols = c(ModifiedPeptide,StrippedSequence,ProteinGroup, ProteinName, GeneNames, species, Gradient),  # columns to keep as identifiers
    names_from = Dilution_ID,        # the column to pivot
    values_from = Avg_Intensity      # values to fill in
  )

write.csv(data13, "45min_Injection_Averaged_Intensity_WO_HumanPeptide.csv", row.names = F)

data14 <- data13[,c(1,2,8:16)]
colnames(data14) <- paste0(colnames(data14),".d")

write.csv(data14, "45min_Injection_Averaged_Intensity_For_MetaLab.csv", row.names = F)



protein.groups <- diann_maxlfq(data2[data2$Q.Value <= 0.01 & data2$PG.Q.Value <= 0.01,], sample.header = "Run",group.header="Protein.Group", id.header = "Precursor.Id", quantity.header = "Precursor.Normalised")

protein.groups <- as.data.frame(protein.groups)
protein.groups$ProteinName <- data2$Protein.Names[match(row.names(protein.groups),data2$Protein.Group)]
protein.groups$GeneNames <- data2$Genes[match(row.names(protein.groups),data2$Protein.Group)]

write.csv(protein.groups[,c(37,38,1:36)], "Filtered_ProteinGroup_Quant_FDR0.01.csv", na="NA",eol = "\n")

## average the intensity for injection replicates
data15 <- protein.groups
data15$ProteinGroup <- row.names(protein.groups)

data16 <- data15 %>%
  pivot_longer(
    cols = -c(ProteinGroup, ProteinName, GeneNames),               
    names_to = "Sample",
    values_to = "Intensity"
  )

data16$Sample <- gsub("NSCtrl","Dil0", data16$Sample)

data17 <- data16 %>%
  mutate(Dilution_Replicate = str_extract(Sample, ".*Dil\\d+"),
         Dilution_ID = str_extract(Sample, ".*Dil\\d+_\\d+"),
         Gradient = str_extract(Sample, ".*(?=_Dil)")
  )


data18 <- data17 %>%
  group_by(Dilution_ID,ProteinGroup,ProteinName,GeneNames,Gradient) %>%
  summarise(
    Avg_Intensity = mean(Intensity, na.rm = TRUE),
    SD_Intensity = sd(Intensity, na.rm = TRUE),
    CV = (SD_Intensity / Avg_Intensity) * 100,
    .groups = "drop"
  )

write.csv(data18, "ProteinGroup_45min_Injection_Averaged_Intensity_CV.csv", row.names = F)

library(tidyr)
library(dplyr)

data19 <- data18 %>%
  pivot_wider(
    id_cols = c(ProteinGroup, ProteinName, GeneNames, Gradient),  # columns to keep as identifiers
    names_from = Dilution_ID,        # the column to pivot
    values_from = Avg_Intensity      # values to fill in
  )

data19$ProteinGroup <- gsub(";", "; ", data19$ProteinGroup)


write.csv(data19, "ProteinGroup_45min_Injection_Averaged_Intensity.csv", row.names = F)
save.image("D:/Feng/DIA_PASEF/20250409_PentaPASEF_Classic_DIA/By_Gradient/45min/Rdata_45min.RData")



