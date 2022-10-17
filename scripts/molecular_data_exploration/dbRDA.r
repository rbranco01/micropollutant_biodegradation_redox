### db-RDA using chemical data

library(tidyverse)
library(ggord)
library(ggtext)
library(glue)
library(lubridate)
library(phia)
library(kableExtra)
library(MASS)

## r db_RDA chemical data

## This analysis uses a phyloseq object sent by Rita on 20220816 and chemical data that should be joined first

# load phyloseq object
phy_rita = readRDS("output_data/psdata_RBRA_sent_by_Rita.rds")

# chemical data
chem_data = readxl::read_xlsx("input_data/RBRA_proj2_Chemical_data.xlsx")

# sample_data
sam_dat = 
  sample_data(phy_rita) %>% data.frame() %>% select(SampleID = sampleid, everything()) %>% tibble

# joined data
chem_data_added = sam_dat %>% inner_join(., chem_data, by = "SampleID")

# NOTE inocula are filtered out because they do not have data on micropollutant removal:
is_innoculum = setdiff(sam_dat$SampleID, chem_data$SampleID)
subset(sam_dat, SampleID %in% is_innoculum)


# overwriting old metadata
phy_rita_redox = subset_samples(phy_rita, !sample_names(phy_rita) %in% is_innoculum)
phy_rita_redox = prune_taxa(taxa_sums(phy_rita_redox)>0, phy_rita_redox)

sam_data = sample_data(chem_data_added)
sample_names(sam_data) = chem_data_added$SampleID
sample_data(phy_rita_redox) = sam_data




### r capscale with partialling out sampletype

# rarefy
min_ss = min(sample_sums(phy_rita_redox))

# rarefy data
set.seed(711)
phy_rita_redox_16095 = rarefy_even_depth(phy_rita_redox, sample.size = min_ss, rngseed = 711, trimOTUs = T)

# total number of reads
tot_read = phy_rita_redox_16095 %>% sample_sums(.) %>% sum()

# create taxa acronyms
df_taxa_names = data.frame(OTU = taxa_names(phy_rita_redox_16095),
                           ASV = paste0("ASV_",seq(1:ntaxa(phy_rita_redox_16095))))

# rename colnames with acronyms
abund = t(otu_table(phy_rita_redox_16095))
cap_data = abund/rowSums(abund)
colnames(cap_data) = df_taxa_names$ASV

rownames(cap_data)
rownames(chem_data_added) = chem_data_added$SampleID

# predictors
chem_pred = 
  chem_data_added %>% 
  dplyr::select(Sample_type, Redox_type, GAB:MET) %>% 
  mutate(MET_ESA = as.numeric(MET_ESA),
         MET_OA = as.numeric(MET_OA)) %>% 
  drop_na(MET_OA) %>% 
  drop_na(MET_ESA) %>% 
  mutate(across(where(is.numeric), ~if_else(.<0, 0, .))) %>% 
  mutate(across(where(is.numeric), ~scale(.)),
         Redox_type = factor(Redox_type, levels = c("aerobic", "nitrate", "iron", "sulphate", "methanogenic")))



# GGally::ggpairs(chem_pred, columns = 2:ncol(chem_pred), 
#                 aes(colour=Sample_type), progress = F)



# remove missing samples (row 23)

cap <- capscale(cap_data[-23,] ~ 
                  GAB + ACK + CBZ + 
                  MCPP + `24_D` + ANP +
                  BAM + BTZ + CLZ_MDP +
                  CLZ + DCB + MET_ESA +
                  MET_OA + MET
                + Condition(Sample_type), # partial analysis now
                chem_pred, 
                dist="jaccard",
                na.action = na.exclude)

anova(cap)
summary(cap)
plot(cap)


# model definitions (all predictor removal efficiencies are added)

full.model <- capscale(cap_data[-23,] ~
                         GAB + ACK + CBZ +
                         MCPP + `24_D` + ANP +
                         BAM + BTZ + CLZ_MDP +
                         CLZ + DCB + MET_ESA +
                         MET_OA + MET,
                       chem_pred[,-c(1:2)],
                       dist="jaccard",
                       na.action = na.exclude)

redox.model <- capscale(cap_data[-23,] ~ 
                          GAB + ACK + CBZ + 
                          MCPP + `24_D` + ANP +
                          BAM + BTZ + CLZ_MDP +
                          CLZ + DCB + MET_ESA +
                          MET_OA + MET +
                          Condition(Sample_type), # partial analysis now
                        chem_pred[,-2], 
                        dist="jaccard",
                        na.action = na.exclude)

# sample_type.model <- capscale(cap_data[-23,] ~ 
#                   GAB + ACK + CBZ + 
#                   MCPP + `24_D` + ANP +
#                   BAM + BTZ + CLZ_MDP +
#                   CLZ + DCB + MET_ESA +
#                   MET_OA + MET +
#                   Condition(Redox_type), # partial analysis now
#                 chem_pred[,-1], 
#                 dist="jaccard",
#                 na.action = na.exclude)

# model selection

# selectedMod_full <- step(full.model)
# selectedMod_st <- step(sample_type.model)
selectedMod_red <- step(redox.model)

smry = summary(full.model)

# anova tables with fdr-corrected q values
# anova_tab_full = anova(selectedMod_full, by = "terms") %>% as_tibble(rownames = "terms")
# anova_tab_full %>% 
#   mutate(R2 = 100*(SumOfSqs/sum(SumOfSqs))) %>% 
#   mutate(fdr_q = p.adjust(`Pr(>F)`, "fdr")) %>% # correct p values for multiple testingfilter(fdr_q < 0.1)
#   filter(fdr_q < 0.1 | is.na(fdr_q))
# 
# anova_tab_st = anova(selectedMod_st, by = "terms") %>% as_tibble(rownames = "terms")
# anova_tab_st %>% 
#   mutate(R2 = 100*(SumOfSqs/sum(SumOfSqs))) %>% 
#   mutate(fdr_q = p.adjust(`Pr(>F)`, "fdr")) %>% # correct p values for multiple testingfilter(fdr_q < 0.1)
#   filter(fdr_q < 0.1 | is.na(fdr_q))

anova_tab_red = anova(selectedMod_red, by = "terms") %>% as_tibble(rownames = "terms")
anova_tab_red = anova_tab_red %>% 
  mutate(R2 = 100*(SumOfSqs/sum(SumOfSqs))) %>% 
  mutate(fdr_q = p.adjust(`Pr(>F)`, "fdr")) %>% # correct p values for multiple testingfilter(fdr_q < 0.1)
  filter(fdr_q < 0.1 | is.na(fdr_q))

print(anova_tab_red)


# plot dbRDA model focussed on redox_conditions (i.e. + Condition(Sample_type))

plot_data = chem_data_added %>% select(SampleID, Sample_type,Redox_type, any_of(smry$biplot %>% rownames()))
df1  <- data.frame(smry$sites[,1:2]) %>% as_tibble(rownames = "SampleID")   # CAP1 and CAP2
plot_data = df1 %>% inner_join(., plot_data, by = "SampleID" )


# plot CAP ordinations
# axis 1 vs 2
ord_red_ax12 = 
  ggord(selectedMod_red,
        xlims = c(-1, 1.25),
        ylims = c(-0.8, 0.6),
        arrow = 0.2,
        axes = c("1","2"),
        grp_in = plot_data$Redox_type,
        ellipse =F, 
        coord_fix = T,
        # ptslab = T, 
        # size = NA,
        addsize = 0,
        add_col = NULL,
        ext = 1.2,
        vec_ext = 0.8,
        var_sub = anova_tab_red$terms[-length(anova_tab_red$terms)]) +
  labs(title = "Partial dbRDA - effect of Redox condition; CAP1 vs CAP3")

# axis 1 vs 3
ord_red_ax13 = 
  ggord(selectedMod_red,
        xlims = c(-1, 1.25),
        ylims = c(-0.9, 0.6),
        arrow = 0.2,
        axes = c("1","3"),
        grp_in = plot_data$Redox_type,
        ellipse =F, 
        coord_fix = T,
        # ptslab = T, 
        # size = NA,
        addsize = 0,
        ext = 1.2,
        vec_ext = 0.8,
        var_sub = anova_tab_red$terms[-length(anova_tab_red$terms)]) +
  labs(title = "Partial dbRDA - effect of Redox condition; CAP1 vs CAP3")


# make ggord, remove points layer and add shapes etc
ord_red_ax12$layers[[1]] <- NULL # to adjust shape as well
ord_red_ax13$layers[[1]] <- NULL # to adjust shape as well

ord_red_ax12$data <- ord_red_ax12$data %>% 
  mutate(newgrp = plot_data$Sample_type, # add extra group for sample_type
         samplecheck = plot_data$SampleID,
         Groups = factor(Groups, levels = c("aerobic", "nitrate", "iron", "sulphate", "methanogenic")),
         newgrp = factor(newgrp, levels = c("soil", "mun", "ditch", "ind")))

ord_red_ax13$data <- ord_red_ax13$data %>% 
  mutate(newgrp = plot_data$Sample_type, # add extra group for sample_type
         samplecheck = plot_data$SampleID,
         Groups = factor(Groups, levels = c("aerobic", "nitrate", "iron", "sulphate", "methanogenic")),
         newgrp = factor(newgrp, levels = c("soil", "mun", "ditch", "ind")))

# shape vals
shapes_red = c(22, 21, 4, 24, 25)
colors_red = c("#C92127", "#F3A484", "#93C4DD", "#0772B2")

# plot ordinations again
plot_ord_red_ax12 = 
  ord_red_ax12 + geom_point(aes(shape = Groups, colour = newgrp, group = newgrp), size=2, stroke=1.5, alpha=1) + 
  scale_shape_manual(name = "",
                     values = shapes_red) +
  scale_color_manual(name = "",
                     values = colors_red,
                     labels = c("Soil","Mun AS","Ditch", "Ind AS")) + 
  guides(fill = guide_legend(override.aes= list(shape = NA))) + 
  theme(legend.title = element_blank())

plot_ord_red_ax13 = 
  ord_red_ax13 + geom_point(aes(shape = Groups, colour = newgrp, group = newgrp), size=2, stroke=1.5, alpha=1) + 
  scale_shape_manual(name = "",
                     values = shapes_red) +
  scale_color_manual(name = "",
                     values = colors_red,
                     labels = c("Soil","Mun AS","Ditch", "Ind AS")) + 
  guides(fill = guide_legend(override.aes= list(shape = NA))) + 
  theme(legend.title = element_blank())

# plot separately and legend
plot_ord_red_ax12 + theme(text = element_text(size = unit(8, "pt")),
                          plot.margin=margin(t=1, b=1, unit="cm"),
                          legend.position = "bottom") + guides(color=guide_legend(nrow=4,byrow=TRUE),
                                                               shape=guide_legend(nrow=5,byrow=TRUE))
plot_ord_red_ax13 + theme(text = element_text(size = unit(8, "pt")),
                          plot.margin=margin(t=1, b=1, unit="cm"))

legend_bar_CAP = 
  get_legend(plot_ord_red_ax12 + theme( legend.position = "bottom",
                                        legend.box = "vertical",
                                        legend.margin=margin(),
                                        legend.key.size = unit(1, 'cm'), #change legend key size
                                        legend.key.height = unit(1, 'cm'), #change legend key height
                                        legend.key.width = unit(1, 'cm'), #change legend key width
                                        legend.title = element_text(size=0), #change legend title font size
                                        legend.text = element_text(size=12)))

# combine plots
CAP_plot = 
  plot_grid(
    plot_ord_red_ax12 + theme(legend.position = "none"), 
    plot_ord_red_ax13 + theme(legend.position = "none"), 
    labels  = c("A", "B"), 
    ncol = 2, 
    nrow = 1, 
    align = "hv")


CAP_AB_plot <- plot_grid(CAP_plot, legend_bar_CAP, ncol=1, nrow = 2, rel_heights = c(5, 1))

# save plot
ggsave(CAP_AB_plot, filename = "figures/RBRA_proj2_CAPscale_Redox_Cond-SampleType.pdf", width = 12, height = 5)


# format anova table
anova_tab_red %>% 
  kableExtra::kbl(digits = 3, 
                  caption = "Constrained ordination using db-RDA (Jaccard distances)<br> Associations between microbial community membership and MP removal efficiencies<br>Terms = MPs with FDR q<0.1 (stepwise selection)") %>% 
  kableExtra::kable_classic()

