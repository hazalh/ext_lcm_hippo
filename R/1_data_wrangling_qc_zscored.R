
# load R files -------------------------------------------------------------
source(here::here("R/library.R"))
source(here::here("R/functions.R"))



# load metadata - all samples ---------------------------------------------------------
metadata <- read_excel("data-raw/ex_hiplcm_metadata.xlsx", col_names = T) %>%
  arrange(area) %>% #ascending order: CA1, CA3, DG
  mutate(expgroup = as.factor(paste(.$group, .$area, sep = "_"))) %>%
  mutate_at(3:12, as.factor) #%>%
  #dplyr::filter(sample_id != "S18") #remove samples from metadata due to high NA

head(metadata)
str(metadata)

save(metadata, file = "data/metadata.RData")


# load metadata - for each dataset ----------------------------------------

metadata_ca1 <- read_excel("data-raw/ex_hiplcm_metadata.xlsx", col_names = T) %>%
  arrange(area) %>% #ascending order: CA1, CA3, DG
  mutate(expgroup = as.factor(paste(.$group, .$area, sep = "_"))) %>%
  mutate_at(3:12, as.factor) %>%
  dplyr::filter(area == "CA1")
  #dplyr::filter(area == "CA3")
  #dplyr::filter(area == "DG")
  #dplyr::filter(sample_id != "S18") #remove samples from metadata due to high NA

head(metadata_ca3)
str(metadata_ca3)

save(metadata_ca1, file = "data/metadata_ca1.RData")
save(metadata_ca3, file = "data/metadata_ca3.RData")
save(metadata_dg, file = "data/metadata_dg.RData")

# load data ---------------------------------------------------------
###### CA1 ######
## 5316 proteins spectronaut >> XXX diann

df_ca1 <- read.delim("data-raw/20250218_151236_ext_ca1_Report.tsv",
                     header = T, check.names = T, sep = "\t") %>%
    dplyr::select("PG.ProteinGroups", "PG.Genes", 6:20)

#df_ca1 <- read.delim("data-raw/newquant_beast/20241126_140011_26112024_ffmd_male_ca1_Report.tsv",
#                 header = T, check.names = T, sep = "\t") %>%
#  dplyr::select("PG.ProteinGroups", "PG.Genes", 5:16)

new_colnames <- sapply(strsplit(colnames(df_ca1), "_"), '[',10)
colnames(df_ca1) <- new_colnames

colnames(df_ca1)[1] <- "protein_names"
colnames(df_ca1)[2] <- "genes"
names(df_ca1)

table(df_ca1$genes == "") #13 doesn't have gene name - ca1
#any(!duplicated(df_ca1$genes))
#remove_rows <- df_ca1[duplicated(df_ca1$genes),] #check protein IDs in uniprot

df_ca1 <- df_ca1 %>%
  mutate(genes = ifelse(genes == "", protein_names, genes)) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate_at(3:17, as.numeric) %>%
  mutate(genes = make.names(genes, unique = TRUE)) %>%
  column_to_rownames(var = "genes") %>%
  dplyr::select(-protein_names) %>%
  log2()


#sample order sorted
sample_order <- metadata_ca1$sample_id
df_ca1 <- df_ca1[ ,sample_order]

#save data
save(df_ca1, file = "data/df_ca1_spec.RData")
write.csv(df_ca1, file = "data/df_ca1_spec.csv", row.names = T, quote = F)



##### CA3 ######
## 5872 proteins spectronaut >> 4508 diann


df_ca3 <- read.delim("data-raw/20250219_081114_ext_ca3_Report.tsv",
                     header = T, check.names = T, sep = "\t") %>%
  dplyr::select("PG.ProteinGroups", "PG.Genes", 6:20)

#df_ca3 <- read.delim("data-raw/newquant_beast/20241127_081037_20240209_Report_ffmd_male_ca3.tsv",
#                 header = T, check.names = T, sep = "\t") %>%
#  dplyr::select("Protein.Group", "Genes",  5:15)

new_colnames <- sapply(strsplit(colnames(df_ca3), "_"), '[',10)
colnames(df_ca3) <- new_colnames

colnames(df_ca3)[1] <- "protein_names"
colnames(df_ca3)[2] <- "genes"
names(df_ca3)

table(df_ca3$genes == "") #19 in ca3
#any(!duplicated(df_ca3$genes))
#remove_rows <- df_ca3[duplicated(df_ca3$genes),] #check protein IDs in uniprot

df_ca3 <- df_ca3 %>%
  mutate(genes = ifelse(genes == "", protein_names, genes)) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  mutate_at(3:17, as.numeric) %>%
  mutate(genes = make.names(genes, unique = TRUE)) %>%
  column_to_rownames(var = "genes") %>%
  dplyr::select(-protein_names) %>%
  log2()

#sample order sorted
sample_order <- metadata_ca3$sample_id
df_ca3 <- df_ca3[ ,sample_order]

#save data
save(df_ca3, file = "data/df_ca3_spec.RData")
write.csv(df_ca3, file = "data/df_ca3_spec.csv", row.names = T, quote = F)



###### DG ######
## 5449 proteins spectronaut >> diann

df_dg <- read.delim("data-raw/20250219_131121_20231012_ext_dg_Report.tsv",
                     header = T, check.names = T, sep = "\t") %>%
  dplyr::select("PG.ProteinGroups", "PG.Genes", 6:20)


#df_dg <- read.delim("data-raw/newquant_beast/20241129_132650_29112024_ffmd_male_dg_Report.tsv",
#                 header = T, check.names = T, sep = "\t") %>%
#    dplyr::select("PG.ProteinGroups", "PG.Genes", 5:13)

new_colnames <- sapply(strsplit(colnames(df_dg), "_"), '[',10)
colnames(df_dg) <- new_colnames

colnames(df_dg)[1] <- "protein_names"
colnames(df_dg)[2] <- "genes"
names(df_dg)

table(df_dg$genes == "") #14 in dg
#any(!duplicated(df_dg$genes))
#remove_rows <- df_dg[duplicated(df_dg$genes),] #check protein IDs in uniprot

df_dg <- df_dg %>%
    mutate(genes = ifelse(genes == "", protein_names, genes)) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    mutate_at(3:11, as.numeric) %>%
    mutate(genes = make.names(genes, unique = TRUE)) %>%
    column_to_rownames(var = "genes") %>%
    dplyr::select(-protein_names) %>%
    log2()


#sample order sorted
sample_order <- metadata_dg$sample_id
df_dg <- df_dg[ ,sample_order]

#save data
save(df_dg, file = "data/df_dg_spec.RData")
write.csv(df_dg, file = "data/df_dg_spec.csv", row.names = T, quote = F)



# NA ----------------------------------------------------------------------
visdat::vis_miss(df_ca1, sort_miss = T)
visdat::vis_miss(df_ca3, sort_miss = T)
visdat::vis_miss(df_dg, sort_miss = T)



# S curve -----------------------------------------------------------------
#### for each dataset separately.


# basic data vis - non-filtered, zscored datasets ----------------------------------------------------------------------

## boxplot
p <- create_boxplot(df_ca1, metadata_ca1)
p + ggtitle("ca1_boxplot")
ggsave("docs/ca1_boxplot.pdf")

p <- create_boxplot(df_ca3, metadata_ca3)
p + ggtitle("ca3_boxplot")
ggsave("docs/ca3_boxplot.pdf")

p <- create_boxplot(df_dg, metadata_dg)
p + ggtitle("dg_boxplot")
ggsave("docs/dg_boxplot.pdf")


##density
pdf("docs/ca1_diann_density.pdf", width = 8, height = 6)
plotDensities(df_ca1, main = "ca1_density")
dev.off()

pdf("docs/ca3_diann_density.pdf", width = 8, height = 6)
plotDensities(df_ca3, main = "ca3_density")
dev.off()

pdf("docs/dg_diann_density.pdf", width = 8, height = 6)
plotDensities(df_dg, main = "dg_density")
dev.off()


######   coefficient of variance ####
df_log10 <- read.delim("data-raw/newquant_cruella/20250114_114258_20240209_Report_ca3.tsv",
                 header = T, check.names = T, sep = "\t") %>%
    #dplyr::select(5:16) %>% #ca1
    dplyr::select(6:16) %>% #ca3
    #dplyr::select(5:13) %>% #dg
    mutate(mean_lfq = rowMeans(., na.rm = T)) %>%
    mutate(mean_cv = apply(., 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100))

ggplot(df_log10,
       aes(x = log2(mean_lfq),
           y = mean_cv)) +
    geom_point(size = 0.5, alpha = 0.5) +
    theme_minimal()

#ggsave("doc/ca1_cv.pdf")
#ggsave("doc/ca3_cv.pdf")
#ggsave("doc/dg_cv.pdf")







# filter each data separately  --------------------------------------------
## min 3 samples per ffmd and lfd
df_ca1_f <- selectGrps(df_ca1, metadata_ca1$group, n=1, 0.50) #5316 proteins >> 5091 diann
df_ca3_f <- selectGrps(df_ca3, metadata_ca3$group, n=1, 0.50) #5872 proteins >> 5522 diann
df_dg_f <- selectGrps(df_dg, metadata_dg$group, n=1, 0.43) #5449 proteins >> 5384 diann

save(df_ca1_f, file = "data/df_ca1_f.RData")
save(df_ca3_f, file = "data/df_ca3_f.RData")
save(df_dg_f, file = "data/df_dg_f.RData")



# z-score all 3 filtered datasets --------------------------------------------------

#df_dg_fmed <- apply(df_dg_f, 2, FUN = median, na.rm = T)
#df_dg_fz2 <- as.data.frame(scale(df_dg_f, center = df_dg_fmed, scale = T))
#df_dg_fz2$protein <- rownames(df_dg_fz2)
#zscoring works correctly for each df.

datasets <- list(ca1 = df_ca1_f,
                 ca3 = df_ca3_f,
                 dg = df_dg_f)

df_fz <- map(datasets, z_score_by_median)

df_fz <- reduce(df_fz, function(x, y) merge(x, y, by = "protein", all = TRUE)) #individually z-scored data and merged together
#6319 proteins - filtered beast >> 4804 filtered diann



# area-specific proteins ----------------------------------------------------------
df_fz <- df_fz %>%
  column_to_rownames(., var = "protein")

sample_order <- metadata$sample_id
df_fz <- df_fz[ ,sample_order]

#save data
save(df_fz, file = "data/df_all_fz.RData")
write_xlsx(df_fz, "data/df_all_fz.xlsx")


#plot
p <- plot_hipp_proteins(df_fz, metadata)

p$plot + ggtitle("df_fz_areaspecific_proteins")
ggsave("docs/df_all_fz_areaspecific_proteins.pdf")



# data vis of filtered and zscored data -----------------------------------

create_boxplot(df_fz, metadata)
ggsave(file = "docs/df_fz_boxplot_allareas.pdf", width=10, height = 8)


pdf("docs/df_fz_density_allareas.pdf", width=10, height = 8)
plotDensities(df_fz) # looks good
dev.off()




# mds plot - after filtering and zscoring ---------------------------------

### CA1 #####
df_fz_ca1 <- df_fz %>%
  select(all_of(colnames(df_fz)[colnames(df_fz) %in% metadata_ca1$sample_id])) %>%
  rownames_to_column(var = "protein") %>%
  rowwise() %>%
  dplyr::filter(sum(!is.na(c_across(-protein))) > 0) %>%
  ungroup() %>%
  column_to_rownames(var = "protein")

#batch corr
ca1_batchcorr <- limma::removeBatchEffect(df_fz_ca1, batch = metadata_ca1$batch)

mds_data <- plotMDS(df_fz_ca1, dim.plot = c(1,2))

mds_data <- data.frame(dim1 = mds_data$x, dim2 = mds_data$y,
                       samples = metadata_ca1$sample_id,
                       area = metadata_ca1$area,
                       group = metadata_ca1$group)

ggplot(mds_data,
       aes(x = dim1, y = dim2, colour = area, shape = group)) +
  geom_point(size = 4) +
  labs(title = "df_ca1_mds", colour = "area") +
  theme_minimal() +
  scale_colour_manual(values = c("CA1" = "#fa7876",
                                 "CA3" = "#edd9a3",
                                 "DG" = "#cc79a7"))

ggsave("docs/df_ca1_fz_mds.pdf")

save(df_fz_ca1, file = "df_fz_ca1.RData")

### CA3 #####
df_fz_ca3 <- df_fz %>%
  select(all_of(colnames(df_fz)[colnames(df_fz) %in% metadata_ca3$sample_id])) %>%
  rownames_to_column(var = "protein") %>%
  rowwise() %>%
  dplyr::filter(sum(!is.na(c_across(-protein))) > 0) %>%
  ungroup() %>%
  column_to_rownames(var = "protein")


#batch corr
ca3_batchcorr <- limma::removeBatchEffect(df_fz_ca3, batch = metadata_ca3$batch)

mds_data <- plotMDS(df_fz_ca3, dim.plot = c(1,2))

mds_data <- data.frame(dim1 = mds_data$x, dim2 = mds_data$y,
                       samples = metadata_ca3$sample_id,
                       area = metadata_ca3$area,
                       group = metadata_ca3$group)

ggplot(mds_data,
       aes(x = dim1, y = dim2, colour = area, shape = group)) +
  geom_point(size = 4) +
  labs(title = "df_ca3_mds", colour = "area") +
  theme_minimal() +
  scale_colour_manual(values = c("CA1" = "#fa7876",
                                 "CA3" = "#edd9a3",
                                 "DG" = "#cc79a7"))

ggsave("docs/df_ca3_fz_mds.pdf")

save(df_fz_ca3, file = "df_fz_ca3.RData")

### DG #####
df_fz_dg <- df_fz %>%
  select(all_of(colnames(df_fz)[colnames(df_fz) %in% metadata_dg$sample_id])) %>%
  rownames_to_column(var = "protein") %>%
  rowwise() %>%
  dplyr::filter(sum(!is.na(c_across(-protein))) > 0) %>%
  ungroup() %>%
  column_to_rownames(var = "protein")

#batch corr
dg_batchcorr <- limma::removeBatchEffect(df_fz_dg, batch = metadata_dg$batch)

mds_data <- plotMDS(df_fz_dg, dim.plot = c(1,2))

mds_data <- data.frame(dim1 = mds_data$x, dim2 = mds_data$y,
                       samples = metadata_dg$sample_id,
                       area = metadata_dg$area,
                       group = metadata_dg$group)

ggplot(mds_data,
       aes(x = dim1, y = dim2, colour = area, shape = group)) +
  geom_point(size = 4) +
  labs(title = "df_dg_mds", colour = "area") +
  theme_minimal() +
  scale_colour_manual(values = c("CA1" = "#fa7876",
                                 "CA3" = "#edd9a3",
                                 "DG" = "#cc79a7"))

ggsave("docs/df_dg_fz_mds.pdf")

save(df_fz_dg, file = "df_fz_dg.RData")



# MDS plots - dim1, dim2---------------------------------------------------------------
#load("data/df_all_fz_diann.RData")

mds_data <- plotMDS(df_fz, dim.plot = c(1,2))

mds_data = data.frame(dim1 = mds_data$x, dim2 = mds_data$y,
                      samples = metadata$sample_id,
                      area = metadata$area,
                      group = metadata$group,
                      termination_date = metadata$termination_date,
                      termination_batch = metadata$termination_batch,
                      sampleprep_date = metadata$sampleprep_date,
                      sampleprep_batch = metadata$sampleprep_batch,
                      cryo_date = metadata$cryo_date,
                      lcm_date = metadata$lcm_date,
                      duration_frozen_tissue = metadata$duration_frozen_slide,
                      duration_frozen_lcmtissue = metadata$duration_frozen_lcmtissue)

create_mds_plots(mds_data, metadata)



# batch correction -------------------------------------------------------
matrix_design_area_diet <- model.matrix(~ 0 + metadata_m$expgroup)

df_fz_bc <- limma::removeBatchEffect(df_fz, batch = metadata_m$batch,
                                     design = matrix_design_area_diet)

mds_data <- plotMDS(df_fz_bc, dim.plot = c(1,2))

mds_data = data.frame(dim1 = mds_data$x, dim2 = mds_data$y,
                      samples = metadata_m$sample_id,
                      area = metadata_m$area,
                      intervention = metadata_m$intervention,
                      batch = metadata_m$batch,
                      sex = metadata_m$sex,
                      cryo_date = metadata_m$cryostat_date,
                      lcm_date = metadata_m$lcm_date,
                      duration_frozen_tissue = metadata_m$freezing__tissue_post_termination,
                      duration_frozen_section = metadata_m$freezing_sections,
                      duration_frozen_postlcm = metadata_m$freezing_post_lcm)

create_mds_plots_bc(mds_data, metadata_m)


## other dims
mds_data <- plotMDS(df_fz, dim.plot = c(4,5))
mds_data <- data.frame(dim4 = mds_data$x, dim5 = mds_data$y,
                       samples = metadata_m$sample_id,
                       area = metadata_m$area,
                       batch = metadata_m$batch,
                       diet = metadata_m$intervention)

ggplot(mds_data,
       aes(x = dim4, y = dim5, colour = area, shape = batch)) +
  geom_point(size = 4) +
  labs(title = "df_dg_mds", colour = "area") +
  theme_minimal() +
  scale_colour_manual(values = c("CA1" = "#fa7876",
                                 "CA3" = "#edd9a3",
                                 "DG" = "#cc79a7"))

ggsave("doc/df_dg_fz_mds_diann.pdf")

# valid counts ------------------------------------------------------------
## use full data, and add the dashes for filter
load("data/df_ca1_diann_crn.RData")
load("data/df_ca3_diann_crn.RData")
load("data/df_dg_diann_crn.RData")


datasets <- list(ca1 = df_ca1,
                 ca3 = df_ca3,
                 dg = df_dg)

df_z <- map(datasets, z_score_by_median)
df_z <- reduce(df_z, function(x, y) merge(x, y, by = "protein", all = TRUE))
# 5427 proteins

df_z <- df_z %>%
  column_to_rownames(., var = "protein")

sample_order <- metadata_m$sample_id
df_z <- df_z[ ,sample_order]

barplot_validcounts(df_z, metadata_m)

ggsave(file = "doc/df_z_validcounts_diann.pdf", width=10, height = 8)


##double-check ids
counts_valid <- colSums(!is.na(df_z))

df_valid <- tibble(
  counts_valid = counts_valid,
  sampleid = metadata_m$sample_id,
  original_sampleid = metadata_m$sample_id_original,
  area = metadata_m$area
) %>%
  mutate(
    sampleid = factor(sampleid, levels = unique(sampleid[order(metadata_m$area)])),
    group = factor(area, levels = unique(area))
  )


# venn diagram ------------------------------------------------------------
## doublecheck total number of proteins by venn diagram

load("data/df_dg_diann_crn.RData")
load("data/df_ca1_diann_crn.RData")
load("data/df_ca3_diann_crn.RData")

protein_list <- list(
    "CA1" = rownames(df_ca1),
    "CA3" = rownames(df_ca3),
    "DG" = rownames(df_dg))

plot(euler(protein_list, shape = "ellipse"), quantities = T)

#dev.off()

venn.plot <- VennDiagram::venn.diagram(
  x = protein_list,
  category.names = c("CA1", "CA3", "DG"),
  filename = NULL,  # Output to plot in the R window
  output = F,
  main = "Venn Diagram of Protein IDs Across 3 Datasets",
  col = "black",
  fill = c("#fa7876", "#edd9a3", "#ca3c97"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.pos = c(0, 0, 0)
)

grid::grid.draw(venn.plot)
##not proportional

#ggVennDiagram::ggVennDiagram(protein_list)
## not proportional

venn_toplot <- calculate_venn_counts(df_ca1, df_ca3, df_dg,
                                     name1 = "ca1",
                                     name2 = "ca3",
                                     name3 = "dg")

## to make a proportional venn diagram
library(eulerr)

pdf("doc/venndiagram_df_z_proteinids_diann.pdf", width = 8, height = 6)
plot(euler(venn_toplot), quantities = T,
     fills = c("#fa7876", "#edd9a3", "#cc79a7"))
dev.off()









