library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)

#read file
young<- Read10X(data.dir = "D:/SKLAB/young HSC data")
aged<- Read10X(data.dir = "D:/SKLAB/aged HSC data")

#create seurat object

young <- CreateSeuratObject(counts = young, project = "young", min.cells = 3, min.features = 200)
aged <- CreateSeuratObject(counts = aged, project = "aged", min.cells = 3, min.features = 200)



#preprocessing

view(young@meta.data)
range(young$nFeature_RNA)
range(young$nCount_RNA)
young[["percent.mt"]] <- PercentageFeatureSet(object = young, pattern = "^mt.")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
young[["percent.mt"]] <- PercentageFeatureSet(young, pattern = "^MT-")


# Visualize QC metrics as a violin plot
VlnPlot(# Visualize QC metrics as a violin plot
VlnPlot(young, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


view(young@meta.data)
range(young$percent.mt)

VlnPlot(young, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

young <- subset(young, subset = nFeature_RNA > 250 & nFeature_RNA < 6000)
VlnPlot(young, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


#data normalization



view(young)

young@assays[["RNA"]]@data@x

young <- NormalizeData(young,
                         normalization.method = "LogNormalize", scale.factor = 10000)

#identification of highly variable features


young <- FindVariableFeatures(young,
                             selection.method = "vst", nfeatures = 2000)

#data scaling

young <- ScaleData(young) #2000 identified variable features

all.genes <- rownames(young)
young <- ScaleData(young, features = all.genes)

young@assays[["RNA"]]@scale.data


young<- RunPCA(young, features = VariableFeatures(object = young))

#cell cycle analysis

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/mus_musculus.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

young <- CellCycleScoring(young, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

young <- RenameIdents(young, "G1" = "Quiescent", "S" = "proliferative", "G2M" = "proliferative" )


cd_genes <- c("Rik",	"Abl1", "Abraxas1", "Actb",	"Actl6a",	"Actl6b",	"Actr2",	"Actr5",	"Actr8",	"Adprs",	"Alkbh1",	"Alkbh2",	"Alkbh3",	"Ankle1",	"Ap5s1",	"Ap5z1",	"Apbb1",	"Apex1",	"Apex2",	"Aplf",	"Aptx",	"Arid1a",	"Arid1b",	"Arid2",	"Ascc1",	"Ascc2",	"Ascc3",	"Asf1a",	"Atm",	"Atr",	"Atrip",	"Atrx", "Atxn7",	"Atxn7l3",	"Aunip",	"Axin2",	"Babam1", "Babam2",	"Bach1",	"Bard1",	"BC055324",	"Bccip",	"Bcl7a",	"Bcl7b",	"Bcl7c",	"Blm",	"Bod1l",	"Brca1",	"Brca2",	"Brcc3",	"Brcc3dc",	"Brd7", "Brd8",	"Brip1",	"Brme1",	"Cbx8",	"Ccdc117",	"Cdc14b",	"Cdc45",	"Cdc5l",	"Cdc7",	"Cdca5",	"Cdk2",	"Cdk7",	"Cdk9",	"Cdkn2d",	"Cebpg",	"Cenps",	"Cenpx",	"Cep164",	"Cetn1",	"Cetn2",	"Cetn4",	"Cgas", "Chaf1a",	"Chaf1b",	"Chchd4",	"Chd1l",	"Chd4",	"Chek1",	"Chek2",	"Chrna4",	"Cinp",	"Clspn",	"Commd1",	"Cul4a",	"Cul4b",	"Cyren",	"Dclre1a",	"Dclre1b",	"Dclre1c",	"Ddb1",	"Ddb2",	"Ddx1",	"Ddx11",	"Dek",	"Dhx9",	"Dmap1",	"Dmc1",	"Dna2",	"Dntt",	"Dot1l",	"Dpf1",	"Dpf2",	"Dpf3", "Dtl",	"Dtx3l",	"Dyrk1b", "Eepd1",	"Egfr",	"Eid3",	"Eme1",	"Eme2",	"Emsy", "Endov"	, "Eny2",	"Ep400",	"Epc1", "Epc2",	"Ercc1",	"Ercc2",	"Ercc3",	"Ercc4",	"Ercc5",	"Ercc6",	"Ercc6l2",	"Ercc8",	"Esco2",	"Etaa1",	"Exd2", "Exo1",	"Exo5",	"Eya1", "Eya2", "Eya3", "Eya4", "Faap100", "Faap20", "Faap24", "Fam111a",	"Fam168a", "Fan1", "Fanca", "Fancb",	"Fancc",	"Fancd2",	"Fancf",	"Fancg",	"Fanci",	"Fancl",	"Fancm",	"Fbh1",	"Fbxo6",	"Fbxw7",	"Fen1",	"Fgf10",	"Fh1",	"Fignl1",	"Fmn2",	"Fmr1",	"Foxm1",	"Fto",	"Fus",	"Fzr1",	"Gen1",	"Ggn",	"Gins2",	"Gins4",	"Gtf2h1",	"Gtf2h2",	"Gtf2h3",	"Gtf2h4",	"Gtf2h5",	"H2ac25",	"H2ax",	"Hdac10",	"Hdac9", "Hdgfl2","Helb",	"Helq",	"Herc2",	"Hinfp",	"Hltf",	"Hmces",	"Hmga1",	"Hmga1b",	"Hmga2",	"Hmgb1",	"Hmgn1",	"Hpf1",	"Hrob",	"Hsf1",	"Hsf2bp",	"Hspa1a",	"Hus1",	"Hus1b",	"Huwe1",	"Ier3",	"Iffo1",	"Ing3",	"Inip",	"Ino80",	"Ino80b",	"Ino80c",	"Ino80d",	"Ino80e",	"Ints3",	"Jmy",	"Kash5",	"Kat2a",	"Kat2b",	"Kat5",	"Kat7",	"Kdm1a",	"Kdm2a",	"Kdm4d",	"Khdc3",	"Kif22",	"Kin",	"Klhl15",	"Kmt5a",	"Kmt5b",	"Kmt5c",	"Lig1",	"Lig3",	"Lig4",	"Macroh2a1",	"Mad2l2",	"Majin",	"Marf1",	"Mbd4",	"Mbtd1",	"Mc1r",	"Mcm2",	"Mcm3",	"Mcm4",	"Mcm5",	"Mcm6",	"Mcm7",	"Mcm8", "Mcm9",	"Mcmdc2",	"Mcrs1",	"Mdc1",	"Meaf6",	"Meiob",	"Meioc",	"Mgme1",	"Mgmt",	"Mlh1",	"Mlh3",	"Mms19",	"Mms22l",	"Mnat1",	"Morc2b",	"Morf4l1", "Morf4l2",	"Mpnd", "Mre11a",	"Mrgbp",	"Mrnip",	"Msh2",	"Msh3",	"Msh4",	"Msh5",	"Msh6",	"Mta1",	"Mus81",	"Mutyh","Nabp1",	"Nabp2",	"Nbn",	"Neil1", "Neil2",	"Neil3",	"Nfrkb",	"Nhej1",	"Nipbl",	"Nono",	"Nop53",	"Npas2",	"Npm1",	"Nscme3l",	"Nsd2",	"Nsmce1",	"Nsmce2",	"Nsmce3",	"Nsmce4a",	"Nthl1",	"Nucks1",	"Nudt16l1",	"Ogg1",	"Ooep",	"Otub1",	"Otub2",	"Palb2",	"Parg",	"Park7",	"Parp1",	"Parp10",	"Parp2",	"Parp3",	"Parp9",	"Parpbp",	"Paxip1",	"Paxx", "Pbrm1",	"Pclaf",	"Pcna",	"Pds5a",	"Pds5b",	"Phf10",	"Pias4",	"Pif1",	"Pml",	"Pms1",	"Pms2",	"Pnkp",	"Pnp",	"Pogz",	"Pola1",	"Polb", "Pold1",	"Pold2",	"old3",	"Pold4",	"Poldip2",	"Pole",	"Pole2",	"Polg2",	"Polh",	"Poli", "Polk",	"Poll",	"Polm",	"Poln",	"Polq",	"Polr2i",	"Pot1a",	"Pot1b",	"Ppp4c",	"Ppp4r2",	"Ppp4r3b",	"Prdm9",	"Primpol",	"Prkcg",	"Prkdc",	"Prmt6",	"Prpf19",	"Psmd14",	"Psme4",	"Pttg1",	"Pwwp3a",	"Rad1",	"Rad17",	"Rad18",	"Rad21",	"Rad21l",	"Rad23a",	"Rad23b",	"Rad50",	"Rad51",	"Rad51ap1",	"Rad51b",	"Rad51c",	"Rad51d",	"Rad52",	"Rad54b",	"Rad54l",	"Rad9a",	"Rad9b",	"Radx",	"Rbbp8",	"Rbx1",	"Rec8",	"Recql",	"Recql4",	"Recql5",	"Rev1",	"Rev3l",	"Rexo4",	"Rfc1",	"Rfc2",	"Rfc3",	"Rfc4",	"Rfc5",	"Rfwd3",	"Rhno1",	"Rif1",	"Rmi1",	"Rmi2",	"Rnaseh2a",	"Rnaseh2b",	"Rnaseh2c",	"Rnf111",	"Rnf138",	"Rnf138rt1",	"Rnf168",	"Rnf169",	"Rnf8",	"Rpa1",	"Rpa2",	"Rpa3",	"Rpain",	"Rps3",	"Rrm1",	"Rrm2b",	"Rtel1",	"Ruvbl1",	"Ruvbl2",	"Samhd1",	"Sem1",	"Senp3",	"Setd2",	"Setmar",	"Setx",	"Sf3b3",	"Sf3b5",	"Sfpq",	"Sfr1",	"Sgf29",	"Shld1",	"Shld2",	"Shld3",	"Shprh",	"Sirt1",	"Sirt6",	"Sirt7",	"Slf1", 	"Slf2",	"Slx1b",	"Slx4",	"Smarca2",	"Smarca4",	"Smarca5",	"Smarcad1",	"Smarcal1",	"Smarcb1",	"Smarcc1",	"Smarcc2",	"Smarcd1",	"Smarcd2",	"Smarcd3",	"Smarce1",	"Smc1a",	"Smc2",	"Smc3",	"Smc4",	"Smc5",	"Smc6",	"Smchd1",	"Smg1",	"Smug1",	"Spata22",	"Spidr",	"Spire1",	"Spire2",	"Spo11",	"Sprtn",	"Ssrp1",	"Stub1",	"Supt16",	"Supt20",	"Supt3",	"Supt7l",	"Suv39h1",	"Swi5",	"Swsap1",	"Sycp1",	"Sycp3",	"Tada1",	"Tada2b",	"Tada3",	"Taf10",	"Taf12",	"Taf2","Taf4",	"Taf5",	"Taf5l",	"Taf6",	"Taf6l",	"Taf7",	"Taf9",	"Taok1",	"Taok3",	"Tdg",	"Tdg-ps",	"Tdp1",	"Tdp2",	"Terb1",	"Terb2",	"Terf2ip",	"Tex12",	"Tex15",	"Tex264",	"Tfpt",	"Ticrr",	"Tigar",	"Timeless",	"Tmem161a",	"Tnks1bp1",	"Tnp1",	"Tonsl",	"Top3a",	"Topbp1",	"Traip",	"Trex1",	"Trex2",	"Trim28",	"Trip12",	"Trip13",	"Trp53",	"Trp53bp1",	"Trpc2",	"Trrap",	"Ttc5",	"Ttf2",	"Twist1",	"Ube2a",	"Ube2b",	"Ube2d3",	"Ube2n",	"Ube2t",	"Ube2v1",	"Ube2v2",	"Ube2w",	"Ubqln4",	"Ubr5",	"Uchl5",	"Ufl1",	"Uhrf1",	"Uimc1",	"Ung",	"Upf1",	"Usp1",	"Usp10",	"Usp22",	"Usp28",	"Usp3",	"Usp45",	"Usp47",	"Usp51",	"Usp7",	"Uvrag",	"Uvssa",	"Vcp",	"Vcpip1",	"Vps72",	"Was",	"Wdhd1",	"Wdr48",	"Wrap53",	"Wrn",	"Wrnip1",	"Xab2", "Xpa",	"Xpc",	"Xrcc1",	"Xrcc2",	"Xrcc3",	"Xrcc4",	"Xrcc5",	"Xrcc6",	"Xrn2",	"Yeats4",	"Yy1",	"Zbtb1",	"
               Zbtb7a",	"Zcwpw1",	"Zfp365",	"Zfp668", "Zfyve26",	"Zmpste24",	"Zmynd8",	"Zranb3",	"Zswim7")

cd_genes <- c("Pcna", "Rpa3","Rfc4", "Rfc5", "Cetn2", "Ddb1", "Neil1", "Neil3", "Xrcc1", "Tdg", "Topbp1", "Palb2", "Mre11a", "Brca1", "Blm")


DoHeatmap(object = young, features = cd_genes)

#deseq

young.markers <- FindAllMarkers(young, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

#merge seurat object

Yonage <- merge(young, y = c(aged),
                add.cell.ids = ls()[1:2],
                project = "Yonage")


library(fgsea)

DefaultAssay(aged) <- "RNA"

aged <- FindMarkers(aged, ident.1 = "aged Proliferative", min.pct = 0.25, logfc.threshold = 0.25)


aged$gene <- rownames(aged)
aged <- aged %>% arrange(desc(avg_log2FC))
write.csv(fgseaaged, file="G:/SKLAB/agedmarker.csv")

fold_changes <- aged$avg_log2FC
names(fold_changes) <- aged$gene
head(fold_changes)

#Load GSEA gene sets: http://www.gsea-msigdb.org/gsea/msigdb/index.jsp

Reactome <- fgsea::gmtPathways("G:/SKLAB/scRNA analysis/fgsea/m2.cp.reactome.v2023.2.Mm.symbols.gmt")
hallmark <- fgsea::gmtPathways("G:/SKLAB/scRNA analysis/fgsea/mh.all.v2023.2.Mm.symbols.gmt")

fgseaaged <- fgsea(pathways = hallmark, stats = fold_changes, eps = 0.0, minSize = 15, maxSize = 500)
fgseaaged <- fgsea(pathways = aged, stats = fold_changes, eps = 0.0, minSize = 15, maxSize = 500)

ggplot(hallmark,aes(x=pathway,y=NES,fill= pval)) + 
  geom_col(position="dodge",width=0.7,linewidth=0.8) +
  coord_flip() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
