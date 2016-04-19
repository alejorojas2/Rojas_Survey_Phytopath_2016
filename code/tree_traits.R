### Tree build for mapping trait data
#Libraries
library(ape)
library(seqinr)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(phangorn)
library(ggplot2)
library(ggtree)

#Convert sequences into alignment (assuming that they are aligned)
ITS.final <- read.fasta("../data/clean/ITS_spp.fasta", as.string = TRUE)
ITS.aln <- as.DNAbin(ape::as.alignment(ITS.final))

#Distance calculation
Oom.d <- dist.dna(ITS.aln)

#UPGMA Tree 
tr <- upgma(Oom.d)
tr <- root(tr, outgroup="Eurychasma_dicksonii")
#plot(tr, cex = 0.7, type = "fan")

trnj <- nj(Oom.d)
trnj <- root(trnj, outgroup="Eurychasma_dicksonii")
#plot(trnj, cex = 0.8, type = "fan")

cls <- list(Aphanomyces=c("Aphanomyces_cladogamus","Aphanomyces_cochlioides"),
            Clade_6=c("Phytophthora_aff_rosacearum","Phytophthora_inundata",
                      "Phytophthora_megasperma","Phytophthora_rosacearum"),
            Clade_7="Phytophthora_sojae",
            Clade_8=c("Phytophthora_drechsleri","Phytophthora_sansomeana"),
            Phytopythium=c("Phytopythium_aff_vexans","Phytopythium_chamaehyphon",
                           "Phytopythium_helicoides","Phytopythium_litorale",
                           "Phytopythium_megacarpum","Pythium_sterilum"),
            Pythiogeton=c("Pythiogeton_zeae"),
            Clade_A=c("Pythium_adhaerens","Pythium_aphanidermatum",
                      "Pythium_chondricola","Pythium_monospermum"),
            Clade_B=c("Pythium_aff_diclinum","Pythium_aff_dictyosporum",
                      "Pythium_aff_dissotocum","Pythium_aff_torulosum","Pythium_angustatum",
                      "Pythium_aristosporum","Pythium_arrhenomanes","Pythium_catenulatum",
                      "Pythium_coloratum","Pythium_conidiophorum","Pythium_contiguanum",
                      "Pythium_inflatum","Pythium_kashmirense","Pythium_lutarium",
                      "Pythium_oopapillum","Pythium_pachycaule","Pythium_periilum",
                      "Pythium_pyrilobum","Pythium_tardicrescens","Pythium_torulosum",
                      "Pythium_vanterpoolii"),
            Clade_D=c("Pythium_acanthicum","Pythium_amasculinum","Pythium_hydnosporum",
                      "Pythium_oligandrum","Pythium_periplocum"),
            Clade_E=c("Pythium_acrogynum","Pythium_aff_hypogynum","Pythium_camurandrum",
                      "Pythium_carolinianum","Pythium_hypogynum","Pythium_longandrum",
                      "Pythium_longisporangium","Pythium_middletonii","Pythium_minus",
                      "Pythium_pleroticum","Pythium_rhizosaccharum","Pythium_rostratifingens"),
            Clade_F=c("Pythium_attrantheridium","Pythium_cryptoirregulare",
                      "Pythium_intermedium","Pythium_irregulare","Pythium_kunmingense",
                      "Pythium_paroecandrum","Pythium_sp_balticum","Pythium_spinosum",
                      "Pythium_sylvaticum","Pythium_terrestris"),
            Clade_G=c("Pythium_aff_iwayamai","Pythium_nagaii"),
            Clade_I=c("Pythium_glomeratum","Pythium_heterothallicum","Pythium_ultimum",
                      "Pythium_ultimum_var_sporangiiferum","Pythium_ultimum_var_ultimum"),
            Clade_J=c("Pythium_acanthophoron","Pythium_aff_perplexum","Pythium_nodosum",
                      "Pythium_nunn","Pythium_orthogonon","Pythium_perplexum")
)

tr <- groupOTU(tr, cls)
trnj <- groupOTU(trnj, cls)

mypal <- colorRampPalette(brewer.pal(8, "Dark2" ))

P <- ggtree(tr, size=0.5, layout = "rectangular", ladderize = TRUE, aes(color=group)) +
  geom_text(aes(label=label), size=1.8, hjust=-0.1) +
  scale_color_manual(values = mypal(18), labels=c("Other",names(cls))) +
  #geom_text(aes(label=node), size=1.5, hjust = -0.2 ) +
  theme(plot.margin=unit(c(1,20,20,1), "mm"),
        legend.position="right") +
  ggplot2::xlim(0,0.8)

trnj$edge.length[165] <- trnj$edge.length[165]/1.9

P2 <- ggtree(trnj, size=0.4, layout = "rectangular", ladderize = TRUE, 
             branch.length = 'branch.length') +
  geom_tiplab(aes(label=label, color=group), size=2) +
  #geom_text(aes(label=node), size=3, hjust=-0.1) +
  scale_color_manual(values = mypal(18), labels=c("Other",names(cls))) +
  guides(colour="none") 

P2.1 <- P2 +
  geom_cladelabel(node = 125, 
                  label = "Clade E", offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 147,"Clade J", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 128,"Clade F", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 124,"Clade G", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node=139,"Clade I", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 92,"Clade B", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 95,"Clade A", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 112,"Clade D", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 159,"Clade 6", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 163,"Clade 8", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 156,"Clade 7", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) +
  geom_cladelabel(node = 153,"Phytopythium", 
                  offset=0.05,  
                  fontsize = 2, barsize = 2, align = T) 

#Read metadata
root.data <- read.table("../data/clean/Root_metadata.txt", sep="\t", header = TRUE, 
                        stringsAsFactors = F, row.names = 1) 
root.data1 <- root.data[c(12:16)]
root.data1 <- root.data1 %>% 
  add_rownames(var = "Species") %>%
  mutate(Sd13.c1 = Sd13.c + 2, Sd20.c1 = Sd20.c + 2) %>%
  select_("Species","ar.c", "ln.c", "wpr.c", "Sd13.c1", "Sd20.c1")


#Recode metadata using p-values
root.data1[is.na(root.data1)] <- "NA"

#Relabel columns
root.data1 <-  dplyr::rename(root.data1, "Root area (cm2)"=ar.c, 
                             "Root length (cm)"=ln.c, 
                             "Weight per root (mg)"=wpr.c,
                             "Seed rot 13ºC"=Sd13.c1,
                             "Seed rot 20ºC"=Sd20.c1)

root.data.hm <- root.data1[,-1]
rownames(root.data.hm) <- root.data1$Species
root.data.hm <- data.frame(root.data.hm)

#Original heatmap
t_p <- gheatmap(P2.1, root.data.hm, offset = -0.1, width = 0.2, colnames=F)
lbl <- get_heatmap_column_position(t_p, by="top")

tree_final_plot <- t_p + geom_text(data=lbl, aes(x, y, label=label), nudge_y = 2.5, 
                nudge_x = 0.01, angle=45, size=2) +
  scale_fill_manual(values=c("#969696","#333366","#969696","#016c59","#cccccc"),
                    name="Pathogenicity", labels=c("Seed\nNon-pathogenic",
                                                   "Pathogenic","Non-pathogenic",
                                                   "Seedling\nPathogenic", 
                                                   "Not determined")) +
  theme(legend.text=element_text(size=5))