#### Zoe Dyson
# example output is Klebs_CG23_spacer_list.csv

library(dplyr)
library(reshape2)
library(tidyr)

# set working direcotry containing crispr detect output
setwd("~/Documents/Data_and_Analyses/Population_structure/CRISPR/CG23_crisprdetect/cutoff2.5")

# list all gff files 
files <- list.files(pattern="gff")

# set up file to contain all crispr data
combined_crispr_data <- NULL

# compile all crispr data
for (file in 1:length(files)){
  crispr_data <- read.table(files[file], header=F)
  # add filename to the crispr data so we know what crispr comes from where
  crispr_data$V10 <- files[file]
  combined_crispr_data <- rbind(combined_crispr_data, crispr_data)
}

# assign headers with informative names to full crispr table
names(combined_crispr_data) <- c("accession_number", "detection_tool",
                                 "annotation", "start_loci", "end_loci",
                                 "length", "orientataion", "data","full_record",
                                 "filename")


# pull out repeat, direct repeat, and binding site data into separate tables
repeat_data <- combined_crispr_data[combined_crispr_data$annotation=="repeat_region",]
direct_repeat_data<- combined_crispr_data[combined_crispr_data$annotation=="direct_repeat",]
binding_site_data<- combined_crispr_data[combined_crispr_data$annotation=="binding_site",]

# convert to dplyr format
repeat_data <- tbl_df(repeat_data)
direct_repeat_data <- tbl_df(direct_repeat_data)
binding_site_data <- tbl_df(binding_site_data)

# Split full record into separate columns and reformat selected data
repeat_data <- repeat_data  %>%
  separate(full_record, c("id", "sequence","dbxref","ontology_term"),sep=";") %>%
  mutate(id=gsub("ID=", "", id)) %>%
  mutate(sequence=gsub("Note=", "", sequence)) %>%
  mutate(dbxref=gsub("Dbxref=", "", dbxref)) %>%
  mutate(ontology_term=gsub("Ontology_term=", "", ontology_term)) %>% 
  select(filename, annotation, start_loci, end_loci, length, orientataion,id,
        sequence,dbxref,ontology_term)

# Split full record into separate columns and reformat selected data
direct_repeat_data <- direct_repeat_data  %>%
  separate(full_record, c("id", "name", "parent","sequence","dbxref","ontology_term"),sep=";") %>%
  mutate(id=gsub("ID=", "", id)) %>%
  mutate(sequence=gsub("Note=", "", sequence)) %>%
  mutate(dbxref=gsub("Dbxref=", "", dbxref)) %>%
  mutate(ontology_term=gsub("Ontology_term=", "", ontology_term)) %>%
  mutate(parent=gsub("Parent=", "", parent)) %>%
  mutate(name=gsub("Name=", "", name)) %>% 
  select(filename, annotation, start_loci,end_loci, length, orientataion,id,
          name,parent,sequence,dbxref,ontology_term)

# Split full record into separate columns and reformat selected data
binding_site_data <- binding_site_data %>%
  separate(full_record, c("id", "name", "parent","sequence","dbxref","ontology_term"),sep=";") %>%
  mutate(id=gsub("ID=", "", id)) %>%
  mutate(sequence=gsub("Note=", "", sequence)) %>%
  mutate(dbxref=gsub("Dbxref=", "", dbxref)) %>%
  mutate(ontology_term=gsub("Ontology_term=", "", ontology_term)) %>%
  mutate(parent=gsub("Parent=", "", parent)) %>%
  mutate(name=gsub("Name=", "", name)) %>%
  select(filename, annotation, start_loci,end_loci, length, orientataion,id,
          name,parent,sequence,dbxref,ontology_term)

# write individual records to file
write.csv(repeat_data, "repeat_data.csv")
write.csv(direct_repeat_data, "direct_repeat_data.csv")
write.csv(binding_site_data, "binding_site_data.csv")

# create repeat fasta files
repeat_fasta_file<-file("repeats.fasta", "w")
for (line in 1:nrow(repeat_data)){
  cat(paste0(">", repeat_data$filename[line], "_",  repeat_data$annotation[line], "_", repeat_data$id[line], "\n"), file=repeat_fasta_file)
  cat(paste0(repeat_data$sequence[line], "\n"), file=repeat_fasta_file)
  }
close(repeat_fasta_file)

# create direct repeat fasta files
direct_repeat_fasta_file<-file("direct_repeats.fasta", "w")
for (line in 1:nrow(direct_repeat_data)){
  cat(paste0(">", direct_repeat_data$filename[line], "_",  direct_repeat_data$annotation[line], "_", direct_repeat_data$id[line], "\n"), file=direct_repeat_fasta_file)
  cat(paste0(direct_repeat_data$sequence[line], "\n"), file=direct_repeat_fasta_file)
}
close(direct_repeat_fasta_file)

# create binding site fasta files
binding_site_fasta_file<-file("binding_sites.fasta", "w")
for (line in 1:nrow(binding_site_data)){
  cat(paste0(">",  binding_site_data$filename[line], "_",  binding_site_data$annotation[line], "_", binding_site_data$id[line], "\n"), file=binding_site_fasta_file)
  cat(paste0(binding_site_data$sequence[line], "\n"), file=binding_site_fasta_file)
}
close(binding_site_fasta_file)


