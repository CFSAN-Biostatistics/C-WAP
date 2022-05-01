library(tidyverse)

m1 = read_tsv("data/gisaid/gisaid-patient-metadata-2020.tsv")
m2 = read_csv("data/gisaid/gisaid-patient-metadata-2021.csv")
m = bind_rows(m1,m2)

pango = tibble()
for(i in Sys.glob("outputs/pangolin/*.csv")) {
    a = read_csv(i, col_types='ccnncnncccccc')
    pango = bind_rows(pango, a)
}

get_epi_id = Vectorize(function(x){strsplit(x,'\\|')[[1]][2]})

pango = mutate(pango, epi_id=get_epi_id(taxon))

p1 = select(m, epi_id="Accession ID", Location, gisaid_lineage=Lineage, collection_date="Collection date", Clade)
p2 = select(pango, taxon, epi_id, pango_lineage=lineage, conflict, ambiguity_score, scorpio_call, scorpio_support, scorpio_conflict, status)

p = full_join(p1, p2, by='epi_id')

filter(p, pango_lineage!="None", status=="passed_qc") %>% write_tsv(path='outputs/pangolin/filtered_metadata.tsv.gz')
p2 = filter(p, pango_lineage!="None", status=="passed_qc", !is.na(scorpio_call)) %>% mutate(k=gsub("_$","",gsub("[\\(\\)/ ]+","_",scorpio_call)))
p2 %>% select(taxon, k) %>% write_tsv(path='outputs/pangolin/lineages_map.tsv.gz')
p2 %>% select(k) %>% unique() %>% arrange(k) %>% write_tsv('outputs/pangolin/lineages_list.txt',col_names = F)
