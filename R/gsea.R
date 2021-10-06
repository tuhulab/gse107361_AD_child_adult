# library
pacman::p_load(stringr, tibble, dplyr, purrr, gprofiler2, tidyr, furrr, BiocParallel, fgsea)
plan(multisession, workers = 6)

# dge result
dge <-
  tibble(path = list.files("data/", pattern = "Adult|Children", full.names = T)) %>%
  mutate(file = map(path, ~readr::read_csv(.x, show_col_types = FALSE)),
         file = map(file, ~rename_with(.x, str_replace, .cols = contains("vs"),
                                      pattern = ".{1,}", replacement = "log2FoldChange")),
         file = map(file, ~select(.x, ENSEMBL, log2FoldChange, AveExpr, padj = adj.P.Val)),
         file = future_map(file, ~ group_by(.x, ENSEMBL) %>% nest %>%
                             mutate(data = map(data, ~.x %>% slice_max(AveExpr, n = 1)))
                           )) %>%
  unnest(cols = "file") %>%
  unnest(cols = "data")
ensembl_query <-
  dge %>% pull(ENSEMBL) %>% unique %>% gconvert() %>% as_tibble() %>%
  filter(!name %>% str_detect("\\d{1,}P$|\\d{1,}P\\d{1,}$|\\.|-AS\\d{1}|-DT"))
dge_gn <-
  dge %>% inner_join(ensembl_query %>% select(target, gene_name = name), by = c("ENSEMBL" = "target")) %>%
  mutate(contrast = path %>% str_remove_all("data/|\\.csv") %>% str_replace("-", "vs")) %>%
  group_by(contrast) %>% nest() %>%
  mutate(data = future_map(data, ~.x %>% group_by(gene_name) %>% nest %>%
                      mutate(data = map(data, ~.x %>% slice_max(AveExpr, n = 1))) %>%
                        ungroup))

dge_result <-
  dge_gn %>% mutate(project = "gse107361") %>%
  rename(dge_result = data) %>%
  select(project, contrast, dge_result)
saveRDS(dge_tibble, "data/dge.rds")

# gs
gs <-
  tibble(gs_name = c("CP:BIOCARTA",
                     "CP:KEGG",
                     "CP:PID",
                     "CP:REACTOME",
                     "CP:WIKIPATHWAYS",
                     "MIR:MIRDB",
                     "TFT:GTRD",
                     "GO:BP",
                     "GO:MF")) %>%
  mutate(gs_list = map(gs_name, function(subcat){
    gs_df <- msigdbr::msigdbr(species = "Homo sapiens",
                              subcategory = subcat) %>%
      group_by(gs_name) %>% tidyr::nest() %>%
      mutate(gene_id = purrr::map(data, ~ .x %>% pull(gene_symbol))) %>% select(-data)
    gs_list <- gs_df$gene_id
    names(gs_list) <- gs_df$gs_name
    return(gs_list)
  }))

# fgsea
gsea_exp <-
  expand.grid(contrast = dge_result$contrast,
              gs_name = gs$gs_name) %>%
  left_join(dge_result, by = "contrast") %>%
  left_join(gs, by = "gs_name") %>%
  mutate(rank = map2(dge_result, gs_name,
                     function(dge_result, gs_name){
                       p_cutoff <- ifelse(gs_name %in% c("MIR:MIRDB",
                                                         "TFT:GTRD"),
                                          1, # .05
                                          1)
                       log2fc_cutoff <- ifelse(gs_name %in% c("MIR:MIRDB",
                                                              "TFT:GTRD"),
                                               0, # 1
                                               0)
                       dge_result %>%
                         mutate(rank = log2FoldChange) %>%
                         filter(padj < p_cutoff, abs(log2FoldChange) > log2fc_cutoff) %>%
                         select(gene_name, rank) %>% deframe}),
         gse_res = map2(rank, gs_list,
                        function(rank, gs_list){
                          gsa_stat <- fgsea::fgsea(pathways = gs_list,
                                                   stats    = rank,
                                                   minSize  = 8,
                                                   maxSize  = 500,
                                                   eps = 0,
                                                   BPPARAM = SnowParam(6)) %>%
                            as_tibble()
                        }))

gsea_res <- gsea_exp %>% select(contrast, gs_name, project, gse_res) %>% unnest()
