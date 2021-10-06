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
saveRDS(dge_result, "data/gse107361_dge.rds")

