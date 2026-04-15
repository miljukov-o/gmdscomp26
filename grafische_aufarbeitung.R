library(tidyverse)

df_list <- list.files(path = "./Datasets/Blinded_data_sets",
                      pattern="\\.csv$")

for (df in df_list) {

  rawdata <- readr::read_delim(paste0("./Datasets/Blinded_data_sets/", df), delim = ",") %>%
    mutate(W = factor(W, levels = c(0,1)))

  cairo_pdf(paste0("./img/", str_sub(df, end = -6), ".pdf"), height = 20, width = 20)
  GGally::ggpairs(rawdata,
                         mapping = ggplot2::aes(color = W))
  dev.off()
}
