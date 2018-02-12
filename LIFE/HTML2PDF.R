auc_max <- 0.94

results3.summary %>%
  filter(method %in% c("logit-simple", "PRS-max")) %>%
  mutate(method = forcats::fct_recode(method, "Our method" = "logit-simple",
                                      "State-of-the-art" = "PRS-max")) %>%
  filter(par.h2 == 0.8, par.dist == "laplace") %>%
  barplot_causal() +
  geom_hline(yintercept = auc_max, linetype = 3, color = "blue") +
  scale_fill_discrete() +
  scale_color_discrete()

ggsave("LIFE/main-AUC-logit.svg", scale = 1/90, width = 844, height = 572)

prefix <- "LIFE/pres"
path <- "../../useR-2017/decktape-1.0.0/"
system(glue::glue(
  "{path}phantomjs",
  " {path}decktape.js", 
  " {prefix}.html {prefix}.pdf", 
  " --load-pause 1000", 
  " -s 1504x1129"
))
