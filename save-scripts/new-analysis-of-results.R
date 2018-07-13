add_AUC <- function(x, pAUC) {
  x %>% 
    mutate(AUC  = map2_dbl(pred, pheno, bigstatsr::AUC),
           pAUC = map2_dbl(pred, pheno, pAUC))
}

add_sens_FDP <- function(x, corr) {
  
  alt <- map2(x$true_set, x$set, ~ {
    
    if (length(.y) == 0) return(c(0, 0))
    
    corr_set <- abs(corr[.x, .y, drop = FALSE])
    
    AltSens <- mean(apply(corr_set, 1, max, na.rm = TRUE))
    AltFDP <- 1 - mean(apply(corr_set, 2, max, na.rm = TRUE))
    c(AltSens, AltFDP)
  }) %>%
    transpose()
  
  x %>%
    mutate(AltSens = unlist(alt[[1]]),
           AltFDP  = unlist(alt[[2]]))
}

results4 <- list.files("results4", full.names = TRUE) %>%
  big_parallelize(., p.FUN = function(files, ind, add_AUC, pAUC, add_sens_FDP, corr) {
    
    library(tidyverse)
    library(Matrix)
    
    files[ind] %>%
      map_dfr(~readRDS(.x)) %>%
      as_tibble() %>%
      mutate(
        par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                            levels = c("30 in HLA", paste(3 * 10^(1:3), "in all")))) %>%
      add_AUC(pAUC) %>%
      add_sens_FDP(corr) %>%
      select(-c(num.simu, true_set, set, pheno, pred))
    
  }, p.combine = bind_rows, ind = seq_along(.), ncores = NCORES, 
  add_AUC = add_AUC, pAUC = pAUC, add_sens_FDP = add_sens_FDP, corr = corr)

print(results4, width = Inf)

####

results4 %>%
  filter(method %in% c("PRS-max", "logit-simple")) %>%
  group_by_at(c(vars(starts_with("par")), "method", "thr.r2")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC), N = n()) %>%
  print(n = Inf) -> results4.summary


c_method_r2 <- function(method, r2) {
  ifelse(is.na(r2), method, paste0(method, "-", r2))
}


barplot_causal <- function(results) {
  
  results %>%
    mutate(method2 = c_method_r2(method, thr.r2)) %>%
    myggplot(aes(par.causal, AUC_mean, fill = method2, color = method2)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_col(position = position_dodge(), alpha = 0.3) +
    geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, ymax = AUC_mean + 2 * AUC_boot),
                  position = position_dodge(width = 0.9), color = "black", width = 0.2) + 
    scale_y_continuous(limits = c(0.5, NA), minor_breaks = 0:20 / 20,
                       oob = scales::rescale_none) +
    labs(x = "Causal SNPs (number and location)", y = "Mean of 100 AUCs",
         fill = "Method", color = "Method")
}

barplot_causal_all(results4.summary)
# ggsave("new-results.pdf", scale = 1/90, width = 1235, height = 917)


####

results4.summary2 <- results4 %>%
  filter(method %in% c("PRS-max", "logit-simple")) %>%
  group_by_at(c(vars(starts_with("par")), "method", "thr.r2")) %>%
  summarise(AltSens_mean = mean(AltSens), AltFDP_mean = mean(AltFDP))


barplot_causal2 <- function(results, h2 = 0.8) {
  
  results %>%
    filter(par.h2 == h2) %>%
    mutate(method2 = c_method_r2(method, thr.r2)) %>%
    myggplot(aes(par.causal, value, fill = method2, color = method2)) +
    geom_col(position = position_dodge(), alpha = 0.3) +
    labs(x = "Causal SNPs (number and location)", y = "Mean of 100 Sensitivity",
         fill = "Method", color = "Method")
}

results4.summary2 %>% 
  gather("metric", "value", c("AltSens_mean", "AltFDP_mean")) %>% 
  barplot_causal2() +
  facet_grid(par.dist ~ metric) +
  scale_y_continuous(limits = c(0, 1))


####

results4 %>%
  filter(method == "logit-simple") %>%
  group_by_at(c(vars(starts_with("par")), "alpha")) %>%
  summarise(n = n()) %>%
  print(n = Inf)
