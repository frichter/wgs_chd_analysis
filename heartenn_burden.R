# Felix Richter
# 7/28/2018
# HeartENN burden tests
##################################

setwd("/Users/felixrichter/Dropbox/PhD/")
options(stringsAsFactors=FALSE)

p = c("magrittr", "extraDistr", "purrr", "dplyr", "ggplot2", "tidyr", "stringi", "readr")
lapply(p, require, character.only = TRUE)

dn_max_score = read_tsv("max_score_hi_qual_dnvs_2018_07_29.txt")


############################################################
# Testing for max score burdens
############################################################

CalculateHeartENNBurden = function(min_diff_score, dn, n_cases_muts, n_ctrl_muts) {
  print(min_diff_score)
  dh_burden = dn %>% 
    filter(max_per_row >= min_diff_score) %>%
    group_by(case_ctrl) %>% 
    summarise(n = n()) %>% 
    ungroup %>% 
    spread(key = case_ctrl, value = n) %>% 
    mutate(
      fisher.p = fisher.test(cbind("Case" = c(Case, n_cases_muts - Case),
                                   "Ctrl" = c(Ctrl, n_ctrl_muts - Ctrl)))$p.value,
      ci_lo = fisher.test(cbind("Case" = c(Case, n_cases_muts - Case),
                                "Ctrl" = c(Ctrl, n_ctrl_muts - Ctrl)))$conf.int[[1]],
      ci_hi = fisher.test(cbind("Case" = c(Case, n_cases_muts - Case),
                                "Ctrl" = c(Ctrl, n_ctrl_muts - Ctrl)))$conf.int[[2]],
      or = (Case/(n_cases_muts-Case))/(Ctrl/(n_ctrl_muts-Ctrl)),
      min_diff_score = min_diff_score)
  return(dh_burden)
}

n_cases_muts = dn_max_score %>% filter(case_ctrl == "Case") %>% nrow
n_ctrl_muts = dn_max_score %>% filter(case_ctrl == "Ctrl") %>% nrow

min_sig_score_list = c(0.01, seq(0.05, 0.3, 0.05))

top_var_enrichment = map_df(min_sig_score_list, CalculateHeartENNBurden, 
                            dn_max_score, n_cases_muts, n_ctrl_muts)


############################################################
# Plot enrichment
############################################################

p = top_var_enrichment %>% 
  mutate(sig = ifelse(fisher.p <= 0.05, "p<0.05", "NS")) %>% 
  ggplot(aes(x = min_diff_score, y = or, ymax = ci_hi, ymin = ci_lo, color = sig)) +
  geom_hline(yintercept = 1, color = "blue") +
  geom_pointrange() +
  scale_color_manual(values = c("black", "red")) +
  xlab("DNV cut-off for minimum difference score") +
  theme_classic()
p
ggsave("case_ctrl_heartenn_burden.png", p, width = 3.75, height = 2)


