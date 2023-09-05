library(dplyr)

here::i_am("analyze/paml_parse.R")

# PARSE FOR LRT CALCULATIONS AND PARAMETER REPORTING

mlc <- readr::read_file(here::here("paml/mlc")) %>%
  stringr::str_split_1(stringr::fixed("NSsites Model ")) %>%
  tibble::enframe() %>%
  filter(substr(value, 1, 1) != "C") %>%
  mutate(name = substr(value, 1, 1)) %>%
  rename(NSSitesModel = name, RawText = value)

# function for AIC from lnL

getAIC <- function(lnL, nParam) {return(-2*(lnL)+2*nParam)}

# Model-null model pairs and their comparisons are manually entered
lrt <- mlc %>%
  mutate(RawLRT = stringr::str_extract(RawText, "lnL(.*)")) %>% # . = any character but new line, * any number of them
  mutate(nParam = stringr::str_extract(RawLRT, "(?<=np:)[:digit:]+")) %>% # (?!...) is a look a head assertion
  mutate(lnL = as.numeric(stringr::str_extract(RawLRT, "[:punct:][:digit:]+[:punct:][:digit:]+"))) %>%
  dplyr::select(!starts_with("Raw")) %>%
  mutate(nParam = as.numeric(nParam)) %>%
  tidyr::pivot_wider(names_from = NSSitesModel, 
                     names_glue = "M{NSSitesModel}{.value}",
                     values_from = c("nParam", "lnL")) %>%
  mutate(M2ChiSq = M2lnL - M1lnL,
         M8ChiSq = M8lnL - M7lnL,
         M3ChiSq = M3lnL - M0lnL,
         M1ChiSq = M1lnL - M0lnL) %>%
  mutate(M2df = M2nParam - M1nParam,
         M8df = M8nParam - M7nParam,
         M3df = M3nParam - M0nParam,
         M1df = M1nParam - M0nParam) %>%
  mutate(M2PValue = pchisq(q = M2ChiSq, df = M2df, lower.tail = FALSE),
         M8PValue = pchisq(q = M8ChiSq, df = M8df, lower.tail = FALSE),
         M3PValue = pchisq(q = M3ChiSq, df = M3df, lower.tail = FALSE),
         M1PValue = pchisq(q = M1ChiSq, df = M1df, lower.tail = FALSE)) %>%
  mutate(M2AIC = getAIC(M2lnL, M2nParam),
         M8AIC = getAIC(M8lnL, M8nParam),
         M7AIC = getAIC(M2lnL, M7nParam),
         M3AIC = getAIC(M3lnL, M3nParam),
         M1AIC = getAIC(M1lnL, M1nParam),
         M0AIC = getAIC(M0lnL, M0nParam)) %>%
  tidyr::pivot_longer(cols = starts_with("M"),
                      names_to = c("NSSitesModel","Stat"),
                      names_pattern = "M(.)(.*)") %>%
  tidyr::pivot_wider(names_from = "Stat",
                     values_from = "value") %>%
  arrange(AIC) # lowest AIC is model 3

params_raw <- mlc %>%
  rowwise() %>%
  mutate(RawParams = stringr::str_extract_all(RawText, "(?s)Detailed output identifying parameters(.*?)dN & dS for each branch")) #(?s) if there are new lines, need this

params <- params_raw %>%
  mutate(Param = stringr::str_extract_all(RawParams, "\\S+\\s+=\\s+[:digit:][:punct:][:digit:]+")) %>% # \\s = whitespace, \\s= non-whitespace, + = one or more
  tidyr::unnest(Param) %>%
  tidyr::separate(col = Param, into = c("Parameter", "Value"), sep = "\\s+=\\s+") %>%
  mutate(Value = as.numeric(Value)) %>%
  dplyr::select(!starts_with("Raw")) %>%
  mutate(Parameter = case_when(
    Parameter == "(dN/dS)" ~ "w",
    Parameter == "(ts/tv)" ~ "K",
    TRUE ~ Parameter
  ))

mle <- params_raw %>%
  rowwise() %>%
  mutate(MLE = stringr::str_extract_all(RawParams, "[:alpha:]:\\s+(.*)")) %>%
  tidyr::unnest(MLE) %>%  
  tidyr::separate(col = MLE, into = c("Parameter", "Value"), sep = ":\\s+") %>%
  mutate(Value = stringr::str_split(Value, "\\s+")) %>%
  tidyr::unnest(Value) %>%
  mutate(Value = as.numeric(Value)) %>%
  group_by(NSSitesModel, Parameter) %>%
  mutate(SiteClass = case_when(
    Parameter == "p" ~ row_number(rev(Value)) - 1,
    Parameter == "w" ~ row_number(Value) - 1,)) %>%
  dplyr::select(!starts_with("Raw")) %>%
  tidyr::unite(Parameter, Parameter, SiteClass)
  
params <- bind_rows(params, mle) %>%
  arrange(NSSitesModel)
  
mlc <- full_join(lrt, params, by = "NSSitesModel", multiple = "all")

readr::write_csv(mlc, here::here("analyze/paml_site_mlc.csv"))

# PARSE BEB RESULTS

rst <- readr::read_file(here::here("paml/rst"))

omega <- rst %>%
  stringr::str_extract_all(stringr::regex("(?<=Bayes Empirical Bayes \\(BEB\\) probabilities for )(.*?)Positively selected sites(.*?\r\n\r\n){2}(.*?)(?=\r\n((\r\n)|\\z))", dotall = TRUE, multiline = TRUE)) %>% #\\z is end of the input
  tibble::enframe() %>%
  tidyr::unnest(value) %>%
  mutate(NSSitesModel = case_when(substr(value, 1, 1) == 3 ~ 2,
                                  substr(value, 1, 2) == 11 ~ 8)) %>%
  mutate(value = stringr::str_split(value, "\r\n\r\nPositively selected sites\r\n\r\n")) %>%
  tidyr::unnest(value)
omega

omega_beb <- omega %>% 
  filter(!stringr::str_detect(value, "^\t")) %>%
  mutate(BEB = stringr::str_split(value, "\r\n")) %>%
  select(NSSitesModel, BEB) %>%
  tidyr::unnest(BEB) %>%
  filter(stringr::str_detect(BEB, "^\\s+")) %>%
  tidyr::separate_wider_delim(cols = BEB, 
                              names = c("Site", "BEB"), 
                              delim = stringr::regex("\\s.\\s{3}")) %>%
  mutate(Site = as.numeric(Site)) %>%
  tidyr::separate_wider_delim(cols = BEB, 
                              names = c("BEB", "w"), 
                              delim = stringr::regex("\\)\\s+")) %>%
  tidyr::separate_wider_delim(cols = w, 
                              names = c("w", "SE"), 
                              delim = stringr::regex("\\s\\+\\-\\s+")) %>%
  dplyr::select(NSSitesModel, Site, w, SE) %>%
  mutate(Site = Site + 31,
         w = as.numeric(w),
         SE = as.numeric(SE))

omega_pos_sel <- omega %>%
  filter(stringr::str_detect(value, "^\t")) %>%
  mutate(PosteriorProb = stringr::str_split(value, "\r\n")) %>%
  select(NSSitesModel, PosteriorProb) %>%
  tidyr::unnest(PosteriorProb) %>%
  filter(stringr::str_detect(PosteriorProb, "^\\s{2,}")) %>%
  tidyr::separate_wider_delim(cols = PosteriorProb, 
                              names = c("Site", "PosteriorProb"), 
                              delim = stringr::regex("(?<=[:digit:]{1,3})\\s[[:alpha:]-]\\s+")) %>%
  mutate(Site = as.numeric(Site) + 31) %>%
  mutate(PosteriorProb = stringr::str_extract(PosteriorProb, "[:digit:].[:digit:]+")) # get first match

omega <- full_join(omega_beb, omega_pos_sel)

readr::write_csv(omega, here::here("analyze/paml_site_omega.csv"))

# PARSE ANCESTRAL RECONSTRUCTION

rst <- readr::read_file(here::here("paml_site_reconstruction/rst"))

arc_rows <- rst %>%
  stringr::str_extract_all(stringr::regex("Marginal reconstruction of ancestral sequences(.*?)Summary of changes along branches", dotall = TRUE, multiline = TRUE)) %>%
  unlist() %>%
  tibble::enframe() %>%
  tidyr::unnest(value) %>%
  # filter(name == 1) %>% 
  pull(value) %>%
  stringr::str_split_1("\r\n") %>%
  tibble::enframe() %>%
  filter(stringr::str_detect(value, "^\\s+[:digit:]+")) %>%
  dplyr::select(!name)

arc_split <- arc_rows %>%
  tidyr::separate_wider_delim(cols = value,
                              names = c("Site", "ARC"),
                              delim = stringr::regex("(?<=^\\s{1,3}[:digit:]{1,3})\\s+")) %>%
  mutate(Site = as.numeric(Site) + 31) %>%
  tidyr::separate_wider_delim(cols = ARC,
                              names = c("Freq", "ARC"),
                              delim = stringr::regex("\\s{2,}(?=[[:alpha:]-]{3})"),
                              too_many = "merge") %>%
  tidyr::separate_wider_delim(cols = ARC,
                               names = c("Tip", "Node"),
                               delim = stringr::regex("\\s:\\s+"))

arc_tips <- arc_split %>%
  rowwise() %>%
  mutate(Tip = stringr::str_split(Tip, stringr::regex("(?<=\\))\\s"))) %>%
  mutate(nTip = length(Tip)) %>%
  mutate(Tip = list(purrr::set_names(Tip, 1:nTip))) %>%
  tidyr::unnest_longer(col = c(Tip),
                       values_to = "Tip",
                       indices_to = "ID") %>%
  mutate(Codon = stringr::str_extract(Tip, stringr::regex("(.*?)(?=\\s)"))) %>%
  mutate(AA = stringr::str_extract(Tip, stringr::regex("(?<=\\()(.)(?=\\))"))) %>%
  dplyr::select(!c(Tip, nTip, Node)) %>%
  mutate(Type = "Tip", .after = ID)
  
arc_nodes <- arc_split %>%
  rowwise() %>%
  mutate(Tip = stringr::str_split(Tip, stringr::regex("(?<=\\))\\s"))) %>%
  mutate(nTip = length(Tip)) %>%
  mutate(Node = stringr::str_extract_all(Node, stringr::regex("(?<=\\()(.{7})(?=\\))"))) %>%
  mutate(nNode = length(Node)) %>%
  mutate(Node = list(purrr::set_names(Node, nTip+(1:nNode)))) %>%
  tidyr::unnest_longer(col = c(Node),
                       values_to = "Node",
                       indices_to = "ID") %>%
  tidyr::separate(Node,
                  into = c("AA", "PosteriorProb"),
                  sep = " ") %>%
  mutate(Type = "Node") %>%
  dplyr::select(!c(Tip, nTip, nNode))

arc <- bind_rows(arc_tips, arc_nodes) %>%
  mutate(ID = as.numeric(ID),
         PosteriorProb = as.numeric(PosteriorProb)) %>%
  relocate(c(ID, Type), .before = Site)

readr::write_csv(arc, here::here("analyze/paml_site_arc.csv"))
