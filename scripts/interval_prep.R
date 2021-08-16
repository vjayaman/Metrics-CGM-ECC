
if (params$intervaltype == "weekly") {
  interval <- "Week"
  interval_list <- metadata$Week %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, all_of(interval))) %>% 
    inner_join(., f2, by = c("Strain" = "isolate")) %>% 
    set_colnames(c("isolate", "ivl", "heightx"))
  
}else if (params$intervaltype == "monthly") {
  interval <- "YearMonth"
  interval_list <- metadata %>% pull(interval) %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, all_of(interval))) %>% 
    inner_join(., f2, by = c("Strain" = "isolate")) %>% 
    set_colnames(c("isolate", "ivl", "heightx"))
  
  
}else if (params$intervaltype == "multiset") {
  # if (interactive()) {
  #   divider <- paste0("Enter date(s) between ", min(metadata$Date), " and ", 
  #                     max(metadata$Date), " to split set into two timepoints ", 
  #                     "\n(YYYY-MM-DD format, separated by commas): ") %>% 
  #     readline() %>% as.character()
  # }else {
  #   cat(paste0("Enter date(s) between ", min(metadata$Date), " and ", 
  #              max(metadata$Date), " to split set into two timepoints ", 
  #              "\n(YYYY-MM-DD format, separated by commas): "))
  #   divider <- readLines("stdin", n = 1) %>% as.character()
  # }
  divider <- readLines("scripts/date.txt")[2]
  setdivider <- strsplit(divider, split = ",") %>% unlist() %>% as.Date(., format = "%Y-%m-%d")
  
  assertion1 <- lapply(setdivider, function(x) nchar(as.character(x)) == 10) %>% unlist()
  assert("Provided input(s) has/have 10 characters", all(assertion1))
  
  assert("Provided input(s) is/are date type, correct format", !any(is.na(setdivider)))
  
  assertion3 <- lapply(setdivider, function(x) x >= min(metadata$Date) & x <= max(metadata$Date)) %>% unlist()
  assert("Provided input(s) is/are between min date and max date", all(assertion3))
  
  sdiv <- c(min(metadata$Date), setdivider, max(metadata$Date))
  fdivs <- tibble(start = sdiv[1:(length(sdiv)-1)], end = sdiv[2:length(sdiv)]) %>% 
    mutate(ivl = paste0("set", 1:nrow(.)))
  
  metadata <- metadata %>% add_column(Multiset = "")
  for (i in 1:nrow(fdivs)) {
    metadata[metadata$Date >= fdivs$start[i] & metadata$Date < fdivs$end[i]]$Multiset <- fdivs$ivl[i]
  }
  metadata[metadata$Date == fdivs$end[i]]$Multiset <- fdivs$ivl[i]
  
  interval <- "Multiset"
  interval_list <- metadata %>% pull(interval) %>% unique() %>% sort()
  
  interval_clusters <- metadata %>% select(c(Strain, Multiset)) %>%
    inner_join(., f2, by = c("Strain" = "isolate")) %>%
    set_colnames(c("isolate", "ivl", "heightx"))
}

