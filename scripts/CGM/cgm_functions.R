
# ---------------------------------------------------------------------------------------------
# CLASS definitions: --------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

Timedata <- R6Class(
  "Timedata", lock_objects = FALSE, 
  public = list(
    name = NULL, raw = NULL, flagged = NULL, cnames = NULL, coded = NULL, 
    melted = NULL, comps = NULL, status = NULL, msg = TRUE, ind_prog = TRUE, 
    
    initialize = function(name, raw, isos, pad_height, pad_cluster, msg, ind_prog) {
      self$name <- name
      self$msg <- msg
      self$ind_prog <- ind_prog
      if (msg) {self$start()}
      self$raw <- raw
      self$coded <- codeIsolates(raw, name, isos, pad_height, pad_cluster)
      self$melted <- meltedIDs(raw, name, pad_height, pad_cluster)
      invisible(self)
    }, 
    start = function() {
      cat(paste0("  Constructing ", toupper(self$name), " data object:\n"))
    }, 
    coded_status = function(nov_code) {
      self$status <- self$coded %>% mutate(status = ifelse(isolate %in% nov_code, "novs", NA))
    }, 
    set_comps = function() {
      self$comps <- self$coded %>% 
        set_colnames(gsub(self$name, "tp", colnames(.))) %>% 
        compsSet(., toupper(self$name), indicate_progress = self$ind_prog)
      invisible(self)
    }, 
    set_cnames = function() {
      self$cnames <- self$coded %>% pull(grep("h", colnames(self$coded), value = TRUE)) %>% unique() %>% sort()
      invisible(self)
    }, 
    flag_clusters = function() {
      self$flagged <- flaggingClusters(self$comps, self$name)
      invisible(self)
    }
  )
)

Heightdata <- R6Class(
  "Heightdata", lock_objects = FALSE, 
  public = list(
    h_before = NULL, h_after = NULL, comps = NULL, changed = tibble(), same = tibble(), 
    tracked = NULL, bef = NULL, aft = tibble(), results = NULL, 
    
    initialize = function(starter, t1_comps, hvals) {
      self$h_after <- self$h_before <- starter
      self$comps <- t1_comps %>% filter(tp1_h == starter) %>% arrange(tp1_h, tp1_cl)
      self$bef <- t1_comps %>% filter(tp1_h == self$h_before) %>% 
        set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
      self$results <- vector(mode = "list", length = length(hvals)) %>% set_names(hvals)
      invisible(self)
    }, 
    clust_tracking = function(t2_comps, t2_cnames, t1_coded, t2_coded, ind_prog) {
      self$changed <- trackClusters(self$comps, t2_comps, t2_cnames, 
                                    t1_coded, t2_coded, ind_prog)
      invisible(self)
    }, 
    post_data = function(t1_comps) {
      self$aft <- t1_comps %>% filter(tp1_h == self$h_after) %>% 
        set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))
      invisible(self)
    }, 
    unchanged = function() {
      self$same <- noChange(self$aft, self$bef, self$tracked)
    }, 
    update_iteration = function() {
      self$tracked <- bind_rows(self$changed, self$same) %>% arrange(tp1_h, tp1_cl)
      self$results[[self$h_after]] <- self$tracked
      invisible(self)
    }, 
    reset_values = function() {
      self$h_before <- self$h_after
      self$h_after <- NULL
      
      self$bef <- self$aft %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
      self$aft <- tibble()
    }
  )
)


# ---------------------------------------------------------------------------------------------
# FORMATTING functions: -----------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

addingType <- function(dfx) {
  df <- dfx %>% add_column(type = NA)
  # type 1: TP1 > 2, TP2 > 2, TP1 == TP2
  inds1 <- which(df$tp1_cl_size > 2 & df$tp2_cl_size > 2 & df$tp1_cl_size == df$tp2_cl_size)
  df$type[inds1] <- "Type1"
  # type 2: TP1 > 2, TP2 > 2, TP2 > TP1
  inds2 <- which(df$tp1_cl_size > 2 & df$tp2_cl_size > 2 & df$tp2_cl_size > df$tp1_cl_size)
  df$type[inds2] <- "Type2"
  # type 3: TP1 < 3, TP2 > 2
  inds3 <- which(df$tp1_cl_size < 3 & df$tp2_cl_size > 2)
  df$type[inds3] <- "Type3"
  # type 4: TP1 < 3, TP2 < 3
  inds4 <- which(df$tp1_cl_size < 3 & df$tp2_cl_size < 3)
  df$type[inds4] <- "Type4"
  assert("No clusters with unassigned type", !any(is.na(df$type)))
  return(df)
}

# Identifying the first and last time each TPX cluster was seen in TPX - flags (for TP1 and TP2 individually)
flaggingClusters <- function(tp_comps, tpx) {
  t_fal <- tp_comps %>% group_by(composition) %>% slice(1, n()) %>% 
    add_column(type = rep(c("first", "last"), nrow(.)/2)) %>% dplyr::ungroup()
  
  t_ff <- t_fal %>% filter(type == "first") %>% select(grep("id", colnames(.)), composition)
  colnames(t_ff) <- c(paste0("first_", tpx, "_flag"), "composition")
  
  t_lf <- t_fal %>% filter(type == "last") %>% select(grep("id", colnames(.)), composition)
  colnames(t_lf) <- c(paste0("last_", tpx, "_flag"), "composition")
  
  full_join(t_ff, t_lf, by = "composition") %>% 
    left_join(tp_comps, ., by = "composition") %>% 
    select(-composition) %>% return()
}

# Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"
timeTaken <- function(pt, sw) {
  if (is.null(sw[["start_time"]]) & is.null(sw[["end_time"]])) {
    paste0("Neither start nor end time were collected") %>% return()
  }else if (is.null(sw[["end_time"]]) & !is.null(sw[["start_time"]])) {
    paste0("End time was not collected.") %>% return()
  }else if (!is.null(sw[["end_time"]]) & is.null(sw[["start_time"]])) {
    paste0("Start time was not collected.") %>% return()
  }else {
    z <- difftime(sw[['end_time']], sw[['start_time']], units = "secs") %>% as.double()
    m <- 60
    h <- m^2
    
    if (z >= h) {
      hrs <- trunc(z/h)
      mins <- trunc(z/m - hrs*m)
      paste0("\nThe ", pt, " process took ", hrs, " hour(s), ", mins, " minute(s), and ", 
             round(z - hrs*h - mins*m), " second(s).") %>% return()
    }else if (z < h & z >= m) {
      mins <- trunc(z/m)
      paste0("\nThe ", pt, " process took ", mins, " minute(s) and ", round(z - mins*m), " second(s).") %>% return()
    }else {
      paste0("\nThe ", pt, " process took ", round(z), " second(s).") %>% return()
    }  
  }
}

# Outputs the same message in two ways, one is directed to standard output and one to a log file
outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}

# Given a dataframe df, two column names c1 and c2 (height and cluster respectively) and a new
# ID prefix tpx (e.g. "tp1"), creates an ID column and adds to df before returning df
newID <- function(df, tpx, c1, c2, ph, pc) {
  newh <- df %>% pull(c1) %>% as.character() %>% as.integer() %>% 
    formatC(., width = max(3, ph), format = "d", flag = "0") %>% paste0("h", .)
  newc <- df %>% pull(c2) %>% as.character() %>% as.integer() %>% 
    formatC(., width = max(3, pc), format = "d", flag = "0") %>% paste0("c", .)
  df %>% add_column(id = paste0(toupper(tpx), "_", newh, "_", newc)) %>% return()
}


# Given the defining filename, read in the data (need the full path from your working directory), 
# indicate to user if file is not found
readBaseData <- function(filename, file_number, delimiter) {
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    read.table(file = filename, stringsAsFactors = FALSE, check.names = FALSE, 
               header = TRUE, sep = delimiter, allowEscapes = TRUE, 
               fileEncoding = checkEncoding(filename)) %>% as_tibble() %>% return()
  }
}

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

# Given a raw time point dataset, the timepoint ID (e.g. "tp1"), and a list of all isolates 
# present at TP1 and TP2, return a dataframe with isolates given numeric codes (e.g. "-1-", "-10-")
codeIsolates <- function(df, tpx, all_iso, ph, pc) {
  hx <- paste0(tpx, "_h")
  cx <- paste0(tpx, "_cl")
  idx <- paste0(tpx, "_id")
  
  df %>% right_join(all_iso, ., by = c("char_isolate" = "isolate")) %>% 
    select(-char_isolate) %>% 
    rename(isolate = num_isolate) %>% as.data.table() %>% 
    melt.data.table(id = "isolate") %>% 
    set_colnames(c("isolate", hx, cx)) %>% 
    factorToInt(., hx) %>% 
    newID(., tpx, hx, cx, ph, pc) %>% 
    mutate(isolate = paste0("-", isolate, "-")) %>% 
    rename_with(., ~gsub("id", idx, .x)) %>% return()
}

meltData <- function(dataset, id_val) {
  melt(dataset, id = id_val) %>% as_tibble() %>% return()
}

# Given a dataframe df, converts all the elements of column c1 from factors to integers
factorToInt <- function(df, c1) {
  df %>% mutate(across(all_of(c1), as.character)) %>% 
    mutate(across(all_of(c1), as.integer)) %>% return()
}

# For padding height and cluster columns with h0..0.., and c0..0.., respectively
padCol <- function(cvals, padval, padchr) {
  ifelse(!is.na(cvals), formatC(cvals, width = padval, format = "d", flag = "0") %>% 
           paste0(padchr, .), NA) %>% return()
}

# Returns a table like the following: 
# A tibble: 7,965 x 4
#   isolate     tp1_h tp1_cl tp1_id       
#   <chr>       <chr>  <dbl> <chr>        
# 1 isolate1    0          4 TP1_h000_c004
# 2 isolate10   0         76 TP1_h000_c076
# 3 isolate100  0         28 TP1_h000_c028
# ...
meltedIDs <- function(df, k, ph, pc) {
  cnames <- paste0(k, c("", "_h", "_cl", "_id"))
  df %>% as.data.table() %>% 
    melt.data.table(id = "isolate") %>% as_tibble() %>% 
    set_colnames(c("isolate", cnames[2:3])) %>% 
    mutate(across(cnames[2], as.character)) %>% 
    newID(., cnames[1], cnames[2], cnames[3], ph, pc) %>% 
    set_colnames(c("isolate", cnames[2:4])) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# Given a TP dataset, we output a set where the composition of each cluster is indicated with a string where 
# each string element represents a single isolate (easier to compare cluster composition later on)
# --------------------------------------------------------------------------------------------------------------
compsSet <- function(tp_coded, tp, indicate_progress) {
  # sort the data and assign general labels
  cb <- tp_coded %>% mutate(num_iso = as.integer(gsub("-", "", isolate)), 
                            composition = isolate, tp_cl_size = isolate) %>% arrange(tp_h, tp_cl, num_iso)
  allheights <- cb$tp_h %>% unique()
  
  # This part returns a data frame of CL ID | CL COMP | CL SIZE, for all heights and clusters at the given TP
  if (indicate_progress) {pb <- txtProgressBar(min = 0, max = length(allheights), initial = 0, style = 3)}
  
  tmp <- lapply(1:length(allheights), function(i) {
    if (indicate_progress) {setTxtProgressBar(pb, i)}
    
    # for height i, we first collect the data from tp_coded, then arrange by isolate
    x <- cb %>% filter(tp_h == allheights[i])
    
    # if there are multiple clusters at a given height:
    if (length(unique(x$tp_id)) > 1) {
      # output a table where each cluster has its isolate composition to the right: iso 1, iso 2, iso 3, ...
      a2 <- aggregate(composition ~ tp_id, data = x, FUN = toString)
      a3 <- aggregate(tp_cl_size ~ tp_id, data = x, FUN = length)
      left_join(a2, a3, by = "tp_id") %>% as_tibble() %>% return()
    }else {
      # if at a height there is only one cluster, we handle it differently
      tibble(tp_id = unique(x$tp_id), composition = paste0(x$isolate, collapse=","), 
             tp_cl_size = length(x$isolate)) %>% return()
    }
  }) %>% bind_rows()
  
  if (indicate_progress) {close(pb)}
  
  # This part binds the height and cluster columns with the cluster composition table
  tpcomps <- tp_coded %>% select(-isolate) %>% unique() %>% 
    left_join(tmp, ., by = "tp_id") %>% select(tp_h, tp_cl, tp_id, composition, tp_cl_size) %>% 
    set_colnames(gsub("tp", tolower(tp), colnames(.))) %>% return()
}


# ---------------------------------------------------------------------------------------------
# TRACKING functions: -------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# Identifying TP1 clusters at height b that did not change from height a --> they inherit the metrics, 
# since the values did not change
# --------------------------------------------------------------------------------------------------------------
noChange <- function(ca, cb, hb) {
  # these are clusters whose composition did not change from cb to ca (height x to height x+1)
  in_both_heights <- inner_join(ca, cb, by = "comp") %>% select(grep("aft", colnames(.)), id_bef)
  
  stayed_the_same <- left_join(in_both_heights, hb, by = c("id_bef"="tp1_id")) %>% 
    select(-grep("tp1", colnames(.)), -id_bef) %>% 
    rename(tp1_h = h_aft, tp1_cl = cl_aft, tp1_cl_size = size_aft, tp1_id = id_aft)
  
  assert("This dataset only has clusters that did not change", !any(is.na(stayed_the_same$tp2_h)))
  return(stayed_the_same)
}

# --------------------------------------------------------------------------------------------------------------
# When using grep doesn't return the needed results (several possible reasons), we track the individual 
# strains in a cluster (to find out which TP2 clusters contain at least these strains)
# --------------------------------------------------------------------------------------------------------------
checkEachIsolate <- function(cluster_i, t2_coded, t2_comps) {
  # no matching TP2 cluster found - the order of the cluster composition is muddling this up
  isolates <- strsplit(cluster_i$composition, split = ",|, ") %>% unlist()
  
  ids <- t2_coded %>% filter(isolate %in% isolates) %>% 
    group_by(tp2_id) %>% summarise(size = n(), .groups = "drop") %>% 
    filter(size == length(isolates)) %>% pull(tp2_id)
  
  df2 <- t2_comps %>% filter(tp2_id %in% ids) %>% select(-composition)
  cluster_i %>% select(-composition) %>% bind_cols(., df2) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# For a set of TP1 singletons, we track them to the TP2 clusters they're found in
# --------------------------------------------------------------------------------------------------------------
trackSingletons <- function(singles, t1_coded, t2_coded, composition_set) {
  isos_in_tp1 <- t1_coded %>% filter(tp1_id %in% singles$tp1_id)
  
  a4 <- t2_coded %>% filter(isolate %in% isos_in_tp1$isolate) %>% 
    left_join(isos_in_tp1, ., by = "isolate") %>% select(-isolate) %>% 
    left_join(., singles, by = "tp1_id")
  
  left_join(a4, composition_set, by = "tp2_id") %>% 
    select(tp1_h, tp1_cl, tp1_id, tp1_cl_size, tp2_h, tp2_cl, tp2_id, tp2_cl_size, num_novs) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# Tracking TP1 clusters to their TP2 equivalents, using grep primarily, with other / backup methods for 
# TP1 singletons and multistrain cases where the first method fails or is not comprehensive enough
# --------------------------------------------------------------------------------------------------------------
# When using grep: if looking for 1, will return 1, 10, 214, etc. So we sandwich with hyphens: -1-,-10-,...
trackClusters <- function(hdata, t2_comps, t2names, t1_coded, t2_coded, indicate_progress) {
  
  t1set <- t2set <- tibble()
  singletons <- hdata %>% filter(tp1_cl_size == 1) %>% select(tp1_id, tp1_cl_size)
  # at least one singleton cluster
  if (nrow(singletons) > 0) {
    t1set <- t2_comps %>% select(tp2_id, tp2_cl_size, num_novs) %>% 
      trackSingletons(singletons, t1_coded, t2_coded, .)
  }
  
  multistrain <- hdata %>% filter(tp1_cl_size > 1)
  if (nrow(multistrain) > 0) { # at least one cluster with size larger than 1
    if (indicate_progress) {tc <- txtProgressBar(min = 0, max = nrow(multistrain), initial = 0, style = 3)}
    
    t2set <- lapply(1:nrow(multistrain), function(i) {
      if (indicate_progress) {setTxtProgressBar(tc, i)}
      results_i <- tibble()
      cluster_i <- multistrain[i,]
      
      if (nchar(cluster_i$composition) < 2000) {
        inds <- grep(cluster_i$composition, t2_comps$composition)
        results_i <- t2_comps[inds,] %>% select(-composition) %>% 
          bind_cols(cluster_i, .) %>% select(-composition)
        
        if ((nrow(results_i) == 0) | (!last(t2names) %in% results_i$tp1_h)) {
          # (1) no matching TP2 cluster found - the order of the cluster composition is muddling this up
          # (2) if the process is being cut off partway - not tracking to all heights (e.g. TP2 height 1634 not found)
          results_i <- checkEachIsolate(cluster_i, t2_coded, t2_comps)
        }
      }else { 
        # (1) too many isolates in the cluster for effective grep use
        results_i <- checkEachIsolate(cluster_i, t2_coded, t2_comps)
      }
      results_i %>% arrange(tp2_h, tp2_cl) %>% return()
    }) %>% bind_rows()
    if (indicate_progress) {close(tc)}
    assert("\nNo NA values\n", all(!is.na(t2set)))
  }
  bind_rows(t1set, t2set) %>% arrange(tp1_h, tp1_cl) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# We introduce the metric columns, once we have collected all the tracking information for the TP1 clusters
# Include: cluster size change, growth rate (using actual sizes), growth rate (using number of novels), 
# number of sneakers, number of novels
# --------------------------------------------------------------------------------------------------------------
oneHeight <- function(df) {
  df %>% 
    mutate(actual_size_change = tp2_cl_size - tp1_cl_size, 
           actual_growth_rate = ((tp2_cl_size - tp1_cl_size) / tp1_cl_size) %>% round(., digits = 3), 
           new_growth = (tp2_cl_size / (tp2_cl_size - num_novs)) %>% round(., digits = 3)) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% return()
}

findingSneakers <- function(novels, q1, q2, matched) {
  compmatches <- matched[!duplicated(matched$tp1_id),]
  
  did_not_chg <- compmatches %>% filter(tp1_cl_size == tp2_cl_size) %>% add_column(add_TP1 = 0)
  chg <- compmatches %>% filter(tp1_cl_size != tp2_cl_size)
  
  if (nrow(chg) > 0) {
    # identifying the number of additional TP1 strains (sneakers) that show up in the TP2 cluster 
    # that each TP1 cluster was first tracked to
    sneakers <- lapply(1:nrow(chg), function(j) {
      kc <- chg %>% slice(j) # key cluster j
      tbl1 <- q1 %>% filter(tp1_id == kc$tp1_id) %>% mutate(status = "tp1_cl_size")
      tbl2 <- q2 %>% filter(tp2_id == kc$tp2_id)
      # number of novels in the TP2 cluster it kc was tracked to
      x2 <- tbl2 %>% filter(status == "novs") %>% nrow()
      # number of sneakers in the TP2 cluster it kc was tracked to
      x3 <- tbl2 %>% filter(!(isolate %in% tbl1$isolate) & is.na(status)) %>% 
        mutate(status = "additional_TP1") %>% nrow()
      
      c2_tally <- tibble(tp1_cl_size = nrow(tbl1), num_novs = x2, add_TP1 = x3, 
                         tp2_cl_size = sum(nrow(tbl1), x2, x3))
      
      assert(paste0("oneHeight(): Novels check failed for ", kc$tp1_id), kc$num_novs == c2_tally$num_novs)
      assert(paste0("Composition not calculated properly for ", kc$tp2_id, 
                    " when tracking ", kc$tp1_id), c2_tally$tp2_cl_size == kc$tp2_cl_size)
      
      left_join(kc, c2_tally, by = c("tp1_cl_size", "tp2_cl_size", "num_novs")) %>% return()
    }) %>% bind_rows()
    
    results <- bind_rows(did_not_chg, sneakers)
  }else {
    results <- did_not_chg
  }
  
  results %>% return()
}