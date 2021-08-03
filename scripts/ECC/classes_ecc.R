processedStrains <- function(base_strains) {
  loc_cols <- intersect(c("Country", "Province", "City"), colnames(base_strains)) %>% sort()
  
  if (length(loc_cols) > 0) {
    strain_data <- base_strains %>% 
      mutate(Date     = as.Date(paste(Year, Month, Day, sep = "-")), 
             Location = do.call(paste, c(base_strains[loc_cols], sep = "_")))
    
    assignments <- strain_data %>% select(Date, Latitude, Longitude, Location) %>% 
      unique() %>% rownames_to_column("dr") %>% as.data.table()
    dr_matches <- left_join(strain_data, assignments, 
                            by = c("Latitude", "Longitude", "Date", "Location")) %>% 
      select(Strain, dr)
  }else {
    strain_data <- base_strains %>% mutate(Date = as.Date(paste(Year, Month, Day, sep = "-")))
    assignments <- strain_data %>% select(Date, Latitude, Longitude) %>% 
      unique() %>% rownames_to_column("dr") %>% as.data.table()
    dr_matches <- left_join(strain_data, assignments, 
                            by = c("Latitude", "Longitude", "Date")) %>% select(Strain, dr)
  }
  
  list("strain_data" = strain_data, 
       "assignments" = assignments, 
       "dr_matches" = dr_matches) %>% return()
}

checkEncoding <- function(fp) {
  readr::guess_encoding(fp) %>% arrange(-confidence) %>% slice(1) %>% pull(encoding) %>% return()
}

Timepoint <- R6Class(
  "Timepoint", lock_objects = FALSE, 
  public = list(
    filepath = NULL, name = NULL, filedata = NULL, proc = NULL, height_list = NULL, 
    initialize = function(fpath, name) {
      self$filepath <- fpath
      self$name <- toupper(name)
      self$readTyping()
      invisible(self)
    }, 
    readTyping = function() {
      self$filedata <- read.table(self$filepath, header = TRUE, sep = "\t", row.names = 1, 
                                  check.names = FALSE, quote = "", stringsAsFactors = FALSE, 
                                  fileEncoding = checkEncoding(self$filepath))
    }, 
    Process = function(hx) {
      self$proc <- self$filedata %>% 
        select(hx$h) %>% 
        set_colnames(paste0(self$name, "_T", colnames(.))) %>% 
        rownames_to_column("Strain") %>% as_tibble() %>% 
        add_column(tpx = 1)
      colnames(self$proc)[colnames(self$proc) == "tpx"] <- self$name
      invisible(self)
    }, 
    listHeights = function(hx) {
      self$height_list <- lapply(1:nrow(hx), function(i) {
        self$filedata[,hx$h[i],drop=FALSE] %>% set_colnames(hx$th[i])
      }) %>% set_names(paste0(self$name, "_", hx$th))
      invisible(self)
    }
  )
)