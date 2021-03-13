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
                                  check.names = FALSE, quote = "", stringsAsFactors = FALSE)
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