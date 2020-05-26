# function to read TargetLynx data output (.txt format) and put it in a tidy format.
# 
get.tl.data    <- function( ) {
  
  require(data.table)
  require(stringr)
  
  # get raw data and meta data file names
  raw.files <- dir(paste0(getwd(),"/TL_data_raw"),pattern = ".txt") 
  # meta.file <- dir(paste0(pd,"/data"),pattern = "^meta_data")
  
  # parse all data to list of  data.table
  dat <- sapply( raw.files, function(fn) { 
    
    # get raw text file
    ql.data.raw <- read.csv(paste0(getwd(),"/TL_data_raw/",fn), sep = "\t", header = FALSE, stringsAsFactors = F)
    
    # Extract metabolite names, metabolite number
    target.string <-  "Compound\\s[0-9]{1,}:\\s{1,}"
    analytes <- str_replace(ql.data.raw$V1[grepl(target.string,ql.data.raw$V1)],target.string,"")
    
    # Tidy up QuanLynx data
    ql.data <- ql.data.raw[,-1]
    resaved <- ifelse(ql.data.raw[2,1] == "", T, F)
    if (resaved == T){
      #For "re-saved quanlynx text files
      n.samples <- as.integer(mean(diff(grep("Compound",ql.data.raw$V1)[-1] ) ) - 4)
      c.names   <- make.names(as.character(ql.data[7, ]))
      r.select  <- as.vector(sapply(seq_along(analytes), function(i) (8 + (i - 1) * (n.samples + 4) ) : (3 + i * (n.samples + 4)  ) ) )
    } else {
      #For untouched quanlynx text files
      n.samples <- as.integer(mean(diff(grep("Compound",ql.data.raw$V1)[-1] ) ) - 2)
      c.names   <- make.names(as.character(ql.data[4, ]))
      r.select  <- as.vector(sapply(seq_along(analytes), function(i) (5 + (i - 1) * (n.samples + 2) ) : (2 + i * (n.samples + 2)  ) ) )
    }
    
    # Set Column names
    ql.data <- ql.data[r.select, ]
    colnames(ql.data)      <- c.names
    colnames(ql.data)[1:2] <- c("Order","File.Name")
    
    # Set data types
    ql.data$Order              <- as.integer(ql.data$Order)
    ql.data$RT                 <- as.numeric(ql.data$RT)
    ql.data$Area               <- as.numeric(ql.data$Area)
    ql.data$Peak.Start.Height  <- as.numeric(ql.data$Peak.Start.Height)
    ql.data$Height             <- as.numeric(ql.data$Height)
    ql.data$S.N                <- as.numeric(ql.data$S.N)
    ql.data$Analyte            <- as.vector(sapply(analytes, function(m)  rep(m,n.samples) ) )
    ql.data$TL_file            <- fn
    
    
    return(ql.data)
    
  },USE.NAMES = F , simplify = F)
  
  # Bind all data sets 
  dat <- do.call(rbind, dat)
  
  # write to csv file
  write.csv(dat,paste0(getwd(),"/TL_data_tidy/tidy_TL_data.csv"))
  
  return(dat)
  
}
