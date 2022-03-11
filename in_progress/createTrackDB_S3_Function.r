#Author: Chao-Jen Wong 
#Modified by Jenny Smith 
#5/10/19 

#Purpose: create a track db for a custom track hub on UCSC genome browser, using BAM(bigwig) file inputs. 


# Example inputs
# bucket <- 'fh-pi-meshinchi-s'
# prefix <- "SR/SMRTSeq/Isoseq3/demultiplexed/"
# region <- 'us-west-2'
# trackName <- "Isoseq3"
# shortLabel <- 'Isoseq3'
# longLabel <- "PacBio Isoseq3 for AML and NBM samples"

function (bucket, 
          prefix,
          trackName,
          pattern = "all_movies.+\\.bam$",
          TrackFileName = "hg38.trackDb.txt", #must be the same one as in the genomes.txt file 
          region='us-west-2',
          url = "https://s3-us-west-2.amazonaws.com",
          concatenate=FALSE,
          col = c(102, 194, 165), 
          shortLabel = NULL, longLabel = NULL) {
  
  require(dplyr)
  require(aws.s3)
  require(aws.signature)
  
  #Set your credentials using ~/.aws/credentials file.
  use_credentials()
  #set your AWS region 
  Sys.setenv("AWS_DEFAULT_REGION" = region)
  #Define type of file for track
  type <- "bam"
  
  if (!bucket_exists(bucket)) 
    stop(bucket, " does not exists.")
  trackName <- gsub("[[:space:]]", "-", trackName)
  if (is.null(shortLabel)) 
    shortLabel <- trackName
  if (is.null(longLabel)) 
    longLabel <- trackName
  message("trackName: ", trackName)
  message("shortLabel: ", shortLabel)
  message("longLabel: ", longLabel)
  
  #files in bucket and prefix to be used, filtered for file type desired by user. 
  fdf <- get_bucket_df(bucket = bucket, prefix=prefix) %>% 
    filter(grepl(pattern, Key))
  
  if (nrow(fdf) == 0) 
    stop("There are no file names in that match the pattern: ",pattern,"in the", bucket)
  urls <- paste(url,bucket,fdf[["Key"]], sep="/")
  
  #Function to create a hubtrack 
  HubTracks <- function(name, col, urls, shortLabel, longLabel) {
    cat("\n")
    cat("track ", name, "\n")
    cat("compositeTrack on\n")
    cat("shortLabel ", shortLabel, "\n")
    cat("longLabel ", longLabel, "\n")
    cat("allButteronPair on\n")
    cat("dragAndDrop on\n")
    cat(paste("type",type,"0 1.0\n"))
    cat("alwaysZero on\n")
    cat("color ", col, "\n")
    cat("\n")
    for (i in 1:length(urls)) {
      trackname <- sub(paste0(".",type), "", basename(urls[i]))
      cat("track ", trackname, "\n")
      cat(paste("type",type,"\n"))
      cat("shortLabel ", trackname, "\n")
      cat("longLabel ", trackname, " reads", "\n")
      cat("parent ", name, "\n")
      cat("bigDataUrl ", urls[i], "\n")
      cat("\n")
    }
  }
  
  #create the trackdb file 
  tFile <- TrackFileName
  message("Creating ", TrackFileName)
  file.create(TrackFileName)
  sink(TrackFileName)
  HubTracks(name = trackName, col, urls, shortLabel, longLabel)
  sink()
  
  #Concatenate to the original TrackDb file if needed. 
  if(concatenate){
    trackDb <- get_object(object=paste0(prefix,TrackFileName),
                          bucket = bucket,  as="parsed")
    trackDb <- strsplit(trackDb,split = "\n")[[1]][-c(1:9)] #remove the header to avoid duplicates
    
    for(line in trackDb){
      cat(paste0(line,"\n"), file= TrackFileName, append = TRUE)
    }
  }
  
  #upload to the bucket 
  put_object(file=TrackFileName, object = paste0(prefix,TrackFileName), bucket = bucket, 
             acl = 'public-read', verbose=TRUE)
  file.remove(TrackFileName)
}

