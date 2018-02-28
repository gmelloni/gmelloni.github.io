# Helper function library
# SHIVA trial retrospective analysis

library(stringr)


# paste wihtout space
`%+%` <- function(a,b) paste0(a,b)

# paste wiht a space in between elements
`%++%` <- function(a,b) paste0(c(a,b), collapse = " ")


################################################################################
#' Convert Amminoacid format
#'
#' This is a custom made helper function that allows to convert AA alteration 
#' position from one format to another. More specifically it is design to recognise
#' the format of the publication (e.g. "Glu454Lys") and to convert it into a 
#' format suitable for the R package CancerPanelSimulator. 
#'
#' @param myvector , a vector of AA positions in the table t3 format from 
#' appending 2
#'
#' @return
#' @export
#'
#' @examples
#' convertAA("Gly454Lys")
convertAA <- function(myString){
  
  #AA table
  aatable <- list( Ala="A"
                  , Ile="I"
                  , Leu="L"
                  , Val="V"
                  , Phe="F"
                  , Trp="W"
                  , Tyr="Y"
                  , Asn="N"
                  , Cys="C"
                  , Gln="Q"
                  , Met="M"
                  , Ser="S"
                  , Thr="T"
                  , Asp="D"
                  , Glu="E"
                  , Arg="R"
                  , His="H"
                  , Lys="K"
                  , Gly="G"
                  , Pro="P"
  )
  
  #readin string
  splitted <- stringr::str_match(myString, "([:alpha:]{3})([:digit:]+)([:alpha:]{3})")
  #check that it is processed
  if (length(splitted) !=4){stop("Error, the input was not correct")}
  if(any(is.na(splitted))){stop("Error, the input was not correct")}
  
  #conver using the list
  partA <- aatable[[splitted[2]]]
  partB <- splitted[3]
  partC <- aatable[[splitted[4]]]
  
  if(is.null(partA) || is.null(partB) || is.null(partC)){
    stop("Error, the AA variant provided is not properly formatted")
  }
  
  #convert the vector
  converted <- paste(c(partA, partB, partC), collapse="")
  
  return(converted)
}

################################################################################
# Simple function to merge an entire list of DF with a common key
mergeThemAll <- function(biglist , by=NULL){
  if(is.null(by)) {
    out <- merge(biglist[[1]] , biglist[[2]])
    for( i in 3:length(biglist)){
      out <- merge(out , biglist[[i]])
    }
    return(out)
  } else {
    out <- merge(biglist[[1]] , biglist[[2]] , by=by)
    for( i in 3:length(biglist)){
      out <- merge(out , biglist[[i]] , by=by)
    }
    return(out)
  }
}

############################################################################
# Another AA converter
# Write helper function to replace 3 chr amminoacid form to 1 chr amminoacid code.
hgvs_long_to_short <- function(x) {
  
  # Amminoacid convertion rules
  # -------------------------------------
  amino <- data.frame(matrix(c(
    "aa3", "aa1",
    "Lys","K", 
    "Asn","N", 
    "Thr","T",
    "Arg","R",
    "Ser","S",
    "Ile","I",
    "Met","M",
    "Gln","Q",
    "His","H",
    "Pro","P",
    "Leu","L",
    "Glu","E",
    "Asp","D",
    "Ala","A",
    "Gly","G",
    "Val","V",
    "Tyr","Y",
    "Cys","C",
    "Trp","W",
    "Phe","F",
    "Ter","*"), ncol=2, byrow = TRUE))
  
  # Convert  3chr amminoacid standard format to 1 chr amminoacid code.
  # -------------------------------------
  for (i in 1:nrow(amino)) {
    x <- gsub(amino[i,1],amino[i,2],x)
  }
  x
}
# extend dplyr
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


## Prepare function for plotting Detection power with Confidenti intervals
# ------------------------------------------------------------------------------
simulate_detectionPower_tbl<- function(df
                                    , img_file = "../Figures/fig4A.svg"
                                    , fill_color = MYPALETTE(10)[2]
                                    , gtitle = "Overall Detection power"){
  
  # Convert to dataframe
  df <- data.frame(matrix(unlist(df)
                          , nrow=length(df)
                          , byrow=T )
                   , stringsAsFactors = FALSE)
  
  # GET STATS
  # ------------------------------------------------------------------------------
  # Calculate for each n° of alteration, the mean
  detectionPower <- apply(df,2,mean)
  
  # Calculate lower CI for bootstrap
  leftCI <- apply(df,2,function(x){left <- quantile(x, 0.025)})
  
  # Calculate upper CI for bootstrap
  rightCI <- apply(df,2,function(x){right <- quantile(x, 0.975)})
  
  # Preapre for column with # alterations
  n_alterations <- 1:length(detectionPower)
  
  # Put them all together in a data frame
  output <- data.frame(n_alterations, detectionPower, leftCI, rightCI)
  
  return(output)
}
  
  
  


## Prepare function for plotting Detection power with Confidenti intervals
# ------------------------------------------------------------------------------
simulate_detectionPower_img <- function(df
                                    , img_file = "../Figures/fig4A.svg"
                                    , fill_color = MYPALETTE(10)[2]
                                    , gtitle = "Overall Detection power"){
  
  # Convert to dataframe
  df <- data.frame(matrix(unlist(df)
                          , nrow=length(df)
                          , byrow=T )
                   , stringsAsFactors = FALSE)
  
  # GET STATS
  # ------------------------------------------------------------------------------
  # Calculate for each n° of alteration, the mean
  detectionPower <- apply(df,2,mean)
  
  # Calculate lower CI for bootstrap
  leftCI <- apply(df,2,function(x){left <- quantile(x, 0.025)})
  
  # Calculate upper CI for bootstrap
  rightCI <- apply(df,2,function(x){right <- quantile(x, 0.975)})
  
  # Preapre for column with # alterations
  n_alterations <- 1:length(detectionPower)
  
  # Put them all together in a data frame
  dtpow <- data.frame(n_alterations, detectionPower, leftCI, rightCI)
  
  # Generate the plot
  image <- dtpow %>%
    ggplot(aes(x=n_alterations, y=detectionPower)) + 
    geom_bar(stat="identity" , position=position_dodge(), fill=fill_color) + 
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(breaks = n_alterations, trans = "reverse") +
    geom_errorbar(aes(ymax=rightCI
                      , ymin=leftCI)
                  , position=position_dodge(0.9)
    ) + 
    ggtitle(gtitle) +
    annotate("text"
             , x=n_alterations
             , y=detectionPower+0.15
             , label= sapply(detectionPower*100
                             , function(x) paste(c(round(x, digits=2), "%")
                                                 , collapse=" ")
             )
    ) +
    coord_flip()
  
  #Export as a svg
  #ggsave(filename=img_file, plot=image, device = "svg")
  #message("file saved in" %++% img_file)
  
  # Show image
  return(image)
}


# Helper functions
# ------------------------
# formula to calculate CI
leftCIprop <- function(p, sample_size){
  se <- sqrt( p*(1 - p)/sample_size)
  left <- p - 1.96*se - 0.5/sample_size
  if (left < 0) left <- 0
  return(left)
}

rightCIprop <- function(p, sample_size){
  se <- sqrt( p*(1 - p)/sample_size)
  right <- p + 1.96*se + 0.5/sample_size
  if (right > 1) right <- 1
  return(right)
}


# function that generates statistics information from bootstrap
get_bootstrap_stats <- function(df){
  
  if(!is.data.frame(df)){stop("input should be a dataframe")}
  if(nrow(df) == 0){stop("input should not be empty")}
  
  # Calculate for each n° of alteration, the mean
  freq <- apply(df,2,mean)
  
  # Calculate lower CI for bootstrap
  leftCI <- apply(df,2,function(x){left <- quantile(x, 0.025)})
  
  # Calculate upper CI for bootstrap
  rightCI <- apply(df,2,function(x){right <- quantile(x, 0.975)})
  
  # Preapre for column with # of raws
  row_counter <- 1:ncol(df)
  
  # Put them all together in a data frame
  output <- data.frame(row_counter, freq, leftCI, rightCI)
  
  return(output)
}

# Data table extension
datatable_withbuttons <- function(df, ...){
  
  DT::datatable(df
                , extensions = 'Buttons'
                , rownames=FALSE
                , options = list(
                  dom = 'lBfrtip'
                  , buttons = c('copy', 'csv', 'excel')
                  
                )
                , ...
  )
}
