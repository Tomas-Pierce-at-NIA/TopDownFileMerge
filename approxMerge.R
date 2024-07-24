
# devtools helps with a couple of tasks
library(devtools)
# must install sqldf from source to work around a compatibility issue
devtools::install_github("ggrothendieck/sqldf")

# sqldf enables the use of SQL to operate on data loaded into R,
# which enables the use of operations that would otherwise be more awkward
library(sqldf)

# readxl is a set of tools related to reading and interpreting excel files using R
library(readxl)

# name of folder containing data
FOLDER <- "T:/TGB/LSS/TGU/Users/Mozhgan/TopPic_analysis/TOMAS/TopPIC"
setwd(FOLDER)

# get the absolute file path of data files by string concatention of folder and file
# seperated by the directory seperator character
# which is always forward slash `/` in the context of R
senescence.filename <- paste(FOLDER, 
                             "Combined_Senescent.xlsx", 
                             sep="/")

quiescent.filename <- paste(FOLDER,
                            "Combined_Quiescent.xlsx",
                            sep="/")

proliferating.filename <- paste(FOLDER,
                                "Combined_Proliferating.xlsx",
                                sep="/")

# read senescent protein sheet into memory
# making a temporary copy of data from the T:\ drive to our system's RAM allows R
# to operate on the data.
# we specify `sheet` explicitly so that if an excel file with multiple sheets
# gets used, we will still get data from the correct sheet
senescence.tibble <- read_excel(path=senescence.filename,
                                   sheet="CombinedSEN_ms2_toppic_proteofo")

# convert senescence data from the `tibble` type to the `data.frame` type.
# `tibble` refers to a more recent way of representing data in R,
# `data.frame` refers to an older way of representing data in R,
# so what can be achieved differs between the two representations,
# sometimes in non-obvious ways.
# In particular, the same function may treat a data.frame or a tibble differently.
# Checking what data representations a function or module expects
# is often a wise idea.
# `data.table` is the more commonly used type, so we default to that generally.
senescence.datatable <- as.data.frame(senescence.tibble)

# read quiescent protein sheet into memory
quiescent.tibble <- read_excel(path=quiescent.filename,
                               sheet="CombinedQSEN_ms2_toppic_proteof")

# convert quiescence data from the `tibble` type to the `data.frame` type
quiescent.datatable <- as.data.frame(quiescent.tibble)

# read proliferating protein sheet into memory
proliferating.tibble <- read_excel(path=proliferating.filename,
                                   sheet="CombinedNonSEN_ms2_toppic_prote")

# convert proliferating data from the `tibble` type to the `data.frame` type
proliferating.datatable <- as.data.frame(proliferating.tibble)

# because spaces in names tend to cause problems, we are going to 
# add new versions of retention time and proteoform mass that do not
# have spaces in each datatable

# instead of implementing these seperately, we are going to create a function
# which will always do this correctly for us

# note that R passes a copy of the data.frame to the function, not the original,
# so that we have to return the modified version and assign it to a name
# in order to make use of it.
addColsUnspaced <- function(dframe) {
  dframe[["retentionTime"]] = dframe[["Retention time"]]
  dframe[["proteoformMass"]] = dframe[["Proteoform mass"]]
  return(dframe)
}

# activity on the right hand side of assignment finishes first,
# then the results are assigned to the name on the left,
# so we can reuse names if we desire.
# it is good practice to only pair 1 name with 1 idea
# using a name called `count` for a total would be confusing
proliferating.datatable <- addColsUnspaced(proliferating.datatable)
senescence.datatable <- addColsUnspaced(senescence.datatable)
quiescent.datatable <- addColsUnspaced(quiescent.datatable)

# now that the setup is done, we can get on with the actual problem
# as stated:
# match each row from senescence table (S), quiescence table (Q), and 
# proliferating table (P) such that
# if a row in S has a mass within 1 dalton and a retention time within 30 s
# of a row in Q and a row in P, combine those rows into a wide row like
# row(S), row(Q), row(P)
# if on the other hand a row is present in S and a row exists in Q but not P
# satisfying these requirements,
# create a row like
# row(S), row(Q), Null
# and for a row present in S with no pairs in Q or P create a wide row like
# row(S), Null, Null
# and do the analogous for each row in Q and P

# step 1: make it easier by sorting on the variables of interest

# the exact order matters less than that the orderings are consistent.
# we use the order() function to sort by column values in R.

# side note: since . is a reserved SQL character, we 
# also take the opportunity to change names to use _,
# in order to avoid a problem
proliferating_sorted <- proliferating.datatable[order(
  proliferating.datatable$retentionTime, 
  proliferating.datatable$proteoformMass),]
senescence_sorted <- senescence.datatable[order(
  senescence.datatable$retentionTime, 
  senescence.datatable$proteoformMass),]
quiescent_sorted <- quiescent.datatable[order(
  quiescent.datatable$retentionTime, 
  quiescent.datatable$proteoformMass),]

# because we want to include cases where some items exist but not all
# in every table, the appropriate operation is a full outer join

results <- sqldf(
  "SELECT P.*, S.*, Q.*
  FROM proliferating_sorted AS P
  FULL OUTER JOIN senescence_sorted AS S
  ON abs(P.proteoformMass - S.proteoformMass) < 1 and abs(P.retentionTime - S.retentionTime) < 30
  FULL OUTER JOIN quiescent_sorted as Q
  ON abs(P.proteoformMass - Q.proteoformMass) < 1 and abs(P.retentionTime - Q.retentionTime) < 30
  "
)

