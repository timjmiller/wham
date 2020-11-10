# Run all WHAM examples

# choose where to save output, otherwise will be saved in working directory
# write.dir <- "choose/where/to/save/output"
if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir)
setwd(write.dir)

library(wham)
wham.dir <- find.package("wham")

# Ex 1
ex1.dir <- file.path(write.dir,"ex1")
if(!dir.exists(ex1.dir)) dir.create(ex1.dir)
setwd(ex1.dir)
source(file.path(wham.dir, "example_scripts", "ex1_SNEMA_yellowtail_flounder.R"))

# Ex 2
ex2.dir <- file.path(write.dir,"ex2")
if(!dir.exists(ex2.dir)) dir.create(ex2.dir)
setwd(ex2.dir)
source(file.path(wham.dir, "example_scripts", "ex2_CPI_recruitment_SNEMA_yellowtail.R"))

# Ex 3
ex3.dir <- file.path(write.dir,"ex3")
if(!dir.exists(ex3.dir)) dir.create(ex3.dir)
setwd(ex3.dir)
source(file.path(wham.dir, "example_scripts", "ex3_projections.R"))

# Ex 4
ex4.dir <- file.path(write.dir,"ex4")
if(!dir.exists(ex4.dir)) dir.create(ex4.dir)
setwd(ex4.dir)
source(file.path(wham.dir, "example_scripts", "ex4_selectivity.R"))

# Ex 5
ex5.dir <- file.path(write.dir,"ex5")
if(!dir.exists(ex5.dir)) dir.create(ex5.dir)
setwd(ex5.dir)
source(file.path(wham.dir, "example_scripts", "ex5_M_GSI.R"))

# Ex 6
ex6.dir <- file.path(write.dir,"ex6")
if(!dir.exists(ex6.dir)) dir.create(ex6.dir)
setwd(ex6.dir)
source(file.path(wham.dir, "example_scripts", "ex6_NAA.R"))
