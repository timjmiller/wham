# Run all WHAM examples

# choose where to save output, otherwise will be saved in working directory
# main.dir <- "choose/where/to/save/output"
if(!exists("main.dir")) main.dir = getwd()
if(!dir.exists(main.dir)) dir.create(main.dir)

library(wham)
wham.dir <- find.package("wham")

# Ex 1
write.dir <- file.path(main.dir,"ex1")
source(file.path(wham.dir, "example_scripts", "ex1_SNEMA_yellowtail_flounder.R"))

# Ex 2
write.dir <- file.path(main.dir,"ex2")
source(file.path(wham.dir, "example_scripts", "ex2_CPI_recruitment_SNEMA_yellowtail.R"))

# Ex 3
write.dir <- file.path(main.dir,"ex3")
source(file.path(wham.dir, "example_scripts", "ex3_projections.R"))

# Ex 4
write.dir <- file.path(main.dir,"ex4")
source(file.path(wham.dir, "example_scripts", "ex4_selectivity.R"))

# Ex 5
write.dir <- file.path(main.dir,"ex5")
source(file.path(wham.dir, "example_scripts", "ex5_M_GSI.R"))

# Ex 6
write.dir <- file.path(main.dir,"ex6")
source(file.path(wham.dir, "example_scripts", "ex6_NAA.R"))
