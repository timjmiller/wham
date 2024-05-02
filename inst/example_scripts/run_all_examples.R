# Run all WHAM examples

# choose where to save output, otherwise will be saved in working directory
# main.dir <- "choose/where/to/save/output"
if(!exists("main.dir")) main.dir = getwd()
if(!dir.exists(main.dir)) dir.create(main.dir)

library(wham)
wham.dir <- find.package("wham")

# Ex 1
write.dir <- file.path(main.dir,"ex1")
source(file.path(wham.dir, "example_scripts", "ex1_basics.R"))

# Ex 2
write.dir <- file.path(main.dir,"ex2")
source(file.path(wham.dir, "example_scripts", "ex2_CPI_recruitment.R"))

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

# Ex 8
write.dir <- file.path(main.dir,"ex8")
source(file.path(wham.dir, "example_scripts", "ex8_compare_asap.R"))

# Ex 9
write.dir <- file.path(main.dir,"ex9")
source(file.path(wham.dir, "example_scripts", "ex9_retro_pred.R"))

# Ex 10
write.dir <- file.path(main.dir,"ex10")
source(file.path(wham.dir, "example_scripts", "ex10_simulation.R"))

# Ex 11
write.dir <- file.path(main.dir,"ex11")
source(file.path(wham.dir, "example_scripts", "ex11_catchability.R"))
