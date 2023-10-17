make_html_figs_tables_fn <- function(od, do.html = TRUE, do.tex = FALSE) {

  wham.dir <- find.package("wham")
  pt = list.files(find.package("wham"), pattern = "figures_tables.Rmd", recursive = T, full.names = T)[1]
  #pt = list.files(find.package("wham"), pattern = "figures_tables_temp.Rmd", recursive = T, full.names = T)[1]
  file.copy(from=pt, to=od, overwrite=TRUE)
  #print(dir(od))
  
  if(do.html) try(rmarkdown::render(file.path(od,"figures_tables.Rmd"), output_format = "html_document", output_file = file.path(od, "wham_figures_tables.html"), 
    quiet = T, envir = new.env()))
  # if(do.html) rmarkdown::render(file.path(od,"figures_tables_temp.Rmd"), output_format = "html_document", output_file = file.path(od, "wham_figures_tables.html"), 
  #   quiet = T, envir = new.env())
  #if(do.tex) rmarkdown::render(file.path(od,"par_tables.Rmd"), output_format = "pdf_document", output_file = file.path(od,"wham_par_tables.pdf"), quiet = T)
  if(do.tex) { #for some reason on windows working outside of the temp directory was causing issues for tinytex::latexmf.
    origdir = getwd()
    setwd(od)
    try(rmarkdown::render("figures_tables.Rmd", output_format = "pdf_document", output_file = file.path(od, "wham_figures_tables.pdf"), quiet = T, envir = new.env()))
    setwd(origdir)
  }

}