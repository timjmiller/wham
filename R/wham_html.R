#' Create HTML file to view output plots in browser
#'
#' Writes a set of HTML files with tabbed navigation between them. Called by
#' \code{\link{plot_wham_output}} if `out.type = 'html'` (default). Opens main file in
#' default browser. Modified from [`r4ss::SS_html`](https://github.com/r4ss/r4ss/blob/master/R/SS_html.R).
#'
#' @param dir.main directory to save html file (\code{\link{plot_wham_output}}
#' makes `.png` plot files in a `plots_png` subdirectory of `dir.main`).
#' @param title Title for HTML page.
#' @param width Width of plots (in pixels).
#' @param openfile Automatically open index.html in default browser?
#'
#' @seealso \code{\link{plot_wham_output}}, `r4ss::SS_html()`
#'
wham_html <- function(dir.main=NULL,
                    title="WHAM Output",
                    width=500,
                    openfile=TRUE){
  # get tabs (= subdirs) and tab labels
  tabs <- list.files(file.path(dir.main, "plots_png"))
  tab.labels <- tabs
  tab.labels[tabs=="diagnostics"] = "Diagnostics"
  tab.labels[tabs=="input_data"] = "Input Data"
  tab.labels[tabs=="misc"] = "Misc"
  tab.labels[tabs=="ref_points"] = "Reference Points"
  tab.labels[tabs=="results"] = "Results"
  tab.labels[tabs=="retro"] = "Retrospective"

  htmlhome <- file.path(dir.main, "wham_output.html")
  htmldir <- file.path(dir.main, "html")
  dir.create(htmldir, showWarnings = FALSE)
  for(t in 0:length(tabs)){
    if(t == 0){
      tab <- "Home"
      htmlfile <- htmlhome
    } else {
      tab <- tabs[t]
      htmlfile <- file.path(htmldir, paste0(tab, ".html"))
    }

    # write HTML head including some CSS stuff about fonts and whatnot
    # source for text below is http://unraveled.com/publications/css_tabs/
    cat('<html><head><title>', title, '</title>\n',
        '    <!-- source for text below is http://unraveled.com/publications/css_tabs/ -->\n',
        '    <!-- CSS Tabs is licensed under Creative Commons Attribution 3.0 - http://creativecommons.org/licenses/by/3.0/ -->\n',
        '    \n',
        '    <style type="text/css">\n',
        '    \n',
        '    body {\n',
        '    font: 100% verdana, arial, sans-serif;\n',
        '    background-color: #fff;\n',
        '    margin: 50px;\n',
        '    }\n',
        '    \n',

        #### this stuff allows scrolling while leaving the tabs in place,
        #### but I'd like to not have to set the height
        ## .container{
        ## }
        ## .panel{
        ## height: 1000px;
        ## overflow: auto;
        ## }

        '    /* begin css tabs */\n',
        '    \n',
        '    ul#tabnav { /* general settings */\n',
        '    text-align: left; /* set to left, right or center */\n',
        '    margin: 1em 0 1em 0; /* set margins as desired */\n',
        '    font: bold 11px verdana, arial, sans-serif; /* set font as desired */\n',
        '    border-bottom: 1px solid #6c6; /* set border COLOR as desired */\n',
        '    list-style-type: none;\n',
        '    padding: 3px 10px 2px 10px; /* THIRD number must change with respect to padding-top (X) below */\n',
        '    }\n',
        '    \n',
        '    ul#tabnav li { /* do not change */\n',
        '    display: inline;\n',
        '    }\n',
        '    \n',
        '    body#tab1 li.tab1, body#tab2 li.tab2, body#tab3 li.tab3, body#tab4 li.tab4 { /* settings for selected tab */\n',
        '    border-bottom: 1px solid #fff; /* set border color to page background color */\n',
        '    background-color: #fff; /* set background color to match above border color */\n',
        '    }\n',
        '    \n',
        '    body#tab1 li.tab1 a, body#tab2 li.tab2 a, body#tab3 li.tab3 a, body#tab4 li.tab4 a { /* settings for selected tab link */\n',
        '    background-color: #fff; /* set selected tab background color as desired */\n',
        '    color: #000; /* set selected tab link color as desired */\n',
        '    position: relative;\n',
        '    top: 1px;\n',
        '    padding-top: 4px; /* must change with respect to padding (X) above and below */\n',
        '    }\n',
        '    \n',
        '    ul#tabnav li a { /* settings for all tab links */\n',
        '    padding: 2px 4px; /* set padding (tab size) as desired; FIRST number must change with respect to padding-top (X) above */\n',
        '    border: 1px solid #6c6; /* set border COLOR as desired; usually matches border color specified in #tabnav */\n',
        '    background-color: #cfc; /* set unselected tab background color as desired */\n',
        '    color: #666; /* set unselected tab link color as desired */\n',
        '    margin-right: 0px; /* set additional spacing between tabs as desired */\n',
        '    text-decoration: none;\n',
        '    border-bottom: none;\n',
        '    }\n',
        '    \n',
        '    ul#tabnav a:hover { /* settings for hover effect */\n',
        '    background: #fff; /* set desired hover color */\n',
        '    }\n',
        '    \n',
        '    /* end css tabs */\n',
        '    \n',
        '    \n',
        '    h2 {\n',
        '    font-size: 20px;\n',
        '    color: #4c994c;\n',
        '    padding-top: 1px;\n',
        '    font-weight: bold;\n',
        '    border-bottom-width: 1px;\n',
        '    border-bottom-style: solid;\n',
        '    border-bottom-color: #6c6;\n',
        '    padding-bottom: 2px;\n',
        '    padding-left: 0px;\n',
        '    }\n',
        '    </style>',
        '</head>\n',
        sep = "", file=htmlfile, append=FALSE)

    # write navigation menu

    #### more stuff related to scroll options
    ## <div class="main">
    ##   <div class="container">
    cat('<!-- Site navigation menu -->\n',
        '  <ul id="tabnav">\n',
        file=htmlfile, append=TRUE)
    for(tt in 0:length(tabs)){
      if(tt == 0){
        cat('    <li class="tab1"><a href=',htmlhome,'>Home</a></li>\n',sep="",
            file=htmlfile, append=TRUE)
      } else {
        cat('    <li class="tab',tt+1,'"><a href=',file.path(htmldir, paste0(tabs[tt], ".html")),'>',tab.labels[tt],'</a></li>\n',sep="",
            file=htmlfile, append=TRUE)
      }
    }
    cat('  </ul>\n', file=htmlfile, append=TRUE)

    # add model summary png on "Home" page
    if(tab=="Home"){
      cat('\n\n<h2><a name=', tab, '>', tab, '</h2>\n', sep="",
          file=htmlfile, append=TRUE)
      cat("<p align=left><a href='", file.path(dir.main,"plots_png","diagnostics","summary_text.png"),
          "'><img src='", file.path(dir.main,"plots_png","diagnostics","summary_text.png"),
          "' border=0 width=", width*1.7, "></a><br>",
          "<i><small>file: <a href='", file.path(dir.main,"plots_png","diagnostics","summary_text.png"),
          "'>", file.path(dir.main,"plots_png","diagnostics","summary_text.png"), "</a></small></i>\n",
          sep="",  file=htmlfile,  append=TRUE)
    } else { # add plots to non-Home pages
      cat('\n\n<h2><a name=', tabs[t], '>', tab.labels[t], '</h2>\n', sep="",
          file=htmlfile, append=TRUE)
      plots <- list.files(file.path(dir.main, "plots_png", tab), full.names=TRUE)
      for(i in 1:length(plots)){
        cat("<p align=left><a href=", plots[i],
            "><img src=", plots[i],
            " border=0 width=", width, "></a><br>",
            "<i><small>file: <a href=", plots[i],
            ">", plots[i], "</a></small></i>\n",
            sep="",  file=htmlfile,  append=TRUE)
      }
    }

    # end html file
    cat("\n\n</body>\n</html>", file=htmlfile, append=TRUE)
  }

  # open HTML file automatically
  if(openfile){
    cat("Opening HTML file in your default web-browser.\n")
    browseURL(htmlhome)
  }
}
