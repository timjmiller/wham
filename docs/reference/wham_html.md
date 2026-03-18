# Create HTML file to view output plots in browser

Writes a set of HTML files with tabbed navigation between them. Called
by
[`plot_wham_output`](https://timjmiller.github.io/wham/reference/plot_wham_output.md)
if \`out.type = 'html'\` (default). Opens main file in default browser.
Modified from
\[\`r4ss::SS_html\`\](https://github.com/r4ss/r4ss/blob/master/R/SS_html.R).

## Usage

``` r
wham_html(dir.main = NULL, title = "WHAM Output", width = 500, openfile = TRUE)
```

## Arguments

- dir.main:

  directory to save html file
  ([`plot_wham_output`](https://timjmiller.github.io/wham/reference/plot_wham_output.md)
  makes \`.png\` plot files in a \`plots_png\` subdirectory of
  \`dir.main\`).

- title:

  Title for HTML page.

- width:

  Width of plots (in pixels).

- openfile:

  Automatically open index.html in default browser?

## See also

[`plot_wham_output`](https://timjmiller.github.io/wham/reference/plot_wham_output.md),
\`r4ss::SS_html()\`
