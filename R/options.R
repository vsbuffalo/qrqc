qrqc_options <- list(
  templates=list(tex=system.file('extdata', 'knitr-report-template.Rnw', package='qrqc'),
                 html=system.file('extdata', 'knitr-report-template.Rhtml', package='qrqc')),
  tex_opts=list(portrait=list(fig.width=8, fig.height=6), landscape=list(fig.width=10, fig.height=7)),
  html_opts=list(fig.width=8, fig.height=5, dpi=100)
  )
