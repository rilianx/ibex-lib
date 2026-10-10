#!/bin/bash
# Version para Kindle (pagina 10.8 x 14.4 cm, 10 pt, tablas ajustadas al ancho):
# genera main-kindle.tex desde main.tex y compila main-kindle.pdf.
cd "$(dirname "$0")"
sed -e 's|^\\usepackage\[margin=2.6cm\]{geometry}|\\usepackage[paperwidth=10.8cm,paperheight=14.4cm,margin=3mm,top=4mm,bottom=6mm,footskip=4mm]{geometry}\n\\usepackage{adjustbox}\\usepackage{etoolbox}\n\\BeforeBeginEnvironment{tabular}{\\begin{adjustbox}{max width=\\linewidth}}\\AfterEndEnvironment{tabular}{\\end{adjustbox}}\n\\sloppy\\emergencystretch=3em|' \
    -e 's|^\\documentclass\[11pt\]{article}|\\documentclass[10pt]{article}|' main.tex > main-kindle.tex
pdflatex -interaction=nonstopmode main-kindle.tex >/dev/null
bibtex main-kindle >/dev/null
pdflatex -interaction=nonstopmode main-kindle.tex >/dev/null
pdflatex -interaction=nonstopmode main-kindle.tex >/dev/null
grep -E 'Output written' main-kindle.log
