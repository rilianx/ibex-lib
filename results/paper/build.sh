#!/bin/sh
# Builds intro-related.pdf and intro-related.txt from the LaTeX source without a TeX
# installation: pandoc (LaTeX -> HTML with MathML, citations via citeproc) and headless
# Chromium (HTML -> PDF). With LaTeX available, `latexmk -pdf intro-related.tex` also works.
set -e
cd "$(dirname "$0")"
# pandoc drops equation numbers and numbers \paragraph: patch a temporary copy.
python3 - <<'PY'
import re
s = open('intro-related.tex').read()
labels = re.findall(r'\\label\{(eq:[^}]*)\}', s)
for k, lab in enumerate(labels, 1):
    s = re.sub(r'\\label\{%s\}(.*?)\\end\{equation\}' % re.escape(lab),
               lambda m: m.group(1) + r'\qquad\qquad (%d)\end{equation}' % k, s, flags=re.S)
    s = s.replace(r'\eqref{%s}' % lab, '(%d)' % k)
s = re.sub(r'\\paragraph\{([^}]*)\}', r'\\noindent\\textbf{\1} ', s)
open('.tmp.tex', 'w').write(s)
PY
pandoc .tmp.tex -s --number-sections --citeproc --bibliography=refs.bib --mathml \
  --metadata link-citations=true --css=style.css --embed-resources -o intro-related.html
pandoc .tmp.tex --number-sections --citeproc --bibliography=refs.bib -t plain --columns=90 \
  -o intro-related.txt 2>/dev/null
node -e "
const {chromium}=require('playwright');
(async()=>{const b=await chromium.launch();const p=await b.newPage();
await p.goto('file://'+process.cwd()+'/intro-related.html');
await p.pdf({path:'intro-related.pdf',format:'A4',margin:{top:'25mm',bottom:'25mm',left:'25mm',right:'25mm'}});
await b.close();})();"
rm -f intro-related.html .tmp.tex
