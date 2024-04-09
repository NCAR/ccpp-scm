# SCM documentation

The SCM technical documentation/users' guide is written in [reStructured Text markup language (RST)](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html). It is designed to be built with the [Sphinx documentation tool](https://www.sphinx-doc.org/en/master/index.html), and is documentation from the `main` branch is configured to build automatically on ReadTheDocs: https://ccpp-scm.readthedocs.io/en/latest/

## Building docs locally

To build the HTML documentation locally (for example, when testing modifications to the documentation), you must have Sphinx [installed on your machine](https://www.sphinx-doc.org/en/master/usage/installation.html), as well as the required python packages listed in the `requirements.txt` file. Once the prerequisite software is installed, you can build the HTML documentation using the following commands:

    python -m sphinx -T -b html -d _build/doctrees -D language=en . html

To build a PDF locally, in addition to Sphinx you must have Perl as well as the Perl-based [`latexmk` utility](https://mg.readthedocs.io/latexmk.html#installation). You can then build a PDF using the following commands:

    python -m sphinx -T -b latex -d _build/doctrees -D language=en . latex
    latexmk -f -pdf -pdflatex="pdflatex" -use-make index.tex

In lieu of running the above commands manually, you can also use `make`. The command `make all` will build both HTML and PDF documentation, while `make html` and `make latex` will build HTML and PDF respectively.

