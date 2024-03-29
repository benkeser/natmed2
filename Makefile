md:
	Rscript -e "install.packages('../natmed2', type = 'source', repos = NULL)" && \
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"

install: 
	Rscript -e "install.packages('../natmed2', type = 'source', repos = NULL)"

site:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"
	Rscript -e "pkgdown::build_site()"

check:
	Rscript -e "devtools::check()"

checkfast:
	Rscript -e "devtools::check(build_args = '--no-build-vignettes')"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

build:
	Rscript -e "devtools::build()"

buildfast:
	Rscript -e "devtools::build(vignettes = FALSE)"

style:
	Rscript -e "styler::style_pkg()"
