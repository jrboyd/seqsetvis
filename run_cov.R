#test by test file
library(seqsetvis)
 for(f in dir("tests/testthat/", pattern = "test", full.names = T)){
   print(f)
   testthat::test_file(f)
 }

devtools::build(vignettes = T)
system("R CMD check ../seqsetvis_0.99.1.tar.gz")

#build properly using Rbuildignore
install.packages()
# devtools::install(build_vignettes = F)

roxygen2::roxygenise()

a = covr::file_coverage(dir("R", full.names = T, pattern = "R$"), dir("tests/testthat/", full.names = T, pattern = "R$"))
covr::report(a)

# covr::codecov()
# covr::codecov(token = "ec9fa5ec-1e23-4e3e-82bb-260ed1ee514a")
