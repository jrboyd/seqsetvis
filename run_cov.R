#test by test file
library(seqsetvis)
 for(f in dir("tests/testthat/", pattern = "test", full.names = T)){
   print(f)
   testthat::test_file(f)
 }

#devtools::build()

#build properly using Rbuildignore
# devtools::install(local = F)


a = covr::file_coverage(dir("R", full.names = T, pattern = "R$"), dir("tests/testthat/", full.names = T, pattern = "R$"))
covr::report(a)

covr::codecov()
covr::codecov(token = "ec9fa5ec-1e23-4e3e-82bb-260ed1ee514a")
