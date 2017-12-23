#test by test file
# for(f in dir("tests/testthat/", pattern = "test", full.names = T)){
#   print(f)
#   testthat::test_file(f)
# }

#build properly using Rbuildignore
# devtools::install(local = F)


a = covr::file_coverage(dir("R", full.names = T, pattern = "R$"), dir("tests/testthat/", full.names = T, pattern = "R$"))
covr::report(a)
