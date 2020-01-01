addDependenciesToDescription <- function(){
  
  # install.packages("renv")
  dependencies <- renv::dependencies(path=getwd(), root = getwd())
  
  uni_dependencies <- unique(dependencies$Package)                   
  uni_dependencies <- uni_dependencies[ uni_dependencies!="mstools"]
  
  for (pkg in uni_dependencies){
    pkg_ver <- packageVersion(pkg)
    cat( sprintf("Working on %s (>= %s)\n", pkg, pkg_ver))
    usethis::use_package(package = pkg, type = "Imports", min_version = pkg_ver)
  }
}