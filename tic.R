do_package_checks()

if (Sys.getenv("DEV_VERSIONS") != "") {
  get_stage("install") %>%
    add_step(step_install_github(c("r-lib/rlang", "tidyverse/dplyr", "tidyverse/tidyr")))
}

if (Sys.getenv("BUILD_PKGDOWN") != "" && ci()$get_branch() == "master") {
  get_stage("before_deploy") %>% 
    add_step(step_setup_ssh()) %>% 
    add_step(step_setup_push_deploy(path = "docs", branch = "gh-pages"))

  get_stage("deploy") %>%
    add_code_step(
      pkgbuild::compile_dll(),
      prepare_call = remotes::install_github("r-lib/pkgbuild")
    ) %>%
    add_code_step(
      pkgdown::build_favicon(),
      prepare_call = install.packages("magick")
    ) %>%
    add_step(step_build_pkgdown(run_dont_run = TRUE)) %>%
    add_code_step(system('echo "fable.tidyverts.org" > docs/CNAME')) %>%
    add_step(step_do_push_deploy(path = "docs"))
 }


