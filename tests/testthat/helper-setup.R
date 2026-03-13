# helper-setup.R — виконується testthat перед кожним тест-файлом
# Завантажує всі модулі проєкту в глобальне середовище

# Знаходимо PROJECT_ROOT (два рівні вище tests/testthat/)
wd <- getwd()
root_candidate <- normalizePath(file.path(wd, "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root_candidate, "R"))) {
  # Можливо, вже в PROJECT_ROOT
  root_candidate <- wd
}

if (dir.exists(file.path(root_candidate, "R"))) {
  old_wd <- getwd()
  setwd(root_candidate)
  if (!exists("CFG")) source("R/config.R")
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
    if (!grepl("config[.]R$", f)) source(f)
  }
  setwd(old_wd)
}
message("[helper-setup.R] Project modules loaded. CFG exists: ", exists("CFG"))
