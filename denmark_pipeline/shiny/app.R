# Run from the project root with launch_config_app(), or shiny::runApp(...).
package_root <- normalizePath("..", winslash = "/")
source(file.path(package_root, "run.R"))
for (part in c("helpers", "demo", "app_ui", "app_server")) source(file.path(package_root, "shiny", paste0(part, ".R")))
defaults <- dp_app_defaults(getOption("denmark.project_root", dirname(package_root)))
options(denmark.app.defaults = defaults)
shiny::shinyApp(dp_app_ui(defaults), dp_app_server)
