using Genie, Logging

Genie.Configuration.config!(
  server_port                     = 8002,
  server_host                     = "127.0.0.1",
  log_level                       = Logging.Info,
  log_to_file                     = false,
  server_handle_static_files      = true,
  path_build                      = "build",
  format_julia_builds             = true,
  format_html_output              = true,
  watch                           = true
)

ENV["JULIA_REVISE"] = "auto"
ENV["JULIA_PARDISO"] = "/Users/edgardovittoria/Applications/panua-pardiso-20230908-mac_arm64/lib"
ENV["PARDISO_LIC_PATH"] = "/Users/edgardovittoria/"
ENV["OMP_NUM_THREADS"] = 1