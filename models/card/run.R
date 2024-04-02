library(plumber)

# Define the keep-alive hook
pr("plumber.R") %>%
  pr_hook("preroute",  function(req) {
    req$headers$`Keep-Alive` <- "timeout=1800, max=2000"
    forward()
  }) %>%
  pr_run(port=80,host="0.0.0.0")

# # for Rscript run in terminal, please set path to "src/" or navigate to this directory
# pr("plumber.R") %>%
#   pr_hook(keep_alive_hook)%>%
#   pr_run(port=80,host="0.0.0.0")
