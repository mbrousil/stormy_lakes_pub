library(targets)
library(future)

plan(future.callr::callr)

tar_make_future(workers = (nbrOfFreeWorkers(background = TRUE)))

# Create a network diagram of the workflow, with a completion timestamp
temp_vis <- tar_visnetwork()

temp_vis$x$main$text <- paste0("Last completed: ", Sys.time())

htmltools::save_html(html = temp_vis,
                     file = "figures/current_visnetwork.html")
