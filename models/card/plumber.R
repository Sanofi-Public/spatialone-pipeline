## ----- api for CARD
library(plumber)
library(jsonlite)
# library(future)

#* @apiTitle CARD (R) model API for Cell Deconvolution
#* @apiDescription Endpoints for accessing CARD(R) model

## ---- filter-logger
#* Log some information about the incoming request
#* @filter logger
function(req) {
  cat(
    as.character(Sys.time()), "-",
    req$REQUEST_METHOD, req$PATH_INFO, "-",
    req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n"
  )
  plumber::forward()
}


## ---- predict (POST)
#* Submit data and get a prediction in return
#* @post /predict
#* @param configs
#* @param atlas_filename_ext
#* @param sample_id
#* @serializer unboxedJSON
#* @return res Success
function(configs, atlas_filename_ext, sample_id, prep_path, results_path, reference_path, prefix, res) {
  # setup
  configs_list <- jsonlite::fromJSON(configs) # conver Json -> list

  # parse
  params <- configs_list[["params"]]

  res$status <- 202

  # 4- prepare to run CARD
  source("CARD_analysis.R")
  print("loaded CARD_analysis ---")

  #prep_path <- paste0("data/prep/", sample_id) #TODO take as param
  #results_path <- paste0("data/results/", sample_id)

  # 5- run CARD (call CARD's run/predict method here)
  deconvolve_with_CARD(
    coord_fname = paste0(prep_path, "tissue_positions_list.csv"),
    h5_fname = paste0(prep_path, "filtered_feature_bc_matrix.h5"),
    spots_fname = paste0(results_path, "spots_df.csv"),
    sc_atlas_fname = paste0(reference_path, params[["atlas_type"]], "_cell_atlas.h5ad"),
    sc_label_key = "cell_type",
    sc_sample_key = params[["sc_sample_key"]],
    out_csv_counts_fname = paste0(results_path, "card_results.csv")
  )
  print("DONE ---")

  res$status <- 200
}

## ---- new endpoint (GET)
#* Get a response from the new endpoint
#* @get /test
function() {
  result <- list(
    message = "This is a test endpoint."
  )
  return(result)
}

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
function(msg = "") {
  list(msg = paste0("The message is: '", msg, "'"))
}
