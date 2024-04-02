## set-up -- load all the libs
# the assumption is that the code will break here if the libraries aren't
# available preventing to 'source' it
# this version outputs counts instead of proportions as earlier
#
library(Matrix)
library(data.table)
library(plyr)
library(rhdf5)
library(anndata)
library(CARD)
library(reticulate)


# h5p <- reticulate::import("hdf5plugin") enable this if things wouldn't work with h5
set.seed(1)

## Define the CARD-based deconvolution function to be used by the pipeline
# references: https://github.com/YingMa0107/CARD, https://pubmed.ncbi.nlm.nih.gov/35501392/
#
# Example invocation: deconvolve_with_CARD(io_prefix='/app/data')
#
deconvolve_with_CARD <- function(# where to look for files
                                 coord_fname = "tissue_positions_list.csv", # tissue positions
                                 h5_fname = "filtered_feature_bc_matrix.h5", # spaceranger output
                                 spots_fname = "spots_df.csv", # spots data for cell count
                                 sc_atlas_fname = "", # the sc atlas to use
                                 sc_label_key = "cell_type",
                                 sc_sample_key = "batch",
                                 out_csv_counts_fname = "card_results.csv") {

  python_anndata <- reticulate::import("anndata", convert = FALSE)
  h5p <- reticulate::import("hdf5plugin", convert = FALSE)
  io_prefix = ""
  ## first -- read and organize the spatial coordinates
  #
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-1: read and organize the spatial coordinates ", file.path(io_prefix, coord_fname)))
  spatial_location <- data.frame(fread(file.path(io_prefix, coord_fname)))
  colnames(spatial_location) <- c(
    "barcode", "in_tissue", "row", "column",
    "pxl_row_in_fullres", "pxl_col_in_fullres"
  )
  spatial_location <- spatial_location[spatial_location$in_tissue == 1, c(1, 3, 4)] # columns hardcoded -- ATTENTION
  colnames(spatial_location) <- c("SpotID", "x", "y")
  rownames(spatial_location) <- spatial_location$SpotID
  spatial_location <- spatial_location[, -1]
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " spatial_location dim ",
               paste(dim(spatial_location), collapse = "x")))

  ## second -- read and organize the spatial expressions
  #
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-2: read and organize the spatial expressions ", file.path(io_prefix, h5_fname)))
  h5 <- rhdf5::h5read(file.path(io_prefix, h5_fname), "matrix")
  spatial_count <- Matrix::sparseMatrix(
    dims = h5$shape,
    i = as.numeric(h5$indices),
    p = as.numeric(h5$indptr),
    x = as.numeric(h5$data),
    index1 = FALSE
  )
  rownames(spatial_count) <- as.data.frame(h5[["features"]])$name # if we read $id -- we will get ensembl ids
  colnames(spatial_count) <- as.character(h5[["barcodes"]])

  # trim to size
  spatial_count <- spatial_count[, which(colnames(spatial_count) %in% rownames(spatial_location))]
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " spatial_count dim ", paste(dim(spatial_count), collapse="x")))
  spatial_location <- spatial_location[match(colnames(spatial_count), rownames(spatial_location)), ]
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " spatial_location dim ", paste(dim(spatial_location), collapse="x")))

  # cleanup a bit
  rm(h5)
  gc()

  ## third -- read the atlas data
  #
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-3: read the atlas data ", file.path(io_prefix, sc_atlas_fname)))

  anndataX_to_r <- function(x, nrow, ncol) {
    Matrix::sparseMatrix(repr = "R",
                         i = as.integer(reticulate::py_to_r(x$indices))+1,
                         p = as.integer(reticulate::py_to_r(x$indptr)),
                         x = as.vector(reticulate::py_to_r(x$data)),
                         dims = c(nrow, ncol)
    )
  }
  transpose_dgRMatrix <- function(inmat) {
    if (class(inmat) != "dgRMatrix") {
      stop("inmat is not of class dgRMatrix")
    }
    out <- new("dgCMatrix",
               i = inmat@j,
               p = inmat@p,
               x = inmat@x,
               Dim = rev(inmat@Dim),
               Dimnames = rev(inmat@Dimnames)
    )
    out
  }

  ad <- python_anndata$read_h5ad(file.path(io_prefix, sc_atlas_fname))

  sc_data <- anndataX_to_r(ad$X, nrow = dim(ad$var)[1], ncol = dim(ad$obs)[1])
  dimnames(sc_data) <- list(rownames(reticulate::py_to_r(ad$var)), rownames(reticulate::py_to_r(ad$obs)))
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " sc_data dim ", paste(dim(sc_data), collapse="x")))
  dim(sc_data)
  # do some sanity checks for the atlas data layout
  if (!all(c(sc_label_key) %in% colnames(reticulate::py_to_r(ad$obs)))) {
    stop(paste0("The CARD deconvolution cant find the column named '",
                sc_label_key, "' in the input file ",
                file.path(io_prefix, sc_atlas_fname)))
  }
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " sc_label_key ", sc_label_key))

  sc_labels <- as.character(reticulate::py_to_r(ad$obs)[[sc_label_key]]) # extract cell types/names
  sc_sample <- rep(1, length(sc_labels))
  if (!(is.na(sc_sample_key))) {
    sc_sample <- as.character(reticulate::py_to_r(ad$obs)[[sc_sample_key]]) # extract batches/samples
  }

  rm(ad)
  gc()

  ## fourth -- align the genes lists between canonical and gene names
  #
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-4: align the genes lists between the observed expressions and atlas data ", file.path(io_prefix, h5_fname)))
  gene_id <- rhdf5::h5read(file.path(io_prefix, h5_fname), "/matrix/features/id")
  gene_name <- rhdf5::h5read(file.path(io_prefix, h5_fname), "/matrix/features/name")
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " gene names: ", length(gene_name), ", unique gene names ", length(unique(gene_name))))

  genes_table <- data.frame(gene_id, gene_name)
  colnames(genes_table) <- c("GENES", "GNAME")

  # remove duplicates from names
  genes_table <- genes_table[which(!duplicated(genes_table$GNAME)), ]
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " dim genes_table: ", paste(dim(genes_table), collapse="x")))

  ## five -- subset the atlas data to the observed gene expressions
  #
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-5: subset the atlas data to the observed gene expressions"))

  # reduce the both genes set to only common names
  intersect_name_set <- base::intersect(rownames(sc_data), genes_table$GNAME)
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " length of intersect_name_set: ", length(intersect_name_set)))

  sc_data <- sc_data[which(rownames(sc_data) %in% intersect_name_set), ]
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " adjusted sc_data dim  ", paste(dim(sc_data), collapse="x")))

  # reduce the atlas -- keep only genes observed in the sample
  genes_table <- genes_table[which(genes_table$GNAME %in% intersect_name_set), ]
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " adjusted genes_table dim  ", paste(dim(genes_table), collapse="x")))

  # rearrange the atlas data
  sc_count <- sc_data[match(rownames(sc_data), genes_table$GNAME), ]
  rm(sc_data)

  sc_meta <- data.frame(cellID = colnames(sc_count), cellType = sc_labels, sampleInfo = sc_sample)
  rownames(sc_meta) <- colnames(sc_count)
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " sc_meta dim  ", paste(dim(sc_meta), collapse="x")))

  ## six -- perform CARD deconvolution
  #
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-6: perform CARD deconv preps"))
  CARD_obj <- CARD::createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5
  )

  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-6: perform CARD deconv"))
  CARD_obj <- CARD::CARD_deconvolution(CARD_object = CARD_obj)
  props <- CARD_obj@Proportion_CARD

  ## seven -- convert proportions to counts and write as a CSV
  #
  print(paste0(format(Sys.time(), "%D,%H:%M:%S"), " CARD-7: write cell counts as a CSV ", file.path(io_prefix, out_csv_counts_fname)))

  # add cell counts as the proportions data first column
  spots_df <- data.frame(fread(file.path(io_prefix, spots_fname)))
  props_df <- as.data.frame(props)
  props_df$barcode<-rownames(props)
  props_df <- merge(spots_df[c("barcode","num_contained_cells")], props_df, by=c("barcode"), all.x=TRUE)
  colnames(props_df)[colnames(props_df) == 'num_contained_cells'] <- 'cell_num'
  colnames(props_df) <- gsub("\\.+", " ", colnames(props_df))

  props_df[is.na(props_df)] <- 0
  rownames(props_df)<-props_df$barcode
  props_df = props_df[ ,colnames(props_df)!="barcode"]

  props <- props_df
  rm(spots_df)


  # define proportions to counts conversion
  props_to_counts <- function(props, cell_num) {
    bins <- length(props) # number of bins
    dist <- rep(0., bins) # initial distribution filled with zeros
    proportion <- props # create a copy
    # proportion <- proportion/sum(proportion)
    for (i in 1:cell_num) { # in the loop figure which bin to increment by 1
      total <- max(sum(dist), 1) # how many items total
      prop <- dist / total # current distribution
      error <- proportion - prop # errors with the desired distribution
      idx <- which.max(error) # which bin has a max error
      dist[idx] <- dist[idx] + 1 # increment it
    }
    dist
  }

  # convert and reshape
  res <- as.matrix(plyr::adply(
    props, 1,
    function(row) {
      props_to_counts(
        as.numeric(row[-1]),
        as.numeric(row[1])
      )
    }
  ))
  colnames(res) <- c(colnames(props), colnames(props)[-1])
  rownames(res) <- rownames(props)
  last_col_to_drop <- length(props[1, ]) # all proportion columns and extra col
  # with count
  res <- res[, -c(1:last_col_to_drop)]

  # produce the output
  data.table::fwrite(data.table::as.data.table(res, keep.rownames = "barcode"),
                     file.path(io_prefix, out_csv_counts_fname),
                     quote = TRUE,
                     sep = ",", col.names = TRUE, row.names = FALSE
  )

  ## teardown -- clean-up the allocated memory
  #
  rm(sc_count, sc_meta, spatial_count, spatial_location, CARD_obj)
  gc()
}
