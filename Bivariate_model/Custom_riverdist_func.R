my.custom.line2network =  function (path = ".", layer, tolerance = 100, reproject = NULL, supplyprojection = NULL) {
  sp <- suppressWarnings(rgdal::readOGR(dsn = path, layer = layer, verbose = F))
  if (class(sp) != "SpatialLinesDataFrame") 
    stop("Specified shapefile is not a linear feature.")
  if (is.na(sp@proj4string@projargs) & !is.null(supplyprojection)) 
    sp@proj4string@projargs <- supplyprojection
  if (is.na(sp@proj4string@projargs)) 
    stop("Shapefile projection information is missing.  Use supplyprojection= to specify a Proj.4 projection to use.  If the input shapefile is in WGS84 geographic (long-lat) coordinates, this will be +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 (in double-quotes).  If so, it must also be reprojected using reproject=.")
  proj4 <- strsplit(sp@proj4string@projargs, split = " ")
  projected <- sp::is.projected(sp)
  if (is.null(reproject) & !projected) 
    stop("Distances can only be computed from a projected coordinate system. Use reproject= to specify a Proj.4 projection to use.")
  if (!is.null(reproject)) {
    sp <- sp::spTransform(sp, sp::CRS(reproject))
    proj4 <- strsplit(sp@proj4string@projargs, split = " ")
  }
  units <- "unknown"
  for (i in 1:length(proj4[[1]])) {
    if (proj4[[1]][i] != "") {
      proj4arg <- strsplit(proj4[[1]][i], split = "=")
      if (proj4arg[[1]][1] == "+units") {
        units <- proj4arg[[1]][2]
        cat("\n", "Units:", proj4arg[[1]][2], "\n")
      }
    }
  }
  if (length(sp@lines) > 1) {
    sp_line <- NA
    sp_seg <- NA
    lines <- list()
    j <- 1
    for (i in 1:length(sp@lines)) {
      for (k in 1:length(sp@lines[i][[1]]@Lines)) {
        lines[[j]] <- sp@lines[i][[1]]@Lines[[k]]@coords
        sp_line[j] <- i
        sp_seg[j] <- k
        j <- j + 1
      }
    }
  }
  if (length(sp@lines) == 1) {
    lines <- sp@lines[1][[1]]@Lines
    length <- length(lines)
    lines.new <- list()
    for (i in 1:length) {
      lines.new[[i]] <- lines[[i]]@coords
    }
    lines <- lines.new
    sp_line <- rep(1, length)
    sp_seg <- 1:length
  }
  length <- length(lines)
  rivID <- 1:length
  lineID <- data.frame(rivID, sp_line, sp_seg)
  connections <- calculateconnections(lines = lines, tolerance = tolerance)
  if (any(connections %in% 5:6)) 
    braided <- TRUE
  lengths <- rep(NA, length)
  for (i in 1:length) {
    lengths[i] <- pdisttot(lines[[i]])
  }
  names <- rep(NA, length)
  mouth.seg <- NA
  mouth.vert <- NA
  mouth <- list(mouth.seg, mouth.vert)
  names(mouth) <- c("mouth.seg", "mouth.vert")
  sequenced <- FALSE
  braided <- NA
  cumuldist <- list()
  for (i in 1:length) {
    xy <- lines[[i]]
    n <- dim(xy)[1]
    cumuldist[[i]] <- c(0, cumsum(sqrt(((xy[1:(n - 1), 1] - 
                                           xy[2:n, 1])^2) + ((xy[1:(n - 1),    2] - xy[2:n, 2])^2))))
  }
  out.names <- c("sp", "lineID", "lines", "connections", "lengths", 
                 "names", "mouth", "sequenced", "tolerance", "units", 
                 "braided", "cumuldist")
  out <- list(sp, lineID, lines, connections, lengths, names, 
              mouth, sequenced, tolerance, units, braided, cumuldist)
  names(out) <- out.names
  class(out) <- "rivernetwork"
  length1 <- length(out$lengths)
  suppressMessages(out <- removeduplicates(out))
  length2 <- length(out$lengths)
  if (length2 < length1) 
    cat("\n", "Removed", length1 - length2, "duplicate segments.", "\n")
  
  # THIS LINE CAUSES ISSUES COMENT IT OUT
  # THIS LINE CAUSES ISSUES COMENT IT OUT
  # suppressMessages(out <- removemicrosegs(out))
  
  length3 <- length(out$lengths)
  if (length3 < length2) 
    cat("\n", "Removed", length2 - length3, "segments with lengths shorter than the connectivity tolerance.", "\n")
  return(out)
}


