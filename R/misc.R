#' @importFrom plyr round_any
round_rand = function(x, accuracy = 1) {
  prob = (x - plyr::round_any(x, accuracy, f = floor)) / accuracy
  return(ifelse(runif(length(x)) < prob,
  							plyr::round_any(x, accuracy, f = ceiling),
  							plyr::round_any(x, accuracy, f = floor)))
	}

