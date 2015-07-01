#' @import tidyr dplyr magrittr
load_micro_output = function(file) {

	n_col <- max(count.fields(file, sep = " "))
	output <- tbl_df(read.table(file ,sep=" ",fill=TRUE, col.names = 1:n_col))
	colnames(output) = c("run", "start", "time", as.character(0:(n_col-4)))
	output[is.na(output)] = 0
	#output = output[, colSums(output[-c(1:3)]) !=0]

	output = gather_(data = output, key_col = "infections", value_col = "population" , gather_cols = names(output)[-c(1:3)], convert=TRUE)

}

#' @import tidyr dplyr magrittr readr
load_micro_output_c = function(file) {

#	n_col <- max(count.fields(file, sep = " "))
#	output <- tbl_df(read.table(file ,sep=" ",fill=TRUE, col.names = 1:n_col))
  output = read_delim('micro.txt', " ", col_names = FALSE)
	colnames(output) = c("run", "start", "time", "control", as.character(0:(ncol(output)-5)))
	#mutate_each_(output, funs(na_to_0), as.character(0:(ncol(output)-5)))
	#output = output[, colSums(output[-c(1:3)]) !=0]
  output = as_data_frame(lapply(output, na_to_0))
	output = gather_(data = output, key_col = "infections", value_col = "population" , gather_cols = names(output)[-c(1:4)], convert=TRUE)

}

na_to_0 = function(x) replace(x, is.na(x), 0)

#' @import dplyr magrittr
restrict_micro_output = function(micro_output) {
  micro_output %>%
		group_by(run, time, start) %>%
		do({
			macro_state = restrict.micro_state(.$population)
			data_frame(N = macro_state[1], P = macro_state[2])
			})
}
