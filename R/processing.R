#' @importFrom tidyr gather
#' @importFrom dply arrange
load_micro_output = function(file) {
	
	n_col <- max(count.fields(file, sep = " "))
	output <- read.table(file ,sep=" ",fill=TRUE, col.names = 1:n_col)
	colnames(output) = c("run", "start", "time", as.character(0:(n_col-4)))
	output[is.na(output)] = 0
	output = output[, colSums(output[-c(1:3)]) !=0]
	
	output %<>%
		gather(infections, population,  -run, -start, -time, convert=TRUE) %>%
		arrange(time, run, start, infections) %>%
		tbl_df
}

restrict_micro_output = function(micro_output) {
  micro_output %>%
		group_by(run, time, start) %>%
		do({
			macro_state = restrict.micro_state(.$population)
			data.frame(N = macro_state[1], P = macro_state[2])
			})
}