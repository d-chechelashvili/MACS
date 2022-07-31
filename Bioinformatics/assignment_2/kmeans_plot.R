
args <- commandArgs(TRUE);
input.file <- as.character(args[1]);
plot.file <- as.character(args[2]);

data.table <- read.delim(file = input.file, header = FALSE, stringsAsFactors = FALSE);
data.table <- cbind(data.table, grepl("Clust", data.table[,3]));
data.table[,3] <- as.numeric(sub("Clust", "", data.table[,3])) + 1;
data.table <- cbind(data.table, (data.table[,4] + 0.4) * 3);
data.table[,4] <- data.table[,4] * -18 + 19;

png(filename = plot.file, height = 1000, width = 1000);
par(mai = c(1,1,0.3,0.3));
plot(
	data.table[,1], data.table[,2],
	xlab = "BRCA1 Gene Expression",
	ylab = "SOX2 Gene Expression",
	col = data.table[,3],
	pch = data.table[,4],
	cex = data.table[,5],
	cex.axis = 1.3,
	cex.lab = 2.5
	);
dev.off();


