otu_data <- data.frame(
  SampleID = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"),
  OTU1 = c(10, 0, 5, 20, 15, 10),
  OTU2 = c(5, 25, 15, 0, 10, 5),
  OTU3 = c(0, 10, 5, 10, 20, 5),
  OTU4 = c(20, 5, 5, 10, 5, 15),
  OTU5 = c(15, 5, 10, 5, 0, 20),
  OTU6 = c(30, 0, 15, 5, 0, 10)
)
write.csv(otu_data, file ="C:/Users/samet/OneDrive/Belgeler/BioParrot/data/otu_table.csv", row.names = FALSE)
library(data.table)
fwrite(otu_data, "otu_table.csv")
