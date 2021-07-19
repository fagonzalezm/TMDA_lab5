require("e1071")
require("RWeka")
require("caret")

######### Se carga la base de datos #########

names <- c("id", "clumpThickness", "uniformityOfCellSize", "uniformityOfCellShape",
           "marginalAdhesion", "singleEpithelialCellSize", "bareNuclei",
           "blandChromatin", "normalNucleoli", "mitoses", "class")
data <- read.table("./breast-cancer-wisconsin.data", sep=",", col.names = names)

set.seed(20)

################################
####### Preprocesamiento #######
################################

####### Missing Values ###########
# - Como existen 16 datos que presentan missing values en la variable barnuclei y el total de datos es 699,
#   se opta por eliminar estos datos.
data.original <- data

data.n = nrow(data)
data.m = ncol(data)

for (row in 1:data.n) {
  for (col in 1:data.m) {
    if (data[row, col] == "?") {
      data[row, col] <- NA
    }
  }
}
data$bareNuclei <- as.integer(data$bareNuclei)
data <- na.omit(data)

# Se quita la variable id
features <- data[,2:11]

# Se transforman los 2 y 4 de la variable class a "benigna" y "cancerosa" según corresponda
features$class <- factor(features$class, levels=c(2,4), labels=c("benigna","cancerosa"))

# Se balancean los datos utilizando un metodo de submuestreo down-sample
table(features$class)
features.balanced <- downSample(x = features[, -ncol(features)],
                         y = features$class)
table(features.balanced$Class)
