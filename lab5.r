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

####### Se quita la variable cell shape ########

features.balanced <- features.balanced[, -3]

####### Seleccion de caracteristicas ###########

ranking<-InfoGainAttributeEval(Class ~ . , data = features.balanced)
print(ranking)

# en base a esto se quita mitosis

features.balanced <- features.balanced[, -8]

ranking<-InfoGainAttributeEval(Class ~ . , data = features.balanced)
print(ranking)


################################
########### SVM ################
################################


formula = features.balanced$Class ~ .
model <- svm(formula, data = features.balanced)
print(model)
summary(model)


x <- subset(features.balanced, select = -Class)
y <- features.balanced$Class
model <- svm(x, y)
# test with train data
pred <- predict(model, x)
table(pred, y)

# compute decision values and probabilities:
pred <- predict(model, x, decision.values = TRUE)
attr(pred, "decision.values")[1:7,]

# visualize (classes by color, SV by crosses):
plot(cmdscale(dist(features.balanced[,-8])), col = as.integer(features.balanced[,8]), pch = c("o","+")[1:478 %in% model$index + 1])


obj <- tune(svm, Class ~ ., data = features.balanced, kernel = "linear", ranges = list(cost = 2^(-1:4)),tunecontrol = tune.control(sampling = "cross", cross = 2 ))
summary(obj)

plot(obj)
summary(obj$best.model)

pred <- predict(obj$best.model, x)
table(pred, y)

obj <- tune(svm, Class ~ ., data = features.balanced, kernel = "radial", ranges = list(gamma = 2^(-2:4), cost = 2^(-1:4), tunecontrol = tune.control(sampling = "cross", cross = 2 )))
summary(obj)
