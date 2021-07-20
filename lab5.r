require("e1071")
require("RWeka")
require("caret")
require("ggplot2")

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

# Se transforman los 2 y 4 de la variable class a "benigna" y "cancerosa" segÃºn corresponda
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

smp_size <- floor(0.70 * nrow(features.balanced))

set.seed(123)
train_ind <- sample(seq_len(nrow(features.balanced)), size = smp_size)

train <- features.balanced[train_ind, ]
test <- features.balanced[-train_ind, ]


################################
########### SVM ################
################################


formula = train$Class ~ .
set.seed(123)
model <- svm(formula, data = train)
print(model)
summary(model)

# test with train data
pred <- predict(model, test)
table(pred, test$Class)

# compute decision values and probabilities:
pred <- predict(model, test, decision.values = TRUE)
attr(pred, "decision.values")[1:7,]

# visualize (classes by color, SV by crosses):
plot(cmdscale(dist(features.balanced[,-4])), col = as.integer(features.balanced[,4]), pch = c("o","+")[1:478 %in% model$index + 1])

############ Lineal ###################


########### solo quitando mitosis y cell shape #################

# set.seed(20)
# obj <- tune(svm, Class ~ ., data = features.balanced, kernel = "linear", ranges = list(cost = 2^(-6:0)),tunecontrol = tune.control(sampling = "cross", cross = 5 ))
# summary(obj)
# plot(obj)
# 
# ########## cost = 0.25, error = 0.02929 ################
# # se quitan mas variables
# 
# ranking<-InfoGainAttributeEval(Class ~ . , data = features.balanced)
# print(ranking)

# se quita clumpThickness, normaNucleoi y marginalAdhesion

features.balanced <- features.balanced[, -1]
features.balanced <- features.balanced[, -6]
features.balanced <- features.balanced[, -2]

# set.seed(20)
# obj <- tune(svm, Class ~ ., data = features.balanced, kernel = "linear", ranges = list(cost = 2^(-6:0)),tunecontrol = tune.control(sampling = "cross", cross = 5 ))
# summary(obj)
# plot(obj)
# summary(obj$best.model)

########## cost = 1, error = 0.04181 ################
# se quitan mas variables

ranking<-InfoGainAttributeEval(Class ~ . , data = features.balanced)
print(ranking)

#se quita singleEpithelial

features.balanced <- features.balanced[, -2]

set.seed(20)
obj <- tune(svm, Class ~ ., data = features.balanced, kernel = "linear", ranges = list(cost = 2^(-6:0)),tunecontrol = tune.control(sampling = "cross", cross = 5 ))
summary(obj)
plot(obj)
summary(obj$best.model)
model <- obj$best.model
plot(cmdscale(dist(features.balanced[,-4])), col = as.integer(features.balanced[,4]), pch = c("o","+")[1:478 %in% model$index + 1])

pred <- predict(obj$best.model, features.balanced)
m <- table(pred, features.balanced$Class)


########## cost = 0.125, error = 0.04184 ################
# con cross = 5
########## cost = 0.25, error = 0.03975 ################


#se quita blandChromatin

# features.balanced <- features.balanced[, -3]
# 
# set.seed(20)
# obj <- tune(svm, Class ~ ., data = features.balanced, kernel = "linear", ranges = list(cost = 2^(-6:0)),tunecontrol = tune.control(sampling = "cross", cross = 2 ))
# summary(obj)
# plot(obj)
# summary(obj$best.model)
# 
# #se quita bareNuceli
# features.balanced <- features.balanced[, -2]
# 
# set.seed(20)
# obj <- tune(svm, Class ~ ., data = features.balanced, kernel = "linear", ranges = list(cost = 2^(-6:0)),tunecontrol = tune.control(sampling = "cross", cross = 2 ))
# summary(obj)
# plot(obj)
# summary(obj$best.model)
# 
# pred <- predict(obj$best.model, features.balanced)
# table(pred, features.balanced$Class)

# 
# plot(obj)
# summary(obj$best.model)
# 
# model <- obj$best.model
# 
# plot(cmdscale(dist(features.balanced[,-8])), col = as.integer(features.balanced[,8]), pch = c("o","+")[1:478 %in% model$index + 1])
# 
# pred <- predict(obj$best.model, features.balanced)
# table(pred, features.balanced$Class)

############ RADIAL ###################

# set.seed(20)
# obj1 <- tune(svm, Class ~ ., data = features.balanced, kernel = "radial", ranges = list(gamma = 2^(-2:4), cost = 2^(-1:4), tunecontrol = tune.control(sampling = "cross", cross = 5 )))
# summary(obj1)
# 
# ggplot(data = obj1$performances, aes(x = cost, y = error, color = as.factor(gamma)))+
#   geom_line() +
#   geom_point() +
#   labs(title = "Error de clasificación vs hiperparámetros C y gamma", color = "gamma") +
#   theme_bw() +
#   theme(legend.position = "bottom")
# summary(obj1$best.model)
# 
# pred <- predict(obj1$best.model, features.balanced)
# m1 <- table(pred, features.balanced$Class)

############################################
# gamma = 0.25, cost = 0.5, error = 0.0376 #
############################################

set.seed(20)
obj2 <- tune(svm, Class ~ ., data = features.balanced, kernel = "radial", ranges = list(gamma = 2^(-2:6), cost = 2^(-2:6), tunecontrol = tune.control(sampling = "cross", cross = 5 )))
summary(obj2)

pred <- predict(obj2$best.model, features.balanced)
m2 <- table(pred, features.balanced$Class)

ggplot(data = obj2$performances, aes(x = cost, y = error, color = as.factor(gamma)))+
  geom_line() +
  geom_point() +
  labs(title = "Error de clasificación vs hiperparámetros C y gamma", color = "gamma") +
  theme_bw() +
  theme(legend.position = "bottom")
summary(obj2$best.model)

model <- obj2$best.model

plot(cmdscale(dist(features.balanced[,-4])), col = as.integer(features.balanced[,4]), pch = c("o","+")[1:478 %in% model$index + 1])

############################################
# gamma = 0.125, cost = 1, error = 0.0376 #
############################################

VP <- m[1]
FP <- m[3]
VN <- m[4]
FN <- m[2]

precision = VP / (VP + FP)
recall = VP / (VP + FN)
calculoF1 <- 2*precision*recall/(precision + recall)


VP2 <- m2[1]
FP2 <- m2[3]
VN2 <- m2[4]
FN2 <- m2[2]

precision2 = VP2 / (VP2 + FP2)
recall2 = VP2 / (VP2 + FN2)
calculoF12 <- 2*precision2*recall2/(precision2 + recall2)