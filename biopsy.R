# Algorytm boosting w klasyfikacji

# IMPORT ZBIORU DANYCH
library(MASS)
library(gbm)
library(tibble)
head(biopsy)
glimpse(biopsy)

# PRZYGOTOWANIE DANYCH
biopsy$class <- ifelse(biopsy$class == "malignant", 1, 0)

biopsy <- na.omit(biopsy)

idx <- sample(seq(1, 2), size = nrow(biopsy), 
              replace = TRUE, prob = c(.8, .2))
train <- (1:nrow(biopsy))[idx == 1]
test <- (1:nrow(biopsy))[idx == 2]

names(biopsy) <- c("ID","thickness","c.size","c.shape","adhesion",
                   "s.size","nuclei","chomatine","nucleoi",
                   "mitoses","class")

# GRADIENT BOOSTING
biopsy.boost <- gbm(class~.-ID, data = biopsy[train,], 
                    distribution = "bernoulli", 
                    n.trees = 1000, 
                    cv.folds = 10)
summary(biopsy.boost, las = 2)

plot(biopsy.boost, 
     as.character(summary(biopsy.boost)$var[1]))
plot(biopsy.boost, 
     as.character(summary(biopsy.boost)$var[2]))
plot(biopsy.boost, 
     as.character(summary(biopsy.boost)$var[3]))

# PREDYKCJA
yhat.boost <- predict(biopsy.boost, 
                      newdata = biopsy[test, ], 
                      n.trees = 1000)
mean((yhat.boost-biopsy[test,]$class)^2)

# MACIERZ BLEDOW
results <- rep(0, length = length(test))
results[yhat.boost > 0.5] = 1
tb <- table(pred = results, 
            actual = biopsy$class[-train])
tb

# ACC
sum(diag(tb))/sum(tb)

# DODATEK
ntree.opt.oob <- gbm.perf(biopsy.boost,
                            method = "OOB", 
                            plot.it = F)
