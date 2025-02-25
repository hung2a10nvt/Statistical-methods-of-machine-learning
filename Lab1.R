library(moments)
library(MASS)
library(fitdistrplus)

sampleSize      <- 150      # общий объём выборки
lambdaExp       <- 1        # интенсивность для экспоненциального распределения
gammaShape      <- 2        # shape-параметр для гамма-распределения
gammaRate       <- 3        # rate-параметр для гамма-распределения
geomProb        <- 0.16     # p-параметр для геометрического распределения
binomN          <- 7        # количество испытаний для биномиального распределения
binomProb       <- 0.01     # вероятность успеха для биномиального распределения

# 1. БИНОМИАЛЬНОЕ РАСПРЕДЕЛЕНИЕ

# 1.1 Генерация выборки
binomSample <- rbinom(sampleSize, size = binomN, prob = binomProb)

# 1.2 Гистограмма + эмпирическая ф. распределения
histBinom <- hist(binomSample,
                  main = "Биномиальное распределение",
                  xlab = "Кол-во успехов")
lines(histBinom$counts ~ histBinom$mids, col = "blue")

ecdfBinom <- ecdf(binomSample) # эмпирическая функция распределения
plot(ecdfBinom, main = "Эмпирическая функция распределения (Binomial)")

# 1.3 Выборочные числовые характеристики
binomMean   <- mean(binomSample)
binomVar    <- var(binomSample)
binomSd     <- sd(binomSample)
binomMedian <- median(binomSample)
binomMode   <- density(binomSample)$x[which.max(density(binomSample)$y)]
binomSkew   <- skewness(binomSample)
binomKurt   <- kurtosis(binomSample)

# 1.4 Теоретические мат.ожидание и дисперсия
theorBinomMean <- binomN * binomProb    # E(X) = n.p
theorBinomVar  <- binomN * binomProb * (1 - binomProb) # Var(X) = n.p(1-p)

# (Пример «переоценки» параметра p по выборке)
binomProbEst <- binomMean / binomN

# 1.6 Проверка гипотезы о виде распределения критерием χ^2 
# (примеры бывают разные, здесь тестируют против нормального с μ=mean, σ=sd)
binomHistCounts <- hist(binomSample, plot = FALSE)$counts
binomHistBreaks <- hist(binomSample, plot = FALSE)$breaks
kBinom          <- length(binomHistCounts)
binomHistBreaks[1]       <- -Inf
binomHistBreaks[kBinom+1] <- Inf
binomTheorCdf   <- pnorm(binomHistBreaks, mean = binomMean, sd = binomSd)
binomProbTheor  <- binomTheorCdf[2:(kBinom+1)] - binomTheorCdf[1:kBinom]
chisq.test(binomHistCounts, p = binomProbTheor)

# 2. ГЕОМЕТРИЧЕСКОЕ РАСПРЕДЕЛЕНИЕ

# 2.1 Генерация выборки
geomSample <- rgeom(sampleSize, prob = geomProb)

# 2.2 fitdistr для оценки параметра
fitdistr(geomSample, "geometric")

# 2.3 Выборочные числовые характеристики
geomMean   <- mean(geomSample)
geomVar    <- var(geomSample)
geomSd     <- sd(geomSample)
geomMedian <- median(geomSample)
densGeom   <- density(geomSample)
geomMode   <- densGeom$x[which.max(densGeom$y)]
geomSkew   <- skewness(geomSample)
geomKurt   <- kurtosis(geomSample)

# 2.4 Теоретические мат.ожидание и дисперсия
theorGeomMean <- (1 - geomProb) / geomProb
theorGeomVar  <- (1 - geomProb) / (geomProb^2)

# 2.5 Гистограмма + эмпирическая функция распределения
histGeom <- hist(geomSample,
                 main = "Геометрическое распределение",
                 xlab = "Кол-во испытаний до первого успеха")
lines(histGeom$counts ~ histGeom$mids, col = "blue")

ecdfGeom <- ecdf(geomSample)
plot(ecdfGeom, main = "Эмпирическая функция распределения (Geometric)")

# 2.6 Проверка гипотезы о виде распределения критерием χ^2
geomHistCounts <- hist(geomSample, plot = FALSE)$counts
geomHistBreaks <- hist(geomSample, plot = FALSE)$breaks
kGeom          <- length(geomHistCounts)
geomHistBreaks[1]       <- -Inf
geomHistBreaks[kGeom+1] <- Inf
geomTheorCdf   <- pgeom(geomHistBreaks, prob = geomProb)
geomProbTheor  <- geomTheorCdf[2:(kGeom+1)] - geomTheorCdf[1:kGeom]
chisq.test(geomHistCounts, p = geomProbTheor)

# 3.1 Генерация выборки
expSample <- rexp(sampleSize, rate = lambdaExp)

# 3.2 Построение гистограммы + полигон и оценка плотности
histExp <- hist(expSample,
                col = "blue",
                main = "Экспоненциальное распределение",
                xlab = "Значения") 
lines(histExp$counts ~ histExp$mids, col = "red")  # полигон частот

densExp <- density(expSample)  # оценить эмпирическую функцию плотности (применимую к непрерывным переменным).
plot(densExp, main = "Оценка плотности (Exponential)")

# 3.3 Выборочные числовые характеристики
expMean   <- mean(expSample)
expVar    <- var(expSample)
expSd     <- sd(expSample) # standard deviation
expMedian <- median(expSample)
expMode   <- densExp$x[which.max(densExp$y)]
expSkew   <- skewness(expSample)
expKurt   <- kurtosis(expSample)

# 3.4 Теоретические мат.ожидание и дисперсия
theorExpMean <- 1 / lambdaExp
theorExpVar  <- 1 / (lambdaExp^2)

# 3.5 Оценка параметров через fitdistr
fitdist(expSample, "exp")
fitdistr(expSample, "exponential")

# 3.6 Проверка гипотезы о виде распределения критерием χ^2
expHistCounts <- hist(expSample, plot = FALSE)$counts
expHistBreaks <- hist(expSample, plot = FALSE)$breaks
kExp          <- length(expHistCounts)
expHistBreaks[1]     <- -Inf
expHistBreaks[kExp+1] <- Inf
expTheorCdf   <- pexp(expHistBreaks, rate = lambdaExp)
expProbTheor  <- expTheorCdf[2:(kExp+1)] - expTheorCdf[1:kExp]
chisq.test(expHistCounts, p = expProbTheor)

# 4. ГАММА-РАСПРЕДЕЛЕНИЕ

# 4.1 Генерация выборки
gammaSample <- rgamma(sampleSize, shape = gammaShape, rate = gammaRate)

# 4.2 Гистограмма + полигон и оценка плотности
histGamma <- hist(gammaSample,
                  col  = "green",
                  main = "Гамма-распределение",
                  xlab = "Значения")
lines(histGamma$counts ~ histGamma$mids, col = "red")
densGamma <- density(gammaSample)
plot(densGamma, main = "Оценка плотности (Gamma)")

# 4.3 Выборочные числовые характеристики
gammaMean   <- mean(gammaSample)
gammaVar    <- var(gammaSample)
gammaSd     <- sd(gammaSample)
gammaMedian <- median(gammaSample)
gammaMode   <- densGamma$x[which.max(densGamma$y)]
gammaSkew   <- skewness(gammaSample)
gammaKurt   <- kurtosis(gammaSample)

# 4.4 Теоретические мат.ожидание и дисперсия
theorGammaMean <- gammaShape / gammaRate
theorGammaVar  <- gammaShape / (gammaRate^2)

# 4.5 Оценка параметров
fitdistr(gammaSample, "gamma")

# 4.6 Проверка гипотезы о виде распределения критерием χ^2
gammaHistCounts <- hist(gammaSample, plot = FALSE)$counts
gammaHistBreaks <- hist(gammaSample, plot = FALSE)$breaks
kGamma          <- length(gammaHistCounts)
gammaHistBreaks[1]       <- -Inf
gammaHistBreaks[kGamma+1] <- Inf
gammaTheorCdf   <- pgamma(gammaHistBreaks, shape = gammaShape, rate = gammaRate)
gammaProbTheor  <- gammaTheorCdf[2:(kGamma+1)] - gammaTheorCdf[1:kGamma]
chisq.test(gammaHistCounts, p = gammaProbTheor)