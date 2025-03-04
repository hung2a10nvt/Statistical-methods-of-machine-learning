# Загрузка необходимых библиотек
library(moments)
library(MASS)
library(fitdistrplus)

set.seed(123)

# Общее значение размера выборки
N <- 150

###############################
# 1. БИНОМИАЛЬНОЕ РАСПРЕДЕЛЕНИЕ
###############################
# Параметры биномиального распределения
binom_n <- 10      # количество испытаний
binom_p <- 0.3     # вероятность успеха

# Генерация выборки
binom_sample <- rbinom(N, size = binom_n, prob = binom_p)

# Построение гистограммы и полигона частот
binom_hist <- hist(binom_sample,
                   main = "Биномиальное распределение",
                   xlab = "Количество успехов",
                   col = "lightblue")
lines(binom_hist$counts ~ binom_hist$mids, col = "blue", lwd = 2)

# Построение эмпирической функции распределения (ECDF)
binom_ecdf <- ecdf(binom_sample)
plot(binom_ecdf, main = "Эмпирическая функция распределения (Биномиальное)",
     xlab = "Количество успехов", ylab = "F(x)")

# Расчет выборочных характеристик
binom_mean    <- mean(binom_sample)
binom_var     <- var(binom_sample)
binom_sd      <- sd(binom_sample)
binom_median  <- median(binom_sample)
# Оценка моды через density (для дискретных данных можно использовать таблицу частот)
binom_mode    <- density(binom_sample)$x[which.max(density(binom_sample)$y)]
binom_skew    <- skewness(binom_sample)
binom_kurt    <- kurtosis(binom_sample)

# Теоретические значения
binom_theor_mean <- binom_n * binom_p
binom_theor_var  <- binom_n * binom_p * (1 - binom_p)

# Оценка параметра (например, p_est = выборочное среднее / число испытаний)
binom_p_est <- binom_mean / binom_n

# Проверка гипотезы о виде распределения с помощью критерия χ²
binom_x_vals     <- 0:binom_n  # возможные значения
binom_obs_counts <- table(factor(binom_sample, levels = binom_x_vals))
binom_exp_probs  <- dbinom(binom_x_vals, size = binom_n, prob = binom_p)
binom_exp_counts <- N * binom_exp_probs
binom_chisq_test <- chisq.test(binom_obs_counts, p = binom_exp_probs, rescale.p = TRUE)
print(binom_chisq_test)

###############################
# 2. ГЕОМЕТРИЧЕСКОЕ РАСПРЕДЕЛЕНИЕ
###############################
# Параметры геометрического распределения
geom_prob <- 0.16  # вероятность успеха

# Генерация выборки (количество неудач до первого успеха)
geom_sample <- rgeom(N, prob = geom_prob)

# Оценка параметра с помощью fitdistr
geom_fit <- fitdistr(geom_sample, "geometric")
print(geom_fit)

# Расчет выборочных характеристик
geom_mean    <- mean(geom_sample)
geom_var     <- var(geom_sample)
geom_sd      <- sd(geom_sample)
geom_median  <- median(geom_sample)
geom_density <- density(geom_sample)
geom_mode    <- geom_density$x[which.max(geom_density$y)]
geom_skew    <- skewness(geom_sample)
geom_kurt    <- kurtosis(geom_sample)

# Теоретические значения
geom_theor_mean <- (1 - geom_prob) / geom_prob
geom_theor_var  <- (1 - geom_prob) / (geom_prob^2)

# Построение гистограммы и эмпирической функции распределения
geom_hist <- hist(geom_sample,
                  main = "Геометрическое распределение",
                  xlab = "Количество неудач до первого успеха",
                  col = "lightgreen")
lines(geom_hist$counts ~ geom_hist$mids, col = "blue", lwd = 2)
geom_ecdf <- ecdf(geom_sample)
plot(geom_ecdf, main = "Эмпирическая функция распределения (Геометрическое)",
     xlab = "Количество неудач", ylab = "F(x)")

# Проверка гипотезы о виде распределения (χ²-тест) для геометрического распределения
geom_hist_counts <- hist(geom_sample, plot = FALSE)$counts
geom_hist_breaks <- hist(geom_sample, plot = FALSE)$breaks
k_geom          <- length(geom_hist_counts)
geom_hist_breaks[1]         <- -Inf
geom_hist_breaks[k_geom+1]    <- Inf
geom_theor_cdf  <- pgeom(geom_hist_breaks, prob = geom_prob)
geom_exp_probs  <- geom_theor_cdf[2:(k_geom+1)] - geom_theor_cdf[1:k_geom]
geom_chisq_test <- chisq.test(geom_hist_counts, p = geom_exp_probs)
print(geom_chisq_test)

###############################
# 3. ЭКСПОНЕНЦИАЛЬНОЕ РАСПРЕДЕЛЕНИЕ
###############################
# Параметры экспоненциального распределения
exp_lambda <- 1  # интенсивность (λ)

# Генерация выборки
exp_sample <- rexp(N, rate = exp_lambda)

# Построение гистограммы, полигона частот и оценка плотности
exp_hist <- hist(exp_sample,
                 col = "lightblue",
                 main = "Экспоненциальное распределение",
                 xlab = "Значения",
                 probability = TRUE)
lines(exp_hist$counts ~ exp_hist$mids, col = "red", lwd = 2)
exp_density <- density(exp_sample)
plot(exp_density, main = "Оценка плотности (Экспоненциальное)",
     xlab = "Значения", ylab = "Плотность")

# Расчет выборочных характеристик
exp_mean    <- mean(exp_sample)
exp_var     <- var(exp_sample)
exp_sd      <- sd(exp_sample)
exp_median  <- median(exp_sample)
exp_mode    <- exp_density$x[which.max(exp_density$y)]
exp_skew    <- skewness(exp_sample)
exp_kurt    <- kurtosis(exp_sample)

# Теоретические значения
exp_theor_mean <- 1 / exp_lambda
exp_theor_var  <- 1 / (exp_lambda^2)

# Оценка параметров через fitdistr и fitdist
exp_fit1 <- fitdistr(exp_sample, "exponential")
exp_fit2 <- fitdist(exp_sample, "exp")
print(exp_fit1)
print(exp_fit2)

# Проверка гипотезы о виде распределения (χ²-тест) для экспоненциального распределения
exp_hist_counts <- hist(exp_sample, plot = FALSE)$counts
exp_hist_breaks <- hist(exp_sample, plot = FALSE)$breaks
k_exp          <- length(exp_hist_counts)
exp_hist_breaks[1]      <- -Inf
exp_hist_breaks[k_exp+1] <- Inf
exp_theor_cdf  <- pexp(exp_hist_breaks, rate = exp_lambda)
exp_exp_probs  <- exp_theor_cdf[2:(k_exp+1)] - exp_theor_cdf[1:k_exp]
exp_chisq_test <- chisq.test(exp_hist_counts, p = exp_exp_probs)
print(exp_chisq_test)

###############################
# 4. ГАММА-РАСПРЕДЕЛЕНИЕ
###############################
# Параметры гамма-распределения
gamma_shape <- 2   # параметр shape
gamma_rate  <- 3   # параметр rate

# Генерация выборки
gamma_sample <- rgamma(N, shape = gamma_shape, rate = gamma_rate)

# Построение гистограммы, полигона частот и оценка плотности
gamma_hist <- hist(gamma_sample,
                   col = "lightgreen",
                   main = "Гамма-распределение",
                   xlab = "Значения",
                   probability = TRUE)
lines(gamma_hist$counts ~ gamma_hist$mids, col = "red", lwd = 2)
gamma_density <- density(gamma_sample)
plot(gamma_density, main = "Оценка плотности (Гамма)",
     xlab = "Значения", ylab = "Плотность")

# Расчет выборочных характеристик
gamma_mean    <- mean(gamma_sample)
gamma_var     <- var(gamma_sample)
gamma_sd      <- sd(gamma_sample)
gamma_median  <- median(gamma_sample)
gamma_mode    <- gamma_density$x[which.max(gamma_density$y)]
gamma_skew    <- skewness(gamma_sample)
gamma_kurt    <- kurtosis(gamma_sample)

# Теоретические значения
gamma_theor_mean <- gamma_shape / gamma_rate
gamma_theor_var  <- gamma_shape / (gamma_rate^2)

# Оценка параметров с помощью fitdistr
gamma_fit <- fitdistr(gamma_sample, "gamma")
print(gamma_fit)

# Проверка гипотезы о виде распределения (χ²-тест) для гамма-распределения
gamma_hist_counts <- hist(gamma_sample, plot = FALSE)$counts
gamma_hist_breaks <- hist(gamma_sample, plot = FALSE)$breaks
k_gamma          <- length(gamma_hist_counts)
gamma_hist_breaks[1]       <- -Inf
gamma_hist_breaks[k_gamma+1] <- Inf
gamma_theor_cdf  <- pgamma(gamma_hist_breaks, shape = gamma_shape, rate = gamma_rate)
gamma_exp_probs  <- gamma_theor_cdf[2:(k_gamma+1)] - gamma_theor_cdf[1:k_gamma]
gamma_chisq_test <- chisq.test(gamma_hist_counts, p = gamma_exp_probs)
print(gamma_chisq_test)
