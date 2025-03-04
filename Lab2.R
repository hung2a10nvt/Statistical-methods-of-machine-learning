# 1. Задание параметров и генерация выборок

set.seed(123)
N_vals <- c(100, 500, 1000)
true_mean <- 2
true_sd   <- 3         # Истинное среднеквадратичное отклонение (sigma)
alpha     <- 0.05      # Уровень значимости

# Критическое значение z для доверительных интервалов при известной дисперсии
z_alpha <- qnorm(1 - alpha/2) 

results <- data.frame(
  N              = integer(),
  sample_mean    = double(),
  left_knownVar  = double(),
  right_knownVar = double(),
  left_unknownVar= double(),
  right_unknownVar=double()
)

# Генерация, построение гистограмм/ qqplot
#      и оценка выборочных характеристик

for (n in N_vals) {
  x <- rnorm(n, mean = true_mean, sd = true_sd)
  
  hist(x,
       main = paste("Гистограмма при N =", n),
       xlab = "Значения выборки",
       freq = FALSE,
       col  = "lightblue",
       border = "gray")
  # Добавим поверх теоретическую плотность
  curve(dnorm(x, mean = true_mean, sd = true_sd),
        col = "red", lwd = 2, add = TRUE)
  
  # qqplot
  qqnorm(x,
         main = paste("QQ-plot при N =", n),
         pch  = 19,
         col  = "blue")
  qqline(x, col = "red", lwd = 2)
  
  # Выборочные оценки
  sample_mean <- mean(x)
  sample_var  <- var(x)
  sample_sd   <- sd(x)
  
  # 4. Доверительный интервал для среднего при известной дисперсии
  
  left_known  <- sample_mean - z_alpha * true_sd / sqrt(n)
  right_known <- sample_mean + z_alpha * true_sd / sqrt(n)
  
  # 5. Доверительный интервал для среднего при неизвестной дисперсии
  t_alpha <- qt(1 - alpha/2, df = n - 1)  # критическое t
  left_unknown  <- sample_mean - t_alpha * sample_sd / sqrt(n)
  right_unknown <- sample_mean + t_alpha * sample_sd / sqrt(n)
  
  # Запишем все результаты в таблицу
  results <- rbind(results,
                   data.frame(N               = n,
                              sample_mean     = sample_mean,
                              left_knownVar   = left_known,
                              right_knownVar  = right_known,
                              left_unknownVar = left_unknown,
                              right_unknownVar= right_unknown))
}

# 6. График зависимости точечных оценок и границ ДИ от объёма выборки

plot(results$N, 
     results$sample_mean,
     type = "b",
     pch  = 19,
     col  = "blue",
     ylim = range(c(results$left_knownVar, results$right_knownVar,
                    results$left_unknownVar, results$right_unknownVar,
                    results$sample_mean, true_mean)),
     xlab = "Объём выборки (N)",
     ylab = "Оценка среднего и границы ДИ",
     main = "Зависимость оценок среднего и ДИ от объёма выборки")

# ДИ при известной дисперсии
lines(results$N, results$left_knownVar,  type = "b", pch = 19, col = "red")
lines(results$N, results$right_knownVar, type = "b", pch = 19, col = "red")

# ДИ при неизвестной дисперсии
lines(results$N, results$left_unknownVar,  type = "b", pch = 19, col = "green")
lines(results$N, results$right_unknownVar, type = "b", pch = 19, col = "green")

# Истинное мат. ожидание
abline(h = true_mean, col = "black", lty = 2, lwd = 2)

legend("topright",
       legend = c("Выборочное среднее",
                  "ДИ (известная σ^2)",
                  "ДИ (неизвестная σ^2)",
                  "Истинное среднее"),
       col    = c("blue", "red", "green", "black"),
       lty    = c(1, 1, 1, 2),
       pch    = c(19, 19, 19, NA),
       lwd    = 2)
