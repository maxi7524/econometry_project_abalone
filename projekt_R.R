# -------- INICJALIZACJA ŚRODOWISKA I ŁADOWANIE PAKIETÓW -------- #
options(encoding = "UTF-8")
Sys.setlocale("LC_ALL", "pl_PL.UTF-8")

# Sprawdzenie i instalacja brakujących pakietów
required_packages <- c(
  "tidyverse", "broom", "car", "performance", 
  "corrplot", "lmtest", "olsrr", "caret", "recipes", "MASS"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(broom)
library(car)
library(performance)
library(corrplot)
library(lmtest)
library(olsrr)
library(caret)
library(recipes)
library(MASS)    # potrzebne do funkcji boxcox()
library(tibble)  # potrzebne do column_to_rownames()

seed <- 1333
set.seed(seed)

# -------- FUNKCJE POMOCNICZE -------- #
vifs_calculation <- function(df, entry_name = NULL) {
  # Obliczanie wartości VIF dla ramki df; entry_name opcjonalnie jako nazwa wiersza
  temp_df <- cbind(df, y_temp = rnorm(nrow(df)))
  vif_results <- car::vif(lm(y_temp ~ ., data = temp_df))
  
  if (!is.null(entry_name)) {
    vif_df <- data.frame(NAME_INDEX = entry_name, t(vif_results)) %>%
      tibble::column_to_rownames("NAME_INDEX")
  } else {
    vif_df <- data.frame(t(vif_results))
  }
  return(vif_df)
}

F_test_permutation <- function(df, y, test_cols, n_perm) {
  # Przeprowadza test F z permutacjami dla zmiennych test_cols w df, wektor y, liczba permutacji n_perm
  results <- list()
  X_reduced <- df %>% select(-all_of(test_cols))
  
  for (test_col in test_cols) {
    X_full <- df %>% select(-all_of(setdiff(test_cols, test_col)))
    
    model_reduced <- lm(y ~ ., data = cbind(X_reduced, y = y))
    model_full    <- lm(y ~ ., data = cbind(X_full, y = y))
    
    f_stat <- anova(model_reduced, model_full)$F[2]
    
    f_stats <- replicate(n_perm, {
      X_perm <- X_full
      X_perm[[test_col]] <- sample(X_full[[test_col]])
      anova(model_reduced, lm(y ~ ., data = cbind(X_perm, y = y)))$F[2]
    })
    
    p_value <- (sum(f_stats >= f_stat) + 1) / (n_perm + 1)
    
    results[[test_col]] <- data.frame(
      `Testowana zmienna` = test_col,
      `Statystyka F`     = f_stat,
      `P value`          = p_value,
      check.names = FALSE
    )
  }
  bind_rows(results) %>% tibble::column_to_rownames("Testowana zmienna")
}

# -------- WCZYTANIE DANYCH -------- #
COL_NAMES <- c(
  'Płeć', 'Długość', 'Średnica', 'Wysokość', 'Cała_waga',
  'Waga_po_obraniu', 'Waga_trzewi', 'Waga_powłoki', 'Pierścienie'
)

if (!file.exists('abalone_data_set.csv')) {
  stop("Brak pliku danych. Proszę zapewnić 'abalone_data_set.csv'")
} else {
  df <- read.csv(
    'abalone_data_set.csv',
    col.names   = COL_NAMES,
    fileEncoding = "UTF-8"
  )
}

# Użycie dplyr::select(), aby uniknąć konfliktów z innymi pakietami
X_reduced <- df %>%
  dplyr::select(-Pierścienie)
y <- df %>%
  dplyr::select(Pierścienie)

# Kodowanie zmiennej Płeć na trzy kolumny binarne
X_reduced <- X_reduced %>%
  mutate(
    Płeć_M = as.numeric(Płeć == "M"),
    Płeć_F = as.numeric(Płeć == "F"),
    Płeć_I = as.numeric(Płeć == "I")
  ) %>%
  dplyr::select(
    Płeć_I, Płeć_M, Płeć_F,
    everything(),
    -Płeć
  )

# -------- PRZETWARZANIE DANYCH DLA WYBRANEJ PŁCI -------- #
SEX <- 'I'
filt <- X_reduced[[paste0('Płeć_', SEX)]] == 1
sex_columns <- c('Płeć_I', 'Płeć_M', 'Płeć_F')

# Wyodrębnienie danych dla danej płci
y_local <- y[filt, , drop = FALSE]
X_old   <- X_reduced[filt, !(names(X_reduced) %in% sex_columns), drop = FALSE]

# Tworzenie nowych cech w X_old
X_new <- X_old %>%
  mutate(
    `Długość/Wysokość`          = Długość / Wysokość,
    `Masa/objętość`             = Cała_waga / (Długość * Średnica * Wysokość),
    `Waga_powłoki/Cała_waga`    = Waga_powłoki / Cała_waga,
    `Waga_po_obraniu/Cała_waga` = Waga_po_obraniu / Cała_waga,
    `Waga_trzewi/Cała_waga`     = Waga_trzewi / Cała_waga,
    `Waga_powłoki/Waga_po_obraniu` = Waga_powłoki / Waga_po_obraniu,
    `Waga_trzewi/Waga_po_obraniu`  = Waga_trzewi / Waga_po_obraniu,
    `Waga_trzewi/Waga_powłoki`     = Waga_trzewi / Waga_powłoki
  )

# Usuwanie wierszy z NA lub Inf
invalid_rows <- rowSums(is.na(X_new) | is.infinite(as.matrix(X_new))) > 0
X_old   <- X_old[!invalid_rows, , drop = FALSE]
X_new   <- X_new[!invalid_rows, , drop = FALSE]
y_local <- y_local[!invalid_rows, , drop = FALSE]

NEW_FEATURE_NAMES <- setdiff(names(X_new), names(X_old))

# -------- ANALIZA MODELU LINIOWEGO (DIAGNOSTYKA) -------- #
model_data <- cbind(X_old, Pierścienie = y_local$Pierścienie)
model_old  <- lm(Pierścienie ~ ., data = model_data)

# Obliczenie reszt i przygotowanie wykresów diagnostycznych
residuals <- residuals(model_old)

par(mfrow = c(1, 2))
plot(fitted(model_old), residuals, main = "Reszty vs dopasowane")
abline(h = 0, col = "red")
qqnorm(residuals, main = "Wykres QQ")
qqline(residuals, col = "red")

cat("Test Shapiro-Wilka p-value:", shapiro.test(residuals)$p.value, "\n")
cat("Test Breusch-Pagan p-value:", bptest(model_old)$p.value, "\n")

# -------- chunk6: TEST F Z PERMUTACJAMI I REDUKCJA MODELU -------- #

# Ustawienie poziomu istotności i progu VIF
significant_level <- 0.05
vif_reduce_value   <- 12

# Dodanie kolumny informującej o odrzuceniu H0 przy danym poziomie istotności
col_name_h0 <- paste0("Odrzucenie H0 dla alpha = ", significant_level)
F_test_before_reduction[[col_name_h0]] <- F_test_before_reduction$`P value` < significant_level

# Wybór zmiennych, dla których nie odrzucono H0 (czyli P value >= 0.05)
reduce_features_F_test <- rownames(
  F_test_before_reduction[
    !F_test_before_reduction[[col_name_h0]],
  ]
)

# Iteracyjna redukcja na podstawie VIF (próg = 12)
reduce_features_vif_test <- c()
vifs_after_F_test      <- list()
index_name            <- "Model rozszerzony przed redukcją vif"

repeat {
  drop_features    <- c(reduce_features_F_test, reduce_features_vif_test)
  current_features <- setdiff(colnames(X_final), drop_features)
  
  local_vifs <- vifs_calculation(
    X_final[, current_features, drop = FALSE],
    index_name
  )
  
  max_vif_column <- names(which.max(local_vifs[1, ]))
  max_vif_value  <-   max(local_vifs[1, ])
  
  vifs_after_F_test[[length(vifs_after_F_test) + 1]] <- local_vifs
  reduce_features_vif_test <- c(reduce_features_vif_test, max_vif_column)
  
  if (max_vif_value < vif_reduce_value) break
}

# Po zebraniu listy zmiennych do usunięcia (z testu F i z iteracyjnych VIF):
# utworzenie nowej macierzy X_final_reduced, z której usuwamy wszystkie te kolumny
to_drop         <- c(reduce_features_F_test, reduce_features_vif_test)
X_final_reduced <- X_final[, setdiff(colnames(X_final), to_drop), drop = FALSE]

# Złączenie wyników VIF w jedną ramkę (wypełnienie brakujących kolumn NA, jeśli zachodzi potrzeba)
vifs_after_F_test <- dplyr::bind_rows(vifs_after_F_test)

# Porównanie R-squared między:
#  1) modelem pełnym na wszystkich zmiennych X_new (bez redukcji VIF i F),
#  2) a modelem końcowym na X_final_reduced (po usunięciu zmiennych)
model_full    <- lm(Pierścienie ~ ., data = cbind(X_new,            Pierścienie = y_local$Pierścienie))
model_reduced <- lm(Pierścienie ~ ., data = cbind(X_final_reduced, Pierścienie = y_local$Pierścienie))

cat("R-squared model pełny:     ", summary(model_full)$r.squared, "\n")
cat("R-squared model po redukcji:", summary(model_reduced)$r.squared, "\n")

# -------- chunk7: REGRESJA Z BOX–COX I USUWANIE OUTLIERÓW -------- #

# Załaduj dodatkowe biblioteki potrzebne do analiz Box–Cox i Cooka
library(MASS)    # do funkcji boxcox()
library(car)     # do funkcji cooks.distance()
library(tibble)  # do funkcji column_to_rownames()

set.seed(seed)

# 1) Podział na zbiór uczący i testowy (75% / 25%)
n_total   <- nrow(X_final_reduced)
train_idx <- sample(seq_len(n_total), size = floor(0.75 * n_total))
test_idx  <- setdiff(seq_len(n_total), train_idx)

train_data <- cbind(
  Pierścienie = y_local$Pierścienie[train_idx],
  X_final_reduced[train_idx, , drop = FALSE]
)
test_data  <- cbind(
  Pierścienie = y_local$Pierścienie[test_idx],
  X_final_reduced[test_idx,  , drop = FALSE]
)

# Przygotowanie nazw metryk i listy wyników
metric_names  <- c(
  "Probe_size",
  "R2_train",
  "R2_test",
  "Lambda_hat",
  "LR_p",
  "Shapiro_train_p",
  "Shapiro_train_bc_p",
  "Shapiro_test_p",
  "Shapiro_test_bc_p"
)
metrics_list <- vector("list", length = 7)
names(metrics_list) <- paste0("iter_", 0:6)

# Indeksy obserwacji w bieżącym zbiorze uczącym (początkowo wszystkie)
current_idx <- seq_len(nrow(train_data))

for (i in seq_len(7)) {
  # 2) Wyodrębnienie podzbioru uczącego po usunięciu ewentualnych outlierów
  train_loc <- train_data[current_idx, , drop = FALSE]
  y_loc     <- train_loc$Pierścienie
  
  # 3) Dopasowanie modelu bez transformacji
  model_raw <- lm(Pierścienie ~ ., data = train_loc)
  
  # 4) Obliczenie R^2 na zbiorze treningowym
  r2_train <- summary(model_raw)$r.squared
  
  # 5) Obliczenie R^2 na zbiorze testowym
  preds_test <- predict(model_raw, newdata = test_data)
  ss_res     <- sum((test_data$Pierścienie - preds_test)^2)
  ss_tot     <- sum((test_data$Pierścienie - mean(test_data$Pierścienie))^2)
  r2_test    <- 1 - ss_res / ss_tot
  
  # 6) Reszty surowe i test Shapiro–Wilk (jeśli >= 3 obserwacje)
  resid_train <- resid(model_raw)
  resid_test  <- test_data$Pierścienie - preds_test
  
  sh_tr_p <- if (length(resid_train) >= 3) shapiro.test(resid_train)$p.value else NA
  sh_te_p <- if (length(resid_test)  >= 3) shapiro.test(resid_test)$p.value  else NA
  
  # 7) Estymacja lambda Box–Cox (MLE) z wykorzystaniem MASS::boxcox()
  bc_profile  <- boxcox(Pierścienie ~ ., data = train_loc, plotit = FALSE)
  lambda_seq  <- bc_profile$x
  loglik_seq  <- bc_profile$y
  lambda_hat  <- lambda_seq[which.max(loglik_seq)]
  
  # 8) Test ilorazu wiarygodności względem lambda = 1
  loglik_hat     <- max(loglik_seq)
  idx_lambda1    <- which.min(abs(lambda_seq - 1))
  loglik_lambda1 <- loglik_seq[idx_lambda1]
  lr_stat        <- 2 * (loglik_hat - loglik_lambda1)
  lr_p           <- 1 - pchisq(lr_stat, df = 1)
  
  # 9) Transformacja y według Box–Cox i ponowne dopasowanie modelu
  if (abs(lambda_hat) < .Machine$double.eps) {
    train_loc$Pierścienie_bc <- log(train_loc$Pierścienie)
    test_y_bc <- log(test_data$Pierścienie)
  } else {
    train_loc$Pierścienie_bc <- (train_loc$Pierścienie^lambda_hat - 1) / lambda_hat
    test_y_bc               <- (test_data$Pierścienie^lambda_hat - 1) / lambda_hat
  }
  # Przygotuj ramkę do regresji na y_bc (usuń oryginalne Pierścienie)
  df_bc             <- train_loc
  df_bc$Pierścienie <- NULL
  
  model_bc <- lm(Pierścienie_bc ~ ., data = df_bc)
  
  # 10) Reszty po transformacji i test Shapiro–Wilk
  resid_train_bc <- resid(model_bc)
  resid_test_bc  <- test_y_bc - predict(model_bc, newdata = test_data)
  
  sh_tr_bc_p <- if (length(resid_train_bc) >= 3) shapiro.test(resid_train_bc)$p.value else NA
  sh_te_bc_p <- if (length(resid_test_bc)  >= 3) shapiro.test(resid_test_bc)$p.value  else NA
  
  # 11) Zapisanie wyników iteracji
  metrics_list[[paste0("iter_", i - 1)]] <- c(
    Probe_size         = nrow(train_loc),
    R2_train           = round(r2_train, 6),
    R2_test            = round(r2_test,  6),
    Lambda_hat         = round(lambda_hat, 6),
    LR_p               = round(lr_p,         6),
    Shapiro_train_p    = round(sh_tr_p,      6),
    Shapiro_train_bc_p = round(sh_tr_bc_p,   6),
    Shapiro_test_p     = round(sh_te_p,      6),
    Shapiro_test_bc_p  = round(sh_te_bc_p,   6)
  )
  
  # 12) Jeżeli to nie jest ostatnia iteracja, usuń obserwacje odstające wg Cooka
  if (i < 7) {
    cooks_d  <- cooks.distance(model_raw)
    thresh   <- 4 / nrow(train_loc)
    keep_idx <- which(cooks_d < thresh)
    current_idx <- current_idx[keep_idx]
  }
}

# 13) Utworzenie ramki wyników
results_df <- as_tibble(do.call(cbind, metrics_list), rownames = "metric")
results_df <- column_to_rownames(results_df, var = "metric")

# 14) Wyświetlenie wyników w konsoli
print(results_df)

# (Opcjonalnie) Zapis do pliku CSV:
# TEGO NIE TRZEBA już jest wywołane 
# write.csv(results_df, file = "boxcox_leverage.csv", row.names = TRUE)

