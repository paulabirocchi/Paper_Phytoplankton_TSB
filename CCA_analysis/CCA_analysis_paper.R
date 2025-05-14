##################
# Adapted R code from Wilson Souza made by Paula Birocchi (paula.birocchi@usp.br)
##################
# LOAD LIBRARIES #
##################
library(vegan)
library(relaimpo)

########################
# DEFINE INPUTS & INFO #
########################

frequencies <- list(
  low = list(
    env_file = "sub_df_env_final_mat.csv",
    phyto_file = "sub_df_phyto_final_mat.csv",
    cca1 = 68,
    cca2 = 22
  ),
  neap_spring = list(
    env_file = "sq_df_env_final_mat.csv",
    phyto_file = "sq_df_phyto_final_mat.csv",
    cca1 = 72.8,
    cca2 = 18.4
  ),
  monthly = list(
    env_file = "month_df_env_final_mat.csv",
    phyto_file = "month_df_phyto_final_mat.csv",
    cca1 = 71.8,
    cca2 = 18.4
  )
)

##########################################
# LOOP THROUGH EACH FREQUENCY SCENARIO  #
##########################################

for (freq in names(frequencies)) {
  
  message("\n\nRunning CCA analysis for frequency: ", freq)
  
  # Extract parameters
  env_file <- frequencies[[freq]]$env_file
  phyto_file <- frequencies[[freq]]$phyto_file
  cca1_pct <- frequencies[[freq]]$cca1
  cca2_pct <- frequencies[[freq]]$cca2
  output_file <- paste0("cca_variance_explained_", freq, "_paper.png")
  
  # Load data
  env_vars <- read.csv(env_file, stringsAsFactors = FALSE)
  phyto_data <- read.csv(phyto_file, stringsAsFactors = FALSE)
  
  rownames(env_vars) <- env_vars[, 1]; env_vars <- env_vars[, -1]
  rownames(phyto_data) <- phyto_data[, 1]; phyto_data <- phyto_data[, -1]
  
  # Standardize abiotic variables
  abiotic_std <- decostand(env_vars, method = "standardize")
  
  # Replace zeros in phytoplankton with NaN
  phyto_data[phyto_data == 0] <- NaN
  complete_rows <- complete.cases(phyto_data)
  phyto_clean <- phyto_data[complete_rows, ]
  abiotic_clean <- abiotic_std[complete_rows, ]
  
  # Normalize to avoid negative values
  abiotic_clean <- abiotic_clean - min(abiotic_clean) + 1
  phyto_clean <- phyto_clean - min(phyto_clean, na.rm = TRUE) + 1
  
  # Remove colinear variable
  abiotic_clean <- abiotic_clean[, -which(names(abiotic_clean) == "Total.Sea.Level..Marinha.")]
  
  # Rename variables
  colnames(abiotic_clean) <- c("River", "Salt", "DO", "Turb", "Chla", "Temp", "CDOM", "Tide", 
                               "Subtidal", "Rain", "SolarRad", "AirTemp", "EWwind", "NSwind", 
                               "EWcurrent", "NScurrent")
  colnames(phyto_clean) <- c("CHAE", "HEMI", "GUISTRI", "RHIPRO", "RHIRO", "GUIDA", "COS")
  
  # Prepare data (remove Coscino and specific variables)
  phyto_nocos <- phyto_clean[, -which(names(phyto_clean) == "COS")]
  coscino_only <- phyto_clean[, which(names(phyto_clean) == "COS")]
  
  abiotic_n1 <- abiotic_clean[, -which(names(abiotic_clean) == "AirTemp")]
  abiotic_n2 <- abiotic_n1[, -which(names(abiotic_n1) == "Chla")]
  abiotic_n3 <- abiotic_n2[, -which(names(abiotic_n2) == "Turb")]
  abiotic_final <- abiotic_n3[, -which(names(abiotic_n3) == "Tide")] # only for HF
  
  # CCA analysis
  cca_result <- cca(phyto_nocos ~ ., data = abiotic_final)
  sp_scores <- scores(cca_result, display = "sp", scaling = 2)
  env_scores <- scores(cca_result, display = "bp", scaling = 2)
  
  # Linear model
  data_comb <- cbind(coscino_only, abiotic_final)
  set.seed(123)
  idx_train <- sample(seq_len(nrow(data_comb)), size = 0.70 * nrow(data_comb))
  train_set <- data_comb[idx_train, ]
  test_set <- data_comb[-idx_train, ]
  lm_model <- lm(coscino_only ~ ., data = train_set)
  lm_summary <- summary(lm_model)
  lm_coeffs <- lm_summary$coefficients[-1, 1]
  lm_importance <- abs(lm_coeffs) / sum(abs(lm_coeffs)) * 100
  var_names <- names(lm_coeffs)
  
  # CCA variable importance
  prop_env <- rowSums(env_scores^2) / sum(rowSums(env_scores^2)) * 100
  
  ########################
  # PLOT AND SAVE FIGURE #
  ########################
  
  png(output_file, width = 2100, height = 2970, res = 300)
  par(mfrow = c(2, 1), mar = c(4, 5, 2, 2))
  
  # CCA Plot
  plot.new()
  plot.window(xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75))
  title(main = paste0("CCA without Coscinodiscus (", freq, ")"),
        xlab = paste0("CCA1 (", cca1_pct, "%)"),
        ylab = paste0("CCA2 (", cca2_pct, "%)"))
  axis(1, at = seq(-0.75, 0.75, by = 0.25))
  axis(2, at = seq(-0.75, 0.75, by = 0.25))
  abline(h = 0, v = 0, lty = 3)
  arrows(0, 0, env_scores[,1], env_scores[,2], length = 0.1, col = "blue")
  text(env_scores[,1], env_scores[,2], labels = rownames(env_scores), col = "blue", cex = 0.8, pos = 3)
  points(sp_scores[,1], sp_scores[,2], pch = 15, col = "red", cex = 0.8)
  text(sp_scores[,1], sp_scores[,2], labels = rownames(sp_scores), col = "red", cex = 0.6, pos = 1)
  
  # Barplot: CCA vs MLR contributions
  variance_data <- rbind(prop_env, lm_importance)
  colnames(variance_data) <- var_names
  bar_positions <- barplot(
    variance_data,
    beside = TRUE,
    col = c("darkorange", "darkblue"),
    main = "Proportion Explained by Environmental Variables",
    ylab = "Proportion explained (%)",
    ylim = c(0, max(variance_data) * 1.1),
    xaxt = "n",
    legend.text = c("CCA", "MLR"),
    args.legend = list(x = "topright", bty = "n")
  )
  text(
    x = colMeans(bar_positions),
    y = par("usr")[3] - max(variance_data) * 0.05,
    labels = var_names, srt = 45, adj = 1, xpd = TRUE, cex = 0.9
  )
  
  dev.off()
  message("Saved plot: ", output_file)
}
