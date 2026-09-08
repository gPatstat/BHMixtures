dom=seq(0,1,length.out=100)
w=dom[2]-dom[1]

nsim <- 50
seed0 <- 03072

dir3 <- "C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/Sim_A"
data_box3 <- numeric()
param_a <- 0.01
for (sim_id in 1:nsim) {
  seed <- sim_id * 100000 + seed0
  fname <- file.path(dir3, paste0("simB_sim ", sim_id, " par ", param_a, " seed ", seed, " .rdata"))
  tryCatch({
    load(fname)
    if (!is.complex(pca_basis)) {
      mat1 <- pca_basis %*% EM_sample$H
      mat2 <- get_basis %*% F_sample$pF$H$coefs
      data_box3 <- rbind(data_box3, c(error_vertex(mat1, mat2), 3))
    }
  }, error = function(e) {
    message("ERROR (3 vertici) sim ", sim_id, ": ", conditionMessage(e))
  })
}
cat("Simulazioni valide a 3 vertici:", nrow(data_box3), "/", nsim, "\n")

## ---- Caso 5 e 8 vertici: Risultati_sim_C ----
dirC <- "C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/Sim_C"
data_box58 <- numeric()
for (m in c(5, 8)) {
  for (sim_id in 1:nsim) {
    seed <- seed0 + sim_id * 100000
    fname <- file.path(dirC, paste0("simC_BHM_sim", sim_id, "par", m, "seed", seed, ".rdata"))
    tryCatch({
      load(fname)
      if (!is.complex(pca_basis)) {
        mat1 <- pca_basis %*% EM_sample$H
        mat2 <- get_basis %*% F_sample$pF$H$coefs
        data_box58 <- rbind(data_box58, c(error_vertex(mat1, mat2), m))
      }
    }, error = function(e) {
      message("ERROR (", m, " vertici) sim ", sim_id, ": ", conditionMessage(e))
    })
  }
}
cat("Simulazioni valide 5/8 vertici:", nrow(data_box58), "\n")
print(table(data_box58[,2]))

## ---- Unione e plot ----
data_box_all <- rbind(data_box3, data_box58)
colnames(data_box_all) <- c("error", "n_vertices")
data_box_all <- as.data.frame(data_box_all)
data_box_all$n_vertices <- factor(data_box_all$n_vertices, levels = c(3,5,8))

saveRDS(data_box_all, "data_box_all.rds")

png("vertex_error_boxplot.png", width = 1600, height = 1200, res = 200)
x11()
boxplot(error ~ n_vertices, data = data_box_all, outline = FALSE,
        xlab = "m", ylab = "L2 error - best vertex permutation",
        main = "Simulation C",
        col = c("#A6CEE3", "#B2DF8A", "#FB9A99"))
dev.off()

cat("Fatto.\n")
