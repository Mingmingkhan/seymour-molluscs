load(file = "results/klb7.metacom.results.RData")
load(file = "results/klb8.metacom.results.RData")
load(file = "results/klb9.metacom.results.RData")

#plot metacom Res 
coh.res.7 <- (klb7.metacom.res.df[2, c(2, 8, 14, 20)])
turn.res.7 <- (klb7.metacom.res.df[2, c(4, 10, 16, 22)])
bc.res.7 <- (klb7.metacom.res.df[1, c(6, 12, 18, 24)])

coh.res.8 <- (klb8.metacom.res.df[2, c(2, 8, 14, 20)])
turn.res.8 <- (klb8.metacom.res.df[2, c(4, 10, 16, 22)])
bc.res.8 <- (klb8.metacom.res.df[1, c(6, 12, 18, 24)])

coh.res.9 <- (klb9.metacom.res.df[2, c(2, 8, 14)])
turn.res.9 <- (klb9.metacom.res.df[2, c(4, 10, 16)])
bc.res.9 <- (klb9.metacom.res.df[1, c(6, 12, 18)])

x.pol <- c(0, 0, 5, 5)
y.pol <- c(-11, 0, 0, -11)

plot(xlim = c(-3, 4), ylim = c(-10, -1), 
     xlab = "Turnover Z-score", ylab = "Coherence Z-score",
     x = -3:2, y = -3:2, col = "white", 
     main = "All Metacommunity Results")

abline(v = 0, lwd = 2)
box()
text(-1.1, -2.0, "Nested Clumped", cex = 1.2)
text(0.8, -2.0, "Clementsian", cex = 1.2, col = "black")

#klb9
points(turn.res.9$raw.turnover_res, coh.res.9$raw.coherence_res, 
       pch = 21, cex = bc.res.9$raw.boundary_clumping_res, 
       bg = "#cce86eff", lwd = 2, col = "black")

#expand
points(turn.res.9$treat1.turnover_res, coh.res.9$treat1.coherence_res, 
       pch = 22, cex = bc.res.9$treat1.boundary_clumping_res, 
       bg = "#cce86eff", lwd = 2, col = "black")

#contract
points(turn.res.9$treat2.turnover_res, coh.res.9$treat2.coherence_res, 
       pch = 25, cex = bc.res.9$treat2.boundary_clumping_res,
       bg = "#cce86eff", lwd = 2)
#text(3.5, -8.0, klb9.metacom.outputs[3])

#klb8
points(turn.res.8$raw.turnover_res, coh.res.8$raw.coherence_res, 
       pch = 21, cex = bc.res.8$raw.boundary_clumping_res, 
       bg = "#a6da5aff", lwd = 2, col = "black")

#left
points(turn.res.8$treat1.turnover_res, coh.res.8$treat1.coherence_res, 
       pch = 0, cex = bc.res.8$treat1.boundary_clumping_res,
       col = "#a6da5aff", lwd = 2)

#right
points(turn.res.8$treat2.turnover_res, coh.res.8$treat2.coherence_res, 
       pch = 2, cex = bc.res.8$treat2.boundary_clumping_res,
       col = "#a6da5aff", lwd = 2)

#contract
points(turn.res.8$treat4.turnover_res, coh.res.8$treat4.coherence_res, 
       pch = 25, cex = bc.res.8$treat4.boundary_clumping_res,
       bg = "#a6da5aff", lwd = 2)

#klb7 
points(turn.res.7$raw.turnover_res, coh.res.7$raw.coherence_res, 
       pch = 21, cex = bc.res.7$raw.boundary_clumping_res, 
       bg = "#00a650ff", lwd = 2, col = "black")

#left
points(turn.res.7$treat1.turnover_res, coh.res.7$treat1.coherence_res, 
       pch = 0, cex = bc.res.7$treat1.boundary_clumping_res,
       col = "#00a650ff", lwd = 2)
#right
points(turn.res.7$treat2.turnover_res, coh.res.7$treat2.coherence_res, 
       pch = 2, cex = bc.res.7$treat2.boundary_clumping_res,
       col = "#00a650ff", lwd = 2)

#contract
points(turn.res.7$treat4.turnover_res, coh.res.7$treat4.coherence_res, 
       pch = 25, cex = bc.res.7$treat4.boundary_clumping_res,
       bg = "#00a650ff", lwd = 2)
#text(-2.5, -2.8, klb7.metacom.outputs[5])

legend("topright", bg = "white", 
       legend = c("raw", "left.shift", 
                  "right.shift", 
                  "contract", 
                  "KLB 7", "KLB 8", 
                  "KLB 9", 
                  "significant structure",
                  "quasi structure"), 
       lwd = 1, pch = c(NA, NA, NA, NA,
                        0, 1, 2, 16, 1), 
       col = c("#007FFF","#6B990F", "#FFC34C",
              "#860086", 
                "black", "black",
               "black","black"), 
       pt.cex = 0.8)
