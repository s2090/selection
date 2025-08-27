# loading libraries and data
library(readxl)
library(metafor)
library(VBsparsePCA)
library(stringr)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
data <- read_excel("data.xlsx", sheet = "data")
data <- data[data$selectiontype != "MISSING DATA",]
data[c(34:55)] <- sapply(data[c(34:55)], as.numeric)

# adding factors and calculating variances
data$model <- paste(data$id, data$uni, data$group, data$fitness.var, sep="_")
data$effect.id <- 1:nrow(data)
data$quadratic.sel <- ifelse(data$quadratic < 0 & data$quadratic.sig == 1, "stabilising", ifelse(data$quadratic > 0 & data$quadratic.sig == 1, "disruptive", "none"))
data$quadratic.sel <- as.factor(data$quadratic.sel)
data$quadratic.sel <- relevel(data$quadratic.sel, ref = "none")
data$traittype <- as.factor(data$traittype)
data$traittype <- relevel(data$traittype, ref = "other")
data$variance <- data$se^2

# significance colours
cols <- c("0" = "black", "1" = "red")
cols <- cols[as.factor(data$optimum.sig)]
cols.quad <- c("0" = "black", "1" = "red")
cols.quad <- cols.quad[as.factor(data$quadratic.sig)]

# loading previously generated RDS files
meta.trait.quadratic <- readRDS("results/meta.trait.quadratic.RDS")
meta.fitness.quadratic <- readRDS("results/meta.fitness.quadratic.RDS")
meta.quadratic <- readRDS("results/meta.quadratic.RDS")
meta.random <- readRDS("results/meta.random.RDS")
samp_dist <- readRDS("results/samp_dist.RDS")
results <- readRDS("results/results.RDS")
diff <- readRDS("results/diff.RDS")
percentages <- readRDS("results/percentages.RDS")



# data preview with a funnel plot (standardised optimum vs inverse sampling variance)
par(mfrow=c(1,2))
plot(data$optimum, 1/data$variance, col=cols, pch=16, xlab="Standardised optimum", ylab="Inverse variance (1/V)", cex = 0.5)
abline(v=mean(data$optimum), lty=2)
abline(v=-3)
abline(v=3)
plot(data$optimum, log(1/data$variance), col=cols, pch=16, xlab="Standardised optimum", ylab="Log of inverse variance (log(1/V))", cex = 0.5)
abline(v=mean(data$optimum), lty=2)
abline(v=-3)
abline(v=3)
par(mfrow=c(1,1))

# verifying high variance for estimates of quadratic selection gradients close to 0
plot(log(abs(data$quadratic)), log(1/data$variance), col=cols.quad, pch=16, xlab="Log of the absolute standardised quadratic selection gradient (log(|γ|))", ylab="Log of inverse variance (log(1/V))", cex = 0.7)



# model without intercept, traittype:quadratic.sel first
meta.trait.quadratic <- rma.mv(yi=optimum, V=variance, mods= ~ 0 + traittype:quadratic.sel, random=list(~1|species, ~1|model, ~1|effect.id), data=data)
summary(meta.trait.quadratic)
saveRDS(meta.trait.quadratic, file = "results/meta.trait.quadratic.RDS")

# model without intercept, fitnesstype:quadratic.sel first
meta.fitness.quadratic <- rma.mv(yi=optimum, V=variance, mods= ~ 0 + fitnesstype:quadratic.sel, random=list(~1|species, ~1|model, ~1|effect.id), data=data, control=list(rel.tol=1e-8))
summary(meta.fitness.quadratic)
saveRDS(meta.fitness.quadratic, file = "results/meta.fitness.quadratic.RDS")

# model without intercept, quadratic.sel first
meta.quadratic <- rma.mv(yi=optimum, V=variance, mods= ~ 0 + quadratic.sel, random=list(~1|species, ~1|model, ~1|effect.id), data=data)
summary(meta.quadratic)
saveRDS(meta.quadratic, file = "results/meta.quadratic.RDS")

# random-effects model without moderators
meta.random <- rma.mv(yi=optimum, V=variance, random=list(~1|species, ~1|model, ~1|effect.id), data=data)
summary(meta.random)
saveRDS(meta.random, file = "results/meta.random.RDS")


# calculating heterogeneity statistics (Nakagawa and Santos, 2012)
v.typical <- (nrow(data) - 1)*sum(1/data$variance)/(sum(1/data$variance)^2-sum((1/data$variance)^2))
tau2 <- sum(meta.random$sigma2)
tau <- sqrt(tau2)
i2.species <- meta.random$sigma2[1]/(tau2 + v.typical) * 100
i2.model <- meta.random$sigma2[2]/(tau2 + v.typical) * 100
i2.effect.id <- meta.random$sigma2[3]/(tau2 + v.typical) * 100
i2.total <- i2.species + i2.model + i2.effect.id
h2 <- (tau2 + v.typical)/v.typical
het <- data.frame(v.typical, i2.total, i2.species, i2.model, i2.effect.id, tau2, tau, h2, r2)
for (i in c("meta.trait.quadratic", "meta.fitness.quadratic", "meta.quadratic")) {
  eval(parse(text = paste("tau2.res <- sum(", i, "$sigma2)", sep="")))
  tau.res <- sqrt(tau2.res)
  eval(parse(text = paste("i2.species.res <- ", i, "$sigma2[1]/(tau2.res + v.typical) * 100", sep="")))
  eval(parse(text = paste("i2.model.res <- ", i, "$sigma2[2]/(tau2.res + v.typical) * 100", sep="")))
  eval(parse(text = paste("i2.effect.id.res <- ", i, "$sigma2[3]/(tau2.res + v.typical) * 100", sep="")))
  i2.total.res <- i2.species.res + i2.model.res + i2.effect.id.res
  h2.res <- (tau2.res + v.typical)/v.typical
  r2 <- max((tau2-tau2.res)/tau2 * 100,0)
  het <- merge(het, data.frame(v.typical, "i2.total" = i2.total.res, "i2.species" = i2.species.res, "i2.model" = i2.model.res, "i2.effect.id" = i2.effect.id.res, "tau2" = tau2.res, "tau" = tau.res, "h2" = h2.res, r2), all = TRUE, sort = FALSE)
}
row.names(het) <- c("random", "traittype:quadratic", "fitnesstype:quadratic", "quadratic selection")



# folded normal transformation (Morissey, 2016)
samp_dist <- list()
results <- data.frame(factor = "a", model = "a", estimate = "a", se = 0, quantile.l = 0, quantile.u = 0, samplesize = 0)
results <- results[-1,]

# trait type:quadratic selection
estimates <- row.names(meta.trait.quadratic$beta)
samp_len <- 0
for (i in 1:13) {
  if (i > 1 & i < 13) {
    if (i %% 3 == 0) {
      i <- (i + 3)/3
    } else if ((i - 1) %% 3 == 0) {
      i <- (i + 14)/3
    } else if ((i - 2) %% 3 == 0) {
      i <- (i + 25)/3
    }
  }
  eval(parse(text = paste("samples <- nrow(subset(data, ", strsplit(estimates[i], "type")[[1]][1], "type == '", strsplit(strsplit(estimates[i], "type")[[1]][2], ":")[[1]][1], "' & quadratic.sel == '", strsplit(strsplit(estimates[i], "quadratic.sel")[[1]][2], ":")[[1]][1], "', drop = TRUE))", sep = "")))
  samp_dist[[samp_len+i]] <- replicate(1000,{
    x <- rnorm(samples, meta.trait.quadratic$beta[i], sqrt(sum(meta.trait.quadratic$sigma2)))
    foldednorm.mean(mean(x), var(x))
  })
  results <- merge(results, data.frame(factor = str_to_sentence(paste(strsplit(estimates[i], "type")[[1]][1], " type: ", strsplit(strsplit(estimates[i], "type")[[1]][2], ":")[[1]][1], " - ", strsplit(strsplit(estimates[i], "quadratic.sel")[[1]][2], "ne")[[1]][1], " selection", sep="")), model = "Trait type:quadratic", estimate = mean(samp_dist[[samp_len+i]]), se = sd(samp_dist[[samp_len+i]]), quantile.l = quantile(samp_dist[[samp_len+i]], c(0.025,0.975))[[1]], quantile.u = quantile(samp_dist[[samp_len+i]], c(0.025,0.975))[[2]], samplesize = samples), all = TRUE, sort = FALSE)
}

# fitness type:quadratic selection
estimates <- row.names(meta.fitness.quadratic$beta)
samp_len <- length(samp_dist)
for (i in 1:9) {
  i <- i + 2 * (i - 1) - floor((i - 1)/3) * 8
  eval(parse(text = paste("samples <- nrow(subset(data, ", strsplit(estimates[i], "type")[[1]][1], "type == '", strsplit(strsplit(estimates[i], "type")[[1]][2], ":")[[1]][1], "' & quadratic.sel == '", strsplit(estimates[i], "quadratic.sel")[[1]][2], "', drop = TRUE))", sep="")))
  samp_dist[[samp_len+i]] <- replicate(1000,{
    x <- rnorm(samples, meta.fitness.quadratic$beta[i], sqrt(sum(meta.fitness.quadratic$sigma2)))
    foldednorm.mean(mean(x), var(x))
  })
  results <- merge(results, data.frame(factor = str_to_sentence(paste(strsplit(estimates[i], "type")[[1]][1], " type: ", strsplit(strsplit(estimates[i], "type")[[1]][2], ":")[[1]][1], " - ", strsplit(strsplit(estimates[i], "quadratic.sel")[[1]][2], "ne")[[1]][1], " selection", sep="")), model = "Fitness type:quadratic", estimate = mean(samp_dist[[samp_len+i]]), se = sd(samp_dist[[samp_len+i]]), quantile.l = quantile(samp_dist[[samp_len+i]], c(0.025,0.975))[[1]], quantile.u = quantile(samp_dist[[samp_len+i]], c(0.025,0.975))[[2]], samplesize = samples), all = TRUE, sort = FALSE)
}

# quadratic selection
estimates <- row.names(meta.quadratic$beta)
samp_len <- length(samp_dist)
for (i in 1:3) {
  eval(parse(text = paste("samples <- nrow(subset(data, quadratic.sel == '", strsplit(estimates[i], "quadratic.sel")[[1]][2], "', drop = TRUE))", sep="")))
  samp_dist[[samp_len+i]] <- replicate(1000,{
    x <- rnorm(samples, meta.quadratic$beta[i], sqrt(sum(meta.quadratic$sigma2)))
    foldednorm.mean(mean(x), var(x))
  })
  results <- merge(results, data.frame(factor = str_to_sentence(paste("quadratic selection - ", strsplit(estimates[i], "quadratic.sel")[[1]][2], sep="")), model = "Quadratic selection", estimate = mean(samp_dist[[samp_len+i]]), se = sd(samp_dist[[samp_len+i]]), quantile.l = quantile(samp_dist[[samp_len+i]], c(0.025,0.975))[[1]], quantile.u = quantile(samp_dist[[samp_len+i]], c(0.025,0.975))[[2]], samplesize = samples), all = TRUE, sort = FALSE)
}

# random-effects model
samples <- nrow(data)
samp_dist[[length(samp_dist)+1]] <- replicate(1000, {
  x <- rnorm(samples, meta.random$beta[1], sqrt(sum(meta.random$sigma2)))
  foldednorm.mean(mean(x), var(x))
})
results <- merge(results, data.frame(factor = "Overall estimate", model = "Random-effects model", estimate = mean(samp_dist[[length(samp_dist)]]), se = sd(samp_dist[[length(samp_dist)]]), quantile.l = quantile(samp_dist[[length(samp_dist)]], c(0.025,0.975))[[1]], quantile.u = quantile(samp_dist[[length(samp_dist)]], c(0.025,0.975))[[2]], samplesize = samples), all = TRUE, sort = FALSE)
results$factor <- str_replace_all(results$factor, "no selection", "no quadratic selection")

# p-values
results$p.value <- vector(mode="numeric", length=length(results$estimate))
results$p.value.full <- vector(mode="numeric", length=length(results$estimate))
for (i in 1:length(results$estimate)) {
  results$p.value.full[i] <- pnorm((mean(samp_dist[[i]]) - 3)/sd(samp_dist[[i]]))
  if (results$p.value.full[i] < 0.00001) {
    results$p.value[i] <- "< 0.00001"
  } else {
    results$p.value[i] <- as.character(format(p.value.full[i], scientific = FALSE))
  }
}

saveRDS(samp_dist, file="results/samp_dist.RDS")
saveRDS(results, file="results/results.RDS")



# forest plot
forest(results$estimate, sei = results$se, ci.lb = results$quantile.l, ci.ub = results$quantile.u, refline = 3, xlim = c(-4.8,5.6), alim = c(0,3), slab = results$model, ilab = data.frame(results$factor, results$p.value, results$samplesize), ilab.lab = c("Moderator level", "p-value", "N"), ilab.pos = c(4,2,4), ilab.xpos = c(-3.2, 3.9, 3.8), header = "Model", xlab = "Absolute standardised optimum", pch=16, shade="zebra")



# qualitative analysis
all <- nrow(data)
directional.positive <- nrow(subset(data, linear.sig == 1 & linear > 0, drop = TRUE))
directional.negative <- nrow(subset(data, linear.sig == 1 & linear < 0, drop = TRUE))
directional <- directional.positive + directional.negative
quadratic.optima <- nrow(subset(data, selectiontype == "stabilising" | selectiontype == "disruptive", drop = TRUE))
stabilising <- nrow(subset(data, selectiontype == "stabilising", drop = TRUE))
disruptive <- nrow(subset(data, selectiontype == "disruptive", drop = TRUE))
no.selection <- nrow(subset(data, selectiontype == "NA" | selectiontype == "not in range", drop = TRUE))
notinrange <- nrow(subset(data, quadratic.sig == 1 & optimum.sig == 0, drop = TRUE))
quadratic <- nrow(subset(data, quadratic.sig == 1, drop = TRUE))
quadratic.negative <- nrow(subset(data, quadratic < 0 & quadratic.sig == 1, drop = TRUE))
quadratic.positive <- nrow(subset(data, quadratic > 0 & quadratic.sig == 1, drop = TRUE))

# testing difference between stabilising and disruptive selection frequency
binom.test(stabilising, n = all, p = 0.5*(quadratic.optima/all), alternative = "two.sided")
# testing difference between quadratic and directional selection frequency
binom.test(quadratic.optima, n = all, p = 0.5*((quadratic.optima + directional)/all), alternative = "two.sided")
# testing difference in results if not checking optima
chisq.test(data.frame(c(stabilising, disruptive, all - quadratic.optima), c(quadratic.negative, quadratic.positive, all - quadratic)))
chisq.test(data.frame(c(quadratic.optima, all - quadratic.optima), c(quadratic, all - quadratic)))

# selection type percentages
percentages <- data.frame(proportion = "stabilising selection", numerator = "stabilising selection", n1 = stabilising, denominator = "all estimates", n2 = all, percentage = stabilising/all * 100)
percentages <- rbind(percentages, c("disruptive selection", "disruptive selection", disruptive, "all estimates", all, disruptive/all * 100))
percentages <- rbind(percentages, c("quadratic selection (with optimum in range)", "quadratic selection (with optimum in range)", quadratic.optima, "all estimates", all, quadratic.optima/all * 100))
percentages <- rbind(percentages, c("directional selection", "directional selection", directional.positive + directional.negative, "all estimates", all, (directional.positive + directional.negative)/all * 100))
percentages <- rbind(percentages, c("positive directional selection", "positive directional selection", directional.positive, "all estimates", all, directional.positive/all * 100))
percentages <- rbind(percentages, c("negative directional selection", "negative directional selection", directional.negative, "all estimates", all, directional.negative/all * 100))
percentages <- rbind(percentages, c("not in range", "not in range", notinrange, "all estimates", all, notinrange/all * 100))
percentages <- rbind(percentages, c("negative quadratic selection", "negative quadratic selection", quadratic.negative, "all estimates", all, quadratic.negative/all * 100))
percentages <- rbind(percentages, c("positive quadratic selection", "positive quadratic selection", quadratic.positive, "all estimates", all, quadratic.positive/all * 100))
percentages <- rbind(percentages, c("quadratic selection", "significant quadratic selection", quadratic, "all estimates", all, quadratic/all * 100))
percentages <- rbind(percentages, c("stabilising selection (of quadratic with optimum in range)", "stabilising selection", stabilising, "quadratic selection with optimum in range", quadratic.optima, stabilising/quadratic.optima * 100))
percentages <- rbind(percentages, c("disruptive selection (of quadratic with optimum in range)", "disruptive selection", disruptive, "quadratic selection with optimum in range", quadratic.optima, disruptive/quadratic.optima * 100))
percentages <- rbind(percentages, c("negative quadratic selection gradients (of quadratic)", "negative quadratic selection gradients", quadratic.negative, "quadratic selection", quadratic, quadratic.negative/quadratic * 100))
percentages <- rbind(percentages, c("positive quadratic selection gradients (of quadratic)", "positive quadratic selection gradients", quadratic.positive, "quadratic selection", quadratic, quadratic.positive/quadratic * 100))
percentages <- rbind(percentages, c("additional stabilising selection (-100%)", "negative quadratic selection", quadratic.negative, "stabilising selection", stabilising, (quadratic.negative/stabilising - 1) * 100))
percentages <- rbind(percentages, c("additional disruptive selection (-100%)", "positive quadratic selection", quadratic.positive, "disruptive selection", disruptive, (quadratic.positive/disruptive - 1) * 100))
percentages <- rbind(percentages, c("additional quadratic selection (-100%)", "quadratic selection", quadratic, "quadratic selection (with optimum in range)", quadratic.optima, (quadratic/quadratic.optima - 1) * 100))
# selection type mismatch with or missed by studies that did not account for nonlinearity, estimate quadratic selection or check for optima in the range of the data
percentages <- rbind(percentages, c("reported selection mismatch (untransformed)", "mismatch in reported selection among untransformed estimates", nrow(subset(data, transformed == 0 & reported != selectiontype.matching & !is.na(reported), drop = TRUE)), "reported untransformed estimates", nrow(subset(data, transformed == 0 & !is.na(reported), drop = TRUE)), nrow(subset(data, transformed == 0 & reported != selectiontype.matching & !is.na(reported), drop = TRUE))/nrow(subset(data, transformed == 0 & !is.na(reported), drop = TRUE)) * 100))
percentages <- rbind(percentages, c("quadratic selection missed", "quadratic selection for studies that only examined linear", nrow(subset(data, quantified == "linear" & (selectiontype == "stabilising" | selectiontype == "disruptive") & !is.na(reported), drop = TRUE)), "quadratic selection where selection was reported", nrow(subset(data, (selectiontype == "stabilising" | selectiontype == "disruptive") & !is.na(reported), drop = TRUE)), nrow(subset(data, quantified == "linear" & (selectiontype == "stabilising" | selectiontype == "disruptive") & !is.na(reported), drop = TRUE))/nrow(subset(data, (selectiontype == "stabilising" | selectiontype == "disruptive") & !is.na(reported), drop = TRUE)) * 100))
# add transformed == 1?
percentages <- rbind(percentages, c("mismatch from unchecked optimum (untransformed)", "reported quadratic selection where optimum was not checked and not in range (untransformed)", nrow(subset(data, checked.optima == 0 & transformed == 0 & (reported == "stabilising" | reported == "disruptive") & optimum.sig == 0, drop = TRUE)), "reported quadratic selection where optimum was not checked (untransformed)", nrow(subset(data, checked.optima == 0 & transformed == 0 & (reported == "stabilising" | reported == "disruptive"), drop = TRUE)), nrow(subset(data, checked.optima == 0 & transformed == 0 & (reported == "stabilising" | reported == "disruptive") & optimum.sig == 0, drop = TRUE))/nrow(subset(data, checked.optima == 0 & transformed == 0 & (reported == "stabilising" | reported == "disruptive"), drop = TRUE)) * 100))
percentages <- rbind(percentages, c("mismatch from unchecked optimum (transformed)", "reported quadratic selection where optimum was not checked and not in range (transformed)", nrow(subset(data, checked.optima == 0 & (transformed == 1 | transformed == "NA") & (reported == "stabilising" | reported == "disruptive") & optimum.sig == 0, drop = TRUE)), "reported quadratic selection where optimum was not checked (transformed)", nrow(subset(data, checked.optima == 0 & transformed == 1 & (reported == "stabilising" | reported == "disruptive"), drop = TRUE)), nrow(subset(data, checked.optima == 0 & transformed == 1 & (reported == "stabilising" | reported == "disruptive") & optimum.sig == 0, drop = TRUE))/nrow(subset(data, checked.optima == 0 & transformed == 1 & (reported == "stabilising" | reported == "disruptive"), drop = TRUE)) * 100))
percentages <- rbind(percentages, c("mismatch from unchecked optimum", "reported quadratic selection where optimum was not checked and not in range", nrow(subset(data, checked.optima == 0 & (reported == "stabilising" | reported == "disruptive") & optimum.sig == 0, drop = TRUE)), "reported quadratic selection where optimum was not checked", nrow(subset(data, checked.optima == 0 & (reported == "stabilising" | reported == "disruptive"), drop = TRUE)), nrow(subset(data, checked.optima == 0 & (reported == "stabilising" | reported == "disruptive") & optimum.sig == 0, drop = TRUE))/nrow(subset(data, checked.optima == 0 & (reported == "stabilising" | reported == "disruptive"), drop = TRUE)) * 100))

saveRDS(percentages, "results/percentages.RDS")
