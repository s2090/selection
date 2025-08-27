# loading libraries
library(readxl)
library(readODS)
library(stringr)
library(coda)
library(retry)
library(MCMCglmm)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# study data loading function
loaddata <- function(id) {
  error = 0
  file <- paste("data", data.id$folder[data.id$id == id], data.id$file[data.id$id == id], sep="/")
  sheet <- data.id$sheet[data.id$id == id]
  if (endsWith(file, ".csv")) {
    if (isTRUE(grepl(";", readLines(file, n = 2)[2]))) {
      dat <- as.data.frame(read.csv2(file))
    } else {
      dat <- as.data.frame(read.csv(file))
    }
  } else if (endsWith(file, ".xlsx") | endsWith(file, ".xls")) {
    dat <- as.data.frame(read_excel(file, sheet = sheet))
  } else if (endsWith(file, ".ods")) {
    dat <- as.data.frame(read_ods(file, sheet = sheet))
  } else if (endsWith(file, ".txt")) {
    if (isTRUE(grepl(";", readLines(file, n = 2)[2])) & isTRUE(grepl(",", readLines(file, n = 2)[2]))) {
      dat <- as.data.frame(read.table(file, header = TRUE, dec = ","))
    } else {
      dat <- as.data.frame(read.table(file, header = TRUE, dec = "."))
    }
  } else {
    error <- 1
    print("Error: File does not end with .csv, .xlsx, .xls, .ods or .txt.")
  }
  if (error != 1) {
    return(dat)
  }
}

# model running function using information from data sheet
runmod <- function(id) {
  dat <- loaddata(id)
  info <- data[data$id == id,]
  
  output <- data.frame(id = 1, study = "a", group = "a", fitness = "a", trait = "a", optimum = 1, sd = 1, optimumtest = 1, samplesize = 1, linear = 1, linear.se = 1, linear.95CI.l = 1, linear.95CI.u = 1, linear.effsamp = 1, linear.pMCMC = 1, quadratic = 1, quadratic.se = 1, quadratic.95CI.l = 1, quadratic.95CI.u = 1, quadratic.effsamp = 1, quadratic.pMCMC = 1)
  output <- output[-1,]
  
  for (q in sapply(list(info$uni), unique)) {
    if (is.na(q)) {
      info.mod <- info[is.na(info$uni),]
    } else {
      info.mod <- info[info$uni == q,]
    }
    if (!all(is.na(sapply(list(info.mod$group), unique)))) {
      groups.list <- list(info.mod$group)
      groups.col.list <- list(info.mod$group.col)
      groups <- as.list(unique(data.frame(groups = groups.list[[1]], groups.col = groups.col.list[[1]]))[,"groups"])
      groups.col <- as.list(unique(data.frame(groups = groups.list[[1]], groups.col = groups.col.list[[1]]))[,"groups.col"])
    } else {
      groups <- as.list("")
      groups.col <- as.list("")
    }
    
    for (j in 1:length(groups)) {
      if (groups[j] != "" & !is.na(groups[j])) {
        info.group <- eval(parse(text = paste("subset(info.mod, group == '", paste(groups[j][[1]], collapse="; "), "', drop = TRUE)", sep="")))
      } else if (!all(!is.na(groups[j]))) {
        info.group <- subset(info.mod, is.na(group), drop = TRUE)
      } else {
        info.group <- info.mod
      }
      dat.group <- dat
      if (all(!is.na(sapply(list(info.group$group), unique)))) {
        if (str_detect(groups[j], "; ")) {
          groups[j] <- as.list(strsplit(groups[j][[1]], "; "))
          groups.col[j] <- as.list(strsplit(groups.col[j][[1]], "; "))
        }
        group <- ""
        for (i in 1:length(groups[j][[1]])) {
          dat.group[,groups.col[j][[1]][i]] <- as.character(dat.group[,groups.col[j][[1]][i]])
          group <- paste(group, "dat.group$`", groups.col[j][[1]][i], "` == '", groups[j][[1]][i], "'", sep="")
          if (i != length(groups[j][[1]])) {
            group <- paste(group, " & ", sep="")
          } else {
            eval(parse(text = paste(" dat.group[", group, ",]", sep="")))
          }
        }
      }
      
      fitness <- sapply(list(na.omit(info.group$fitness.var)), unique)
      for (k in 1:length(fitness)) {
        info.group.fit <- eval(parse(text = paste("subset(info.group, fitness.var == '", fitness[k], "', drop = TRUE)", sep="")))
        dat.group[,fitness[k]] <- as.numeric(unlist(dat.group[,fitness[k]]))
        family <- info.group.fit$family[1]
        dat.group.fit <- dat.group[!is.na(dat.group[,fitness[k]]),]
        if (!("none" %in% info.group.fit$link[1])) {
          if (!is.na(info.group.fit$link[1])) {
            eval(parse(text = paste("dat.group.fit[,fitness[k]] <- ", str_replace_all(info.group.fit$link[1], "x", "(as.numeric(unlist(dat.group.fit[,fitness[k]])))"), sep="")))
          } else if (family == "gaussian") {
            dat.group.fit[,fitness[k]] <- dat.group.fit[,fitness[k]]/mean(dat.group.fit[,fitness[k]])
          }
        }
        traits.base <- sapply(list(na.omit(info.group.fit$trait.var)), unique)
        traits <- ""
        for (i in traits.base) {
          trait <- as.list(strsplit(i, ":"))[[1]]
          for (l in trait) {
            dat.group.fit[,l] <- as.numeric(dat.group.fit[,l])
            stand <- info.group.fit$standardisation[info.group.fit$trait.var == i][1]
            if (!is.na(stand)) {
              if (stand != "none") {
                for (m in sapply(list(na.omit(dat.group.fit[,stand])), unique)) {
                  dat.group.fit[,l][dat.group.fit[,stand] == m] <- as.data.frame(scale(unlist(dat.group.fit[,l][dat.group.fit[,stand] == m])))[[1]]
                }
              }
            } else {
              stand <- ""
              dat.group.fit[[l]] <- as.data.frame(scale(unlist(dat.group.fit[[l]])))[[1]]
            }
            dat.group.fit <- dat.group.fit[!is.na(dat.group.fit[,l]),]
            if (length(trait) > 1) {
              if (l == trait[length(trait)]) {
                for (m in 1:length(trait)) {
                  if (m == length(trait)) {
                    traits <- paste(traits, "`", trait[m], "` + ", sep="")
                  } else {
                    traits <- paste(traits, "`", trait[m], "` * ", sep="")
                  }
                }
                for (m in 1:length(trait)) {
                  if (m == length(trait)) {
                    traits <- paste(traits, "I(`", trait[m], "`^2) + ", sep="")
                  } else {
                    traits <- paste(traits, "I(`", trait[m], "`^2) * ", sep="")
                  }
                }
              }
            } else {
              traits <- paste(traits, "`", l, "` + I(`", l, "`^2) + ", sep="")
            }
          }
        }
        
        fixed <- ""
        if (!is.na(info.group.fit$fixed[1])) {
          for (i in strsplit(info.group.fit$fixed[1], " \\+ ")[[1]]) {
            for (m in 1:length(strsplit(i, " \\* ")[[1]])) {
              n <- strsplit(strsplit(i, " \\* ")[[1]][m], "as.factor\\(")[[1]]
              if (n[[1]] == "") {
                n <- substr(n[[2]], 1, nchar(n[[2]]) - 1)
                dat.group.fit[,n] <- as.factor(unlist(dat.group.fit[,n]))
              } else {
                n <- n[[1]]
                dat.group.fit[,n] <- as.numeric(unlist(dat.group.fit[,n]))
              }
              dat.group.fit <- dat.group.fit[!is.na(dat.group.fit[,n]),]
              if (length(strsplit(i, " \\* ")[[1]]) > 1 & m > 1) {
                fixed <- paste(fixed, " * `", n, "`", sep="")
              } else {
                fixed <- paste(fixed, " + `", n, "`", sep="")
              }
            }
          }
        }
        
        if (!is.na(info.group.fit$random[1])) {
          p <- 1
          for (i in 1:length(strsplit(info.group.fit$random[1], " \\+ ")[[1]])) {
            for (l in 1:length(strsplit(strsplit(info.group.fit$random[1], " \\+ ")[[1]][i], "\\/")[[1]])) {
              n <- strsplit(strsplit(strsplit(info.group.fit$random[1], " \\+ ")[[1]][i], "\\/")[[1]][l], "as.factor\\(")[[1]]
              if (n[[1]] == "") {
                n <- substr(n[[2]], 1, nchar(n[[2]]) - 1)
                dat.group.fit[,n] <- as.factor(unlist(dat.group.fit[,n]))
              } else {
                n <- n[[1]]
                dat.group.fit[,n] <- as.numeric(unlist(dat.group.fit[,n]))
              }
              dat.group.fit <- dat.group.fit[!is.na(dat.group.fit[,n]),]
              if (length(strsplit(strsplit(info.group.fit$random[1], " \\+ ")[[1]][i], "\\/")[[1]]) > 1) {
                if (l == 1) {
                  eval(parse(text = paste("dat.group.fit['random", i, "'] <- dat.group.fit[,n]", sep="")))
                } else {
                  eval(parse(text = paste("dat.group.fit$random", i, " <- paste0(dat.group.fit$random", i, ", dat.group.fit[,n])", sep="")))
                }
                if (l == length(strsplit(strsplit(info.group.fit$random[1], " \\+ ")[[1]][i], "\\/")[[1]])) {
                  n <- paste("random", i, sep="")
                }
              }
            }
            if (p == 1) {
              random <- paste(", random = ~ `", n, "`", sep="")
              p <- 0
            } else {
              random <- paste(random, " + `", n, "`", sep="")
            }
          }
        } else {
          random <- ""
        }
        
        if (nrow(dat.group.fit) > 0) {
          model <- paste("`", fitness[k], "` ~ ", substr(traits, 1, nchar(traits) - 3), fixed, random, sep="")
          mcmc <- paste("MCMCglmm(", model, ", data = dat.group.fit, family = '", family, "', verbose = FALSE)", sep="")
          cat(paste("\n\nid: ", id, sep=""))
          if (groups[j][[1]][1] != "" & !is.na(groups[j][[1]][1])) {
            cat(paste(", ", "group: ", sep=""))
            cat(paste(as.character(groups[j][[1]]), collapse=" & "))
          }
          cat(paste("\n\n", mcmc, sep=""))
          retry(mod <- eval(parse(text = mcmc)), when = "Mixed model equations singular", silent = FALSE)
          sum <- summary(mod)
          print(sum)
          samplesize <- nrow(dat.group.fit)
          cat(paste("\nSample size: ", samplesize, sep=""))
          
          square <- sqrt(length(traits.base))
          if (((trunc(square)) < (square)) && ((square) < (trunc(square) + 0.5))) {
            par(mfrow = c(trunc(square), trunc(square) + 1), oma = c(1, 1, 2, 1), mar = c(2, 4, 3, 1))
          } else {
            par(mfrow = c(ceiling(square), ceiling(square)), oma = c(1, 1, 2, 1), mar = c(2, 4, 3, 1))
          }
          for (i in 1:length(traits.base)) {
            if (traits.base[i] %in% colnames(mod$Sol) & paste("I(", traits.base[i], "^2)", sep="") %in% colnames(mod$Sol)) {
              ticks = 1
            } else if (paste("`", traits.base[i], "`", sep="") %in% colnames(mod$Sol) & paste("I(`", traits.base[i], "`^2)", sep="") %in% colnames(mod$Sol)) {
              ticks = 2
            } else {
              ticks = 0
            }
            if (ticks > 0) {
              trait <- as.list(strsplit(traits.base[i], ":"))[[1]]
              if (length(trait) > 1) {
                multi <- ""
                for (m in 1:length(trait)) {
                  if (ticks == 1) {
                    if (m == length(trait)) {
                      multi <- paste(multi, "I(", trait[m], "^2)", sep="")
                    } else {
                      multi <- paste(multi, "I(", trait[m], "^2):", sep="")
                    }
                  } else {
                    if (m == length(trait)) {
                      multi <- paste(multi, "I(`", trait[m], "`^2)", sep="")
                    } else {
                      multi <- paste(multi, "I(`", trait[m], "`^2):", sep="")
                    }
                  }
                }
                if (ticks == 1) {
                  optimum <- (-1 * mod$Sol[,traits.base[i]]) / (2 * mod$Sol[,multi])
                } else {
                  optimum <- (-1 * mod$Sol[,paste("`", traits.base[i], "`", sep="")]) / (2 * mod$Sol[,multi])
                }
                median <- median(optimum)
                sd <- sd(optimum)
                within <- mean(-3 < optimum & 3 > optimum)
                if (ticks == 1) {
                  output <- merge(output, data.frame(id = info.group.fit$id[i], study = info.group.fit$authors[i], group = paste(groups[j][[1]], collapse="; "), fitness = info.group.fit$fitness[i], trait = traits.base[i], optimum = median, sd = sd, optimumtest = within, samplesize = samplesize, linear = sum$solutions[traits.base[i],1], linear.se = sd(mod$Sol[,traits.base[i]]), linear.95CI.l = sum$solutions[traits.base[i],2], linear.95CI.u = sum$solutions[traits.base[i],3], linear.effsamp = sum$solutions[traits.base[i],4], linear.pMCMC = sum$solutions[traits.base[i],5], quadratic = sum$solutions[multi,1], quadratic.se = sd(mod$Sol[,multi]), quadratic.95CI.l = sum$solutions[multi,2], quadratic.95CI.u = sum$solutions[multi,3], quadratic.effsamp = sum$solutions[multi,4], quadratic.pMCMC = sum$solutions[multi,5]), all = TRUE, sort = FALSE)
                } else {
                  output <- merge(output, data.frame(id = info.group.fit$id[i], study = info.group.fit$authors[i], group = paste(groups[j][[1]], collapse="; "), fitness = info.group.fit$fitness[i], trait = traits.base[i], optimum = median, sd = sd, optimumtest = within, samplesize = samplesize, linear = sum$solutions[traits.base[i],1], linear.se = sd(mod$Sol[,paste("`", traits.base[i], "`", sep="")]), linear.95CI.l = sum$solutions[traits.base[i],2], linear.95CI.u = sum$solutions[traits.base[i],3], linear.effsamp = sum$solutions[traits.base[i],4], linear.pMCMC = sum$solutions[traits.base[i],5], quadratic = sum$solutions[multi,1], quadratic.se = sd(mod$Sol[,multi]), quadratic.95CI.l = sum$solutions[multi,2], quadratic.95CI.u = sum$solutions[multi,3], quadratic.effsamp = sum$solutions[multi,4], quadratic.pMCMC = sum$solutions[multi,5]), all = TRUE, sort = FALSE)
                }
              } else {
                if (ticks == 1) {
                  optimum <- (-1 * mod$Sol[,traits.base[i]]) / (2 * mod$Sol[,paste("I(", traits.base[i], "^2)", sep="")])
                } else {
                  optimum <- (-1 * mod$Sol[,paste("`", traits.base[i], "`", sep="")]) / (2 * mod$Sol[,paste("I(`", traits.base[i], "`^2)", sep="")])
                }
                median <- median(optimum)
                sd <- sd(optimum)
                within <- mean(-3 < optimum & 3 > optimum)
                if (ticks == 1) {
                  output <- merge(output, data.frame(id = info.group.fit$id[i], study = info.group.fit$authors[i], group = paste(groups[j][[1]], collapse="; "), fitness = info.group.fit$fitness[i], trait = traits.base[i], optimum = median, sd = sd, optimumtest = within, samplesize = samplesize, linear = sum$solutions[traits.base[i],1], linear.se = sd(mod$Sol[,traits.base[i]]), linear.95CI.l = sum$solutions[traits.base[i],2], linear.95CI.u = sum$solutions[traits.base[i],3], linear.effsamp = sum$solutions[traits.base[i],4], linear.pMCMC = sum$solutions[traits.base[i],5], quadratic = sum$solutions[paste("I(", traits.base[i], "^2)", sep=""),1], quadratic.se = sd(mod$Sol[,paste("I(", traits.base[i], "^2)", sep="")]), quadratic.95CI.l = sum$solutions[paste("I(", traits.base[i], "^2)", sep=""),2], quadratic.95CI.u = sum$solutions[paste("I(", traits.base[i], "^2)", sep=""),3], quadratic.effsamp = sum$solutions[paste("I(", traits.base[i], "^2)", sep=""),4], quadratic.pMCMC = sum$solutions[paste("I(", traits.base[i], "^2)", sep=""),5]), all = TRUE, sort = FALSE)
                } else {
                  output <- merge(output, data.frame(id = info.group.fit$id[i], study = info.group.fit$authors[i], group = paste(groups[j][[1]], collapse="; "), fitness = info.group.fit$fitness[i], trait = traits.base[i], optimum = median, sd = sd, optimumtest = within, samplesize = samplesize, linear = sum$solutions[paste("`", traits.base[i], "`", sep=""),1], linear.se = sd(mod$Sol[,paste("`", traits.base[i], "`", sep="")]), linear.95CI.l = sum$solutions[paste("`", traits.base[i], "`", sep=""),2], linear.95CI.u = sum$solutions[paste("`", traits.base[i], "`", sep=""),3], linear.effsamp = sum$solutions[paste("`", traits.base[i], "`", sep=""),4], linear.pMCMC = sum$solutions[paste("`", traits.base[i], "`", sep=""),5], quadratic = sum$solutions[paste("I(`", traits.base[i], "`^2)", sep=""),1], quadratic.se = sd(mod$Sol[,paste("I(`", traits.base[i], "`^2)", sep="")]), quadratic.95CI.l = sum$solutions[paste("I(`", traits.base[i], "`^2)", sep=""),2], quadratic.95CI.u = sum$solutions[paste("I(`", traits.base[i], "`^2)", sep=""),3], quadratic.effsamp = sum$solutions[paste("I(`", traits.base[i], "`^2)", sep=""),4], quadratic.pMCMC = sum$solutions[paste("I(`", traits.base[i], "`^2)", sep=""),5]), all = TRUE, sort = FALSE)
                }
              }
              cat(paste("\n\nTrait: ", info.group.fit$trait[i], "\nOptimum median: ", median, "\nOptimum SD: ", sd, "\nWithin 3 SD: ", within, sep=""))
              plot(density(optimum), xlim = c(-(sd-median),(sd+median)), main = NA)
              title(main = paste(str_to_sentence(info.group.fit$trait[i]), " optimum distribution", sep=""), line = 1, font.main = 1)
            } else {
              output <- merge(output, data.frame(id = info.group.fit$id[i], study = info.group.fit$authors[i], group = paste(groups[j][[1]], collapse="; "), fitness = info.group.fit$fitness[i], trait = traits.base[i], optimum = NA, sd = NA, optimumtest = NA, samplesize = NA, linear = NA, linear.se = NA, linear.95CI.l = NA, linear.95CI.u = NA, linear.effsamp = NA, linear.pMCMC = NA, quadratic = NA, quadratic.se = NA, quadratic.95CI.l = NA, quadratic.95CI.u = NA, quadratic.effsamp = NA, quadratic.pMCMC = NA), all = TRUE, sort = FALSE)
              cat(paste("\n\n", traits.base[i], " has been removed from the model.", sep=""))
            }
          }
        } else {
          for (i in 1:length(traits.base)) {
            output <- merge(output, data.frame(id = info.group.fit$id[i], study = info.group.fit$authors[i], group = paste(groups[j][[1]], collapse="; "), fitness = info.group.fit$fitness[i], trait = traits.base[i], optimum = NA, sd = NA, optimumtest = NA, samplesize = NA, linear = NA, linear.se = NA, linear.95CI.l = NA, linear.95CI.u = NA, linear.effsamp = NA, linear.pMCMC = NA, quadratic = NA, quadratic.se = NA, quadratic.95CI.l = NA, quadratic.95CI.u = NA, quadratic.effsamp = NA, quadratic.pMCMC = NA), all = TRUE, sort = FALSE)
          }
          cat("\n\n")
          cat(paste(groups[j][[1]], collapse=" & "))
          cat(paste(" has no complete data points.", sep=""))
        }
        if (groups[j][[1]][1] == "" | is.na(groups[j][[1]][1])) {
          title <- data.id$group[data.id$id == id]
        } else {
          title <- paste(data.id$group[data.id$id == id], as.character(groups[j][[1]][1]), sep=" - ")
          for (m in 1:length(groups[j][[1]])) {
            if (m != 1) {
              title <- paste(title, as.character(groups[j][[1]][m]), sep=" & ")
            }
          }
        }
        if (!is.na(q) && length(sapply(list(info$uni), unique)) > 1) {
          title <- paste(title, " (", q, ")", sep="")
        }
        title <- paste(title, info.group.fit$fitness[i], sep=" - ")
        mtext(title, side = 3, line = -0.2, font = 2, outer = TRUE)
        saveRDS(mod, paste("results/", id, " ", title, " - ", "mod.rds", sep=""))
      }
    }
  }
  write.csv(output, paste("results/", id, " ", data.id$group[data.id$id == id], ".csv", sep=""))
  return(output)
}

# loading data sheet. id sheet includes relative locations and filenames of study data
data <- read_excel("data.xlsx", sheet = "data")
data.id <- read_excel("data.xlsx", sheet = "id")
data <- data[!(is.na(data$id)),]

# running analysis for all studies
results <- data.frame(id = 1, study = "a", group = "a", fitness = "a", trait = "a", optimum = 1, sd = 1, optimumtest = 1, samplesize = 1, linear = 1, linear.se = 1, linear.95CI.l = 1, linear.95CI.u = 1, linear.effsamp = 1, linear.pMCMC = 1, quadratic = 1, quadratic.se = 1, quadratic.95CI.l = 1, quadratic.95CI.u = 1, quadratic.effsamp = 1, quadratic.pMCMC = 1)
results <- results[-1,]
for (id in data.id$id) {
  results <- merge(results, runmod(id), all = TRUE, sort = FALSE)
}
