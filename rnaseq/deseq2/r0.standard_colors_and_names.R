# ======================================================================
# VISUAL - establish standard colors and labels for all figures
# ======================================================================

# CUSTOMIZE - these are Zajac values

# display NAMES
displayNames = list()
displayNames$cell <- subset <- list(naive="naive", "unstim"="unstim", "IL2-"="minus", "IL2+"="plus")
displayNames$cell <- type <- list(CD8="cd8", "effector"="eff", "memory"="mem")
displayNames$cond <- abbrev <- list(CD8.naive="cd8.naive",
                                    "effector unstim"="eff.ustim", "effector IL2-"="eff.minus", "effector IL2+"="eff.plus",
                                    "memory unstim"="mem.ustim", "memory IL2-"="mem.minus", "memory IL2+"="mem.plus"
                                    )

# display COLORS
displayColors = list()

# GFP+=>blue, GFP-=>red, unstim=>gray, naive=>white
displayColors$cell <- subset <-        c("gray90", "gray40",  "blue", "blue", "red",    "red")
names(displayColors$cell <- subset) <- c("naive",  "unstim", "plus", "IL2+", "minus", "IL2-")

displayColors$cell <- type <-        c("gray90", "gray90", "deeppink","deeppink",  "purple","purple")
names(displayColors$cell <- type) <- c("cd8",    "CD8",    "eff",     "effector",  "mem",   "memory")


displayColors$anno <- type <- subset <- list(
                                  cell <- type   = displayColors$cell <- type,
                                  cell <- subset   = displayColors$cell <- subset
                                  )
# add aliases
displayColors$anno <- type <- subset$'Cell Type' <- displayColors$anno <- type <- subset$cell <- type
displayColors$anno <- type <- subset$'Subset' <- displayColors$anno <- type <- subset$cell <- subset

#
