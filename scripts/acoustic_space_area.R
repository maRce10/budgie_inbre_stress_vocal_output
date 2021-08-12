# need to load these packages first
library(MASS)  
library(adehabitatHR) 
library(sp)
library(pbapply)  

# each level must have at least 5 data points
# ... additional parameters to be passed to kernelUD() for kernel area estimation

acoustic_space_area <- function(X, dimensions = NULL, group, cl = 1, type = "mcp", pb = TRUE, distance.matrix = FALSE, outliers = 0.95, rarefaction = FALSE, iteracions = 30, proportional = TRUE, ...) {
    
    if (is.null(dimensions) & !distance.matrix)
        stop("'dimensions' must be supplied when working with tabular (i.e. non-distance matrix) data")

    if (distance.matrix & length(group) < nrow(X))
        stop("'group' must have the same length than the number of rows in 'X' when the latter is a distance matrix")
    
        
    # reset progress bar when exiting
    on.exit(pbapply::pboptions(type = .Options$pboptions$type))
    
    # set progress bar
    pbapply::pboptions(type = ifelse(pb, "timer", "none"))
    
    # convert to distance
    if (distance.matrix){

        # mds for 2 dimensions
        mds <- isoMDS(X, y = cmdscale(dst_mt, k = 2), k = 2, maxit = 5, trace = FALSE, tol = 1e-3, p = 2)

        # put in a data frame
        X <- data.frame(group, mds$points)
        names(X) <- c("group", "MDS1", "MDS2")
        
        dimensions <- c("MDS1", "MDS2")
        }    
    
    # split in a list of vectors
    X_l <- split(x = X, f = X[, group])
    
    # for rarefaction estimation
    min.n <- min(sapply(X_l, nrow))
    
    # stop if too small sample sizes
    if (min.n < 6 & rarefaction & type == "kernel")
        stop("There is at least one group with less than 6 observations which is the minimum needed for kernel area estimation")

    if (min.n < 5 & rarefaction & type == "mcp")
        stop("There is at least one group with less than 5 observations which is the minimum needed for 'mcp' area estimation")
    
   # function to calculate areas
    area_fun <- function(W, tp, ol) {
        
        # get area
        # kernel area
        if (tp == "kernel"){
            area <- if (length(W) >= 6) kernel.area(kernelUD(W, extent = 1.5, ...), percent = ol * 100)[[1]] else NA
        } 
        
        # mcp area
        if (tp == "mcp"){
            area <- if (length(W) >= 5) mcp(xy = W, percent = ol * 100)$area else area <- NA 
        } 
        
        return(area)
        }
            
   # function to calculated rarefacted areas
    raref_fun <- function(W, tp, ol, min.n){
        
        Z <- W[sample(1:length(W), min.n),]
        
        area <- area_fun(Z, tp, ol)
    
        return(area)    
    }
    
    # calculate all areas
    areas_l <- pblapply(X_l, cl = cl, function(Y, dims = dimensions, typ = type, mn.n = min.n, outlrs = outliers, reps = iteracions) {
        
        # subset for each individual 
        W <- Y[, dims]
        
        # add dimensions as coordinates to calculate acoustic space area
        coordinates(W) <-  as.formula(paste("~ ", dims[1], "+", dims[2]))
      
        area <- if (rarefaction) mean(replicate(n = reps, expr = raref_fun(W, typ, outlrs, mn.n)), na.rm =  TRUE) else area_fun(W, typ, outlrs)
        
        # put in a data frame
        out_df <- data.frame(group = Y[1, group], n = nrow(Y), area = area)
        
        # add stress if distance matrix was supplied
        if (distance.matrix) 
            out_df$stress <- mds$stress
    
        return( out_df)
    }) 
    
    # put in a data frame
    acous.areas <- do.call(rbind, areas_l)
    
    # rename rows
    row.names(acous.areas) <- 1:nrow(acous.areas)
    
    # make it proportional to largest area
    if (proportional)
        acous.areas$area <- acous.areas$area / max(acous.areas$area, na.rm = TRUE)
        
    # rename columns
    if (rarefaction)
        acous.areas$rarefation.n <- min.n
    
    return(acous.areas)
}


acoustic_space_overlap <- function(X, dimensions = NULL, group, cl = 1, type = "density", pb = TRUE, distance.matrix = FALSE, outliers = 0.95, mean.overlap = TRUE, ...) {
    
    if (is.null(dimensions) & !distance.matrix)
        stop("'dimensions' must be supplied when working with tabular (i.e. non-distance matrix) data")
    
    if (distance.matrix & length(group) < nrow(X))
        stop("'group' must have the same length than the number of rows in 'X' when the latter is a distance matrix")
    
    
    if (type == "distance") mean.overlap <- TRUE
    
    # reset progress bar when exiting
    on.exit(pbapply::pboptions(type = .Options$pboptions$type))
    
    # set progress bar
    pbapply::pboptions(type = ifelse(pb, "timer", "none"))
    
    # convert to distance
    if (distance.matrix){
        
        # mds for 2 dimensions
        mds <- isoMDS(X, y = cmdscale(dst_mt, k = 2), k = 2, maxit = 5, trace = FALSE, tol = 1e-3, p = 2)
        
        # put in a data frame
        X <- data.frame(group, mds$points)
        names(X) <- c("group", "MDS1", "MDS2")
        
        dimensions <- c("MDS1", "MDS2")
    }    
    
    # force dimensions to range 0-1
    X[, dimensions[1]] <- (X[, dimensions[1]] - min( X[, dimensions[1]])) / max((X[, dimensions[1]] - min( X[, dimensions[1]])))
    X[, dimensions[2]] <- (X[, dimensions[2]] - min( X[, dimensions[2]])) / max((X[, dimensions[2]] - min( X[, dimensions[2]])))
    
    # split in a list of vectors
    X_l <- split(x = X, f = X[, group])
    
    # stop if too small sample sizes
    if (min(sapply(X_l, nrow)) < 2)
        stop("There is at least one group with less than 2 observations which is the minimum needed for overlap estimation")
    
    # get densities
    if (type == "density"){
            
      total_coors <- as.ppp(as.matrix(X[, dimensions]), c(range(X[, dimensions[1]]), range(X[, dimensions[2]])))
      total_space <-  raster(density.ppp(total_coors))
       
    densities <- lapply(X_l, function(Y, dims = dimensions){
        
        coors <- as.ppp(as.matrix(Y[, dimensions]), c(range(Y[, dimensions[1]]), range(Y[, dimensions[2]])))
        
        # raster_dens <- raster(as.matrix.data.frame(Y[, dimensions]))
        raster_dens <- raster(density.ppp(coors))
        
        raster_dens <- crop(extend(raster_dens, total_space), total_space)
        
        # mask based on total acoustic space
        # raster_dens <- mask(raster_dens, total_space)
        
        # convert to range 0-1
        values(raster_dens) <- (values(raster_dens) - min(values(raster_dens), na.rm = TRUE)) 
        values(raster_dens) <- values(raster_dens)/ max(values(raster_dens), na.rm = TRUE)
        
        # keep 95% interval
        values(raster_dens)[values(raster_dens) < 1 - outliers] <- NA 
     
        # same resolution and extend as total space
        raster_dens <- resample(raster_dens, total_space)
        extent(raster_dens)  <- c(0,1, 0, 1)
        
        return(raster_dens)
        })
    
    names(densities) <- names(X_l)
        
    }
    
    # function to calculate areas
    ovlp_fun <- function(W, Z, tp) {
        
        # get area
        # kernel area
        if (tp == "density"){
            
          # convert NAs to 0
          values(W)[is.na(values(W))] <- 0
          values(Z)[is.na(values(Z))] <- 0

            # 1 vs 2
            wgt_ovlp.WvZ <- sum((values(W) * values(Z)) / values(Z), na.rm = TRUE) / sum(values(Z), na.rm = TRUE)
            
            # 2 vs 1
            wgt_ovlp.ZvW <- sum((values(Z) * values(W)) / values(W), na.rm = TRUE) / sum(values(W), na.rm = TRUE) 
            
            # convert to 1 if higher than 1
            if (wgt_ovlp.WvZ > 1) wgt_ovlp.WvZ <- 1
            if (wgt_ovlp.ZvW > 1) wgt_ovlp.ZvW <- 1
            
            
            
      # put results in am matrix     
      output <- matrix(c(wgt_ovlp.WvZ, wgt_ovlp.ZvW), nrow = 1)     
        } 
        
        # distance among points
        if (tp == "distance"){
            
            U <- rbind(W, Z)
            
            dists <- as.matrix(dist(U[ , dimensions]))
            
            dist.1v2 <- dists[U[, group] == W[, group][1], U[, group] == Z[, group][1]]
            
            output <- matrix(mean(dist.1v2), nrow = 1)
        } 
        
        return(output)
    }
    
    # function to calculated rarefacted areas
    raref_fun <- function(W, Z, tp, min.n){
        
        Q <- W[sample(1:length(W), min.n),]
        P <- Z[sample(1:length(Z), min.n),]
        ovlp <- ovlp_fun(Z, tp)
        
        return(ovlp)    
    }
    
    # get all combinations to get pairwise overlaps
    group_combs <- t(combn(unique(X[, group]), 2))
    
    # calculate all areas
    ovlps_l <- pblapply(1:nrow(group_combs), cl = cl, function(i, dens = densities, gc = group_combs, dims = dimensions, typ = type, mn.n = min.n, reps = iteracions) {
        
        if (type == "density"){
            W <- densities[[which(names(densities) == gc[i, 1])]]
            Z <- densities[[which(names(densities) == gc[i, 2])]]
        } else {
            W <- X_l[[which(names(X_l) == gc[i, 1])]]
            Z <- X_l[[which(names(X_l) == gc[i, 2])]]
        }
        
        overlaps <- ovlp_fun(W, Z, typ)
        
        # put in a data frame
        
        out_df <- if (type == "density" & !mean.overlap) 
            data.frame(group.1 = gc[i, 1], group.2 = gc[i, 1], overlap.1v2 = overlaps[1], overlap.2v1 = overlaps[2]) else
         data.frame(group.1 = gc[i, 1], group.2 = gc[i, 1], overlap = mean(overlaps)) 
        return(out_df)
    }) 
    
    # put in a data frame
    acous.ovlp <- do.call(rbind, ovlps_l)
    
    #create a similarity matrix with the max xcorr
    mat <- matrix(nrow = length(X_l), ncol = length(X_l))
    mat[] <- if (type == "distance") 0 else 1
    colnames(mat) <- rownames(mat) <- names(X_l)
    
    # put data into a matrix
    if (mean.overlap){
        # add max correlations
        mat[lower.tri(mat, diag=FALSE)] <- acous.ovlp$overlap
        mat <- t(mat)
        mat[lower.tri(mat, diag=FALSE)] <- acous.ovlp$overlap
    
        out_mat <- mat
        } else {
        mat2 <- mat
        
        # mat 1 vs 2
        mat[lower.tri(mat, diag=FALSE)] <- acous.ovlp$overlap.1v2
        mat <- t(mat)
        # mat[lower.tri(mat, diag=FALSE)] <- acous.ovlp$overlap.1v2
        mat[lower.tri(mat, diag=FALSE)] <- acous.ovlp$overlap.2v1
        
        # mat 2 vs 1
        # mat2[lower.tri(mat2, diag=FALSE)] <- acous.ovlp$overlap.2v1
        # mat2 <- t(mat2)
        # mat2[lower.tri(mat2, diag=FALSE)] <- acous.ovlp$overlap.2v1
        # 
        # out_mat <- list(overlap.1v2 = mat, overlap.2v1 = mat2)
        out_mat <- mat
    }
    return(out_mat)
}

#to get the binary matrices
binary_matrix <- function(labels) {
    # mat= original matrix 
    # labels = labels of rows/columns belonging to the same level
    #create new mat
    x <- matrix(nrow = length(labels), ncol = length(labels))  
    
    # with only NA
    x[!is.na(x)] <- NA
    
    # fill it out with 0s and 1s 
    for(i in 1:ncol(x))
        for(j in 1:length(x[,i]))
            if (labels[j] == labels[i]) x[j, i] <- 0 else x[j, i] <- 1
    
    # rename cols and rows (rows cannot be duplicated)
    rownames(x) <- paste0(labels, " (row ", 1:nrow(x), ")")
    colnames(x) <- labels
    return(x)      
}

# acu_dens_ovlp <- function(Y, Z. acu.sp.sp, mds_prox, usr)
# {   
#     Ypre <- Y[Y$period == "pre", ]
#     Ypost <- Y[Y$period == "post", ]
#     
#     if(nrow(Ypre) > 0) 
#     {
#         pp.coords.pre <- as.ppp(as.matrix(Ypre[, c("MDS1", "MDS2")]), usr)
#         r.coords.pre <- raster(density.ppp(pp.coords.pre))
#         
#         msk.pre <- mask(r.coords.pre, acu.sp.sp)
#         
#         # convert to range 0-1
#         values(msk.pre) <- (values(msk.pre) - min(values(msk.pre), na.rm = TRUE)) 
#         values(msk.pre) <- values(msk.pre)/ max(values(msk.pre), na.rm = TRUE)
#         
#         # keep 95% interval
#         values(msk.pre)[values(msk.pre) < 0.05] <- NA
#     }
#     
#     if(nrow(Ypost) > 0)  
#     {  
#         pp.coords.post <- as.ppp(as.matrix(Ypost[, c("MDS1", "MDS2")]), usr)
#         r.coords.post <- raster(density.ppp(pp.coords.post))
#         
#         msk.post <- mask(r.coords.post, acu.sp.sp)
#         
#         # convert to range 0-1
#         values(msk.post) <- (values(msk.post) - min(values(msk.post), na.rm = TRUE)) 
#         values(msk.post) <- values(msk.post) / max(values(msk.post), na.rm = TRUE)
#         
#         # keep 95% interval
#         values(msk.post)[values(msk.post) < 0.05] <- NA
#     }
#     
#     # test method
#     # msk.pre <- msk.post
#     # values(msk.pre)[sample(1:length(values(msk.pre)), size = (length(values(msk.pre)))*0.75)] <- 0
#     
#     
#     try(values(msk.pre)[!is.na(values(msk.post)) & is.na(values(msk.pre))] <- 0, silent = TRUE)
#     
#     
#     wgt.ovlp <- try_na(sum(values(msk.post) * (values(msk.post) - abs(values(msk.post) - values(msk.pre))), na.rm = TRUE) / sum(values(msk.post) * values(msk.post), na.rm = TRUE))
#     return(wgt.ovlp)
# }

# plot_acu_space <- function(Y, x, mds_prox)
# {
#     Ypre <- Y[Y$period == "pre", ]
#     Ypost <- Y[Y$period == "post", ]
#     Gpost <- mds_prox[mds_prox$period == "post" & mds_prox$Cage == Y$Cage[1] & mds_prox$Bird.Name != bnm[x], ]
#     # W <- prepost[prepost$POSTpair == Y$POSTpair[1] & prepost$treatment == "Pre", ]
#     Y <- rbind(Y, Gpost)
#     
#     #select calls closest to mean acoustic space
#     if(nrow(Ypre) > 0)
#     {
#         mt <- rbind(data.frame(MDS1 = mean(Ypre[,"MDS1"]),MDS2 = mean(Ypre[,"MDS2"])), Ypre[,c("MDS1", "MDS2")])
#         
#         fspre <- Ypre$sound.files[which.min(as.matrix(dist(mt, upper = T))[-1,1])]
#         Zpre <- est[est$sound.files == fspre,, drop = FALSE][1,]
#         predur <- Zpre$end - Zpre$start
#     } else predur <- 0
#     
#     if(nrow(Ypost) > 0)
#     {    
#         mt2 <- rbind(data.frame(MDS1 = mean(Ypost[,"MDS1"]),MDS2 = mean(Ypost[,"MDS2"])), Ypost[,c("MDS1", "MDS2")])
#         fspost <- Ypost$sound.files[which.min(as.matrix(dist(mt2, upper = T))[-1, 1])]
#         Zpost <- est[est$sound.files == fspost, , drop = FALSE][1,]
#         postdur <- Zpost$end - Zpost$start
#     } else postdur <- 0
#     
#     if(nrow(Gpost) > 0)
#     {
#         mt3 <- rbind(data.frame(MDS1 = mean(Gpost[,"MDS1"]),MDS2 = mean(Gpost[,"MDS2"])), Gpost[,c("MDS1", "MDS2")])
#         fspost2 <- Gpost$sound.files[which.min(as.matrix(dist(mt3, upper = T))[-1, 1])]
#         Zpost2 <- est[est$sound.files == fspost2, , drop = FALSE][1, , drop = FALSE]
#         postdur2 <- Zpost2$end - Zpost2$start
#     }  else postdur2 <- 0
#     
#     #and closest to group centroid
#     # mean(Y$MDS1[Y$indiv != x])
#     # G <- prepost[prepost$indiv != x,]
#     #  fsgroup <- as.character(G$file.sel[which.min(dist(G[,2:3], y = data.frame(mean(G[,2]), mean(G[,3]))))])
#     # Gcent <- sel.data[sel.data$sf.sel == fsgroup,]
#     # Gdur <- Gcent$end - Gcent$start
#     
#     #adjust all to the same duration
#     maxdur <- max(postdur, predur, postdur2)
#     if(maxdur > predur & nrow(Ypre) > 0) 
#     {
#         Zpre$start <- (Zpre$start - (maxdur-predur)/2)
#         Zpre$end <- (Zpre$end + (maxdur-predur)/2)
#     } 
#     
#     if(maxdur > postdur & nrow(Ypost) > 0) 
#     {
#         Zpost$start <- (Zpost$start - (maxdur-postdur)/2)
#         Zpost$end <- (Zpost$end + (maxdur-postdur)/2)
#     } 
#     
#     if(maxdur > postdur2 & nrow(Zpost2) > 0) 
#     {
#         Zpost2$start <- (Zpost2$start - (maxdur-postdur)/2)
#         Zpost2$end <- (Zpost2$end + (maxdur-postdur)/2)
#     } 
#     
#     if(nrow(Ypre) > 0)   
#         wpre <- read_wave(est, index = which(est$sound.files == Zpre$sound.files), from = Zpre$start, to = Zpre$end)
#     
#     if(nrow(Ypost) > 0)
#         wpost <- read_wave(est, index = which(est$sound.files == Zpost$sound.files), from = Zpost$start, to = Zpost$end)
#     
#     if(nrow(Zpost2) > 0)
#         wpost2 <- read_wave(est, index = which(est$sound.files == Zpost2$sound.files), from = Zpost2$start, to = Zpost2$end)
#     
#     # tiff(filename = file.path("/media/twright/E88673588673266A/Budgie analysis/Paper 3 social context FoxP/Acoustics space change with spectros and group centroid", paste0(x, ".tiff")), units = "px", width = 630, height = 480)
#     
#     nf <- layout(matrix(c(9, 9, 9, 9, 1, 2, 5, 8, 1, 3, 6, 8, 1, 4, 7, 8), nr = 4, byrow=TRUE), widths=c(0.2, 1.5, 0.11, 3), heights=c(0.2, 1, 1, 1.4))
#     par(mar = c(0, 0, 0, 0))
#     plot(1, col = "white", frame.plot = F)
#     text(x = 1, y = 1.05, "Frequency (kHz)", srt = 90, cex = 2) 
#     par(mar = c(0, 2, 0, 0))
#     
#     if(nrow(Ypre) > 0)  
#         spectro(wpre, f = wpre@samp.rate, wl = wl, ovlp = ovlp, scale = F, osc = F, flim = c(1.1, 5.9), palette = reverse.gray.colors.1, grid = F, axisX = F, tlab = "", cexaxis = 2) else
#             plot(1, col = "white", frame.plot = F, axes = F)
#     
#     
#     if(nrow(Ypost) > 0)  
#         spectro(wpost, f = wpost@samp.rate, wl = wl, ovlp = ovlp, scale = F, osc = F, flim = c(1.1, 5.9), palette = reverse.gray.colors.1, grid = F,axisX = F, tlab = "", cexaxis = 2) else
#             plot(1, col = "white", frame.plot = F, axes = F)
#     
#     
#     if(nrow(Zpost2) > 0)  
#     {
#         par(mar = c(4, 2, 0, 0))  
#         spectro(wpost2, f = wpost2@samp.rate, wl = wl, ovlp = ovlp, scale = F, osc = F, flim = c(1.1, 5.9), palette = reverse.gray.colors.1, grid = F, cexaxis = 2, cexlab = 2)
#     } else  plot(1, col = "white", frame.plot = F, axes = F)
#     
#     alpha <- 0.5
#     colpre <- adjustcolor( "red2", alpha.f = alpha)  
#     colpost <- adjustcolor( "cyan4", alpha.f = alpha)  
#     colgrp <- adjustcolor( "gray", alpha.f = alpha)  
#     
#     # Y$col <- as.character(Y$period)
#     #   Y$col[Y$period ==  "pre"] <- cols
#     #   Y$col[Y$period ==  "post"] <- col1 
#     #     Y$col[Y$period ==  "post"] <- col3 
#     
#     par(mar = c(0, 0, 0, 0))
#     plot(1, col = "white", frame.plot = F, axes = F)
#     rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = colpre)
#     
#     par(mar = c(0, 0, 0, 0))
#     plot(1, col = "white", frame.plot = F, axes = F)
#     rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = colpost)
#     
#     par(mar = c(4, 0, 0, 0))
#     plot(1, col = "white", frame.plot = F, axes = F, xlab = "")
#     rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = colgrp)
#     
#     par(mar = c(4, 6, 0, 0))  
#     basecex <- 2.3
#     
#     #fix plot if points are overlapped by legend
#     # a <- Y[Y$MDS1 > max(Y$MDS1) - (diff(range(Y$MDS1)) * 0.2) & Y$MDS2 > max(Y$MDS2) - (diff(range(Y$MDS2)) * 0.2),]
#     # xlim <- c(min(Y$MDS1) - (diff(range(Y$MDS1)) * 0.035), max(Y$MDS1) + (diff(range(Y$MDS1)) * 0.35))
#     # ylim <- c(min(Y$MDS2) - (diff(range(Y$MDS2)) * 0.035), max(Y$MDS2) + (diff(range(Y$MDS2)) * 0.35))
#     
#     xlim <- range(acou.space[, 1])
#     xlim[2] <-  xlim[2] * 1.15
#     
#     ylim <- range(acou.space[, 2])
#     ylim[2] <-  ylim[2] * 1.45
#     
#     # if(nrow(a) > 0)  
#     plot(x = Y$MDS1, y = Y$MDS2, col = "white", pch = 20, cex = basecex,  xlab = "MDS 1", ylab= "MDS 2", xlim = xlim, ylim = ylim, cex.lab = 2) #else
#     # plot(x = Y$MDS1, y = Y$MDS2, col = "white", pch = 20, cex = basecex,  xlab = "MDS 1", ylab= "MDS 2", cex.lab = 2)
#     
#     #plot groups Post calls
#     # points(x = Y$MDS1[Y$indiv != x], y = Y$MDS2[Y$indiv != x], col = col3, pch = 20, cex = basecex * 1.4)
#     
#     polygon(x = acou.space[,1], y = acou.space[,2], col = adjustcolor(cols[4], 0.2), border = adjustcolor(cols[4], 0.2), lwd =4)
#     # plot(acou.space, add = T)
#     #plot pre and post indiv calls
#     if(nrow(Zpost2) > 0)  
#         points(x = Gpost$MDS1, y = Gpost$MDS2, col = colgrp, pch = 20, cex = basecex* 1.4)
#     
#     if(nrow(Ypre) > 0)  
#         points(x = Ypre$MDS1, y = Ypre$MDS2, col = colpre, pch = 20, cex = basecex* 1.4)
#     
#     if(nrow(Ypost) > 0)  
#         points(x = Ypost$MDS1, y = Ypost$MDS2, col = colpost, pch = 20, cex = basecex* 1.4)
#     
#     if(nrow(Gpost) > 0) 
#     {
#         points(x = Y$MDS1[Y$sound.files == fspost2], y = Y$MDS2[Y$sound.files == fspost2], col = "white", pch = 20, cex = basecex)
#         points(x = Y$MDS1[Y$sound.files == fspost2], y = Y$MDS2[Y$sound.files == fspost2], col = colgrp, pch = 20, cex = basecex)}    
#     
#     
#     if(nrow(Ypre) > 0) 
#     {points(x = Y$MDS1[Y$sound.files == fspre], y = Y$MDS2[Y$sound.files == fspre], col = "black", pch = 20, cex = basecex * 1.4)
#         points(x = Y$MDS1[Y$sound.files == fspre], y = Y$MDS2[Y$sound.files == fspre], col = "white", pch = 20, cex = basecex)  
#         points(x = Y$MDS1[Y$sound.files == fspre], y = Y$MDS2[Y$sound.files == fspre], col = colpre, pch = 20, cex = basecex)
#     }
#     
#     if(nrow(Ypost) > 0) 
#     {   points(x = Y$MDS1[Y$sound.files == fspost], y = Y$MDS2[Y$sound.files == fspost], col = "black", pch = 20, cex = basecex * 1.4)
#         points(x = Y$MDS1[Y$sound.files == fspost], y = Y$MDS2[Y$sound.files == fspost], col = "white", pch = 20, cex = basecex)
#         points(x = Y$MDS1[Y$sound.files == fspost], y = Y$MDS2[Y$sound.files == fspost], col = colpost, pch = 20, cex = basecex)
#     }
#     
#     if(nrow(Gpost) > 0) 
#     {
#         points(x = Y$MDS1[Y$sound.files == fspost2], y = Y$MDS2[Y$sound.files == fspost2], col = "black", pch = 20, cex = basecex)
#     }    
#     
#     #   #plot groups Post calls centroid
#     #   mex <- mean(Y$MDS1[Y$indiv != x])
#     #   mey <- mean(Y$MDS2[Y$indiv != x])
#     #  points(x = mex, y = mey, col = "black", pch = 20, cex = basecex * 1.4)
#     # points(x = mex, y = mey, col = "white", pch = 20, cex = basecex)
#     #   points(x = mex, y = mey, col = col3, pch = 20, cex = basecex)
#     #   
#     lims <- par("usr")  
#     
#     legend(lims[2]-((lims[2] - lims[1]) * 0.55), lims[4]-((lims[4] - lims[3]) * 0.002), lty=c(1,1), col=c(colpre, colpost, colgrp, "black"), legend = c("week 1", "week 3", "Cage mate", "Cage mate centroid"), pch = 20, lwd = 0, pt.cex= basecex * 1.5, cex = basecex * 0.5)
#     
#     par(mar = c(0, 0, 0, 0))
#     plot(1, col = "white", frame.plot = F, axes = F)
#     text(1, labels = paste0(bnm[x], " (", Y$Treatment[1], ")"), cex = 2)
#     # dev.off()
# }


# example of acoustic space overlap
# 
# n <- 5000
# W <- data.frame(PC1 = runif(n), PC2 = runif(n), ID = rep(c("a"), n))
# # W <- data.frame(PC1 = rnorm(n), PC2 = rnorm(n), ID = rep(c("a"), n))
# W$PC1 <- W$PC1 + abs(min(W$PC1))
# W$PC1 <- W$PC1 / max(W$PC1)
# W$PC2 <- W$PC2 + abs(min(W$PC2))
# W$PC2 <- W$PC2 / max(W$PC2)
# 
# W2 <- W
# W2$ID <- "b"
# W2 <- W2[W2$PC1 <= 0.20, ]
# W2$PC2 <- sample(W2$PC2)
# W2$PC1 <- W2$PC1 + 0.2
# 
# 
# W3 <- W
# W3$ID <- "c"
# W3 <- W3[W3$PC1 >= 0.40, ]
# W3$PC2 <- sample(W3$PC2)
# W3$PC1 <- W3$PC1 + 0.2
# # W <- W[W$PC1 >= 0.25, ]
# 
# 
# W4 <- rbind(
#   W,
#   W2, 
#   W2[W2$PC1 >= 0.1 & W2$PC1 <= 0.4 & W2$PC2 >= 0.2 & W2$PC2 <= 0.8, ], 
#   W[W$PC1 >= 0.1 & W$PC1 <= 0.9 & W$PC2 >= 0.2 & W$PC2 <= 0.8, ],
#   W3
#   )
# 
# plot(W4$PC1[W4$ID == "a"], W4$PC2[W4$ID == "a"], col = "red", xlim  = c(-0.1, 1.2))
# points(W4$PC1[W4$ID == "b"], W4$PC2[W4$ID == "b"], col = "blue", pch = 3, cex = 0.5)
# points(W4$PC1[W4$ID == "c"], W4$PC2[W4$ID == "c"], col = "black", pch = 3, cex = 0.5)
# 
# 
# # run with density
# acoustic_space_overlap(X = W4, dimensions = c("PC1", "PC2"), group = "ID", cl = 10, pb = TRUE, distance.matrix = FALSE, outliers = 0.95, mean.overlap = FALSE) 

# acoustic_space_overlap(X = W3, dimensions = c("PC1", "PC2"), group = "ID", cl = 10, pb = TRUE, distance.matrix = FALSE, outliers = 0.95, mean.overlap = FALSE, type = "distance") 


