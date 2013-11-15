plot.glarma <- function(x, which = c(1L:6L), fits = 1L:3L,
                        ask = prod(par("mfcol")) <
                              length(which) && dev.interactive(),
                        lwdObs = 1, lwdFixed = 1, lwdGLARMA = 1,
                        colObs = "black", colFixed = "blue", colGLARMA = "red",
                        ltyObs = 2, ltyFixed = 1, ltyGLARMA = 1,
                        pchObs = 1, legend = TRUE, residPlotType = "h",
                        bins = 10, line = TRUE, colLine = "red",
                        colHist = "royal blue", lwdLine = 2,
                        colPIT1 = "red", colPIT2 = "black",
                        ltyPIT1 = 1, ltyPIT2 = 2, typePIT = "l", ...)
{
    ## ask = prod(par("mfcol")) < length(which) && dev.interactive()
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    showFits <- rep(FALSE, 3)
    showFits[fits] <- TRUE

    if (x$type == "Poi" | x$type == "NegBin") {
        if (show[1L] & any(showFits == TRUE)) {
            obs.logic <- showFits[1]
            fixed.logic <- showFits[2]
            glarma.logic <- showFits[3]
            fits <- list(obs = x$y, fixed = exp(x$eta), glarma = x$mu)
            legendNames <- names(fits)
            titleNames <- c("Observed", "Fixed", "GLARMA")
            ltyAll <- c(ltyObs, ltyFixed, ltyGLARMA)
            lwdAll <- c(lwdObs, lwdFixed, lwdGLARMA)
            colAll <- c(colObs, colFixed, colGLARMA)
            yRange <- c(diff(range(fits$obs)), diff(range(fits$fixed)),
                        diff(range(fits$glarma)))
            yRange[!showFits] <- -Inf
            yLim <- which.max(yRange)
        }
        if (any(show[2L:4L] == TRUE)) {
            residuals <- x$residuals
            caption <- c(paste("ACF of", x$residType, "Residuals"),
                         paste(x$residType, "Residuals"),
                         paste("Normal Q-Q Plot of ",x$residType, "Residuals"))
        }
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
        if (show[1L] & any(showFits == TRUE)) {
            dev.hold()

            ts.plot(fits[[yLim]], ylab = "Counts", xlab = "Time",
                    col = NA, main = paste(titleNames[showFits],
                    collapse = " vs "), ...)
            if (obs.logic)
                lines(fits$obs, lwd = lwdObs, lty = ltyObs, col = colObs)

            if (fixed.logic)
                lines(fits$fixed, lwd = lwdFixed, lty = ltyFixed,
                      col = colFixed)

            if (glarma.logic)
                lines(fits$glarma, lwd = lwdGLARMA, lty = ltyGLARMA,
                      col = colGLARMA)

            if (legend & any(showFits == TRUE)) {
                par(xpd = NA)
                mfrow <- par("mfrow")
                grap.param <- legend("top", legend = legendNames[showFits],
                                     lty = ltyAll[showFits],
                                     ncol = 3,
                                     cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                                     bty = "n", plot = FALSE)
                legend(grap.param$rect$left,
                       grap.param$rect$top + grap.param$rect$h,
                       legend = legendNames[showFits], col = colAll[showFits],
                       lwd = lwdAll[showFits], lty = ltyAll[showFits],
                       ncol = 3,
                       cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                       bty = "n",
                       text.font = 4)
                par(xpd = FALSE)
            }
            dev.flush()
        }
        if (show[2L]) {
            dev.hold()
            acf(residuals, main = caption[1], ...)
            dev.flush()
        }
        if (show[3L]) {
            dev.hold()
            plot.ts(residuals, ylab = "Residuals", type = residPlotType,
                    main = caption[2], ...)
            dev.flush()
        }
        if (show[4L]) {
            dev.hold()
            qqnorm(residuals, ylab = "Residuals", main = caption[3], ...)
            abline(0, 1, lty = 2)
            dev.flush()
        }
        if (show[5L]) {
          dev.hold()
          histPIT(x, bins = bins, line = line, colLine = colLine,
                  colHist = colHist, lwdLine = lwdLine, ...)
          dev.flush()
        }
        if (show[6L]) {
          dev.hold()
          qqPIT(x, bins = bins, col1 = colPIT1, col2 = colPIT2,
                lty1 = ltyPIT1, lty2 = ltyPIT2, type = typePIT, ...)
          dev.flush()
        }
    }

    if (x$type == "Bin") {
        if (show[1L] & any(showFits = TRUE)) {
            obs.logic <- showFits[1]
            fixed.logic <- showFits[2]
            glarma.logic <- showFits[3]
            observed <- x$y[, 1]/apply(x$y, 1, sum)
            fits <- list(fixed = 1/(1 + exp(-x$eta)), glarma = 1/(1 + exp(-x$W)))
            legendNames <- c("obs", names(fits)[1], names(fits)[2])
            titleNames <- c("Observed", "Fixed", "GLARMA")
            pchAll <- c(pchObs, NA, NA)
            ltyAll <- c(ltyObs, ltyFixed, ltyGLARMA)
            lwdAll <- c(NA, lwdFixed, lwdGLARMA)
            colAll <- c(colObs, colFixed, colGLARMA)
        }
        if (any(show[2L:4L] == TRUE)) {
            residuals <- x$residuals
            caption <- c(paste("ACF of", x$residType, "Residuals"),
                         paste(x$residType, "Residuals"),
                         paste("Normal Q-Q Plot of ",x$residType, "Residuals"))
        }
        if (ask) {
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))
        }
        if (show[1L] & any(showFits == TRUE)) {
            dev.hold()
            plot(1:length(observed), observed, ylab = "Counts", xlab = "Time",
                 col = NA, main = paste(titleNames[showFits], collapse = " vs "), ...)
            if (obs.logic)
                points(observed, pch = pchObs, col = colObs)

            if (fixed.logic)
                lines(fits$fixed, lwd = lwdFixed, lty = ltyFixed,
                      col = colFixed)

            if (glarma.logic)
                lines(fits$glarma, lwd = lwdGLARMA, lty = ltyGLARMA,
                      col = colGLARMA)

            if (legend & any(showFits == TRUE)) {
                par(xpd = NA)
                mfrow <- par("mfrow")
                grap.param <- legend("top", legend = legendNames[showFits],
                                     pch = pchAll[showFits],
                                     lty = ltyAll[showFits], ncol = 3,
                                     cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                                     text.font = 4, plot = FALSE)

                legend(grap.param$rect$left, grap.param$rect$top + grap.param$rect$h,
                       legend = legendNames[showFits], pch = pchAll[showFits],
                       col = colAll[showFits], lwd = lwdAll[showFits],
                       lty = ltyAll[showFits], ncol = 3,
                       cex = 0.7 - (mfrow[1] - 1)/10 - (mfrow[2] - 1)/10,
                       bty = "n", text.font = 4)
                par(xpd = FALSE)
            }
            dev.flush()
        }
        if (show[2L]) {
            dev.hold()
            acf(residuals, main = caption[1], ...)
            dev.flush()
        }
        if (show[3L]) {
            dev.hold()
            plot.ts(residuals, ylab = "Residuals", type = residPlotType,
                    main = caption[2], ...)
            dev.flush()
        }
        if (show[4L]) {
            dev.hold()
            qqnorm(residuals, ylab = "Residuals", main = caption[3], ...)
            abline(0, 1, lty = 2)
            dev.flush()
        }
        if (show[5L]) {
           dev.hold()
           histPIT(x, bins = bins, line = line, colLine = colLine,
                   colHist = colHist, lwdLine = lwdLine, ...)
           dev.flush()
        }
        if (show[6L]) {
           dev.hold()
           qqPIT(x, bins = bins, col1 = colPIT1, col2 = colPIT2,
                 lty1 = ltyPIT1, lty2 = ltyPIT2, type = typePIT, ...)
           dev.flush()
        }
    }
    invisible()
}
