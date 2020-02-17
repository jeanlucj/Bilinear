#' AMMIplot function
#'
#' This function will make linear and winner plots from a bilinear() object like Figure 2 and Figure 3 of Gauch & Zobel (1997).
#'
#' @param bilinearObject object from output of a call to biliner()
#' @param plots character vector of maximum length 2 containing the names of the plots to be produced. Possible arguments are "linear" and "winner". If both are specified in a single call, both plots will be plotted to the same device
#' @param color character vector of length 2 containing colors for the plots. The first specifies the genotype color and the second the environment color. If only one color is specified, only genotypes will be colored
#' @param PC integer. Principal component to plot. The default is 1
#' @param f numeric. Scale parameter in (0, 1) for exponent on eigenvalues for weighting the genotype scores. Environment scores are weighted 1 - f. Default is 0.5
#' @param ... Additional arguments.
#' @details
#' Many arguments can be passed on to plot() through (...)
#'
#' @examples
#'
#' data(soyMeanMat)
#' AMMIfit <- bilinear(x = soyMeanMat)
#' AMMIplot(AMMIfit)
#' AMMIplot(AMMIfit, "winner")
#' AMMIplot(AMMIfit, "winner", color = "hotpink")
#' AMMIplot(AMMIfit, c("linear", "winner"), color = c("hotpink", "darkorchid"))
#'
#' @keywords AMMI
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline par plot points text
#' @export

AMMIplot <- function(bilinearObject, plots = "linear", color = c("darkgreen", "darkblue"),  PC = 1, f = 0.5, ...){

	argumentChange <- function(defaultArgs, userArgs){
		userArgs <- list(...)
		defaultArgs[names(userArgs)] <- userArgs
		return(defaultArgs)
	}

	lineIntersect <- function(s1, i1, s2, i2){
		c( (i2 - i1) / (s1 - s2), (s1 * i2 - s2 * i1) / (s1 - s2))
	}

	Escores <- bilinearObject$scores$Escores
	Gscores <- bilinearObject$scores$Gscores

	Lambda <- bilinearObject$svdE$d
	mu <- bilinearObject$mu
	Geffect <- bilinearObject$Geffect
	Eeffect <- bilinearObject$Eeffect

	I <- length(Geffect)
	J <- length(Eeffect)

	M <- length(Lambda) - 1
	if(PC >= M) stop("Must specify 'PC' less than ", M)
	nmPC <- paste0("PC", PC)

	Kstar <- bilinearObject$sigPC
	model <- bilinearObject$model

	percExpl <- round(Lambda/sum(Lambda)*100,1)

	if(!all(plots %in% c("linear", "winner"))) stop("'plots' must be a character vector with arguments 'linear' and/or 'winner'")

	nplots <- length(plots)
	if(nplots > 1) par(mfrow = c(1, nplots))

	Gintercept <- mu + Geffect
	Erange <- range(Escores[,nmPC])

	x <- seq(Erange[1], Erange[2], length.out = 1000)
	nominalLines <- sweep(x %*% t(Gscores[,nmPC]), 2, Gintercept, "+")
	whichWin <- apply(nominalLines, 1, which.max)
	winnerPos <- names(Geffect)[whichWin]
	whichWin <- unique(whichWin)
	winner <- unique(winnerPos)

	if("linear" %in% plots){
		labelx <- tapply(x, factor(winnerPos, levels = winner), mean)
		labely <- NULL

		for (i in 1:length(winner)){
			labely <- c(labely, Gintercept[winner[i]] + Gscores[winner[i],nmPC] * labelx[i])
		}

		nominal <- bilinearObject$DF$Y - Eeffect[bilinearObject$DF$E]
		labely <- labely + diff(range(nominal)) * 0.05

		xlimits <- range(Escores[,nmPC]) * 1.1
		ylimits <- range(nominal - mu) * 1.1 + mu

		linecol <- rep("#000000", I)
		names(linecol) <- names(Geffect)
		linecol[winner] <- color[1]
		winnerMin <- winner[which.min(Gscores[winner, nmPC])]
		winnerMax <- winner[which.max(Gscores[winner, nmPC])]
		linecol[winnerMin] <- 2
		linecol[winnerMax] <- 4

		linewd <- rep(0.8, I)
		linewd[names(Geffect) %in% winner] <- 2

		linPlotArgs <- argumentChange(list(main = "Linear AMMI plot", ylim = ylimits, xlim = xlimits, yaxs="i", pch = 17,
							ylab = "Nominal", xlab = paste0("Environment PC", PC), cex = 1.5), list(...))

		do.call(plot, c(list(x = Escores[,nmPC], y = rep(ylimits[1], dim(Escores)[1])), linPlotArgs))

		# linLineArgs <- list(col = linecol[i], lwd = linewd[i])

		for (i in (1:length(Gintercept))[-whichWin]){
		  abline(Gintercept[i], Gscores[i, nmPC], col = linecol[i], lwd = linewd[i])
		}
		for (i in whichWin){
		  abline(Gintercept[i], Gscores[i, nmPC], col = linecol[i], lwd = linewd[i])
		}
		text(labelx, labely, labels = winner, col = linecol[winner])
	}
	points(Escores[, nmPC], nominal[bilinearObject$DF$G == winnerMin], pch=16, cex=0.6, col=2)
	points(Escores[, nmPC], nominal[bilinearObject$DF$G == winnerMax], pch=16, cex=0.6, col=4)

	if("winner" %in% plots){
		winnerInt <- list()
		for (i in 2:length(winner)){
			winnerInt[[i-1]] <- lineIntersect(Gscores[[winner[i-1],nmPC]], Gintercept[[winner[i-1]]], Gscores[[winner[i],nmPC]], Gintercept[[winner[i]]])
		}

		xlimits <- range(Eeffect) * 1.1 + mu
		ylimits <- range(Escores[, nmPC]) * 1.1

		winPLotArgs <- argumentChange(list(xlim = xlimits, ylim = ylimits, main = "Winner AMMI plot",
										   ylab = paste0("Environment PC", PC, "Score"), xlab = "Environmental mean",
										   pch = 16, col = color[2]), list(...))
		do.call(plot, c(list(x = mu + Eeffect, y = Escores[,nmPC]), winPLotArgs))
		abline(h = sapply(winnerInt, function(x) x[1]))
		text(mu + Eeffect * 1.1, Escores[,nmPC] * 1.1, labels = rownames(Escores), col = color[2])

		winPointArgs <- argumentChange(list(pch = 16, col = color[1]), list(...))
		do.call(points, c(list(x = mu + Geffect [names(Geffect) %in% winner], y = Gscores[names(Geffect) %in% winner, nmPC]), winPointArgs))
		text(mu + Geffect [names(Geffect) %in% winner] * 1.1, Gscores[names(Geffect) %in% winner,nmPC] * 1.1, labels = rownames(Gscores)[names(Geffect) %in% winner], col = color[1])
	}
	if(nplots > 1) par(mfrow = c(1, 1))
}


