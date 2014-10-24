
get.pve  <- function(model){
	anova.out <- anova(model);

	pve <- anova.out$Deviance[-1] / sum(anova.out$Deviance[-1], na.rm = TRUE);
	out.df <- data.frame(
		predictor = rownames(anova.out)[-1],
		pve = pve
		);
	return(out.df);
	}
