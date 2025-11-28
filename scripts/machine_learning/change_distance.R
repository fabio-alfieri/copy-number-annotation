change_distance <- function(df) {
  distances <- c("dist.to.closest.OG", "dist.to.closest.TSG", "dist.to.closest.FGS",
                 "distance.to.centromere", "distance.to.telomere", "Ess.distance_pancancer")
  df[,distances] <- -1*df[,distances]
  return(df)
}