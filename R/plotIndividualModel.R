#' @export
plotIndividualModel <- function(modelList, num) {

  miniList <- lapply(modelList, "[[", num)

  time <- attributes(modelList)[['time']]
  cond <- attributes(modelList)[['cond']]
  loopCols <- attributes(modelList)[['loopCols']]
  mod$m$lhs()
  mod$m$fitted()

  plotDf <- data.frame(originalval = mod$m$lhs(),
                       predictedval = mod$m$fitted(),
                       time = time[1:7],
                       condition = cond[1:7])

  ggplot(data = plotDf) +
    geom_line(aes(x = time, y = predictedval)) +
    geom_point(aes(x = time, y = originalval), color = 'red') +
    theme_classic()


}
