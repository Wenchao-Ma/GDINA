#' Tatsuoka's fraction subtraction data
#'
#' Fraction Subtraction data (Tatsuoka, 2002) consists of responses of 536 examinees to 20 items measuring 8 attributes.
#'
#' @format A list of responses and Q-matrix with components:
#' \describe{
#' \item{\code{dat}}{responses of 536 examinees to 20 items}
#' \item{\code{Q}}{The \eqn{20 \times 8} Q-matrix}
#' }
#' @examples
#' \dontrun{
#' mod1 <- GDINA(frac20$dat,frac20$Q,model="DINA")
#' mod1
#' summary(mod1)
#' # Higher order model
#' mod2 <- GDINA(frac20$dat,frac20$Q,model="DINA",att.dist="higher.order")
#' mod2
#' anova(mod1,mod2)
#' }
#' @references
#' Tatsuoka, C. (2002). Data analytic methods for latent partially ordered classification models. \emph{Journal of the Royal Statistical Society, Series C, Applied Statistics, 51}, 337-350.
"frac20"
