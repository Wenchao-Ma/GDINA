# Calculating Item parameter location.
#
# @param Q A required \eqn{J \times K} binary Q-matrix. J represents test
#    length and K represents the number of attributes of this test.
#    Entry 1 at row j and column k represents the \eqn{k^{th}} attribute
#    is measured by item \eqn{j}, and 0 means item \eqn{j} does not
#    measure attribute \eqn{k}.
# @seealso \code{\link{alpha}}.
# @return An \eqn{J \times L} item location matrix. In the G-DINA model,
#    if an item j measures \eqn{Kj} attributes, L latent classes can
#    be combined into \eqn{2^{Kj}} latent groups. This matrix gives
#    which latent group each of L latent classes belongs to for item j.
#    The attribute pattern of \eqn{2^{Kj}} latent groups can be obtained
#    via \code{\link{alpha}} function.
#
eta.loc <- function(Q) {
  K <- ncol(Q)
  J <- nrow(Q) #

  pattern <- alpha(K, T, Q)

  L <- no_LC(Q)  # The number of latent groups


  par.loc <- matrix(0, J, L)

  for (j in 1:J) {
    # for each item
    loc <- which(Q[j, ] >= 1) #which attributes are required
    if (length(loc) == 1) {
      # if one attribute is required only
      reduced_pattern <-
        (pattern[, loc] >= Q[j, loc]) * 1 #reduced attributes
      par.loc[j, ] <- reduced_pattern + 1
    } else{
      #2 or more attribute
      reduced_pattern <- (t(pattern[, loc]) >= c(Q[j, loc])) * 1
      patternj <- t(alpha(length(loc)))
      par.loc[j, ] <-
        apply(reduced_pattern, 2, function(x)
          which.max(colSums(x == patternj)))
    }

  }
  ##it is a J x L matrix
  return (par.loc)
}
