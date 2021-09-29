#' Estimate the parameters in CCM
#'
#' The function for estimating the parameters in CCM.
#'
#' @usage EstimateCCM = function(profiles, phytree, ip=0.1, pen=0.5, ... )
#'
#' @param profiles a matrix containing the profiles. Columns are profiles and rows are species.
#' @param phytree a phylogenetic tree.
#' @param ip the initial values for optimizer. For a better convergence, a good set of starting values could be the estimates by a large tuning parameter.
#' @param pen the tuning parameter \eqn{\lambda}. Default value is 0.5. Value of 0 means no regularization.
#' @param ... control parameters available to optimizer `nlminb` such as `trace`, `rel.tol`, .... .Example See `?nlminb`, 'control' argument .
#' @return a list with following elements
#' \itemize{
#' \item alpha: estimated intrinsic rates.
#' \item B: estimated association matrix.
#' \item nlm.par: all estimated parameters in the order as `c(alpha, diag(B), B[upper.tri(B)])`.
#' \item nlm.converge: convergence message. See `?nlminb`.
#' \item nlm.hessian: estimated Hessian matrix used for calculating standard errors.
#' }
#' @references
#' Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
#' @references
#' David M. Gay (1990), Usage summary for selected optimization routines. Computing Science Technical Report 153, AT&T Bell Laboratories, Murray Hill.
#'
#' @examples
#' set.seed(123)
#' # generate a ranom 200-tip tree.
#' # Larger tree contains more information and tends to give better MLEs.
#' t <- rtree(200)
#' # setting arbatrunderlying parameters for a pair
#' n <- 2
#' alpha <- c(0.1, 0.1)
#' B <- matrix(0, n, n)
#' B[1,2] <- B[2,1] <- 1
#' diag(B) <- c(-0.3, 0.3)
#'
#' # using the same set of parameter to simulate 20 pairs
#' # and estimate with CCM.
#'
#' trueP = c(alpha, diag(B), B[upper.tri(B)]) # true parameters
#' estP = matrix(NA, nrow=20, ncol=length(trueP))
#' for (i in 1:20){
#'     simDF <- SimulateProfiles(t, alpha, B)
#'     aE <- EstimateCCM(profiles = simDF, phytree=t)
#'     estP[i,] = c(aE$alpha, diag(aE$B), aE$B[upper.tri(aE$B)])
#'     # or estP[i,] = aE$par same order as above
#'     print(i)
#' }
#'
#' # plot the estimates
#' boxplot(estP)
#' points(1:length(trueP), trueP, pch=8, col="red")
#'
#'
#' @export

EstimateCCM <- function(profiles, phytree, ip=0.1, pen=0.5,  ...){
    ### Main function for parameter estimation based on `ace()` from `ape` package.
    if (nrow(profiles) != length(phytree$tip.label)){
        stop("the profile matrix size doesn't match the number of tips on the tree.")
    }

    if (ncol(profiles)>7){
        warning("It may take much longer to estimate a large community. Consider to split it into small subcommunities.")
    }
    profiles = profiles[phytree$tip.label,]
    n = ncol(profiles)
    lvls = apply(expand.grid(replicate(n, c(0,1), simplify = F)),1,paste0, collapse="")
    x <- apply(profiles, 1, paste0, collapse="")
    x = factor(x, levels = lvls)
    nl <- nlevels(x)
    x <- as.integer(x)

    # extract information from tree
    phy=phytree
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    liks <- matrix(0, nb.tip + nb.node, nl)
    TIPS <- 1:nb.tip
    liks[cbind(TIPS, x)] <- 1
    phy <- reorder(phy, "postorder")
    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]
    EL <- phy$edge.length

    ## construct likelihood function and estimate parameters
    dev <- function(p) {
        if (any(is.nan(p)) || any(is.infinite(p)))
            return(1e+50)
        comp <- numeric(nb.tip + nb.node)
        C = matrix(0, nrow=n, ncol=n)

        # assign initial values
        p0s = p[1:n]
        gldifs = p[(n+1):(2*n)]
        C[upper.tri(C, diag=F)] = p[(2*n+1):length(p)]
        diag(C) = gldifs
        C[lower.tri(C)] = t(C)[lower.tri(C)]
        Q = ConstructQunsym(C, p0s)
        diag(Q) <- -rowSums(Q)
        decompo <- eigen(Q)
        lambda <- decompo$values
        GAMMA <- decompo$vectors
        invGAMMA <- solve(GAMMA)
        for (i in seq(from = 1, by = 2, length.out = nb.node)) {
            j <- i + 1L
            anc <- e1[i]
            des1 <- e2[i]
            des2 <- e2[j]
            v.l <- GAMMA %*% diag(exp(lambda * EL[i])) %*%
                    invGAMMA %*% liks[des1, ]
            v.r <- GAMMA %*% diag(exp(lambda * EL[j])) %*%
                    invGAMMA %*% liks[des2, ]
            v <- v.l * v.r
            comp[anc] <- sum(v)
            liks[anc, ] <- v/comp[anc]
        }
        dev <- -1 * sum(log(comp[-TIPS]))
        penalty = pen*sum( p^2)
        dev=dev+penalty
        if (is.na(dev))
            Inf
        else dev

    }

    # estimate the rates
    np=n*(n+1)/2 + n
    obj <- list()
    if (length(ip) == 1){
        ip = rep(ip, length.out = np)
    }
    iter.history <- capture.output(
        out <- nlminb(ip, function(p) dev(p),
                      control = list(...)), split=T)
    obj$loglik <- -out$objective/2
    obj$rates <- out$par
    out.nlm <- try(nlm(function(p) dev(p), p = obj$rates, iterlim = 1, stepmax = 0, hessian = TRUE), silent = TRUE)
    Cestimate = matrix(0, nrow=n, ncol=n)
    Cestimate[upper.tri(Cestimate,diag=F)] = obj$rates[(2*n+1):length(obj$rates)]
    Cestimate[lower.tri(Cestimate)] = t(Cestimate)[lower.tri(Cestimate)]
    diag(Cestimate) = obj$rates[(n+1):(2*n)]

    return(list(alpha = obj$rates[1:n], B = Cestimate, nlm.par=out$par, nlm.converge = out$convergence, nlm.hessian=out.nlm$hessian))

}

# construct Q matrix
ConstructQunsym <- function(C, c0s){
    n = nrow(C)
    mutdif = diag(C)
    diag(C) = 0
    statesM = expand.grid(replicate(n, c(0,1), simplify = F))
    Q = matrix(0, nrow=nrow(statesM), ncol=nrow(statesM))
    rownames(Q) <- colnames(Q)<- apply(statesM, 1, paste0, collapse=",")

    ### convert binary statesM to {-1, 1}
    for (i in 1:nrow(statesM)) {
        statei = as.numeric(statesM[i, ])

        repi =  matrix(rep(statei,length(statei)), nrow=length(statei), byrow = T)

        diag(repi) = 1 -  diag(repi)
        rowlable = paste0(statei, collapse=",")
        collable = apply(repi,1, paste0, collapse=",")

        diag(repi) = statei
        # same states:
        sm = (repi == statei)
        dm = (!sm)

        qs = (c0s - ifelse(statei==0,-1,1)*mutdif)  - diag(sm %*% C) + diag(dm %*% C)
        qs = exp(qs)
        Q[rowlable, collable] = qs
    }

    return(Q)
}


