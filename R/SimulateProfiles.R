#' @title  Simulate profiles
#'
#' @description  Simulate binary profiles based on the Community Coevolution Model
#' using the input tree and user-defined parameters
#'
#' @name SimulateProfiles
#' @usage SimulateProfiles(phytree, alpha, B, root=F)
#' @param phytree a phylogenetic tree.
#' @param alpha a vector of the intrinsic rates (\eqn{\alpha}) in CCM. The length of the vector decides how many profiles to simulate.
#' @param B a symmetric association matrix (\eqn{\Beta}) in CCM. If the input matrix is not symmetric,
#' the upper triangle of the matrix will be used.
#' @param root a vector of states at root node. `root=F` as default (random states will be assigned at root).
#' @return a matrix containing simulated profiles.
#' @examples
#' set.seed(123)
#' # generate a random 50-tip tree
#' t <- rtree(50)
#' # simulate 5 profiles
#' n <- 5
#' # assigning parameter values.
#' # The parameters should be in a reasonable scale otherwise the simulated profiles may have all 0s or 1s.
#' alpha =runif(n, -0.5, 2) # assign n random intrinsic rates.
#' association matrix B
#' B <- matrix(0,n,n)
#' B[1:3,1:3] <- 1 # first 3 genes are correlated. The other 2 have no interaction with others
#' diag(B) <- runif(n, -0.5,0.5) # assign half the difference between gain and loss rates for each gene
#' # simulate 5 random profiles
#' simDF <- SimulateProfiles(t, alpha, B)
#' # plot the profiles
#' d <- TreeToDend(t) # convert the tree to dendrogram
#' ProfilePlot(simDF, d)
#' @export
### Simulate n profiles
SimulateProfiles = function(phytree, alpha, B, root=F){
    if (nrow(B)!=ncol(B)){
        stop("the association matrix B should be a square matrix")
    }
    if (nrow(B)!=ncol(B) | length(alpha)!=nrow(B)){
        stop("the dimemsion of the asscoiation matrix and the number of intrinsic rates don't match")
    }

    if (!isSymmetric(B)){
        B[lower.tri(B)] = t(B)[lower.tri(t(B))]
    }

    # tree information
    phy=phytree
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    Nedges <- dim(phy$edge)[1]
    ROOT <- ntips + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    el <- phy$edge.length

    n = length(alpha)
    N <- dim(phy$edge)[1]
    if (is.numeric(root)){
        S = root
    } else {
        S = sample(c(0,1), n, replace=T)
    }

    ## generate the profiles
    xprofiles = matrix(NA, nrow=N+1, ncol=n)
    xprofiles[ROOT,] = S
    for (i in N:1) {

        Sanc = xprofiles[anc[i], ]
        evolbranch = SimulateOneBranch(S=Sanc, l = el[i], alpha, B)
        xprofiles[des[i],] = evolbranch$Snew

    }
    simProfiles = xprofiles[1:ntips,]
    rownames(simProfiles)=phy$tip.label
    colnames(simProfiles)=paste0("c", 1:n)
    return(simProfiles)

}


UpdateRates = function(S, alpha, B){
    mutdif = diag(B)
    diag(B) = 0
    if(!isSymmetric(B)){
        B[lower.tri(B)] = t(B)[lower.tri(t(B))]
    }
    S = ifelse(S==0,-1,1)
    newRates = sapply(seq_along(S), function(x)   exp(alpha[x]- S[x]*mutdif[x] -sum(B[x, S == S[x]]) + sum(B[x, S!=S[x]]) )  )
    return(newRates)

}

### simulation on one branch

SimulateOneBranch = function(S, l, alpha, B){
    n = length(S)
    lambdas = UpdateRates(S, alpha, B)
    history = S
    changeTime = c()
    repeat{
        et = rexp(n, lambdas)
        minet = min(et) * 1
        if (l < minet){
            break
        }
        l = l- minet
        changeT = which(et <= minet) ## <= ==
        changeTime = c(changeTime, minet)
        S[changeT] = as.numeric(!S[changeT])
        history=rbind(history, S)
        lambdas = UpdateRates(S, alpha, B)

    }
    return(list(Snew = S, lambdasNew = lambdas, history=history, changeTime=changeTime))

}


