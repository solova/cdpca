if(!require('matlib')) {
  install.packages('matlib')
  library('matlib')
}

CDpca <- function(data, class = NULL,  P, Q, SDPinitial = FALSE,  tol=10^(-5), maxit, r, cdpcaplot= TRUE) {
  # CDpca performs a clustering and disjoint principal components analysis
  # on the given numeric data matrix and returns a list of results
  #
  # Args:
  #  data: data frame (numeric).
  #  class: vector (numeric) or 0, if classes of objects are unknown.
  #  fixAtt: vector (numeric) or 0, for a selection of attributes.
  #  SDPinitial: 1 or 0, for random initialization of U and V
  #  nnloads: 1 or 0, for nonnegative loadings
  #  cdpcaplot: integer: 1 or 0, if no plot is displayed.
  #  Q: integer, number of clusters of variables.
  #  P: integer, number of clusters of objects.
  #  tol: real number, small positive tolerance.
  #  maxit: integer, maximum of iterations.
  #  r: number of runs of the cdpca algoritm for the final solution.
  #
  # Returns:
  #  iter: iterations used in the best loop for computing the best solution
  #  loop: best loop number
  #  timebestloop: computation time on the best loop
  #  timeallloops: computation time for all loops
  #  Y: the component score matrix
  #  Ybar: the object centroids matrix in the reduced space
  #  A: the component loading matrix
  #  U: the partition of objects
  #  V: the partition of variables
  #  F: function to maximize
  #  bcdev: between cluster deviance
  #  bcdevTotal: between cluster deviance over the total variability
  #  tableclass: cdpca classification
  #  pseudocm: pseudo confusion matrix of the real and cdpca classifications
  #  Enorm: error norm for the obtained cdpca model

  #CDPCA is applied on normalized data (mean zero and unit variance)

  #############################################
  #             Function RandMat              #
  #############################################

  RandMat <- function(dim1, dim2) {
    # Generates a random binary and row stochastic matrix.
    #
    # Args:
    #  dim1: number of rows.
    #  dim2: number of columns.
    #
    # Returns:
    #  The random matrix (dim1xdim2) with only one nonzero element per row.
    #
    U0 <- matrix(0, dim1, dim2)
    U0[1:dim2, 1:dim2] <- diag(dim2)
    for (j in (dim2+1):dim1){
      p <- sample(1:dim2, 1)
      U0[j, p] <- 1
    }
    U0
  }  # end RandMat function

  ######################################
  #
  # SDP initialization
  #
  ######################################
  SDPinitialization <- function(dim1,dim2,D,data){
    # Generates a binary and row stochastic matrix
    # using SDP-based approach
    #
    # Args:
    #  dim1: number of rows.
    #  dim2: number of columns.
    #
    # Returns:
    #  The random matrix (dim1xdim2) with only one nonzero element per row.
    #
    # approximate solution for clustering
    em <- cbind(rep(1,dim1))
    matIem <- diag(dim1)-(1/dim1)*em%*%t(em)
    if (dim1 == dim(data)[1]){
      data <- D
    }else{
      data <- t(D)
    }
    So <- matIem%*%data%*%t(data)%*%matIem
    matF <- matrix(svd(So)$v[, 1:(dim2-1)], ncol = dim2-1)
    matFFT <- matF%*%t(matF)
    # approximate solution to the SDP-based model for clustering
    Z.bar <- matFFT + 1/dim1 *em%*%t(em)

    # refine solutions
    cent <- Z.bar%*%data
    centunique <- unique(cent)
    t <- dim(centunique)[1]

    # initial kmeans step
    indcentroids <- sort(sample(1:t, dim2, replace=FALSE), decreasing=FALSE)
    centroids <- centunique[indcentroids, ]
    solkmeans <- kmeans(data, centroids, iter.max = 100000)
    M <- matrix(0, dim1, dim2)
    for (i in 1:dim1){
      M[i, solkmeans$cluster[i]] <- 1
    }
    M
  }



  data.sd <- data.frame(scale(data, center = TRUE, scale = TRUE))

  I <- dim(data)[1]  					# number of objects I
  J <- dim(data)[2]  					# number of variables J
  Xs <- as.matrix(data.sd*sqrt(I/(I-1)))  	# matrix of normalized data (with variance divided by I)
  # matriz X (I x J) (obs x vars), numeric matrix argument
  tfbest <- matrix(0, r, 1)  # computational time at each loop

  # Run CDPCA ALS algorithm r times
  for (loop in 1:r) {
    t1 <- proc.time()
    # Initialization
    iter <- 0

    # SDP initialization or not?
    if (SDPinitial){
      U <- SDPinitialization(I, P, Xs, data)
      V <- SDPinitialization(J, Q, Xs, data)
    }else{
      U <- RandMat(I, P)
      V <- RandMat(J, Q)
    }

    #Set of the variables

    X.bar <- diag(1/colSums(U))%*%t(U)%*%Xs  	# (PxJ) object centroid matrix. It identifies the P centroids in the variable space
    X.group <- U%*%X.bar  				# (IxJ) matrix. Each object in the data matrix is replaced by its centroid.
    vJ <- matrix(1:J, ncol = 1)
    zJ <- matrix(0, J, 1)
    zQ <- matrix(0, 1, Q)

    # matrix A
    A <- matrix(0, J, Q)
    for (j in 1: Q) {
      Jxx <- which(V[, j] == 1)
      len.Jxx <- length(Jxx)
      if (sum(V[, j])>1) {
        if (I >= len.Jxx) {
          # CASE 1: I >= J (more observations than variables)
          S.group <- t(X.group[, Jxx])%*%X.group[, Jxx]
          A[Jxx, j] <- powerMethod(S.group)$vector
        }else{
          # CASE 2: I < J ( more variables than observations)
          SS.group <- X.group[, Jxx]%*%t(X.group[, Jxx])
          PMinitial <- powerMethod(SS.group)
          A[Jxx, j] <- t(X.group[, Jxx])%*%PMinitial$vector/sqrt(PMinitial$value)
        }  # end if
      }else{
        A[Jxx, j] <- 1
      }  # end if
    }  # end for

    A0 <- A
    Y.bar <- X.bar%*%A  				# (PxQ) object centroid matrix. It identifies the P centroids in the reduced space of the Q principal components.
    F0 <- sum(diag(t(U%*%Y.bar)%*%U%*%Y.bar))   # Objective function to maximize.
    Fmax <- F0
    conv <- 2*tol
    while (conv > tol) {
      iter <- iter+1
      Y <- Xs%*%A 					 # (IxQ) component score matrix. It identifies the I objects in the reduced space of the Q principal components.
      U <- matrix(0, I, P)
      # update U
      for (i in 1:I) {
        dist <- rep(0, P)
        for (p in 1:P) {
          dist[p] <- sum((Y[i, ]-Y.bar[p, ])^2)
        } # end for
        min.dist <- which.min(dist)
        U[i, min.dist] <- 1
      }  # end for
      su <- colSums(U)
      while (sum(su == 0) > 0) {
        ind.max <- which.max(su)  # p2
        su.max <- max(su)  # m2
        ind.min <- which.min(su)  # p1
        su.min <- min(su)  # m1
        ind.nzU <- which(U[, ind.max] == 1)
        ind.sel <- ind.nzU[1:floor(su.max)/2]
        U[ind.sel, ind.min] <- 1
        U[ind.sel, ind.max] <- 0
        su <- colSums(U)
      }  # end while

      # Given U and A compute X.bar
      X.bar <- diag(1/colSums(U))%*%t(U)%*%Xs
      Y.bar <- X.bar%*%A
      X.group <- U%*%X.bar

      # Update V and A
      for (j in 1:J) {

        posmax <- which(V[j, ] == 1)
        for (g in 1:Q) {
          V[j, ] <- diag(Q)[g, ]
          ####
          if (g!=posmax) {  ### if 1
            for (gg in c(g,posmax)) {
              xx <- V[, gg]
              Jxx <- which(xx == 1)
              len.Jxx <- length(Jxx)
              A[, gg] <- zJ
              if (sum(xx) > 1) {
                if (I >= len.Jxx) {
                  # CASE 1: I >= J (more observations than variables)
                  S.group <- t(X.group[, Jxx])%*%X.group[, Jxx]
                  A[Jxx, gg] <- matrix(powerMethod(S.group)$vector, nrow = len.Jxx)
                }else{
                  # CASE 2: I < J (more variables than observations)
                  SS.group <- X.group[, Jxx]%*%t(X.group[, Jxx])
                  PMgeneral <- powerMethod(SS.group)
                  A[Jxx, gg] <- matrix((t(X.group[, Jxx])%*%PMgeneral$vector/sqrt(PMgeneral$value))[, 1], nrow = len.Jxx)
                }  # end if
              }else{
                if (sum(xx)==1) {
                  A[Jxx, gg] <- 1
                }  # end if
              }  # end if
            }  # end for
          }  # end ### if 1
          ####
          Y.bar <- X.bar%*%A
          F <- sum(diag(t(U%*%Y.bar)%*%U%*%Y.bar))
          if (F > Fmax) {
            Fmax <- F
            posmax <- g
            A0 <- A
          }else{
            A <- A0
          }  # end if
        }  # end for
        V[j, ]=diag(Q)[posmax, ]
      }  # end for


      Y <- Xs%*%A
      Y.bar <- X.bar%*%A
      F <- sum(diag(t(U%*%Y.bar)%*%U%*%Y.bar))
      conv <- F-F0
      if (conv > tol) {
        F0 <- F
        A0 <- A
      }else{
        break
      }  # end if
      if (iter == maxit) {
        print("Maximum iterations reached.")
        break
      }  # end if

    }  # end while

    # Computation time for each loop
    t2 <- proc.time()
    tfinal <- t2-t1
    tfbest[loop] <- tfinal[3]

    #Results to be observed in each run of CDpca
    # BetClusDevTotal <- F/(I*J)*100  			# between cluster deviance (two different formula for the same thing)
    BetClusDev <- (F/sum(diag(t(Y)%*%Y)))*100  		# between cluster deviance
    tabperloop <- data.frame(cbind(loop, iter, tfinal[3], BetClusDev, F, conv))
    rownames(tabperloop) <- c(" ")
    colnames(tabperloop) <- c("Loop", "Iter", "Loop time", "Between cluster deviance(%):", "F", "Convergence")
    print(tabperloop)  # Loop display
    if (loop == 1) {
      Vcdpca <- V
      Ucdpca <- U
      X.barcdpca <- diag(1/colSums(Ucdpca))%*%t(Ucdpca)%*%Xs
      Acdpca <- A
      Ycdpca <- Xs%*%Acdpca
      Y.barcdpca <- X.barcdpca%*%Acdpca
      Fcdpca <- F
      loopcdpca <- 1
      itercdpca <- iter
      convcdpca <- conv
    }  # end if
    if (F > Fcdpca) {
      Vcdpca <- V
      Ucdpca <- U
      X.barcdpca <- diag(1/colSums(Ucdpca))%*%t(Ucdpca)%*%Xs
      Acdpca <- A
      Ycdpca <- Xs%*%Acdpca
      Y.barcdpca <- X.barcdpca%*%Acdpca
      Fcdpca <- F
      loopcdpca <- loop
      itercdpca <- iter
      convcdpca <- conv
    }  # end if
  }  # end for loop

  # Computation time for all loops
  tftotal=sum(tfbest)

  # maximum between cluster deviance
  BetClusDevTotal <- Fcdpca/(I*J)*100  # (sum(diag(var(Y)))*(I-1))*100
  BetClusDev <- (sum(diag(t(Ucdpca%*%Y.barcdpca)%*%(Ucdpca%*%Y.barcdpca)))/sum(diag(t(Ycdpca)%*%Ycdpca)))*100

  # error in cdpca model
  Enormcdpca <- 1/I*norm(Xs-Ucdpca%*%Y.barcdpca%*%t(Acdpca), "F")

  # Table Real vs CDPCA classification
  classcdpca <- Ucdpca%*%as.matrix(1:ncol(Ucdpca))
  if (!is.null(class)){
    class <- data.frame(class)
    maxclass <- max(class)
  }

  # Variability of Ycdpca and in decreasing order
  varYcdpca <- var(Ycdpca)
  d <- round(diag(varYcdpca)*100/J, 2)
  dorder <- d[order(d, decreasing = TRUE)]
  tablevar <- data.frame(1:Q, d)
  colnames(tablevar) <- c("Dim", "Var (%)")


  # Presentation of the matrices associated to the CDPCA component sorted by their variances.
  Yorder <- t(t(rbind(d, Ycdpca))[order(d, decreasing = TRUE), ])[-1, ]
  Ybarorder <- t(t(rbind(d, Y.barcdpca))[order(d, decreasing=TRUE),])[-1, ]
  Aorder <- t(t(rbind(d, Acdpca))[order(d, decreasing = TRUE), ])[-1, ]
  Vorder <- t(t(rbind(d, Vcdpca))[order(d, decreasing = TRUE), ])[-1, ]
  #
  # We can check the model using these column sort matrices and observe that
  # Ucdpca%*%Y.barcdpca%*%t(Acdpca) = Ucdpca%*%Ybarorder%*%t(Aorder)

  if ( is.null(class) ) {
    realclasstrue <- 2
    pseudocm <- NULL
    tabrealcdpca <- data.frame(classcdpca)
    colnames(tabrealcdpca) <- c("CDPCA Class")
  }else{
    realclasstrue <- 3
    tabrealcdpca <- data.frame(class, classcdpca)
    colnames(tabrealcdpca) <- c("Real Class", "CDPCA Class")
    # pseudo-confusion matrix
    pseudocm <- table(tabrealcdpca)
  }  # end if


  # PLOTS: CDPCA classification
  if (cdpcaplot) {
    displaygraphics <- par(no.readonly = TRUE)
    par(mfrow = c(1, realclasstrue))
    if (realclasstrue == 3) {
      # plot1
      matplot(Yorder[, 1], Yorder[, 2], xlab = "Dim 1", ylab = "Dim 2", type = "n", main = "Real classification")
      for (i in 1:maxclass) { points(Yorder[, 1][class == i], Yorder[, 2][class == i], pch = 1, col = i+1) }
      # plot2
      matplot(Yorder[, 1], Yorder[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"), type = "n", main = "CDPCA classification")
      for (i in 1:P) { points(Yorder[, 1][classcdpca == i], Yorder[, 2][classcdpca == i], pch = 2, col = i+1) }
      points(Ybarorder[, 1], Ybarorder[, 2],pch = 15)  # introduce the centroids into the plot
      # plot3 (legend)
      matplot(Y.barcdpca, type = "n", axes = FALSE, xlab = "", ylab = "")
      legend("topleft", pch = 1, col=c(2:(maxclass+1)), legend = c(1:maxclass), ncol = maxclass)
      legend("bottomleft", pch = 2, col=c(2:(P+1)), legend = 1:P, ncol = P)
    }else{
      # plot2
      matplot(Yorder[, 1], Yorder[, 2], xlab = paste(" Dim 1  (", (dorder[1]), " % )"), ylab = paste(" Dim 2  (", (dorder[2]), " % )"),type = "n", main = "CDPCA classification")
      for (i in 1:P) { points(Yorder[, 1][classcdpca == i], Yorder[, 2][classcdpca == i], pch = 2, col = i+1)}
      points(Ybarorder[, 1],Ybarorder[, 2], pch = 15)  # introduce the centroids into the plot
      # plot3 (legend)
      matplot(Y.barcdpca, type = "n", axes = FALSE, xlab = "", ylab = "")
      legend("bottomleft", pch = 2, col=c(2:(P+1)), legend = 1:P, ncol = P)
    }	 # end if
  } #end if cdpcaplot


  # OUTPUT
  list(timeallloops = tftotal,
       timebestloop = tfbest[loopcdpca],
       loop = loopcdpca,
       iter = itercdpca,
       bcdevTotal = BetClusDevTotal ,
       bcdev = BetClusDev ,
       V = Vorder,
       U = Ucdpca,
       A = Aorder,
       F = Fcdpca,
       Enorm = Enormcdpca,
       tableclass = tabrealcdpca,
       pseudocm = pseudocm,
       Y = Yorder,
       Ybar = Ybarorder,
       dorder = dorder,
       Xscale = Xs)

}  # end CDpca function
