GrayResize <- function (imgList) {
  matList <- list()
  imgNum <- min(length(imgList), all)
  
  for(i in 1:imgNum) {
    if(i%%1000==0){
      cat("Reading images",i,"\n")
    }
    magic <- readBin(imgList[i], what = 0L, n = 8, size = 1L, signed = FALSE)
    judgejpg <- isTRUE(all.equal(magic[1:2], c(0xFF, 0xD8)))
    judgepng <- isTRUE(all.equal(magic[1:8], c(0x89,0x50,0x4E,0x47,0x0D, 0x0A, 0x1A, 0x0A)))
    
    if(judgejpg || judgepng){
      img <- readImage(imgList[i])
      grayImg <- channel(img,"gray")
      grayImg <- resize(grayImg,96,96)
      gray_Mat <- as.matrix(grayImg)
      gray_Mat <- round(gray_Mat,digits=2)
      
      matList[[i]] <- gray_Mat
    }
  }
  return(matList)
}


patch_extraction <- function(imgList,pNum,unit) {
  count <- 0
  patch_vec_list <- list()
  img_num <- length(imgList)
  
  num <- floor(pNum/img_num)
  
  for(i in 1:img_num) {
    for(j in 1:num) {
      posRow <- sample(1:(96-unit),1)
      posCol <- sample(1:(96-unit),1)
      
      patch <- imgList[[i]][posRow:(posRow+unit-1),posCol:(posCol+unit-1)]
      
      patch_vec <- as.vector(patch)
      count <- count + 1
      
      patch_vec_list[[count]] <- patch_vec
    }
  }
  patch_mat <- matrix(unlist(patch_vec_list),nrow = unit*unit,byrow = FALSE)
  return(patch_mat)
}


pre_process <- function(Mat) {
  sd <- sqrt(apply(Mat, 2, var))
  sub_mean <- t(t(Mat) - colMeans(Mat))
  sd[sd==0]<-1
  normMat <- sweep(sub_mean,2,sd,FUN = "/")
  normMat<-round(normMat, digits=2)
  return(normMat)
}

ZCAwhite <- function(x,epsilon) {
  row <- nrow(x)
  col <- ncol(x)
  sigma <- x %*% t(x) / col
  duv <- svd(sigma)
  xPCA <- diag(1./sqrt((duv$d)+epsilon)) %*% t(duv$u) %*% x
  xZCA <- duv$u %*% xPCA
  xZCA<-round(xZCA,digits=2)
  return(xZCA)
}

amapkmeans<-function (x, centers, iter.max = 10, nstart = 1, method = "euclidean") 
{
  dokmeans <- function() {
    Z <- .C("kmeans_Lloyd2", as.double(x), as.integer(m), 
            as.integer(ncol(x)), centers = as.double(centers), 
            as.integer(k), c1 = integer(m), iter = as.integer(iter.max), 
            nc = integer(k), wss = double(k), method = as.integer(method), 
            PACKAGE = "amap")
    if (Z$iter > iter.max) 
      warning("did not converge in ", iter.max, " iterations", 
              call. = FALSE)
    if (any(Z$nc == 0)) 
      warning("empty cluster: try a better set of initial centers", 
              call. = FALSE)
    Z
  }
  
  METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
               "binary", "pearson", "correlation", "spearman", "kendall", 
               "abspearson", "abscorrelation")
  
  method <- pmatch(method, METHODS)
  
  if (is.na(method)) 
    stop("invalid distance method")
  if (method == -1) 
    stop("ambiguous distance method")
  
  if (class(x) == "exprSet") {
    library(Biobase)
    x <- Biobase::exprs(x)
  }
  
  x <- as.matrix(x)
  m <- nrow(x)
  if (missing(centers)) 
    stop("'centers' must be a number or a matrix")
  
  if (length(centers) == 1) {
    k <- centers
    if (nstart == 1) 
      centers <- x[sample(1:m, k), , drop = FALSE]
    if (nstart >= 2 || any(duplicated(centers))) {
      cn <- unique(x)
      mm <- nrow(cn)
      if (mm < k) 
        stop("more cluster centers than distinct data points.")
      centers <- cn[sample(1:mm, k), , drop = FALSE]
    }
  }
  else {
    centers <- as.matrix(centers)
    if (any(duplicated(centers))) 
      stop("initial centers are not distinct")
    cn <- NULL
    k <- nrow(centers)
    if (m < k) 
      stop("more cluster centers than data points")
  }
  if (iter.max < 1) 
    stop("'iter.max' must be positive")
  
  if (ncol(x) != ncol(centers)) 
    stop("must have same number of columns in 'x' and 'centers'")
  
  Z <- .C("kmeans_Lloyd2", as.double(x), as.integer(m), as.integer(ncol(x)), 
          centers = as.double(centers), as.integer(k), c1 = integer(m), 
          iter = as.integer(iter.max), nc = integer(k), wss = double(k), 
          method = as.integer(method), PACKAGE = "amap")
  
  if (Z$iter > iter.max) 
    warning("did not converge in ", iter.max, " iterations", 
            call. = FALSE)
  if (any(Z$nc == 0)) 
    warning("empty cluster: try a better set of initial centers", 
            call. = FALSE)
  if (nstart >= 2 && !is.null(cn)) {
    best <- sum(Z$wss)
    for (i in 2:nstart) {
      centers <- cn[sample(1:mm, k), , drop = FALSE]
      ZZ <- dokmeans()
      if ((z <- sum(ZZ$wss)) < best) {
        Z <- ZZ
        best <- z
      }
    }
  }
  centers <- matrix(Z$centers, k)
  dimnames(centers) <- list(1:k, dimnames(x)[[2]])
  cluster <- Z$c1
  if (!is.null(rn <- rownames(x))) 
    names(cluster) <- rn
  out <- list(cluster = cluster, centers = centers, withinss = Z$wss, 
              size = Z$nc)
  class(out) <- "kmeans"
  out
}

matsplitter<-function(M, u, s) {
  cv<-matrix(NA,nrow=unit*unit,ncol=441,byrow=FALSE)
  rn<-(96-u)/s+1     #21
  group<-floor(rn/3)  #7
  count<-0
  for(r in 0:2){
    for(c in 0:2){
      r_off<-r*group*s+1
      c_off<-c*group*s+1
      for(x in 0:(group-1)){
        for(y in 0:(group-1)){
          x_off<-r_off+x*s
          y_off<-c_off+y*s
          mat<-M[x_off:(x_off+u-1),y_off:(y_off+u-1)]
          vec<-as.vector(mat)
          count<-count+1
          cv[,count]<-vec
        }
      }
    }
  }
  return(cv)
} 


euclidean_dist<-function(pMat,centroids){
  G<-t(pMat)%*%centroids
  Ss<-diag(t(pMat)%*%pMat)
  Rr<-diag(t(centroids)%*%centroids)
  n<-ncol(pMat)
  S<-matrix(Ss,nrow=n,ncol=k,byrow=FALSE)
  R<-matrix(Rr,nrow=n,ncol=k,byrow=TRUE)
  
  dis<-sqrt(S+R-2*G)
  return(dis)
}


feature_Generator <- function(inputList,centroids) {
  num <- length(inputList)
  fMat <- matrix(NA,nrow=9*k,ncol=num,byrow=FALSE)
  count<-0
  for(i in 1:num) {
    if(i%%100==0){
      cat("Converting images",i,"\n")
    }
    
    temp<-inputList[[i]]
    number<-ncol(imageData(temp))
    
    if(!is.null(number)){
      pMat <- matsplitter(temp,unit,s)
      
      pMat <- pre_process(pMat)
      dis<-euclidean_dist(pMat,centroids)
      udis<-rowMeans(dis)
      dif_dis<-dis-udis
      dif_dis[dif_dis<0]<-0
      fvList <- matrix(,nrow=k,ncol=9,byrow=FALSE)
      for(j in 0:8){
        row<-j*49+1
        temp<-dif_dis[row:(row+48),]
        local_sum<-colMeans(temp)
        fvList[,(j+1)] <- local_sum
      }
      count<-count+1
      vec<-as.vector(fvList)
      fMat[,count]<-vec
    }
  }
  fMat<-fMat[,1:count]
  return(fMat)
}

euclidean_dist_kmean<-function(pMat,centroids){
  G<-pMat%*%t(centroids)
  Ss<-diag(pMat%*%t(pMat))
  Rr<-diag(centroids%*%t(centroids))
  n<-nrow(pMat)
  k<-nrow(centroids)
  S<-matrix(Ss,nrow=n,ncol=k,byrow=FALSE)
  R<-matrix(Rr,nrow=n,ncol=k,byrow=TRUE)
  
  dis<-sqrt(S+R-2*G)
  return(dis)
}


mykmean<-function(x,centers,distFun,nIters)
{
  x<-as.matrix(x)
  kcenters<-x[sample(nrow(x),centers),]
  
  for(i in 1:nIters){
    if(distFun=="euclidean"){
      distsToCenters <- euclidean_dist_kmean(x,kcenters)
    }
    else if(distFun=="correlation"){
      temp<-cor(t(x),t(kcenters))
      temp[is.na(temp)]<-0
      distsToCenters <-1-temp
    }
    else{
      cat("Wrong distFun:euclidean or correlation")
    }
    
    clusters<-apply(distsToCenters,1,which.min)
    kcenters<-apply(x,2,tapply,clusters,mean)
    
  }
  return(kcenters)
}

mymodel<-function(testpath,model,kcentroids,offset,num,label)
{
  setwd(testpath)
  names <- list.files(pattern="*.jpg")
  mat<-list()
  for(i in offset:(offset+num-1)) {
    
    if(i%%1000==0){
      cat("Reading images",i,"\n")
    }
    magic <- readBin(names[i], what = 0L, n = 8, size = 1L, signed = FALSE)
    judgejpg <- isTRUE(all.equal(magic[1:2], c(0xFF, 0xD8)))
    judgepng <- isTRUE(all.equal(magic[1:8], c(0x89,0x50,0x4E,0x47,0x0D, 0x0A, 0x1A, 0x0A)))
    
    if(judgejpg || judgepng){
      img <- readImage(names[i])
      grayImg <- channel(img,"gray")
      grayImg <- resize(grayImg,96,96)
      gray_Mat <- as.matrix(grayImg)
      gray_Mat <- round(gray_Mat,digits=2)
      mat[[i-offset+1]] <- gray_Mat
    }
  }
  
  mat_feature <- feature_Generator(mat,kcentroids)
  test<-rbind(mat_feature,label)
  test_x<-test[-nrow(test),]
  test_y<-test[nrow(test),]
  test_x<-t(test_x)
  
  test_pred<-predict(model,test_x)
  acc<-sum(test_pred==test_y)/length(test_pred)
  return(acc)
}

mymodel_single<-function(testpath,imagename,model,kcentroids)
{
  setwd(testpath)
  
  mat<-list()
    magic <- readBin(imagename, what = 0L, n = 8, size = 1L, signed = FALSE)
    judgejpg <- isTRUE(all.equal(magic[1:2], c(0xFF, 0xD8)))
    judgepng <- isTRUE(all.equal(magic[1:8], c(0x89,0x50,0x4E,0x47,0x0D, 0x0A, 0x1A, 0x0A)))
    
    if(judgejpg || judgepng){
      img <- readImage(imagename)
      grayImg <- channel(img,"gray")
      grayImg <- resize(grayImg,96,96)
      gray_Mat <- as.matrix(grayImg)
      gray_Mat <- round(gray_Mat,digits=2)
      mat[[1]] <- gray_Mat
    }
    else
      return(0)
  
  
  mat_feature <- feature_Generator(mat,kcentroids)
  test_x<-mat_feature
  
  test_x<-as.matrix(test_x,nrow=length(test_x),ncol=1,byrow=T)
  
  test_pred<-predict(model,t(test_x))
  if(test_pred==1)
    return(+1)
  else
    return(-1)
  
  return(test_pred)
}
