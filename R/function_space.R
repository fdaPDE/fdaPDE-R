.FunctionSpaceCtr <- setRefClass(
  Class = "FunctionSpaceObject",
  fields = c(
    mesh  = "ANY", 
    fe_order="integer"
  )
)

# constructor
setGeneric("FunctionSpace", function(mesh,fe_order) standardGeneric("FunctionSpace"))
setMethod("FunctionSpace", signature = c(mesh="ANY", fe_order="numeric"),
          function(mesh,fe_order){
            if(fe_order == 1){
              return(.FunctionSpaceCtr(mesh=mesh, fe_order=1L))
            }else if(fe_order == 2){
              return(.FunctionSpaceCtr(mesh=mesh, fe_order=2L))
            }
                                                  
})

setMethod("FunctionSpace", signature = c(mesh="ANY", fe_order="missing"),
          function(mesh){
              return(.FunctionSpaceCtr(mesh=mesh, fe_order=1L))
})

## finite element function
.FunctionCtr <- setRefClass(
  Class = "FunctionObject",
  fields = c(
    FunctionSpace  = "ANY", 
    coeff = "matrix",
    pde = "ANY"
  ),
  methods = list(
    eval_at = function(X) {
      M = dim(FunctionSpace$mesh$nodes())[2]
      if(is.vector(X)) {
        pde$eval(FunctionSpace$mesh$data, coeff, matrix(X, nrow=1,ncol=M))
      } else {           
        if(dim(X)[2] != M) {
          stop(paste("matrix of evaluation points should be an N x", M, "matrix"))
        }
        pde$eval(FunctionSpace$mesh$data, coeff, as.matrix(X))
      }
    }
  )
)
## constructor
Function <- function(FunctionSpace) {
  coeff = matrix(ncol = 1, nrow = 0)
  .FunctionCtr(coeff = coeff, FunctionSpace = FunctionSpace)
}

## FunctionObject plot overload
setMethod("plot", signature=c(x="FunctionObject"), function(x, ...){
  times <- x$FunctionSpace$mesh$times
  is_parabolic <- FALSE
  if(length(times)!=0) is_parabolic <- TRUE
  if(!is_parabolic){
  plot_data <- data.frame(X=x$FunctionSpace$mesh$nodes()[,1], 
                          Y=x$FunctionSpace$mesh$nodes()[,2],
                          coeff=x$coeff[1:nrow(x$FunctionSpace$mesh$nodes())])
  I=x$FunctionSpace$mesh$elements()[,1]
  J=x$FunctionSpace$mesh$elements()[,2]
  K=x$FunctionSpace$mesh$elements()[,3]
  fig<- plot_ly(plot_data, x=~X, y=~Y, z=~coeff,
          i = I, j = J, k = K,
          intensity=~coeff,color = ~coeff, type="mesh3d", 
          colorbar=list(title=""), ...) %>%
      layout(scene = list(
        aspectmode = "data", 
        xaxis = list(
          title = '', showgrid = F, zeroline = F, showticklabels = F),
        yaxis = list(
          title = '', showgrid = F, zeroline = F, showticklabels = F),
        zaxis = list(
          title = '', showgrid = F, zeroline = F, showticklabels = F)),
        camera = list(
          eye = list(x = 1.25, y = -1.25, z = 1.25))) %>%
      colorbar(len = 1, title="")
  }else{
    plot_data <- data.frame(X=rep(x$FunctionSpace$mesh$nodes()[,1], times=length(times)), 
                            Y=rep(x$FunctionSpace$mesh$nodes()[,2], times=length(times)),
                            coeff=as.vector(x$coeff[1:nrow(x$FunctionSpace$mesh$nodes()),]),
                            times = rep(times, each=nrow(x$FunctionSpace$mesh$nodes())))
    limits = c(min(x$coeff), max(x$coeff))
    I=x$FunctionSpace$mesh$elements()[,1]
    J=x$FunctionSpace$mesh$elements()[,2] 
    K=x$FunctionSpace$mesh$elements()[,3]
    fig<- plot_ly(plot_data, x=~X, y=~Y, z=~coeff, frame=~times,
                  i = I, j = J, k = K, cmin = limits[1], cmax=limits[2],
                  intensity=~coeff,color = ~coeff, type="mesh3d",
                  colorbar=list(title=""), ...) %>%
        layout(scene = list(
          aspectmode = "data", 
          xaxis = list(
            title = '', showgrid = F, zeroline = F, showticklabels = F),
          yaxis = list(
            title = '', showgrid = F, zeroline = F, showticklabels = F),
          zaxis = list(
            title = '', showgrid = F, zeroline = F, showticklabels = F)),
          camera = list(
            eye = list(x = 1.25, y = -1.25, z = 1.25))) %>%
        colorbar(len = 1, title="") %>%
        animation_opts(frame=5) %>% 
        animation_slider(currentvalue = list(prefix ="t = "))
  }
  fig 
}
)

## FunctionObject contour overload
setMethod("contour", signature=c(x="FunctionObject"), function(x, ...){
  times <- x$FunctionSpace$mesh$times
  is_parabolic <- FALSE
  if(length(times)!=0) is_parabolic <- TRUE
  
  if(!is_parabolic){
  plot_data <- data.frame(X=x$FunctionSpace$mesh$nodes()[,1], 
                          Y=x$FunctionSpace$mesh$nodes()[,2],
                          coeff=x$coeff[1:nrow(x$FunctionSpace$mesh$nodes())])
  I=x$FunctionSpace$mesh$elements()[,1]
  J=x$FunctionSpace$mesh$elements()[,2] 
  K=x$FunctionSpace$mesh$elements()[,3]
  fig <- plot_ly(plot_data, type="contour", x=~X, y=~Y, z=~coeff, 
                 i = I, j = J, k = K,
                 intensity=~coeff, color = ~coeff,
                 contours=list(showlabels = TRUE),
                 colorbar=list(title=""), ...) %>%
    layout(xaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F),
           yaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F))
  }else{
    plot_data <- data.frame(X=rep(x$FunctionSpace$mesh$nodes()[,1], times=length(times)), 
                            Y=rep(x$FunctionSpace$mesh$nodes()[,2], times=length(times)),
                            coeff=as.vector(x$coeff[1:nrow(x$FunctionSpace$mesh$nodes()),]),
                            times = rep(times, each=nrow(x$FunctionSpace$mesh$nodes())))
    limits = c(min(x$coeff), max(x$coeff))
    I=x$FunctionSpace$mesh$elements()[,1]
    J=x$FunctionSpace$mesh$elements()[,2] 
    K=x$FunctionSpace$mesh$elements()[,3]
    fig <- plot_ly(plot_data, type="contour", x=~X, y=~Y, z=~coeff, frame=~times, 
                   i = I, j = J, k = K, zmin=limits[1], zmax=limits[2],
                   intensity=~coeff, color = ~coeff,
                   contours=list(showlabels = TRUE),
                   colorbar=list(title=""), ...) %>%
      layout(xaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F),
             yaxis = list(title = "", showgrid=F, zeroline=F, ticks="", showticklabels=F)) %>%
      animation_opts(frame=5) %>% 
      animation_slider(currentvalue = list(prefix ="t = "))
  }
  fig
}
)