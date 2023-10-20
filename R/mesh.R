## Mesh Class Definition
.MeshCtr <- setRefClass(
  Class = "MeshObject",
  fields = c(
    data  = "ANY",
    m = "integer",   # local dim
    n = "integer",   # embedding dim
    times = "vector" # for parabolic problems
  ),
  methods = c(
    nodes = function(){
      data$nodes()
    },
    elements = function(){
      data$elements()
    },
    boundary = function(){
      data$boundary()
    },
    neighbors = function(){
      data$neighbors()
    },
    times = function(){
      return(times)
    }
  )
)

setGeneric("Mesh", function(domain) standardGeneric("Mesh"))
setMethod("Mesh", signature = c(domain="list"),
          function(domain){
            domain$elements <- domain$elements - 1
            storage.mode(domain$elements) <- "integer"
            m <- ncol(domain$elements) - 1
            n <- ncol(domain$nodes) 
            if(m == 2 & n == 2)
              .MeshCtr(data=new(Mesh_2D, domain), m=as.integer(m),n=as.integer(n), 
                       times=vector(mode="double"))                                                  
            else
              stop("wrong input argument provided.")
})

setGeneric("%X%", function(op1, op2) standardGeneric("%X%"))
setMethod("%X%", signature=c(op1="MeshObject", op2="numeric"),
          function(op1, op2){
            if(op2[1] > op2[length(op2)])
              stop("Error! First time instant is greater than last time instant.")
            op1$times <- times
            op1          
})

## Mesh - auxiliary methods
unroll_edges_aux <- function(Mesh){
  mesh <- Mesh$data
  edges <- matrix(nrow=3*nrow(mesh$elements()), ncol=2)
  for(i in 1:nrow(mesh$elements())){
    edges[(3*(i-1) + 1),]   = mesh$elements()[i,c(1,2)] + 1
    edges[(3*(i-1) + 2),]   = mesh$elements()[i,c(2,3)] + 1
    edges[(3*(i-1) + 3),]   = mesh$elements()[i,c(3,1)] + 1
  }
  edges
}

setGeneric("unroll_edges", function(Mesh) standardGeneric("unroll_edges"))
setMethod("unroll_edges", "MeshObject", function(Mesh){
  unroll_edges_aux(Mesh)
})

plot_mesh_aux <- function(x, ...){
  edges <- unroll_edges(x)
  plot_ly(...) %>% 
    add_markers(x = x$data$nodes()[,1],
                y = x$data$nodes()[,2],
                color = I('black'), size = I(1),
                hoverinfo = 'text',
                text = paste('</br><b> Coordinates:', round(x$data$nodes()[,1],2),
                             round(x$data$nodes()[,2],2)),
                showlegend = T,
                visible = T) %>%
    add_segments(x = x$data$nodes()[edges[,1],1],
                 y = x$data$nodes()[edges[,1],2],
                 xend = x$data$nodes()[edges[,2],1],
                 yend = x$data$nodes()[edges[,2],2], 
                 color = I('black'), size = I(1),
                 showlegend = F) %>%
    layout(
      xaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F
      ),
      yaxis = list(
        title = '',
        showgrid = F,
        zeroline = F,
        showticklabels = F
      ))
}

setMethod("plot", "MeshObject", function(x, ...){
  plot_mesh_aux(x, ...)  
})

