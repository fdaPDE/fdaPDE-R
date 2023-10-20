# Pde Class wraps C++ R_PDE class
.PdeCtr <- setRefClass(
  Class = "PdeObject",
  fields = c(
   # L  = "ANY",                         # DiffOpObject
    is_dirichletBC_set = "logical",      
    times =       "vector",              # vector of time instants (parabolic problem only)
    is_initialCondition_set = "logical", # initial condition       (parabolic problem only)
    pde_ = "ANY" ,                       # Wraps of C++ class
    is_parabolic = "logical",
    is_init = "logical"                  
  ),
  methods = c(
    solve = function(){
      if(!is_dirichletBC_set){
        warning("dirichletBC not setted. Assuming homogeneus boudary condition.")
        if(!is_parabolic){
          dirichletBC_ = as.matrix(rep(0,times=nrow(pde_$get_dofs_coordinates())))
        }else 
          dirichletBC_ = matrix(0, nrow=nrow(pde_$get_dofs_coordinates()), ncol=length(times))
      pde_$set_dirichlet_bc(dirichletBC_)
      is_dirichletBC_set <<- TRUE
      }
      if(is_parabolic & (!is_initialCondition_set)){
        stop("initialCondition must be provided.")
      }
      if(!is_init) pde_$init()
      pde_$solve()
    },
    solution = function(){
      pde_$solution()
    },
    get_dofs_coordinates = function(){
      pde_$get_dofs_coordinates()
    },
    set_dirichletBC = function(dirichletBC){
      if(!is_parabolic){
        dirichletBC_ <- as.matrix(dirichletBC(pde_$get_dofs_coordinates()))
      }else{
        
        dirichletBC_ <- dirichletBC(pde_$get_dofs_coordinates(),times)
        pde_$set_dirichlet_bc(dirichletBC_)
      }
      is_dirichletBC_set <<- TRUE
      pde_$set_dirichlet_bc(dirichletBC_)
    },
    set_initialCondition = function(initialCondtion){
      if(!is_parabolic)
        stop("Cannot set initial condition for elliptic problem.")
      is_initialCondition_set <<- TRUE
      pde_$set_initial_condition(initialCondition(pde_$get_dofs_coordinates()))
    },
    get_mass = function(){
      pde_$get_mass()
    },
    get_stiff = function(){
      pde_$get_stiff()
    }
  )
)


setGeneric("Pde", function(L,u,dirichletBC, initialCondition) standardGeneric("Pde"))
setMethod("Pde", signature=c(L="DiffOpObject", u="ANY", dirichletBC="ANY", 
                             initialCondition="missing"),
          function(L,u,dirichletBC){
            D = L$f$FunctionSpace$mesh$data ## C++ R_Mesh class
            
            is_parabolic = FALSE
            if( length(L$f$FunctionSpace$mesh$times) != 0 ) is_parabolic = TRUE
            times <- L$f$FunctionSpace$mesh$times
            
            ## set pde type
            pde_type <- 0
            pde_parameters <- NULL
            if ("diffusion" %in% names(L$params) & !is.matrix(L$params$diffusion)){
              ## specialized implementation for laplace operator
              pde_type <- 1L 
              
              pde_parameters$diffusion <- 0.0
              
            } else {
              ## general diffusion-transport-reaction problem, constant coefficients
              pde_type <- 2L
              
              pde_parameters$diffusion <- matrix(0, nrow = 2, ncol = 2)
              
            }
            
            ## prepare pde_parameters list
            pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
            pde_parameters$reaction  <- 0.0
            pde_parameters$time <- 0L
            for(i in 1:length(L$tokens)) {
              pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
            }
            
            fe_order <- L$f$FunctionSpace$fe_order
            ## define Rcpp module
            pde_ <- NULL
            if (fe_order == 1) { ## linear finite elements
              pde_ <- new(PDE_2D_ORDER_1, D, pde_type, pde_parameters, L$f)
            }
            if (fe_order == 2) { ## quadratic finite elements
              pde_ <- new(PDE_2D_ORDER_2, D, pde_type, pde_parameters, L$f)
            }
            
            L$f$pde = pde_
            
            quad_nodes <- as.matrix(pde_$get_quadrature_nodes())
            ## evaluate forcing term on quadrature nodes
            pde_$set_forcing(as.matrix(u(quad_nodes)))
            
            ## initialize solver 
            pde_$init()
            is_init = TRUE
            
            ## set Dirichlet BC 
            dirichletBC_ <- as.matrix(dirichletBC(pde_$get_dofs_coordinates()))
            pde_$set_dirichlet_bc(dirichletBC_)
            is_dirichletBC_set = TRUE
            
            ## return
            .PdeCtr(times = times,
                    is_dirichletBC_set = is_dirichletBC_set,
                    is_initialCondition_set = FALSE,
                    pde_ = pde_,                  # Wraps of C++ class
                    is_parabolic = is_parabolic,
                    is_init=is_init)
          })

setMethod("Pde", signature=c(L="DiffOpObject", u="ANY", dirichletBC="ANY", 
                             initialCondition="ANY"),
          function(L,u,dirichletBC, initialCondition){
            D = L$f$FunctionSpace$mesh$data ## C++ R_Mesh class
            
            is_parabolic = FALSE
            if( length(L$f$FunctionSpace$mesh$times) != 0 ) is_parabolic = TRUE
            times <- L$f$FunctionSpace$mesh$times
            
            ## set pde type
            pde_type <- 0
            pde_parameters <- NULL
            if ("diffusion" %in% names(L$params) & !is.matrix(L$params$diffusion)){
              ## specialized implementation for laplace operator
              pde_type <- 3L 
              
              pde_parameters$diffusion <- 0.0
              
            } else {
              ## general diffusion-transport-reaction problem, constant coefficients
              pde_type <- 4L
              
              pde_parameters$diffusion <- matrix(0, nrow = 2, ncol = 2)
              
            }
            
            ## prepare pde_parameters list
            pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
            pde_parameters$reaction  <- 0.0
            pde_parameters$time <- 0L
            for(i in 1:length(L$tokens)) {
              pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
            }
            
            fe_order <- L$f$FunctionSpace$fe_order
            ## define Rcpp module
            pde_ <- NULL
            if (fe_order == 1) { ## linear finite elements
              pde_ <- new(PDE_2D_ORDER_1, D, pde_type, pde_parameters, L$f)
            }
            if (fe_order == 2) { ## quadratic finite elements
              pde_ <- new(PDE_2D_ORDER_2, D, pde_type, pde_parameters, L$f)
            }
            
            L$f$pde = pde_
            
            ## evaluate forcing term on quadrature nodes
            quad_nodes <- as.matrix(pde_$get_quadrature_nodes())
            pde_$set_forcing(u(quad_nodes,times))
            
            ## initialize solver 
            pde_$init()
            is_init = TRUE
            
            ## set Dirichlet BC
            dirichletBC_ <- dirichletBC(pde_$get_dofs_coordinates(), times)
            pde_$set_dirichlet_bc(dirichletBC_)
            is_dirichletBC_set = TRUE
            
            ## set initial condition
            pde_$set_initial_condition(initialCondition(pde_$get_dofs_coordinates()))
            is_initialCondition_set = TRUE
            
            ## return
            .PdeCtr(times = times,
                    is_dirichletBC_set = is_dirichletBC_set,
                    is_initialCondition_set = is_initialCondition_set,
                    pde_ = pde_,                  # Wraps of C++ class
                    is_parabolic = is_parabolic,
                    is_init=is_init)
          })

setMethod("Pde", signature=c(L="DiffOpObject", u="ANY", dirichletBC="missing", 
                             initialCondition="missing"),
          function(L,u,dirichletBC){
            D = L$f$FunctionSpace$mesh$data ## C++ R_Mesh class
            
            is_parabolic = FALSE
            if( length(L$f$FunctionSpace$mesh$times) != 0 ) is_parabolic = TRUE
            times <- L$f$FunctionSpace$mesh$times
            
            ## set pde type
            pde_type <- 0
            pde_parameters <- NULL
            if ("diffusion" %in% names(L$params) & !is.matrix(L$params$diffusion)){
              ## specialized implementation for laplace operator
              pde_type <- 1L 
              
              pde_parameters$diffusion <- 0.0
              
            } else {
              ## general diffusion-transport-reaction problem, constant coefficients
              pde_type <- 2L
              
              pde_parameters$diffusion <- matrix(0, nrow = 2, ncol = 2)
              
            }
            
            if(is_parabolic){
              if(pde_type == 1L){
                pde_type = 3L
              }else{
                pde_type = 4L
              }
            }
            
            ## prepare pde_parameters list
            pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
            pde_parameters$reaction  <- 0.0
            pde_parameters$time <- 0L
            for(i in 1:length(L$tokens)) {
              pde_parameters[[paste(L$tokens[i])]] <- L$params[[paste(L$tokens[i])]]
            }
            
            fe_order <- L$f$FunctionSpace$fe_order
            ## define Rcpp module
            pde_ <- NULL
            if (fe_order == 1) { ## linear finite elements
              pde_ <- new(PDE_2D_ORDER_1, D, pde_type, pde_parameters, L$f)
            }
            if (fe_order == 2) { ## quadratic finite elements
              pde_ <- new(PDE_2D_ORDER_2, D, pde_type, pde_parameters, L$f)
            }
            
            L$f$pde = pde_
            
            quad_nodes <- as.matrix(pde_$get_quadrature_nodes())
            ## evaluate forcing term on quadrature nodes
            if(!is_parabolic){
              pde_$set_forcing(as.matrix(u(quad_nodes)))
            }else{
              pde_$set_forcing(u(quad_nodes, times))
            }
            ## initialize solver 
            pde_$init()
            is_init = TRUE
            
            ## return
            .PdeCtr(times = times,
                    is_dirichletBC_set = FALSE,
                    is_initialCondition_set = FALSE,
                    pde_ = pde_,                  # Wraps of C++ class
                    is_parabolic = is_parabolic,
                    is_init=is_init)
})

