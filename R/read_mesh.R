# Reads .mesh file (FreeFem++) 
read_mesh <-function(filename){
  
  # Reading file
  file <- readLines(filename)
  
  format.version <- strtoi(strsplit(file[1],split=" ")[[1]][2])
  
  if(format.version == 1){
    
    #2D 
    embedding_dimension <- strtoi(file[4])
  }else if(format.version == 2){
    
    # 2.5D or 3D
    embedding_dimension <- strtoi(strsplit(file[3],split=" ")[[1]][2])
    
  }
  
  # setting local dimension
  local_dimension <- as.integer(ifelse( identical(which(file=="Tetrahedra"),integer(0)), 2, 3))
  
  # Nodes and boundary markers 
  idx.nodes <- which(file == "Vertices")
  num_nodes <- strtoi(file[idx.nodes+1])
  nodes <- matrix(nrow=num_nodes, ncol=embedding_dimension)
  boundary <- matrix(nrow=num_nodes, ncol=1)
  for( i in (idx.nodes+2):(idx.nodes+1+strtoi(file[idx.nodes+1]))){
    line <- strsplit(file[i], " ")[[1]]
    nodes[i-(idx.nodes+1),] <- as.numeric(line)[1:embedding_dimension] # coordinates
    boundary[i-(idx.nodes+1)] <- as.numeric(line)[(embedding_dimension+1)] # boundary markers -> label
  }
  
  # Tringles
  idx.triangles <- which(file == "Triangles")
  num_triangles <- strtoi(file[idx.triangles+1]) 
  triangles <- matrix(nrow=num_triangles, ncol=3)
  for( i in (idx.triangles+2):(idx.triangles+1+strtoi(file[idx.triangles+1]))){
    line <- strsplit(file[i], " ")[[1]]
    triangles[i-(idx.triangles+1),] <- as.integer(line)[1:3]
  }
  
  if(local_dimension == 3 & embedding_dimension ==3){
    # Tetrahedra
    idx.tetrahedra <- which(file == "Tetrahedra")
    num_tetrahedra <- strtoi(file[idx.tetrahedra+1]) 
    tetrahedra <- matrix(nrow=num_tetrahedra, ncol=(local_dimension+1))
    for( i in (idx.tetrahedra+2):(idx.tetrahedra+1+strtoi(file[idx.tetrahedra+1]))){
      line <- strsplit(file[i], " ")[[1]]
      tetrahedra[i-(idx.tetrahedra+1),] <- as.integer(line)[1:(local_dimension+1)]
    }
  }
  
  boundary <- matrix((boundary != 0))
  storage.mode(boundary) <- "integer"
  storage.mode(triangles)<- "integer"
  
  if(local_dimension == 2){
    return(list(nodes=nodes, 
                elements=triangles,
                boundary= boundary))
  }else{
    storage.mode(tetrahedra) <- "integer"
    return(list(nodes=nodes, 
                elements=tetrahedra,
                boundary= boundary))
  }
}
