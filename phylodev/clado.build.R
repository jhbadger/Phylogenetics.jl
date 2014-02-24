boopy <- function (tp) 
{
    add.internal <- function() {
        edge[j, 1] <<- current.node
        node <<- node + 1
        edge[j, 2] <<- current.node <<- node
        index[node] <<- j
        j <<- j + 1
    }
    add.terminal <- function() {
        edge[j, 1] <<- current.node
        edge[j, 2] <<- tip
        index[tip] <<- j
        tip.label[tip] <<- tpc[k]
        k <<- k + 1
        tip <<- tip + 1
        j <<- j + 1
    }
    go.down <- function() {
        l <- index[current.node]
        node.label[current.node - nb.tip] <<- tpc[k]
        k <<- k + 1
        current.node <<- edge[l, 1]
    }
    if (!length(grep(",", tp))) {
        obj <- list(edge = matrix(c(2, 1), 1, 2), Nnode = 1)
        tp <- unlist(strsplit(tp, "[\\(\\);]"))
        obj$tip.label <- tp[2]
        if (tp[3] != "") 
            obj$node.label <- tp[3]
        class(obj) <- "phylo"
        return(obj)
    }
    
    # Split the string
    tsp <- unlist(strsplit(tp, NULL))
    tp <- gsub(")", ")NA", tp)
    tp <- gsub(" ", "", tp)
    
    tpc <- unlist(strsplit(tp, "[\\(\\),;]"))
    tpc <- tpc[tpc != ""]
    skeleton <- tsp[tsp == "(" | tsp == ")" | tsp == "," | tsp == 
        ";"]
        

    nsk <- length(skeleton)
    

    nb.node <- length(skeleton[skeleton == ")"])
    nb.tip <- length(skeleton[skeleton == ","]) + 1
    nb.edge <- nb.node + nb.tip

    node.label <- character(nb.node)
    tip.label <- character(nb.tip)
    edge <- matrix(NA, nb.edge, 2)
    

    current.node <- node <- nb.tip + 1
    edge[nb.edge, 1] <- 0
    edge[nb.edge, 2] <- node

    index <- numeric(nb.edge + 1)
    index[node] <- nb.edge
    j <- k <- tip <- 1
    # Up to this point, R and Julia code seems to act the same for this function.


    for (i in 2:nsk) {
        if (skeleton[i] == "(") 
            add.internal()
        if (skeleton[i] == ",") {
            if (skeleton[i - 1] != ")") 
                add.terminal()
        }
        if (skeleton[i] == ")") {
            if (skeleton[i - 1] == ",") {
                add.terminal()
                go.down()
            }
            if (skeleton[i - 1] == ")") 
                go.down()
        }
    }
    
    edge <- edge[-nb.edge, ]
    obj <- list(edge = edge, tip.label = tip.label, Nnode = nb.node, 
        node.label = node.label)
    obj$node.label <- if (all(obj$node.label == "NA")) 
        NULL
    else gsub("^NA", "", obj$node.label)
    class(obj) <- "phylo"
    attr(obj, "order") <- "cladewise"
    obj
}