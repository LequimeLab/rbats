#' @import ape
#' @import TreeTools
NULL

#' Shuffler
#'
#' Shuffle the states.
#' @param xml the xml dictionary
#' @param userinput the state attribute that should be shuffled
#' @return the new shuffled xml dictionary
#' @export
shuffle <- function(xml, userinput){
  all_states <- c()
  # Get all the states in a vector
  for(taxon in xml$states){
    all_states <- c(all_states, taxon[[userinput]])
  }
  # Randomize the order
  random_states <- sample(all_states)
  # Put the randomized states back in the xml dictionary
  for(i in 1:length(xml$states)){
    xml$states[[i]][[userinput]] <- random_states[i]
  }
  return(xml)
}


#' Process the xml file
#'
#' Processes the xml file to create the xml dictionary
#' @param inputfile the xml file
#' @return the xml dictionary
#' @export
process_xml = function(input_file) {
  output_list <- list(name="List containing all taxons and their attributes", taxons=list())
  all_attributes = c()
  all_names <- c()
  attval <- list()
  attvalout <- list()
  inattr <- FALSE
  
  if(file.exists(input_file)){
    # Open the file connection
    con = file(input_file, "r")
    on.exit(close(con))
    
    while(TRUE){
      line = readLines(con, n = 1)
      if(grepl("<taxon id=\"", line, fixed=TRUE)){
        # Split the line, slice the needed characters, collapse the split string and remove the whitespaces
        name <- gsub("[[:space:]]", "", paste(strsplit(line, "")[[1]][14:(length(strsplit(line,"")[[1]])-2)], collapse=" "))
        all_names <- c(all_names, name)
      } else if(grepl("<attr name=\"", line, fixed=TRUE)){
        attribute <- toupper(gsub("[[:space:]]", "", paste(strsplit(line, "")[[1]][16:(length(strsplit(line,"")[[1]])-2)], collapse=" ")))
        all_attributes <- c(all_attributes, attribute)
        inattr <- TRUE
      } else if(inattr){
        value <- gsub("[[:space:]]", "", paste(strsplit(line, "")[[1]][5:length(strsplit(line,"")[[1]])], collapse=" "))
        attval[[attribute]] <- value
        inattr <- FALSE
      } else if(grepl("</taxon>", line, fixed=TRUE)){
        attvalout[[name]] <- attval
        attval <- list()
      } else if(grepl("</taxa>", line, fixed=TRUE) | length(line) == 0){
        # Stop reading the file after the taxa part ends
        break
      }
    }
  }
  output_list$names <- all_names
  output_list$attributes <- unique(all_attributes)
  output_list$states <- attvalout
  return(output_list)
}


#' Get the total distance
#'
#' Get the total distance of all the nodes
#' @param node the node
#' @return the total distance
#' @export
get_total_distance <- function(node){
  totaldist <- 0
  if(length(node$daughters) > 0){
    for(daughter in node$daughters){
      # Recursion
      newdist <- get_total_distance(daughter)
      totaldist <- totaldist + newdist
    }
    totaldist <- totaldist + node$distance
  } else {
    # If the node has no daughters, return its distance
    totaldist <- node$distance
  }
  return(totaldist)
}


#' Get the internal distance
#'
#' Get the total distance between all nodes
#' @param node the node
#' @return the internal distance
#' @export
get_internal_distance <- function(node){
  totaldist <- 0
  if(length(node$daughters) > 0){
    for(daughter in node$daughters){
      # Recursion to g
      newdist <- get_internal_distance(daughter)
      totaldist <- totaldist + newdist
    }
    totaldist <- totaldist + node$distance
  } # If the node has no daughters, return 0
  return(totaldist)
}


#' Make the tree
#'
#' Creates the tree by making all the nodes and connecting them
#' @param tree the string of the tree
#' @param pos the position in the string
#' @param statesdict the xml dictionary of the states
#' @param state_attribute the attribute that determines the state
#' @param PS_method the method of how the Parsimony Score will be calculated
#' @return the tree
#' @export
make_tree <- function(tree, statesdict, state_attribute, PS_method) {
  # Split the tree so its a vector of all the characters
  treesplit <- strsplit(tree, "")[[1]]
  thisnode <- list(name="The tree", distance=0, ismono=FALSE, PS=0, daughters=list())
  
  i <- 1
  while(i<=length(treesplit)){
    char <- treesplit[i]
    if(char=="("){
      if(treesplit[i + 1] == "("){
        # If the next character is also a "(", make an internal node
        newnode <- make_internal_node(treesplit, i + 1, statesdict, state_attribute, PS_method)
        thisnode$daughters <- append(thisnode$daughters, list(newnode))
        i <- newnode$position - 1
      } else {
        # If not, make a terminal node of the next character
        newnode <- make_terminal_node(i + 1, treesplit, statesdict, state_attribute)
        thisnode$daughters <- append(thisnode$daughters, list(newnode))
        i <- newnode$position - 1
      }
    } else if(char == ",") {
      if(treesplit[i+1] == "("){
        # If the next character is also a "(", make an internal node
        newnode <- make_internal_node(treesplit, i + 1, statesdict, state_attribute, PS_method)
        thisnode$daughters <- append(thisnode$daughters, list(newnode))
        i <- newnode$position - 1
      } else {
        # If not, make a terminal node of the next character
        newnode <- make_terminal_node(i + 1, treesplit, statesdict, state_attribute)
        thisnode$daughters <- append(thisnode$daughters, list(newnode))
        i <- newnode$position - 1
      }
    } 
    i <- i + 1
  }
  thisnode <- finish_node(thisnode, PS_method)
  thisnode <- association_index(thisnode)
  return(thisnode)
}


#' Make an internal node
#'
#' Creates an internal node by making the intternal nodeds and terminal nodes in it
#' @param tree the split string of the tree
#' @param pos the position in the string
#' @param statesdict the xml dictionary of the states
#' @param state_attribute the attribute that determines the state
#' @param PS_method the method of how the Parsimony Score will be calculated
#' @return the internal node
#' @export
make_internal_node <- function(treesplit, pos, statesdict, state_attribute, PS_method) {
  thisnode <- list(daughters=list(), ismono=FALSE, PS=0)
  closeme <- FALSE
  
  i <- pos + 1
  while(!closeme){
    char <- treesplit[i]
    if(char=="("){
      # Recursion: If the character is a "(" make an internal node in this internal node
      newnode <- make_internal_node(treesplit, i, statesdict, state_attribute, PS_method)
      lastnode <- newnode
      thisnode$daughters <- append(thisnode$daughters, list(newnode))
      i <- newnode$position - 1
    } else if(char==",") {
      if(treesplit[i + 1] == "("){
        # Recursion: If the character is a "(" make an internal node in this internal node
        newnode <- make_internal_node(treesplit, i + 1, statesdict, state_attribute, PS_method)
        lastnode <- newnode
        thisnode$daughters <- append(thisnode$daughters, list(newnode))
        i <- newnode$position - 1
      } else {
        # If the character is a ",", add another terminal node to this internal node
        newpos <- i + 1
        newnode <- make_terminal_node(newpos, treesplit, statesdict, state_attribute)
        thisnode$daughters <- append(thisnode$daughters, list(newnode))
        i <- newnode$position - 1
      }
    } else if(char == ")") {
      # Close this node
      closeme <- TRUE
      thisnode$position <- i
      if(i < length(treesplit) & treesplit[i+1] == ":"){
        # If this internal node has a distance
        i <- i + 1
        distbuffer <- c()
        while(i < length(treesplit)){
          char <- treesplit[i]
          if(char=="," | char == ")"){
            break
          } else {
            distbuffer <- c(distbuffer, char)
          }
          i <- i + 1
        }
        # Formatting
        distance <- paste(distbuffer, collapse = " ")
        distance <- gsub("[[:space:]]", "", distance)
        distance <- gsub(":", "", distance)
        thisnode$distance <- as.double(distance)
        thisnode$position <- i
      }
      thisnode <- finish_node(thisnode, PS_method)
      
      return(thisnode)
    } else {
      # If the character is a number or letter, make an internal node
      newpos <- i
      newnode <- make_terminal_node(newpos, treesplit, statesdict, state_attribute)
      thisnode$daughters <- append(thisnode$daughters, list(newnode))
      i <- newnode$position - 1
    }
    i <- i + 1
    # Insurance that the while loop will stop
    if(i > length(treesplit)){
      break
      closeme <- FALSE
    }
  }
}


#' Make a terminal node
#'
#' Creates a terminal node with the name and distance, and uses the xml dictionary to add the state
#' @param newpos the new position in the string
#' @param treesplit the split string of the tree
#' @param statesdict the xml dictionary of the states
#' @param state_attribute the attribute that determines the state
#' @return the terminal node
#' @export
make_terminal_node <- function(newpos, treesplit, statesdict, state_attribute) {
  thisnode <- list(name = "empty", position=newpos, daughters=list())
  
  isdist <- FALSE
  newchar <- treesplit[newpos+1]
  if(newchar == "," | newchar == ")"){
    # If there is no distance
    thisnode$position <- newpos + 1
  } else {
    buffer <- c()
    i <- newpos
    while(i<=length(treesplit)){
      # Loop through the string
      char <- treesplit[i]
      if(char == ":"){
        # If the character is a ":", the previous characters are the name and the next few characters are the distance
        nodename <- paste(buffer, collapse = " ")
        nodename <- gsub("[[:space:]]", "", nodename)
        nodename <- gsub("\\(", "", nodename)
        thisnode$name <- nodename
        buffer <- c()
        isdist <- TRUE
      } else if(char == ")" | char == ","){
        # If the character is a ")" or a "," the terminal node should be closed
        if(isdist){
          # If there was a distance, add the distance
          distance <- paste(buffer, collapse = " ")
          distance <- gsub("[[:space:]]", "", distance)
          thisnode$distance <- as.double(distance)
          break
        } else {
          # If not, the previous characters were the name. End the internal node
          thisnode$name <- buffer
          break
        }
      } else {
        # These are characters of either the name or the distance, so they are saved to be added later
        if(char == " "){
          char <- "_"
        }
        buffer <- c(buffer, char)
      }
      i <- i + 1
    }
    # Set the position
    thisnode$position <- i
  }
  
  # Set the state of the terminal node using the node name (=taxon name) and the state attribute.
  thisnode$state <- statesdict$states[[thisnode$name]][[state_attribute]]
  return(thisnode)
}


#' Finish internal node
#'
#' Finishes the internal node by setting the state, the parsimony score, whether or not the node is monophyletic, the association index,
#'     and calculating the monophyletic clade if necessary.
#' @param thisnode the internal node that should be finished
#' @param PS_method the method of how the Parsimony Score will be calculated
#' @return the finished internal node
#' @export
finish_node <- function(thisnode, PS_method){
  has_nm <- FALSE
  has_m <- FALSE
  # Check if thisnode has a monophyletic daughter, a non-monophyletic daughter, neiter or both.
  for(daughter in thisnode$daughters){
    if(!is.null(daughter$ismono)){
      if(daughter$ismono == FALSE){
        has_nm = TRUE
      } else if(daughter$ismono == TRUE){
        has_m = TRUE
      }
    }
  }
  d1state <- thisnode$daughters[[1]]$state
  d2state <- thisnode$daughters[[2]]$state
  
  if(has_nm == TRUE){
    # If thisnode has a non-monophyletic daughter
    thisnode$ismono <- FALSE
    
    if(PS_method == "fitch"){
      # Calculate the Parsimony Score using Fitch's algorithm
      if(length(intersect(d1state, d2state)) > 0){
        # If they have one or more matching states
        thisnode$state <- intersect(d1state, d2state)
      } else {
        thisnode$state <- c(d1state, d2state)
        thisnode$PS <- 1
      }
    } else if(PS_method == "legacy"){
      # Calculate the Parsimony Score using Java BaTS' method
      if(all(d1state %in% d2state)){
        thisnode$state <- d1state
      } else {
        thisnode$PS <- 1
        if(length(intersect(d1state, d2state)) == 0){
          # If they have no states in common
          thisnode$state <- c(d1state, d2state)
        } else {
          thisnode$state <- d1state[d1state %in% d2state]
        }
      }
    }
    # If one of the daughters is monophyletic, calculate the monophyletic clade for that node's state
    if(has_m == TRUE){
      for(daughter in thisnode$daughters){
        if(daughter$ismono == TRUE){
          thisnode$monoweights[[daughter$state]] <- count_leaves(daughter)
        }
      }
    }
  } else {
    if(d1state == d2state){
      # If they have the same states
      thisnode$ismono <- TRUE
      thisnode$state <- d1state
    } else {
      # If they have differing states
      thisnode$ismono <- FALSE
      thisnode$state <- c(d1state, d2state)
      thisnode$PS <- 1
      if(has_m == TRUE){
        # If one or both of the daughters is monophyletic, calculate the monophyletic clade for that node's state
        for(daughter in thisnode$daughters){
          thisnode$monoweights[[daughter$state]] <- count_leaves(daughter)
        }
      }
    }
  }
  thisnode <- association_index(thisnode)
  return(thisnode)
}


#' Get streaks
#'
#' Make a list with every monophyletic clade for each state
#' @param node the node of which the monophyletic clades are calculated
#' @param statelist a list with a vector for every state
#' @return the statelist with all the monophyletic clades
#' @export
topmono <- function(node, statelist){
  if(length(node$daughters)>0){
    for(daughter in node$daughters){
      statelist <- topmono(daughter, statelist)
      if(!is.null(daughter$monoweights)){
        for(name in names(daughter$monoweights)){
          statelist[[name]] <- c(statelist[[name]], daughter$monoweights[[name]])
        }
      }
    }
  }
  return(statelist)
}



#' Get MC scores
#'
#' Get the largest monophyletic clade for each state
#' @param treelist a tree object
#' @param possible_states a vector of all possible states
#' @return a list with the MC score for each monophyletic clade
#' @export
highest_mono <- function(treelist, possible_states){
  statelist <- list()
  for(state in possible_states){
    statelist[[state]] <- c(1)
  }
  
  statelist <- topmono(treelist, statelist)
  
  for(name in names(statelist)){
    statelist[[name]] <- max(statelist[[name]])
  }
  return(statelist)
}


#' Get possible states
#'
#' Get tall the possible states for the state attribute
#' @param xmldict the xml dictionary
#' @param userinput the attribute the user marked as the state attribute
#' @return a vector with all possible states for that attribute
#' @export
get_possible_states <- function(xmldict, userinput){
  statelist <- c()
  for(taxon in xmldict$states){
    statelist <- c(statelist, taxon[[userinput]])
  }
  return(unique(statelist))
}


#' Get parsimony score
#'
#' Get the parsimony score following Fitch's algorithm
#' @param node the node that's checked for parsimony
#' @return the parsimony score
#' @export
get_parsimony <- function(node){
  total_parsimony <- 0
  if(length(node$daughters) > 0){
    for(daughter in node$daughters){
      total_parsimony <- total_parsimony + get_parsimony(daughter)
    }
    total_parsimony <- total_parsimony + node$PS
  } # If the node has no daughters, return 0
  return(total_parsimony)
}



#' Get phylogenetic distance
#'
#' Get the phylogenetic distance following Faith's algorithm. This distance is calculated like so:
#'  - 1 x the distance of all terminal nodes
#'  - 1 x the distance of all monophyletic internal nodes
#'  - 2 x the distance of all non-monophyletic internal nodes
#' @param node the node of which the phylogenetic distance is calculated
#' @return the phylogenetic distance
#' @export
get_phylogenetic_distance <- function(node){
  totaldist <- 0
  if(length(node$daughters) > 0){
    for(daughter in node$daughters){
      newdist <- get_phylogenetic_distance(daughter)
      totaldist <- totaldist + newdist
    }
    totaldist <- totaldist + node$distance
    if(!node$ismono){
      totaldist <- totaldist + node$distance
    }
  } else {
    totaldist <- totaldist + node$distance
  }
  return(totaldist)
}


#' Get daughter states
#'
#' Get all the states of all the terminal nodes under this node
#' @param node the node whose terminal daughter states should be obtained
#' @return a vector with all states
#' @export
get_daughter_states <- function(node){
  all_states <- c()
  if(length(node$daughters)>0){
    for(daughter in node$daughters){
      all_states <- c(all_states, get_daughter_states(daughter))
    }
  } else {
    all_states <- c(all_states, node$state)
  }
  return(all_states)
}


#' Get leaf count
#'
#' Count all the leaves under this node
#' @param node the node whose leaves should counted
#' @return the total amount of leaves
#' @export
count_leaves <- function(node){
  total_leaves <- 0
  if(length(node$daughters)>0){
    for(daughter in node$daughters){
      total_leaves <- total_leaves + count_leaves(daughter)
    }
  } else {
    total_leaves <- 1
  }
  return(total_leaves)
}


#' Calculate AI
#'
#' Calculate the association index of one single node
#' @param node the node whose AI should be calculated
#' @return the association index
#' @export
association_index <- function(thisnode){
  daughter_states <- get_daughter_states(thisnode)
  daughter_count <- count_leaves(thisnode)
  frequency <- max(table(daughter_states)) / daughter_count
  
  AI <- (1-frequency)/(2^(daughter_count - 1))
  
  thisnode$AI <- AI
  return(thisnode)
}


#' Calculate AI
#'
#' Calculate the association index of the whole tree
#' @param node the tree whose AI should be calculated
#' @return the association index
#' @export
calculate_AI <- function(node){
  total_ai <- 0
  if(length(node$daughters)>0){
    for(daughter in node$daughters){
      total_ai <- total_ai + calculate_AI(daughter)
    }
    total_ai <- total_ai + node$AI
  }
  
  return(total_ai)
}


#' Calculate UniFrac
#'
#' Calculate the UniFrac of the whole tree
#' @param node the tree whose AI should be calculated
#' @return the association index
#' @export
get_UniFrac <- function(node){
  totaldist <- 0
  if(length(node$daughters) > 0){
    for(daughter in node$daughters){
      newdist <- get_UniFrac(daughter)
      totaldist <- totaldist + newdist
    }
    # Count double if the node is monophyletic
    if(node$ismono == TRUE){
      totaldist <- totaldist + node$distance
    }
  }
  if(!is.null(node$name)){
    if(node$name == "The tree"){
      # For the last loop, devide the total UniFrac by the internal distance of the tree
      totaldist <- totaldist / get_internal_distance(node)
    }
  }
  return(totaldist)
}



#' Calculate NTI and NRI
#'
#' Calculate both the Net Relatedness Index and the Nearest Taxa Index.
#' @param utree the tree in the ape phylo format
#' @param possible_states a vector with all possible states
#' @param userinput the state attribute
#' @param xmldict the xml dictionary
#' @return the association index
#' @export
NTI_NRI <- function(utree, possible_states, userinput, xmldict){
  # Create the distance matrix
  distmatrix <- data.frame(ape::dist.nodes(utree))
  node_names <- utree$tip.label
  
  distmatrix[utree$Nnode+2 : ncol(distmatrix)] <- list(NULL)
  distmatrix <- distmatrix[0:utree$Nnode+1,]
  
  colnames(distmatrix) <- node_names
  rownames(distmatrix) <- node_names
  
  total_NTI <- 0
  total_NRI <- 0
  for(state in possible_states){
    # Split a distance matrix for each state
    names_by_state <- c()
    for(name in utree$tip.label){
      if(xmldict$states[[name]][[userinput]] == state){
        names_by_state <- c(names_by_state, name)
      }
    }
    statematrix <- distmatrix[names_by_state,names_by_state]
    for(col in statematrix){
      total_NTI <- total_NTI + min(col[col>0])
      total_NRI <- total_NRI + sum(col)
    }
  }
  NTI <- total_NTI / length(possible_states)
  NRI <- total_NRI / length(possible_states) / 2
  return(list(NTI=NTI,NRI=NRI))
}



#' Analyse a single tree file
#'
#' Calculate all statistics of a tree file
#' @param treefile the path to the tree file
#' @param xmlfile the path to the .xml file
#' @param reps how many shuffled trees should be made for each normal tree
#' @param userinput the state attribute
#' @param PS_method the method of how the Parsimony Score will be calculated
#' @return a data frame containing the output statistics
#' @export
bats <- function(treefile, xmlfile, reps=1, userinput=NULL, PS_method="legacy"){
  start.time <- Sys.time()
  print("Reading the tree file...")
  apetrees <- c(ape::read.nexus(treefile))
  xmldict <- process_xml(xmlfile)
  
  # Get possible states
  state_attribute <- set_state_attribute(xmldict, userinput)
  possible_states <- get_possible_states(xmldict, state_attribute)
  
  shuffled_xmldicts <- c()
  for(rep in 1:reps){
    shuffled_xmldicts <- c(shuffled_xmldicts, list(shuffle(xmldict, state_attribute)))
  }
  
  print_start(treefile, apetrees, possible_states, state_attribute, reps)
  
  all_normal_stats <- list()
  all_shuffled_stats <- list()
  
  # Shuffled trees
  for(rep in 1:reps){
    shuffled_stats <- list()
    shuffled_xmldict <- shuffled_xmldicts[rep][[1]]
    for(tree in apetrees){
      # Get Trees
      tree <- ape::as.phylo(tree)
      utree <- tree
      tree <- TreeTools::NewickTree(tree)
      
      shuffled_treelist <- make_tree(tree, shuffled_xmldict, state_attribute, PS_method)
      shuffled_stats <- calculate_all_stats(shuffled_treelist, utree, possible_states, state_attribute, shuffled_xmldict, shuffled_stats)
      shuffled_stats <- lapply(shuffled_stats, sum)
    }
    shuffled_stats <- lapply(shuffled_stats, function(x){return(x/length(apetrees))})
    if(length(all_shuffled_stats) == 0){
      all_shuffled_stats <- shuffled_stats
    } else {
      stat_names <- names(shuffled_stats)
      all_shuffled_stats <- lapply(names(shuffled_stats), function(x){append(all_shuffled_stats[[x]], shuffled_stats[[x]])})
      names(all_shuffled_stats) <- stat_names
    }
  }
  
  # Normal trees
  for(tree in apetrees){
    # Get Trees
    tree <- ape::as.phylo(tree)
    utree <- tree
    tree <- TreeTools::NewickTree(tree)
    
    treelist <- make_tree(tree, xmldict, state_attribute, PS_method)
    all_normal_stats <- calculate_all_stats(treelist, utree, possible_states, state_attribute, xmldict, all_normal_stats)
  }
  
  output <- get_output(all_normal_stats, all_shuffled_stats)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(output)
}


#' Function to be added for different formats
#' Read the trees in the correct format
#read_trees <- function(treefile){
#  if(file.exists(treefile)){
#    # Open the file connection
#    con = file(treefile, "r")
#    on.exit(close(con))
#    
#    header = readLines(con, n = 1)
#    
#    if(header == "#NEXUS"){
#      apetrees <- ape::read.nexus(treefile)
#    } else if(header == "#MRBAYES"){
#      apetrees <- treeio::read.mrbayes(treefile)
#    }
#  }
#  return(apetrees)
#}
    
    
#' Set the state attribute
#' 
#' Calculate the UniFrac of the whole tree
#' @param node the tree whose AI should be calculated
#' @return the association index
#' @export
set_state_attribute <- function(xmldict, userinput){
  warn <- FALSE
  
  if(!is.null(userinput)){
    # If the user did specify a state attribute, set it to that attribute
    state_attribute <- toupper(userinput)
  } else if(length(xmldict$attributes) == 1) {
    # if the user did not specify a state attribute and there is only one attribute, set it to that attribute
    state_attribute <- xmldict$attributes
  } else {
    # If the user did not specify a state attribute and there are multiple attributes
    warn <- TRUE
    state_attribute <- smart_pick(xmldict)
  }
  
  statelist <- get_possible_states(xmldict, state_attribute)
  
  if(warn == TRUE){
    warning("There are multiple attributes and none was specified to be the state attribute.\"", state_attribute,
            "\" was automatically chosen to be the state attribute. If this attribute is not the state attribute,
            please specify which attribute it is.\n", "The found attributes are: \n", paste(xmldict$attributes, collapse=", "))
  }
  
  if(length(unique(statelist)) == 1){
    warning("All taxons have the same state for this attribute, this means that this attribute is not fit to be the state attribute")
  } else if(length(unique(statelist)) == length(xmldict$states)){
    warning("This attribute has only unique states, this means that this attribute is not fit to be the state attribute")
  }
  
  return(state_attribute)
}


#' Calculate all the statistics
#'
#' Calculate all the statistics for the given tree and add them to the list of calculated stats
#' @param treelist the tree
#' @param utree the raw tree
#' @param possible_states a vector of all possible states
#' @param state_attribute the state attribute
#' @param xmldict the xml dictionary
#' @param all_stats the list containing all the stats of all previous trees
#' @return the list containing all the stats of all previous trees with the stats of this tree
#' @export
calculate_all_stats <- function(treelist, utree, possible_states, state_attribute, xmldict, all_stats){
  NTI_and_NRI <- NTI_NRI(utree, possible_states, state_attribute, xmldict)
  
  # Call the function for each statistic, and add them to the all_stats list
  all_stats$Total_distance <- c(all_stats$Total_distance, get_total_distance(treelist))
  all_stats$Internal_distance <- c(all_stats$Internal_distance, get_internal_distance(treelist))
  all_stats$Association_index <- c(all_stats$Association_index, calculate_AI(treelist))
  all_stats$Parsimony_score <- c(all_stats$Parsimony_score, get_parsimony(treelist))
  all_stats$UniFrac_score <- c(all_stats$UniFrac_score, get_UniFrac(treelist))
  all_stats$Nearest_taxa_index <- c(all_stats$Nearest_taxa_index, NTI_and_NRI$NTI)
  all_stats$Net_relatedness_index <- c(all_stats$Net_relatedness_index, NTI_and_NRI$NRI)
  all_stats$Phylogenetic_distance <- c(all_stats$Phylogenetic_distance, get_phylogenetic_distance(treelist))
  # The monophylteic clade is calculated for each state
  for(state in possible_states){
    stat_name <- paste0("Monophyletic_clade_", state, collapse="")
    all_stats[[stat_name]] <- c(all_stats[[stat_name]], highest_mono(treelist, possible_states)[[state]])
  }
  return(all_stats)
}


#' Get the output data frame
#'
#' Calculate the means, upper and lower confidence intervals for both the normal and randomized trees. Then calculate the significance.
#' @param all_stats the list containing all the values of all the stats for all the normal trees
#' @param shuffled_stats the list containing all the values of all the stats for all the randomized trees
#' @return Undecided
#' @export
get_output <- function(all_stats, shuffled_stats){
  output_list <- list()
  
  # Sort the stat lists
  all_stats <- lapply(all_stats, sort)
  shuffled_stats <- lapply(shuffled_stats, sort)
  
  medians <- lapply(all_stats, median)
  # Calculate the mean, UCI and LCI and add it to the dataframe
  output_frame <- t(data.frame(lapply(all_stats, mean)))
  output_frame <- cbind(output_frame, t(data.frame(lapply(all_stats, function(x) {return(x[round(length(x)*0.95)])}))))
  output_frame <- cbind(output_frame, t(data.frame(lapply(all_stats, function(x) {return(x[max(c(round(length(x)*0.05), 1))])}))))
  output_frame <- cbind(output_frame, t(data.frame(lapply(shuffled_stats, mean))))
  output_frame <- cbind(output_frame, t(data.frame(lapply(shuffled_stats, function(x) {return(x[round(length(x)*0.95)])}))))
  output_frame <- cbind(output_frame, t(data.frame(lapply(shuffled_stats, function(x) {return(x[max(c(round(length(x)*0.05), 1))])}))))
  # Calculate the significance
  for(stat in names(medians)){
    if("Monophyletic" %in% strsplit(stat, "_")[[1]]){
      count <- length(shuffled_stats[[stat]][shuffled_stats[[stat]]<medians[[stat]]])
    } else {
      count <- length(shuffled_stats[[stat]][shuffled_stats[[stat]]>medians[[stat]]])
    }
    significance <- 1 - (count / length(shuffled_stats[[stat]]))
    output_list$significance[[stat]] <- significance
  }
  
  output_frame <- cbind(output_frame, t(data.frame(output_list$significance)))
  
  colnames(output_frame) <- c("Normal_mean", "Normal_Upper_CI", "Normal_Lower_CI", "Random_mean", "Random_Upper_CI", "Random_Lower_CI", "Significance")
  
  output_frame <- data.frame(output_frame)
  
  return(output_frame)
}


#' Starting message
#'
#' Print a message before starting the analyzing process to inform the user of the input that is used
#' @param treefile the string of the path to the .trees file
#' @param apetrees the raw read trees from the .trees file
#' @param possible_states a vector of all states for the given state attribute
#' @param state_attribute the state attribute
#' @param reps the number indicating the amount of randomized trees that should be made
#' @export
print_start <- function(treefile, apetrees, possible_states, state_attribute, reps){
  if(length(apetrees) == 1){
    filename <- strsplit(treefile, "/")[[1]][length(strsplit(treefile, "/")[[1]])]
    cat("\nAnalysing 1 tree from file ", filename, " with ", reps, " shuffled trees.\n")
    cat("Using state attribute \"", state_attribute, "\" with ", length(possible_states), " states. ", "(", paste(possible_states, collapse=", "), ")\n")
  } else if(length(apetrees) > 1){
    filename <- strsplit(treefile, "/")[[1]][length(strsplit(treefile, "/")[[1]])]
    cat("\nAnalysing", length(apetrees) , "trees from file ", filename, " with ", reps * length(apetrees), " shuffled trees.\n")
    cat("Using state attribute \"", state_attribute, "\" with ", length(possible_states), " states. ", "(", paste(possible_states, collapse=", "), ")\n")
  }
}


#' Pick a good state attribute
#'
#' Pick an attribute that can be used as a state attribute. Returns either a good state attribute or the first attribute.
#' @param xmldict the xml dictionary
#' @return a good or the first attribute
#' @export
smart_pick <- function(xmldict){
  for(attribute in xmldict$attributes){
    states <- c()
    # Get all the states for this attribute
    for(taxon in xmldict$states){
      states <- c(states, taxon[[attribute]])
    }
    # Check if the attribute can be used as a state attribute
    if(!length(unique(states))==1 && !length(unique(states))==length(states) && !min(table(states))==1){
      return(attribute)
    }
  }
  return(xmldict$attributes[[1]])
}









