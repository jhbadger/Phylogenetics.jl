## Functions and utilities for reading in trees from newick strings and outputting arrays of either Clado or Phylo ##

# readtree - the primary function for reading in trees in the available formats.
# Read in a set of trees from a file
# Returns an array of trees
function readtree(filepath::ASCIIString, format="nwk")
	instream = open(expanduser(filepath))
	instring = readall(instream)
	close(instream)
	trees = split(instring, ';')
	trees = [replace(i, r"(\r|\n|\s)", "") for i in trees]
	trees = trees[bool([length(t) > 0 for t in trees])]
	output_trees = Array(Phylogeny, length(trees))
	for i in 1:length(trees)
		if search(trees[i], ":") == 0:-1
			output_trees[i] = cladobuild(trees[i])
		elseif search(trees[i], ":") != 0:-1
			output_trees[i] = treebuild(trees[i])
		end
	end
	return output_trees
end


# Used to build Clado structs from newick strings during operation of the
# readtree function.
function cladobuild(tp::ASCIIString)
	# AddInternal - used to add an internal node to the tree.
	function AddInternal(edge, currentNode, node, index, j)
		edge[j, 1] = currentNode 			# Assigns the currentNode to the edge array, at the j'th row and the first column.
		node += 1 							# Increaces the node variable by 1.
		edge[j, 2] = currentNode = node 	# Assigns to the edge array, at the j'th row and 2nd column, the increaced node value. 
		index[node] = j 					# Assigns to the index array at the node'th element, the j variables.
		j += 1 								# Increaces the j variable.
		return currentNode, node, j
	end
	# AddTerminal - used to add a terminal node to the tree.
	function AddTerminal(edge, currentNode, tip, index, tipLabel, tpc, k, j)
		edge[j, 1] = currentNode 	# Make the j'th row and first column of edge array the current node.
		edge[j, 2] = tip 			# Make the j'th row and second column of the edge array the tip variable.
		index[tip] = j 				# Make the tip'th element of the index array j.
		tipLabel[tip] = tpc[k] 		# Make the tip'th element of the tip labels array the kth element of tpc.
		k += 1 						# Increace k, tip, and j.
		tip += 1
		j += 1
		return currentNode, tip, k, j
	end
	# GoDown - used to go down...
	function GoDown(index, currentNode, nodeLabel, nbTip, tpc, k, edge)
		l = index[currentNode] 								# Set variable l to the value of the currentNode'th index value. 
		nodeLabel[currentNode - nbTip] = tpc[k] 			# Set the nodeLabel to tpc[k] to get the node label.
		k += 1 												# Increase k by 1.
		currentNode = edge[l, 1] 							# Set current node to edge[l, 1], basically moving back to the previous node.
		return k, currentNode
	end
	tp = "$tp;"												# Add a semi colon to the tree string fed in.
	if ismatch(r"^[^\(]+\(", tp)							# Detect if the tree read in has a name.
		cutoff = length(match(r"^[^\(]+\(", tp).match)		# Get the character at which the tree starts.
		treeName = chop(match(r"^[^\(]+\(", tp).match)		# Get the tree name.
		treeName = treeName[length(treeName)] == ' ' ? treeName[1:length(treeName)-1] : treeName 	# Make sure the name does not have trailing spaces.
		tp = tp[cutoff:length(tp)]
	else treeName = ""
	end
	if search(tp, ",") == 0:-1 		# If the tree contains no commas, just build the tree it obviously is and skip the more complex build process.
		edge = Array(Int, 2,2)
		edge[1,1:2] = [2,1]
		edge[2,1:2] = [1,2]
		tp = split(tp, r"[\\(\\):;]")
		edgeLength = tp[3]
		Nnode = 1
		tipLabel = tp[2]
		if tp[4] != ""
			nodeLabel = tp[4]
		end
		phyloobject = Phylo(edge, Nnode, tipLabel, edgeLength, nodeLabel, -1.0)
		return phyloobject
	end
	tsp = split(tp, "")		# Split the input string up into an array.
	tp = replace(tp, r"\s", "")		# Make sure to remove any spaces.
	tp = replace(tp, ")", ")NA")	# Make any closing brackets change from ")" to ")NA" so as to handle missing clade names/labels.
	tp = replace(tp, "(", "rem(")	# Make any opening brackets cange from "(" to "rem("
	tpc = split(tp, r"[\\(\\),;]")	# Get the characters from the input tree and put them in an array.
	tpc = tpc[1:length(tpc)-1]		
	tpc = tpc[tpc .!= "rem"]		# Get rid of the characters in the array that say "rem" or "remove".
	skeleton = tsp[bool([i == "(" || i == ")" || i == "," || i == ";" for i in tsp])]
	nsk = length(skeleton)				# Length of the nexus skeleton.
	nbNode = sum(skeleton .== ")")		# Get the number of internal nodes as represented by closing brackets in the skeleton: ")". 
	nbTip = sum(skeleton .== ",") + 1 	# Get the number of tips in the tree, given by the number of commas "," + 1.
	nbEdge = nbNode + nbTip 			# Get the number of edges for tree is number of tips and nodes.
	nodeLabel = ["" for i in 1:nbNode]	# Generate an array to contain node labels and also to contain tip labels.
	tipLabel = ["" for i in 1:nbTip]
	edge = Array(Int, nbEdge,2)			# Make an array called edge to contain the relationships between nodes.
	currentNode = node = nbTip + 1 		# Start variables that track position in the tree, it starts a the number of tips + 1. (Recall ape tree structure, nodes are all numbers higher than the number of tips.)
	edge[nbEdge, 1] = 0					# Fill in the root to the edge matrix and it's relationship to the starting node.
	edge[nbEdge, 2] = node
	index = [0 for i in 1:nbEdge+1]		# Create an array called index, it's job is to track which row of the edge matrix a given node is represented by.
	index[node] = nbEdge 				# Make the node'th element of the index array, the nbEdge variable - the number of edges.
	j = k = tip = 1 					# Make j = k = tip = 1.
	for i in 2:nsk 						# Start a loop for all in 2 to the length of the tree skeleton to loop over the skeleton symbols.
		if skeleton[i] == "("			# If the current skeleton symbol = "(", use the function AddInternal.
			currentNode, node, j = AddInternal(edge, currentNode, node, index, j)
		end
		if skeleton[i] == "," && skeleton[i-1] != ")" 		# If the current skeleton symbol = "," and the previous is not ")", use the function AddTerminal.
			currentNode, tip, k, j = AddTerminal(edge, currentNode, tip, index, tipLabel, tpc, k, j)
		end
		if skeleton[i] == ")"
			if skeleton[i - 1] == ","
				currentNode, tip, k, j = AddTerminal(edge, currentNode, tip, index, tipLabel, tpc, k, j)
				k, currentNode = GoDown(index, currentNode, nodeLabel, nbTip, tpc, k, edge)
			end
			if skeleton[i - 1] == ")"
				k, currentNode = GoDown(index, currentNode, nodeLabel, nbTip, tpc, k, edge)
			end
		end
	end
	edge = edge[1:nbEdge-1, 1:2] 		# Take the final unnessecery row in the edge array off.
        for i in 1:length(nodeLabel)
            nodeLabel[i] = replace(nodeLabel[i],r"^NA","")  	# Make sure any node-labels missing are properly dealt with.
        end
	phyloobject = Clado(treeName, edge, nbNode, tipLabel, nodeLabel)
	return phyloobject
end



# Sub function for creation of a Phylo structs from newick format strings.
# Used to build Phylo structs from newick strings during operation of the
# readtree function.
function treebuild(tp::ASCIIString)
	function AddInternal(edge, j, currentNode, node, index)
		edge[j, 1] = currentNode
		node += 1
        edge[j, 2] = currentNode = node
        index[node] = j
        j += 1
        return node, currentNode, j
	end
	function AddTerminal(edge, j, currentNode, tip, index, tpc, k, tipLabel, edgeLength)
		edge[j, 1] = currentNode
        edge[j, 2] = tip
        index[tip] = j
        X = split(tpc[k], ":")
        tipLabel[tip] = X[1]
        edgeLength[j] = X[2]
        k += 1
        tip +=1
        j += 1
        return k, tip, j
	end
	function GoDown(index, currentNode, tpc, k, nodeLabel, nbTip, edgeLength)
		l = index[currentNode]
        X = split(tpc[k], ":")
        nodeLabel[currentNode - nbTip] = X[1]
        edgeLength[l] = X[2]
        k += 1
        currentNode = edge[l, 1]
        return currentNode, k
	end
	tp = "$tp;"
	if ismatch(r"^[^\(]+\(", tp)
		cutoff = length(match(r"^[^\(]+\(", tp).match)
		treeName = chop(match(r"^[^\(]+\(", tp).match)
		treeName = treeName[length(treeName)] == ' ' ? treeName[1:length(treeName)-1] : treeName
		tp = tp[cutoff:length(tp)]
	else treeName = ""
	end
	if search(tp, ",") == 0:-1
		edge = Array(Int, 2,2)
		edge[1,1:2] = [2,1]
		edge[2,1:2] = [1,2]
		tp = split(tp, r"[\\(\\):;]")
		edgeLength = tp[3]
		Nnode = 1
		tipLabel = tp[2]
		if tp[4] != ""
			nodeLabel = tp[4]
		end
		phyloobject = Phylo(edge, Nnode, tipLabel, edgeLength, nodeLabel, -1.0)
		return phyloobject
	end
	tsp = split(tp, "")
	if ismatch(r"\)[^\(:\)\r\n]*;", tp) # Match the root if there is no colon.
		m = match(r"\)[^\(:\)\r\n]*;", tp)
		mstring = m.match[1:length(m.match)-1]
		newString = "$mstring:NA;"
		st1 = tp[1:m.offset-1]
		tp = "$st1$newString"
	end
	tp = replace(tp, ")", ")NA")
	tp = replace(tp, r"\s", "")
	tp = replace(tp, "(", "rem(")
	tpc = split(tp, r"[\\(\\),;]")
	tpc = tpc[1:length(tpc)-1]
	tpc = tpc[tpc .!= "rem"]
    skeleton = tsp[bool([i == "(" || i == ")" || i == "," || i == ";" for i in tsp])]
    nsk = length(skeleton)
    nbNode = sum(skeleton .== ")")
    nbTip = sum(skeleton .== ",") + 1
    nbEdge = nbNode + nbTip
    nodeLabel = ["" for i in 1:nbNode]
	tipLabel = ["" for i in 1:nbTip]
	edgeLength = ["" for i in 1:nbEdge]
    edge = Array(Int, nbEdge, 2)
    currentNode = node = nbTip + 1
    edge[nbEdge, 1] = 0
    edge[nbEdge, 2] = node
    index = [0 for i in 1:(nbEdge + 1)]
    index[node] = nbEdge
    j = k = tip = 1
    for i in 2:nsk
        if skeleton[i] == "("
            node, currentNode, j = AddInternal(edge, j, currentNode, node, index)
        end
        if skeleton[i] == ","
            if skeleton[i - 1] != ")"
                k, tip, j = AddTerminal(edge, j, currentNode, tip, index, tpc, k, tipLabel, edgeLength)
            end
        end
        if skeleton[i] == ")"
            if skeleton[i - 1] == ","
                k, tip, j = AddTerminal(edge, j, currentNode, tip, index, tpc, k, tipLabel, edgeLength)
                currentNode, k = GoDown(index, currentNode, tpc, k, nodeLabel, nbTip, edgeLength)
            end
            if skeleton[i - 1] == ")"
                currentNode, k = GoDown(index, currentNode, tpc, k, nodeLabel, nbTip, edgeLength)
            end
        end
    end
    edge = edge[1:nbEdge-1, 1:2]
    rootEdge = edgeLength[nbEdge]
    if rootEdge != "NA" # Resolve whether there is a rootedge for the tree.
        rootEdge = float64(rootEdge)
    else
        rootEdge = -1.0
    end
    edgeLength = float64([i == "" ? -1.0 : float64(i) for i in edgeLength[1:nbEdge-1]])
    for i in 1:length(nodeLabel)
        nodeLabel[i] = replace(nodeLabel[i],r"^NA","")
    end
    # nodeLabel = [replace(i, r"^NA", "") for i in nodeLabel]
    phyloobject = Phylo(treeName, edge, nbNode, tipLabel, edgeLength, nodeLabel, rootEdge)
    return phyloobject
end



# Function method for writing a single Phylogeny type to file.
function treewrite(tree::Phylogeny,
                   file::ASCIIString = "output.nwk",
                   append::Bool = false,
                   treeNames::Bool = true)
	output = newick(tree, treeNames)
	if append == true
		outstream = open(file, "a")
	else
		outstream = open(file, "w")
	end
	println(outstream, output)
	close(outstream)
end


# Function method for writing multiple trees to file. In order to do this they must be as an array.
function treewrite(tree::Array{Phylogeny},
                   file::ASCIIString = "output.nwk",
                   append::Bool = false,
                   treeNames::Bool = false)
	outputarray = Array(ASCIIString, length(tree))
	for i in 1:length(tree)
		outputarray[i] = newick(tree[i], treeNames)
	end
	if apppend == true
		outstream = open(file, "a")
	else
		outstream = open(file, "w")
	end
	for i in outputarray
		println(outstream, i)
	end
	close(outstream)
end



# function used by the function which write newick trees to check labels.
function checklabels(labels)
	labels = replace(labels, r"^[[:space:]\\(]+", "")
	labels = replace(labels, r"[[:space:]\\)]+$", "")
	labels = replace(labels, r"[[:space:]]", "_")
	labels = replace(labels, r"[,:;]", "")
	labels = replace(labels, r"[\\(\\)]", "-")
	labels = replace(labels, r"_{2,}", "_")
	labels = replace(labels, r"-{2,}", "-")
	return labels
end



# Function used by the functions which write newick trees.
function cp(x, k, STRING)
	STRING[k] = x
	k += 1
	return k, STRING
end



# Returns the internal node which is the root of the tree.
function getroot(children, parents)
	r = parents[[find(in(children, i)) == [1] ? false : true for i in parents]][1]
	return r
end



# This function returns an array of arrays, showing which children
# nodes have.
function getkids(phy::Phylogeny)
	N = length(phy.tipLabel)
	kids = Array(Array{Int64}, N + phy.Nnode)
	for i in 1:N + phy.Nnode
		logic = bool([n == i for n in phy.edge[1:size(phy.edge,1), 1]])
		kids[i] = phy.edge[logic,2]
	end
	return kids
end


# Construct a newick string from a Cladogram
function newick(phy::Clado, name::Bool)
	function addInternal(i, k, STRING, N, nodelab, tiplab, ind)
		k, STRING = cp("(", k, STRING)
		desc = kids[i]
		for j in desc
			if j > N
				STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab)
			else
				k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			end
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
		k, STRING = cp(")", k, STRING)
		if i > N
			k, STRING = cp(nodelab[i - N], k, STRING)
		end
		return STRING, k
	end
	function addTerminal(i, k, STRING, tiplab, children)
		i = i[1]
		k, STRING = cp(tiplab[children[i]], k, STRING)
		return k, STRING
	end
	if name == true
		prefix = phy.name
	else
		prefix = ""
	end
	nodelab = [i != "" ? checklabels(i) : "" for i in phy.nodeLabel]
	tiplab = [i != "" ? checklabels(i) : i for i in phy.tipLabel]
	children = phy.edge[1:size(phy.edge,1), 2]
	parents = phy.edge[1:size(phy.edge,1), 1]
	N = length(phy.tipLabel)
	kids = getkids(phy)
	LS = (4 * N) + 5
	LS = LS + N # if there are nodelabels.
	ind = [findin(children, i) for i in 1:maximum(phy.edge)]
	STRING = ["" for i in 1:LS]
	k = 1
	k, STRING = cp(prefix, k, STRING)
	k, STRING = cp("(", k, STRING)
	desc = kids[getroot(children, parents)]
	for j in desc
		if j > N
			STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab, ind)
		else
			k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
	end
	k, STRING = cp(")", k, STRING)
	k, STRING = cp(nodelab[1], k, STRING)
	k, STRING = cp(";", k, STRING)
	outstring = ""
	for i in STRING
		outstring = "$outstring$i"
	end
	if name == true && prefix != ""
		namebit = chop(match(r"^[^\(]+\(", outstring).match)
		replace(outstring, namebit, "$namebit ")
	end
	return outstring
end


# Function that creates
function newick(phy::Phylo, name::Bool)
	function addInternal(i, k, STRING, N, nodelab, tiplab, ind)
		k, STRING = cp("(", k, STRING)
		desc = kids[i]
		for j in desc
			if j > N
				STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab)
			else
				k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			end
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
		k, STRING = cp(")", k, STRING)
		if i > N
			k, STRING = cp(nodelab[i - N], k, STRING)
		end
		k, STRING = cp(":", k, STRING)
		edgechar = phy.edgeLength[ind[i]][1]
		k, STRING = cp("$edgechar", k, STRING)
		return STRING, k
	end
	function addTerminal(i, k, STRING, tiplab, children)
		i = i[1]
		k, STRING = cp(tiplab[children[i]], k, STRING)
		k, STRING = cp(":", k, STRING)
		edgechar = phy.edgeLength[i]
		k, STRING = cp("$edgechar", k, STRING)
		return k, STRING
	end
	if name == true
		prefix = phy.name
	else
		prefix = ""
	end
	brl = phy.edgeLength
	nodelab = [i != "" ? checklabels(i) : "" for i in phy.nodeLabel]
	tiplab = [i != "" ? checklabels(i) : i for i in phy.tipLabel]
	children = phy.edge[1:size(phy.edge,1), 2]
	parents = phy.edge[1:size(phy.edge,1), 1]
	N = length(phy.tipLabel)
	kids = getkids(phy)
	LS = (4 * N) + 5
	LS = LS + N
	LS = LS +(4 * N)
	ind = [findin(children, i) for i in 1:maximum(phy.edge)]
	STRING = ["" for i in 1:LS]
	k = 1
	k, STRING = cp(prefix, k, STRING)
	k, STRING = cp("(", k, STRING)
	desc = kids[getroot(children, parents)]
	for j in desc
		if j > N
			STRING, k = addInternal(j, k, STRING, N, nodelab, tiplab, ind)
		else
			k, STRING = addTerminal(ind[j], k, STRING, tiplab, children)
			if j != desc[length(desc)]
				k, STRING = cp(",", k, STRING)
			end
		end
	end
	k, STRING = cp(")", k, STRING)
	k, STRING = cp(nodelab[1], k, STRING)
	k, STRING = cp(";", k, STRING)
	outstring = ""
	for i in STRING
		outstring = "$outstring$i"
	end
	if name == true && prefix != ""
		namebit = chop(match(r"^[^\(]+\(", outstring).match)
		replace(outstring, namebit, "$namebit ")
	end
	return outstring
end


# Macro for making a tree from a string.
macro tr_str(s)
	s = replace(s, r"(\r|\n|\s)", "")
	if search(s, ";") != 0:-1
		s = s[1:length(s)-1]
	end
	if search(s, ":") == 0:-1
		tree = cladobuild(s)
		return tree
	elseif search(s, ":") != 0:-1
		tree = treebuild(s)
		return tree
	end
end


# Functions that creates trees from PhyloXML files:
function recursiveBuild(xmlclade, cladeArray, currentClade, parentClade::Int)
	# Update the node tracker.
	currentClade.nodeIndex += 1
	current = currentClade.nodeIndex # Initialize a local variable called current, taken from the currentClade variable to keep as the variable to pass to furthur recursive calls as the parent index.
	# Get name of clade element.
	name = ""
	# Get and process all additional data.... TODO
	# Process taxonomy...
	taxonomy = Taxonomy(xmlclade)
	sequences = Sequences(xmlclade)
	children = get_elements_by_tagname(xmlclade, "clade")
	# Build the clade element.
	cladeArray[currentClade.nodeIndex] = PhyXClade(name, taxonomy, sequences, parentClade)
	for i in children
		recursiveBuild(i, cladeArray, currentClade, current)
	end
end

function phyXMLbuild(xmltree)
	treestring = replace(string(xmltree), r"(\r|\n|\t)", "")
	treestring = replace(treestring, r"(\s{2,})", "")
	tstable = split(treestring, "><")
	startclade = 0
	endclade = 0
	for i in tstable
		if i == "clade"
			startclade += 1
		end
		if i == "/clade"
			endclade += 1
		end
	end
	startclade != endclade ? println("Warning! There are unequal numbers of clade begins and clade ends in this tree") : println("Current tree has $startclade clade nodes")
	# Ok the number of nodes has been established.

	Clade = Array(TestClade, startclade)		# Make an array to contain the Clade elements.
	BackTrack = zeros(Int, startclade)			# Make an array which tracks the parent of a Clade, to allow backtracking.
	Current = nodeTracker(0)                    # Start the nodetracker type.
	BackTrack[1] = 0
	XML = get_elements_by_tagname(xmltree, "clade")[1]
	recursiveBuild(XML, Clade, Current, 0)
end
