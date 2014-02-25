using LightXML


filepath = "/Users/axolotlfan9250/Desktop/phylodev/phyxml2"
treedoc = parse_file(filepath)

phylogenies = get_elements_by_tagname(root(treedoc), "phylogeny")

xmltree = phylogenies[1]

type nodeTracker
	nodeIndex::Int
end

type TestClade
	name::ASCIIString
	parent::Int
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



function recursiveBuild(xmlclade, cladeArray, currentClade, parentClade::Int)
	# Update the node tracker.
	currentClade.nodeIndex += 1
	current = currentClade.nodeIndex # initialize a loca variable called current, taken from the currentClade variable to keep as the variable to pass to furthur recursive calls as the parent index.
	# Get name of clade element.
	name = ""
	# Get and process all additional data.... TODO
	children = get_elements_by_tagname(xmlclade, "clade")
	# Build the clade element.
	cladeArray[currentClade.nodeIndex] = TestClade(name, parentClade)
	for i in children
		recursiveBuild(i, cladeArray, currentClade, current)
	end
end






