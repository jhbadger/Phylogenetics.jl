using LightXML


filepath = "/Users/wardb/Desktop/phylodev/phyxml2"
treedoc = parse_file(filepath)

phylogenies = get_elements_by_tagname(root(treedoc), "phylogeny")

xmltree = phylogenies[1]

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
	Clade = []			# Make an array to contain the Clade elements.
	Track = []			# Make an array which tracks the parent of a Clade, to allow backtracking.
	ChildrenDone =[] 	# Make an array to mark the number of children that have been moved over and completed.
	Current = 1
	Current = cladegen(Current, Track)




end

function cladegen(current, track, childrendone, clade)
	if current == 1
		root = true
	else
		root = false
	end

	clade[current] = 

end

type TestClade
	name::ASCIIString
	parent::Int
	children::Array{Int}
end




