using LightXML


filepath = "/Users/axolotlfan9250/Desktop/phylodev/phyxml2"
treedoc = parse_file(filepath)

phylogenies = get_elements_by_tagname(root(treedoc), "phylogeny")

xmltree = phylogenies[1]

type nodeTracker
	nodeIndex::Int
end

type ID
	provider::String
	identifier::String
end

type Taxonomy
	IDs::Array{ID}
	Code::String
	ScientificName::String
end

type TestClade
	name::ASCIIString
	taxonomy::Taxonomy
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
	current = currentClade.nodeIndex # Initialize a local variable called current, taken from the currentClade variable to keep as the variable to pass to furthur recursive calls as the parent index.
	# Get name of clade element.
	name = ""
	# Get and process all additional data.... TODO
	# Process taxonomy...
	taxonomy = Taxonomy(xmlclade)
	children = get_elements_by_tagname(xmlclade, "clade")
	# Build the clade element.
	cladeArray[currentClade.nodeIndex] = TestClade(name, taxonomy, parentClade)
	for i in children
		recursiveBuild(i, cladeArray, currentClade, current)
	end
end


function Taxonomy(xml::XMLElement)
	taxxml = get_elements_by_tagname(xml, "taxonomy")
	if !isempty(taxxml)
		idxml = get_elements_by_tagname(taxxml[1], "id")
		if !isempty(idxml)
			idarray = [ID(attribute(i, "provider"; required=false), content(i)) for i in idxml]
		else
			idarray = Array(ID, 0)
		end
		code = get_elements_by_tagname(taxxml[1], "code")
		if !isempty(code)
			codeval = content(code[1])
		else
			codeval = ""
		end
		sciname = get_elements_by_tagname(taxxml[1], "scientific_name")
		if !isempty(sciname)
			name = content(sciname[1])
		else
			name = ""
		end
		return Taxonomy(idarray, codeval, name)
	else
		return Taxonomy(Array(ID, 0), "", "")
	end
end









	




