# Abstract type definition for Phylogenetic trees.
abstract Phylogeny

# The type definition for a Phylogenetic Tree with branch lengths.
immutable Phylo <: Phylogeny
	name::String
	edge::Array{Int,2}
	Nnode::Int
	tipLabel::Array{String}
	edgeLength::Array{Float64}
	nodeLabel::Array{String}
	rootEdge::Float64
	Phylo(name, edge,
              Nnode, tipLabel,
              edgeLength, nodeLabel,
              rootEdge) = new(name, edge,
                              Nnode, tipLabel,
                              edgeLength, nodeLabel, rootEdge)
end

# The type definition for a Phylogenetic Tree without branch lengths.
immutable Clado <: Phylogeny
	name::String
	edge::Array{Int,2}
	Nnode::Int
	tipLabel::Array{String}
	nodeLabel::Array{String}
	Clado(name,
              edge,
              Nnode,
              tipLabel,
              nodeLabel) = new(name, edge, Nnode, tipLabel, nodeLabel)
end

# Type definition for a small simple representation of a tree.
immutable ReducedTopology <: Phylogeny
	name::String
	indiesArray::Array{Int}

	function ReducedTopology(phy::Phylogeny)
		children = phy.edge[1:size(phy.edge,1), 2]
		parents = phy.edge[1:size(phy.edge,1), 1]
		IArray = [length(findin(children, i)) != 0 ? parents[findin(children, i)[1]] : -1 for i in 1:max(phy.edge)]
		new(phy.name, IArray)
	end
end

function ReducedTopology(phy::Array{Phylogeny})
	outarray = Array(ReducedTopology, length(phy))
	for i in 1:length(phy)
		outarray[i] = ReducedTopology(phy[i])
	end
	return outarray
end

# Definitions relating to the PhyXClade class, it's construction and structures based on it.
# Nodetracker is nessecery in the building of PhyXClade based trees because an incremental values was needed that could be passed by sharing.
type nodeTracker
	nodeIndex::Int
end

# ID contains ID information, both the code and the provider of the ID. It is nexted in the Taxonomy class.
type ID
	provider::String
	identifier::String
end

# The Taxonomy class contains the Taxonomy information for a clade. Tis includes IDs provided by various providers, codes, and Scientific Names.
type Taxonomy
	IDs::Array{ID}
	Code::String
	ScientificName::String
end

type Accession
	Source::String
	AccessionNumber::String
end

type Sequence
	Symbol::String
	Accession::Accession
	Name::String
	MolecularSequence::String
	Annotations::Array{String}
end

name(x::PhyXElement) = return x.name 
taxonomy(x::PhyXElement) = return x.taxonomy
sequences(x::PhyXElement) = return x.sequences
is_root(x::PhyXElement) = return x.root
is_tip(x::PhyXElement) = return x.tip
parent(x::PhyXElement) = return x.parent

type PhyXTree
	name::String
	phyxNodes::Array{PhyXElement}
end

tipLabels(x::PhyXTree) = map(name, x.phyxNodes[map(is_tip(x.phyxNodes))])
cladeLabels(x::PhyXTree) = map(name, x.phyxNodes[map(!is_tip(x.phyxNodes))])
elementLabels(x::PhyXTree) = map(name, x.phyxNodes)
tips(x::PhyXTree) = return findin(map(is_tip, x.phyxNodes), true)
clades(x::PhyXTree) = return findin(map(is_tip, x.phyxNodes), false)

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


function Sequences(xml::XMLElement)
	seqxml = get_elements_by_tagname(xml, "sequence")
	if !isempty(seqxml)
		outsequences = Array(Sequence, length(seqxml))
		for n in 1:length(seqxml)
			symbolxml = get_elements_by_tagname(seqxml[1], "symbol")
			if !isempty(symbolxml)
				symbol = content(symbolxml[1])
			else
				symbol = ""
			end
			accessionxml = get_elements_by_tagname(seqxml[1], "accession")
			if !isempty(accessionxml)
				accession = Accession(attribute(accessionxml[1], "source", required=false),content(accessionxml[1]))
			else
				accession = Accession("","")
			end
			name = get_elements_by_tagname(seqxml[1], "name")
			if !isempty(name)
				name = content(name[1])
			else
				name = ""
			end
			molseqxml = get_elements_by_tagname(seqxml[1], "mol_seq")
			if !isempty(molseqxml)
				sequence = content(molseqxml[1])
			else
				sequence = ""
			end
			annotationsxml = get_elements_by_tagname(seqxml[1], "annotation")
			if !isempty(annotationsxml)
				annotations = [attribute(i, "ref", required=false) for i in annotationsxml]
			else
				annotations = Array(String, 0)
			end
			outsequences[n] = Sequence(symbol, accession, name, sequence, annotations)
		end
		return outsequences
	else
		return Array(Sequence, 0)
	end
end


# Macro to define the type of PhyXElement at compile time, depending on the information it will have to contain.
macro PhyXPopulate(attributes...)
	code = :(begin end)
	push!(code.args, :($(symbol("Label"))::String))
	push!(code.args, :($(symbol("Root"))::Bool))
	push!(code.args, :($(symbol("Tip"))::Bool))
	push!(code.args, :($(symbol("Parent"))::Int))
	for i in 1:length(attributes)
		push!(code.args, :($(symbol("attr$(attributes[i])"))::$(attributes[i])))
	end
	eval(:(type PhyXElement
		$code
	end))
end

