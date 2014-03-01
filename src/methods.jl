function hash{T<:Phylogeny}(p::T)
    h = 0
    for field in names(T)
        h = bitmix(hash(getfield(p,field)),h)
    end
    return h
end

function isequal{T<:Phylogeny}(p1::T,p2::T)
    for field in names(T)
        if !isequal(getfield(p1,field),getfield(p2,field))
            return false
        end
    end
    return true
end