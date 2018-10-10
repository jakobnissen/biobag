module t

using BioSequences

function composition(x::BioSequences.EachKmerIterator{T,S}) where {T,S}
    if x.step == 1
        return composition_nostep(T, x.seq)
    else
        return composition_step(x)
    end
end

############## Everything below here is internals ##########################

# Array lookups are faster to determine the two-bit value of a nucleotide
# than using trailing_zeros - this shaves off about 15% time.
const kmerof = [0x04, 0x00, 0x01, 0x04,
                0x02, 0x04, 0x04, 0x04,
                0x03, 0x04, 0x04, 0x04,
                0x04, 0x04, 0x04, 0x04]

# We define this function by itself, so that `composition` type-specializes for
# K > 8 for Dict and Array, respectively.
increment!(collection::Array, T, index) = @inbounds collection[index] += 1

function increment!(collection::Dict, T, index)
    key = reinterpret(T, index-1)
    # Defining "value" and using that in the setindex! instead of using += is faster
    get!(collection, key, 0)
    collection[key] = value + 1
    return nothing
end

function convert_to_composition(counts, T)
    if counts isa Dict
        return Composition{T}(counts)
    else
        dictionary = Dict{T,Int}()
        for i in eachindex(counts)
            @inbounds c = counts[i]
            if c > 0
                dictionary[reinterpret(T, i-1)] = c
            end
        end
        return Composition{T}(dictionary)
    end
end

function composition_nostep(::Type{T}, x::BioSequence{A}) where {T<:Kmer, A<:NucAlphs}
    BioSequences.checkkmer(T)
    kmer = typemin(UInt64)

    # When the unfilled bases reaches 0, we have observed K unambiguous nucleotides.
    # and we can increment the Kmercounter. Reset unfilled if we see an
    # ambiguous nucleotide.
    unfilled = BioSequences.kmersize(T) - 1
    mask = UInt64(1) << (2 * BioSequences.kmersize(T)) - UInt64(1)

    if BioSequences.kmersize(T) ≤ 8
        counts = zeros(Int, 4^BioSequences.kmersize(T))
    else
        counts = Dict{T, Int}()
    end

    for i in 1:length(x)
        nuc = BioSequences.inbounds_getindex(x, i)
        @inbounds value = kmerof[reinterpret(UInt8, nuc) + 1]
        if value == 0x04 # ambiguous nucleotide
            unfilled = BioSequences.kmersize(T) - 1
        elseif unfilled == 0
            kmer |= value
            increment!(counts, T, kmer & mask + 1)
            kmer <<= 2
        else
            kmer |= value
            unfilled -= 1
            kmer <<= 2
        end
    end
    return convert_to_composition(counts, T)
end

function comp_kmers(x::BioSequences.OneStepKmerIterator{T}) where {T}
    counts = zeros(Int, 4^BioSequences.kmersize(T))
    for (pos, kmer) in x
        @inbounds counts[reinterpret(Int, kmer) + 1] += 1
    end
    convert_to_composition(counts, T)
end

function composition_nostep2(::Type{T}, x::BioSequence{A}) where {T<:Kmer, A<:NucAlphs}
    BioSequences.checkkmer(T)
    kmer = typemin(UInt64)

    # When the unfilled bases reaches 0, we have observed K unambiguous nucleotides.
    # and we can increment the Kmercounter. Reset unfilled if we see an
    # ambiguous nucleotide.
    unfilled = BioSequences.kmersize(T) - 1
    mask = UInt64(1) << (2 * BioSequences.kmersize(T)) - UInt64(1)

    if BioSequences.kmersize(T) ≤ 8
        counts = zeros(Int, 4^BioSequences.kmersize(T))
    else
        counts = Dict{T, Int}()
    end

    for i in 1:length(x)
        nuc = BioSequences.inbounds_getindex(x, i)
        if isambiguous(nuc)
            unfilled = BioSequences.kmersize(T) - 1
        elseif unfilled == 0
            kmer |= trailing_zeros(nuc)
            increment!(counts, T, kmer & mask + 1)
            kmer <<= 2
        else
            kmer |= trailing_zeros(nuc)
            unfilled -= 1
            kmer <<= 2
        end
    end
    return convert_to_composition(counts, T)
end

function composition_step(x::BioSequences.EachKmerIterator{T}) where {T}
    BioSequences.checkkmer(T)
    kmer = typemin(UInt64)

    # When the unfilled bases reaches 0, we have observed K unambiguous nucleotides.
    # and we can increment the Kmercounter. Reset unfilled if we see an
    # ambiguous nucleotide.
    unfilled = BioSequences.kmersize(T) - 1
    mask = UInt64(1) << (2 * BioSequences.kmersize(T)) - UInt64(1)
    jumpsize = max(1, x.step - BioSequences.kmersize(T) + 1)

    if BioSequences.kmersize(T) ≤ 8
        counts = zeros(Int, 4^BioSequences.kmersize(T))
    else
        counts = Dict{T, Int}()
    end

    frame = 0
    i = x.start
    while i <= length(x.seq)
        nuc = BioSequences.inbounds_getindex(x.seq, i)
        @inbounds value = kmerof[reinterpret(UInt8, nuc) + 1]
        if value == 0x04 # ambiguous nucleotide
            unfilled = BioSequences.kmersize(T) - 1
            continue
        elseif unfilled == 0
            kmer |= value
            if frame == 0
                increment!(counts, T, kmer & mask + 1)
            end
            kmer <<= 2
        else
            unfilled -= jumpsize
            kmer |= value
            kmer <<= 2
        end
        i += jumpsize
        unfilled = jumpsize - 1
        frame += jumpsize
        if frame == x.step
            frame = 0
        end
    end
    return convert_to_composition(counts, T)
end

function extract_kmer_impl(seq, from, k)
    kmer::UInt64 = 0
    isok = true
    for i in 1:k
        nt = inbounds_getindex(seq, from+i-1)
        kmer = kmer << 2 | trailing_zeros(nt)
        isok &= iscertain(nt)
    end
    return kmer, isok
end

####################### Benchmarking ########################
# The new implementation is slower

# The new implementation is faster in all cases, but the difference
# is only appreciable for long sequences, especially at low K and higher stepsizes
# In the best cases, the new implementation is ~6x faster on my laptop.
# In more typical cases, it's 3-4 times faster.

using BenchmarkTools

function benchmark()
    rand_100_4bit = randdnaseq(100)
    rand_100_2bit = BioSequence{DNAAlphabet{2}}(rand_100_4bit)
    rand_10k_4bit = randdnaseq(10000)
    rand_10k_2bit = BioSequence{DNAAlphabet{2}}(rand_10k_4bit)
    rand_1m_4bit = randdnaseq(1000000)
    rand_1m_2bit = BioSequence{DNAAlphabet{2}}(rand_1m_4bit)


    # rand_100_4bit, rand_100_2bit,
    for sequence in (rand_10k_4bit, rand_10k_2bit,
                     rand_1m_4bit, rand_1m_2bit)
        for k in [4, 20] # One which fits in array, one which doesn't
            for step in [1,] # k is worst case, highvalue (e.g. 25) is best case, 1 is usual case
                kmerit = each(Kmer{DNA, k}, sequence, step)
                bits = typeof(sequence).parameters[1].parameters[1]
                println(length(sequence), " nt, ", bits, " bit, k = ", k, ", step = ", step)

                print("Old: ")
                x = BioSequences.composition(kmerit)
                @btime BioSequences.composition($kmerit);

                print("New: ")
                x = composition(kmerit)
                @btime composition($kmerit);
                println("")
            end
        end
    end
end

end # module
