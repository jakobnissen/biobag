module t

using BioSequences

composition(x::BioSequences.EachKmerIterator{T,S}) where {T,S} = composition_internal(x, Val(x.step == 1))

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
    value = get!(collection, key, 0)
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

# Onestep value dispatch means the operations which keep track of the frame are
# compiled away if they're not needed. This takes about 20% off the time.
function composition_internal(x::BioSequences.EachKmerIterator{T, S}, onestep::T2) where {T, S, T2<:Union{Val{true},Val{false}}}
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

    frame = 1
    for i in x.start:length(x.seq)
        if typeof(onestep) == Val{false}
            if frame == x.step
                frame = 1
            else
                frame += 1
            end
        end

        kmer <<= 2
        @inbounds nuc = x.seq[i]
        @inbounds value = kmerof[reinterpret(UInt8, nuc) + 1]
        if value == 0x04 # ambiguous nucleotide
            unfilled = BioSequences.kmersize(T) - 1
            continue
        end
        if unfilled == 0
            if typeof(onestep) == Val{true}
                kmer |= value
                increment!(counts, T, kmer & mask + 1)
            else
                if frame == 1
                    kmer |= value
                    increment!(counts, T, kmer & mask + 1)
                end
            end
        else
            unfilled -= 1
        end
    end
    return convert_to_composition(counts, T)
end

####################### Benchmarking ########################
# The new implementation is faster in all cases, but the difference
# is only appreciable for long sequences, especially at low K.
# In the best cases, the new implementation is ~4x faster on my laptop.

using BenchmarkTools

function benchmark()
    rand_100_4bit = randdnaseq(100)
    rand_100_2bit = BioSequence{DNAAlphabet{2}}(rand_100_4bit)
    rand_10k_4bit = randdnaseq(10000)
    rand_10k_2bit = BioSequence{DNAAlphabet{2}}(rand_10k_4bit)
    rand_1m_4bit = randdnaseq(1000000)
    rand_1m_2bit = BioSequence{DNAAlphabet{2}}(rand_1m_4bit)

    for sequence in (rand_100_4bit, rand_100_2bit,
                     rand_10k_4bit, rand_10k_2bit,
                     rand_1m_4bit, rand_1m_2bit)
        for k in [4, 20]
            bits = typeof(sequence).parameters[1].parameters[1]
            println(length(sequence), " nt, ", bits, " bit, k = ", k)

            print("Old: ")
            x = BioSequences.composition(each(Kmer{DNA,k}, sequence))
            @btime BioSequences.composition(each(Kmer{DNA,$k}, $sequence));

            print("New: ")
            x = composition_new(each(Kmer{DNA,k}, sequence))
            @btime composition_new(each(Kmer{DNA,$k}, $sequence));
            println("")
        end
    end
end

end # module
