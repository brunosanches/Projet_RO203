# This file contains methods to generate a data set of instances (i.e., sudoku grids)
include("io.jl")
include("resolution.jl")

using Random

"""
Generate an n*n grid with a given density

Argument
- n: size of the grid
- density: percentage in [0, 1] of repeated values in the grid
"""
function generateInstance(n::Int64, density::Float64)

    # TODO
    println("In file generation.jl, in method generateInstance(), TODO: generate an instance")
    board = zeros(Int64, n, n)

    function recursiveBoardConstructor(board::Array{Int64, 2}, posi::Int, posj::Int)
        n = size(board,1)
        if posi > n || posj > n
            return true
        end

        unavailable_values = union(Set(board[posi,:]), Set(board[:, posj]))
        available_values = setdiff(1:n, unavailable_values)
        if isempty(available_values)
            return false
        end

        shuffle!(available_values)
        for value in available_values
            board[posi, posj] = value
            if posj == n
                next_posj = 1
                next_posi = posi+1
            else
                next_posj = posj+1
                next_posi = posi
            end
            if recursiveBoardConstructor(board, next_posi, next_posj) == true
                return true
            end
        end
        board[posi, posj] = 0
        return false
    end

    recursiveBoardConstructor(board, 1, 1)

    valuesToMask = Int64(round(n*n*density))
    maskedValues = 0
    maskGrid = zeros(Int64, n, n)

    function recursiveMasking(maskGrid::Array{Int64,2}, masked::Int64, 
                              toMask::Int64, positionsAvailable::Set{Tuple{Int64,Int64}})
        # Check if marked enough elements or if there are any other element that can be marked
        if masked == toMask
            return true
        elseif isempty(positionsAvailable)
            return false
        end

        # Take a random available position
        n = size(maskGrid, 1)
        pos = rand(positionsAvailable)
        
        # Try masking it and then calling recursion if it does not violate connexitivity
        maskGrid[pos...] = 1
        delete!(positionsAvailable, pos)
        if checkConnectivity(maskGrid) == 1
            # Check for values that are in the set, take in count the values in border
            existingValues = intersect(positionsAvailable, [pos .+ (0, 1), pos .+ (0, -1), pos .+ (1, 0), pos .+ (-1, 0)])
            setdiff!(positionsAvailable, existingValues)
            if recursiveMasking(maskGrid, masked+1, toMask, positionsAvailable) == true
                return true
            end
            union!(positionsAvailable, existingValues)
        end
        maskGrid[pos...] = 0
        
        if recursiveMasking(maskGrid, masked, toMask, positionsAvailable) == true
            return true
        end
        push!(positionsAvailable, pos)
        return false
    end

    positionsAvailable = Set(collect(Iterators.product(ntuple(i -> 1:n, 2)...)))

    foundInstance = recursiveMasking(maskGrid, maskedValues, valuesToMask, positionsAvailable)

    if foundInstance == true
        maskedPositions = Set(findall(x -> x == 1, maskGrid))
        for pos in maskedPositions
            valPositions = setdiff(union(collect(Iterators.product(pos[1], 1:n)),
                                         collect(Iterators.product(1:n, pos[2]))), Tuple(pos))
            val = board[rand(valPositions)...]
            board[pos] = val
        end
    end
    
    return foundInstance, board
end 

"""
Generate all the instances

Remark: a grid is generated only if the corresponding output file does not already exist
"""
function generateDataSet()

    # TODO
    println("In file generation.jl, in method generateDataSet(), TODO: generate an instance")
    sizes = [4, 6, 8, 10, 12, 14]
    densities = [0.05, 0.1, 0.2, 0.3]
    for n in sizes
        size_instances = zeros(Int64, length(densities), n, n)
        for i in 1:length(densities)
            filepath = "../data/instance_t" * string(n) * "_d" * string(densities[i]) * ".txt"
            if !isfile(filepath)
                solveable, instance = generateInstance(n, densities[i])
                while any([instance == x for x in size_instances])
                    solveable, instance = generateInstance(n, densities[i])
                    if solveable == false
                        break
                    end
                end
                if solveable == false
                    println("Size: " * string(n) * " density: " * string(densities[i]) * " could not be generated")
                    continue
                end
                size_instances[i, :, :] = instance

                file = open(filepath, "w")
                for lin in 1:n
                    for col in 1:(n-1)
                        print(file, string(instance[lin,col]) * ",")
                    end
                    println(file, instance[lin,n])
                end
                close(file)
            end

        end
    end
    
end



