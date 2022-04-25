# This file contains methods to generate a data set of instances
include("io.jl")
using Random

"""
Generate an n*n grid with a given density

Argument
- n: size of the grid
"""
function generateInstance(n::Int64)

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

    instance = zeros(Int64, 4, n)
    # Construct each perspective
    # Perspective 1
    for j in 1:n
        maxelem = 0
        for i in 1:n
            if board[i,j] > maxelem
                instance[1, j] += 1
                maxelem = board[i,j]
            end
        end
    end

    # Perspective 2
    for i in 1:n
        maxelem = 0
        for j in n:-1:1
            if board[i,j] > maxelem
                instance[2, i] += 1
                maxelem = board[i,j]
            end
        end
    end

    # Perspective 3
    for j in 1:n
        maxelem = 0
        for i in n:-1:1
            if board[i,j] > maxelem
                instance[3, j] += 1
                maxelem = board[i,j]
            end
        end
    end

    # Perspective 4
    for i in 1:n
        maxelem = 0
        for j in 1:n
            if board[i,j] > maxelem
                instance[4, i] += 1
                maxelem = board[i,j]
            end
        end
    end

    return instance
end 

"""
Generate all the instances

Remark: a grid is generated only if the corresponding output file does not already exist
"""
function generateDataSet()
    quantity_per_size = 5
    sizes = [3, 5, 7, 9, 11]
    for n in sizes
        size_instances = zeros(Int64, quantity_per_size, 4, n)
        for i in 1:quantity_per_size
            filepath = "../data/instance_t" * string(n) * "_i" * string(i) * ".txt"
            if !isfile(filepath)
                instance = generateInstance(n)
                while any([instance == x for x in size_instances])
                    instance = generateInstance(n)
                end

                size_instances[i, :, :] = instance

                file = open(filepath, "w")
                for lin in 1:4
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
