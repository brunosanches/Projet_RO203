# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX

include("generation.jl")
include("io.jl")

TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function cplexSolve(limits::Array{Int64,2})

    # Create the model
    m = Model(with_optimizer(CPLEX.Optimizer))

    n = size(limits,2)

    #Defining variables
    @variable(m, x[1:n, 1:n, 1:n], Bin) # Variable 1 si tour i,j = k 0 sinon
    @variable(m, v[1:n, 1:n, 1:4], Bin) # Variable 1 si tour i,j est visible en visant de k
    @variable(m, g_l[1:n, 1:n, 1:n], Bin) # 1 si tour est plus grande qu'autre dans la meme ligne
    @variable(m, g_c[1:n, 1:n, 1:n], Bin) # 1 si tour est plus grande qu'autre dans la meme colonne
    @variable(m, tour[1:n, 1:n] >= 0, Int) # valeur d'hauteur de la tour

    #Defining contraints
    @constraint(m, [i in 1:n, j in 1:n], sum(x[i,j,k] for k in 1:n)==1) # 1 valeur par tour
    @constraint(m, [i in 1:n, k in 1:n], sum(x[i,j,k] for j in 1:n)==1) # 1 valeur par ligne
    @constraint(m, [j in 1:n, k in 1:n], sum(x[i,j,k] for i in 1:n)==1) # 1 valeur par colonne

    ## Constraintes de visualisation (combien de tours je peux voir)
    @constraint(m, [k in [1,3], j in 1:n], sum(v[i,j,k] for i in 1:n)==limits[k,j])
    @constraint(m, [k in [2,4], i in 1:n], sum(v[i,j,k] for j in 1:n)==limits[k,i])

    ## Premiere tour toujours visible
    @constraint(m, [j in 1:n], v[1, j, 1] == 1)
    @constraint(m, [i in 1:n], v[i, n, 2] == 1)
    @constraint(m, [j in 1:n], v[n, j, 3] == 1)
    @constraint(m, [i in 1:n], v[i, 1, 4] == 1)

    ## Si tour a n de hauteur, visible après toutes les perspectives
    @constraint(m, [i in 1:n, j in 1:n, k in 1:4], v[i,j,k] >= x[i,j,n])

    ## Prendre le valeur de la tour
    @constraint(m, [i in 1:n, j in 1:n], sum(k*x[i, j, k] for k in 1:n) == tour[i,j])

    ## Impose les bonnes valeurs pour les variables g_l et g_c
    @constraint(m, [i in 1:n, j in 1:n, j2 in 1:n], tour[i, j2] - tour[i,j] <= (1-g_l[i,j,j2])*n)
    @constraint(m, [i in 1:n, j in 1:n, j2 in 1:n], tour[i, j2] - tour[i,j] >= -g_l[i,j,j2]*n)

    @constraint(m, [i in 1:n, j in 1:n, i2 in 1:n], tour[i2, j] - tour[i,j] <= (1-g_c[i,j,i2])*n)
    @constraint(m, [i in 1:n, j in 1:n, i2 in 1:n], tour[i2, j] - tour[i,j] >= -g_c[i,j,i2]*n)

    # Perspective 1
    @constraint(m, [i in 2:n, j in 1:n], sum(g_c[i, j, i2] for i2 in 1:(i-1)) >= (i-1)*v[i,j,1])
    @constraint(m, [i in 2:n, j in 1:n], sum(g_c[i, j, i2] for i2 in 1:(i-1)) <= i-2 + v[i,j,1])

    # Perspective 2
    @constraint(m, [i in 1:n, j in 1:(n-1)], sum(g_l[i, j, j2] for j2 in (j+1):n) >= (n-j)*v[i,j,2])
    @constraint(m, [i in 1:n, j in 1:(n-1)], sum(g_l[i, j, j2] for j2 in (j+1):n) <= (n-j)-1 + v[i,j,2])

    # Perspective 3
    @constraint(m, [i in 1:(n-1), j in 1:n], sum(g_c[i, j, i2] for i2 in (i+1):n) >= (n-i)*v[i,j,3])
    @constraint(m, [i in 1:(n-1), j in 1:n], sum(g_c[i, j, i2] for i2 in (i+1):n) <= (n-i)-1 + v[i,j,3])

    # Perspective 4
    @constraint(m, [i in 1:n, j in 2:n], sum(g_l[i, j, j2] for j2 in 1:(j-1)) >= (j-1)*v[i,j,4])
    @constraint(m, [i in 1:n, j in 2:n], sum(g_l[i, j, j2] for j2 in 1:(j-1)) <= j-2 + v[i,j,4])

    @objective(m, Max, tour[1,1])
    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(m)


    if JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT
        vtour = Array{Int64}(floor.(JuMP.value.(tour)))
        displaySolution(limits, vtour)
    end
    # Return:
    # 1 - true if an optimum is found
    # 2 - the resolution time
    return JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, vtour
    
end

function insertValue!(possibilities::Array{Int64, 3}, grid::Array{Int64, 2}, 
                                val::Int64, i::Int64, j::Int64)
    n = size(possibilities,1)
    possibilities[i, j, 1:n .!= val] .= 0
    grid[i,j] = val
    possibilities[i, 1:n .!= j, val] .= 0
    possibilities[1:n .!= i, j, val] .= 0
end

function insertValue(possibilities::Array{Int64, 3}, grid::Array{Int64, 2},
                        val::Int64, i::Int64, j::Int64)
    n = size(possibilities, 1)
    newPossibilities = deepcopy(possibilities)
    newPossibilities[i, j, 1:n .!= val] .= 0
    grid[i,j] = val
    newPossibilities[i, 1:n .!= j, val] .= 0
    newPossibilities[1:n .!= i, j, val] .= 0

    return newPossibilities
end

function findUniqueCandidate(possibilities::Array{Int64, 3}, grid::Array{Int64})
    n = size(possibilities, 1)
    # Check if on the line there is one value that is present in only one place
    for i in 1:n
        if any(grid[i,j] == 0 for j in 1:n)
            v = [sum(possibilities[i,j,k] for j in 1:n if grid[i,j] == 0) == 1 for k in 1:n]
            if any(v)
                k = findfirst(x-> x, v)
                j = findfirst(x->possibilities[i,x,k] ==1 && grid[i,x] == 0, 1:n)
                return true, [k], (i,j)
            end
        end
    end
    # Check if on the column there is one value that is present in only one place
    for j in 1:n
        if any(grid[i,j] == 0 for i in 1:n)
            v = [sum(possibilities[i,j,k] for i in 1:n if grid[i,j] == 0) == 1 for k in 1:n]
            if any(v)
                k = findfirst(x-> x, v)
                i = findfirst(x->possibilities[x,j,k] ==1 && grid[x,j] == 0, 1:n)
                return true, [k], (i,j)
            end
        end
    end

    return false, [], (0,0)

end

function findLeastPossibilities(possibilities::Array{Int64, 3}, grid::Array{Int64, 2})
    n = size(possibilities, 1)
    minPossisibilities = n+1
    pos = (0, 0)
    for i in 1:n
        for j in 1:n
            if sum(possibilities[i,j,k] for k in 1:n) < minPossisibilities &&
                grid[i,j] == 0
                minPossisibilities = sum(possibilities[i,j,k] for k in 1:n)
                pos = (i, j)
            end
        end
    end
    return minPossisibilities, pos
end

function findBestPosition(possibilities::Array{Int64, 3}, grid::Array{Int64, 2})
    n = size(possibilities, 1)
    vals = []
    found, vals, position = findUniqueCandidate(possibilities, grid)
    if found
        return vals, position
    end

    nPossibilities, position = findLeastPossibilities(possibilities, grid)
    if position != (0, 0)
        vals = [k*possibilities[position..., k] for k in 1:n if possibilities[position..., k] != 0]
    end

    return vals, position
end

function checkPerspectives(grid::Array{Int64, 2}, limits::Array{Int64, 2}, lin::Int64, col::Int64)
    n = size(grid, 1)

    if col != 0
        # Check first perspective
        countV = 0
        maxH = 0
        for i in 1:n
            if grid[i,col] > maxH
                countV+=1
                maxH = grid[i,col]
            end
        end
        if countV != limits[1, col]
            return false
        end

        # Check second perspective
        countV = 0
        maxH = 0
        for i in n:-1:1
            if grid[i,col] > maxH
                countV+=1
                maxH = grid[i,col]
            end
        end
        if countV != limits[3, col]
            return false
        end
    end

    if lin != 0
        # Check second perspective
        countV = 0
        maxH = 0
        for j in n:-1:1
            if grid[lin,j] > maxH
                countV+=1
                maxH = grid[lin,j]
            end
        end
        if countV != limits[2, lin]
            return false
        end

        # Check fourth perspective
        countV = 0
        maxH = 0
        for j in 1:n
            if grid[lin,j] > maxH
                countV+=1
                maxH = grid[lin,j]
            end
        end
        if countV != limits[4, lin]
            return false
        end
    end
    return true
end


function checkSolution(grid::Array{Int64, 2}, limits::Array{Int64, 2})
    n = size(grid, 1)

    # Check first perspective
    for j in 1:n
        if checkPerspectives(grid, limits, 0, j) == false
            return false
        end
    end

    # Check second perspective
    for i in 1:n
        if checkPerspectives(grid, limits, i, 0) == false
            return false
        end
    end

    return true
end

"""
Heuristically solve an instance
"""
function heuristicSolve(limits::Array{Int64,2})

    n = size(limits,2)

    #Initialize solution and possible values arrays
    possibilities = ones(Int64, n, n, n)
    solution = zeros(Int64, n, n)

    # trivial-cases, limit[1,j] = 1 or == n

    # First perspective
    for j in 1:n
        if limits[1,j] == 1
            insertValue!(possibilities, solution, n, 1, j)
        elseif limits[1,j] == n
            for i in 1:n
                insertValue!(possibilities, solution, i, i, j)
            end
            possibilities[1, j, n] = 0
        else
            possibilities[1, j, n] = 0
        end
    end

    # Second Perspective
    for i in 1:n
        if limits[2, i] == 1
            insertValue!(possibilities, solution, n, i, n)
        elseif limits[2, i] == n
            for j in 1:n
                val = n - j + 1
                insertValue!(possibilities, solution, val, i, j)
            end
            possibilities[i, n, n] = 0
        else
            possibilities[i, n, n] = 0
        end
    end

    # Third perspective
    for j in 1:n
        if limits[3, j] == 1
            insertValue!(possibilities, solution, n, n, j)
        elseif limits[3, j] == n
            for i in 1:n
                val = n - i + 1
                insertValue!(possibilities, solution, val, i, j)
            end
            possibilities[n, j, n] = 0
        else
            possibilities[n, j, n] = 0
        end
    end

    # Fourth perspective
    for i in 1:n
        if limits[4, i] == 1
            insertValue!(possibilities, solution, n, i, 1)
        elseif limits[4, i] == n
            for j in 1:n
                insertValue!(possibilities, solution, j, i, j)
            end
            possibilities[i, 1, n] = 0
        else
            possibilities[i, 1, n] = 0
        end
    end

    function recursiveSolution(possibilities, grid, limits)
        vals, position = findBestPosition(possibilities, grid)
        if position == (0,0)
            return checkSolution(grid, limits)
        elseif length(vals) == 0
            return false
        end

        # If there is only one possibility left
        if length(vals) == 1
            val = vals[1]
            newPossibilities = insertValue(possibilities, grid, val, position...)

            if all(grid[position[1], :] .!= 0) && checkPerspectives(grid, limits, position[1], 0) == false ||
                all(grid[:, position[2]] .!= 0) && checkPerspectives(grid, limits, 0, position[2]) == false
                grid[position...] = 0
                return false
            end

            if recursiveSolution(newPossibilities, grid, limits) == true
                return true
            else
                grid[position...] = 0
                return false
            end
        else
            for val in vals
                newPossibilities = insertValue(possibilities, grid, val, position...)

                if all(grid[position[1], :] .!= 0) && checkPerspectives(grid, limits, position[1], 0) == false ||
                    all(grid[:, position[2]] .!= 0) && checkPerspectives(grid, limits, 0, position[2]) == false
                    grid[position...] = 0
                    continue
                end

                if recursiveSolution(newPossibilities, grid, limits) == true
                    return true
                else
                    grid[position...] = 0
                end
            end
            return false
        end

        return false
    end

    isSolution = recursiveSolution(possibilities, solution, limits)

    return isSolution, solution    
end 

"""
Solve all the instances contained in "../data" through CPLEX and heuristics

The results are written in "../res/cplex" and "../res/heuristic"

Remark: If an instance has previously been solved (either by cplex or the heuristic) it will not be solved again
"""
function solveDataSet()

    dataFolder = "../data/"
    resFolder = "../res/"

    # Array which contains the name of the resolution methods
    #resolutionMethod = ["cplex"]
    resolutionMethod = ["cplex", "heuristique"]

    # Array which contains the result folder of each resolution method
    resolutionFolder = resFolder .* resolutionMethod

    # Create each result folder if it does not exist
    for folder in resolutionFolder
        if !isdir(folder)
            mkdir(folder)
        end
    end
            
    global isOptimal = false
    global solveTime = -1

    # For each instance
    # (for each file in folder dataFolder which ends by ".txt")
    for file in filter(x->occursin(".txt", x), readdir(dataFolder))  
        
        println("-- Resolution of ", file)
        instance = readInputFile(dataFolder * file)
        
        # For each resolution method
        for methodId in 1:size(resolutionMethod, 1)
            
            outputFile = resolutionFolder[methodId] * "/" * file

            # If the instance has not already been solved by this method
            if !isfile(outputFile)
                
                fout = open(outputFile, "w")  

                resolutionTime = -1
                isOptimal = false
                
                # If the method is cplex
                if resolutionMethod[methodId] == "cplex"
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime, solution = cplexSolve(instance)
                    
                    # If a solution is found, write it
                    if isOptimal
                        writeSolution(fout, solution)
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Start a chronometer 
                    startingTime = time()
                    
                    # While the grid is not solved and less than 100 seconds are elapsed
                    while !isOptimal && resolutionTime < 100
                        
                        # Solve it and get the results
                        isOptimal, solution = heuristicSolve(instance)

                        # Stop the chronometer
                        resolutionTime = time() - startingTime
                        
                    end

                    # Write the solution (if any)
                    if isOptimal
                        writeSolution(fout, solution)
                    end 
                end

                println(fout, "solveTime = ", resolutionTime) 
                println(fout, "isOptimal = ", isOptimal)
                
                # TODO Necessary ?
                #println("In file resolution.jl, in method solveDataSet(), TODO: write the solution in fout") 
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end

# A = readInputFile("/home/bruno/Projet_RO203/jeu1/data/instance_t3_i5.txt")
# heuristicSolve(A)