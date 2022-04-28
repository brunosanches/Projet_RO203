# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX

include("generation.jl")

TOL = 0.00001

"""
Check if the solution gave by CPLEX is connected
Input: maskGrid - Int Array where element i,j is masked if it equals 1
"""
function checkConnectivity(maskGrid::Array{Int64})
    visitedElements = zeros(Int64, size(maskGrid))

    n = size(maskGrid, 1)

    # Get the first element that is not masked
    firstElement = findfirst(x -> x == 0, maskGrid)

    function recursiveVisiting(maskGrid::Array{Int64, 2}, visited::Array{Int64, 2}, pos::CartesianIndex)
        n = size(maskGrid, 1)
        if !checkbounds(Bool, visited, pos)
            return
        elseif visited[pos] == 1
            return
        elseif maskGrid[pos] == 1
            return
        end

        visited[pos] = 1
        
        for (i,j) in [(1, 0), (-1, 0), (0, 1), (0, -1)]
            newPos = pos + CartesianIndex(i, j)
            recursiveVisiting(maskGrid, visited, newPos)
        end
    end

    recursiveVisiting(maskGrid, visitedElements, firstElement)
    # Return true if all elements are visited or masked
    return all(visitedElements .| maskGrid .== 1)
end

"""
Solve an instance with CPLEX
"""
function cplexSolve(grid::Array{Int64,2})

    # Create the model
    m = Model(with_optimizer(CPLEX.Optimizer))

    # TODO
    println("In file resolution.jl, in method cplexSolve(), TODO: fix input and output, define the model")
    n = size(grid,1)
    vX = zeros(n, n)

    # Correspond aux données du problème en forme binaire
    valeur = zeros(Int64, n, n, n)
    for i in 1:n
        for j in 1:n
            valeur[i,j,grid[i,j]] = 1
        end
    end

    # Variables
    @variable(m, x[1:n, 1:n], Bin) # 1 si la case est masqué

    # Constraintes
    # Constraintes de non répetition
    @constraint(m, [i in 1:n, k in 1:n], sum(valeur[i,j,k]*(1-x[i,j]) for j in 1:n) <= 1)
    @constraint(m, [j in 1:n, k in 1:n], sum(valeur[i,j,k]*(1-x[i,j]) for i in 1:n) <= 1)

    #@constraint(m, [i in 2:(n-1), j in 2:(n-1)], x[i,j] + x[i, j+1] + x[i, j-1] + x[i+1,j] + x[i-1,j] <= 1)

    @constraint(m, [i in 2:n, j in 1:n], x[i,j] + x[i-1, j] <= 1)
    @constraint(m, [i in 1:(n-1), j in 1:n], x[i,j] + x[i+1, j] <= 1)
    @constraint(m, [i in 1:n, j in 2:n], x[i,j] + x[i, j-1] <= 1)
    @constraint(m, [i in 1:n, j in 1:(n-1)], x[i,j] + x[i, j+1] <= 1)

    @objective(m, Min, sum(x[i,j] for i in 1:n for j in 1:n))
    # Start a chronometer
    start = time()

    # Solve the model
    optimize!(m)

    if JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT
        vX = Array{Int64}(floor.(JuMP.value.(x)))

        while !checkConnectivity(vX)
            g = findall(x -> x == 1, vX)
            @constraint(m, sum(x[idx] for idx in g) <= length(g)-1)

            optimize!(m)
            if JuMP.primal_status(m) != JuMP.MathOptInterface.FEASIBLE_POINT
                break
            end
            vX = Array{Int64}(floor.(JuMP.value.(x)))
        end

        displaySolution(grid, vX)
    end

    
    # Return:
    # 1 - true if an optimum is found
    # 2 - the resolution time
    # 3 - Masked tiles
    return JuMP.primal_status(m) == JuMP.MathOptInterface.FEASIBLE_POINT, time() - start, vX
    
end

"""
Heuristically solve an instance
"""
function heuristicSolve()

    # TODO
    println("In file resolution.jl, in method heuristicSolve(), TODO: fix input and output, define the model")
    
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
    resolutionMethod = ["cplex"]
    #resolutionMethod = ["cplex", "heuristique"]

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
        readInputFile(dataFolder * file)

        # TODO
        println("In file resolution.jl, in method solveDataSet(), TODO: read value returned by readInputFile()")
        
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
                    
                    # TODO 
                    println("In file resolution.jl, in method solveDataSet(), TODO: fix cplexSolve() arguments and returned values")
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime = cplexSolve()
                    
                    # If a solution is found, write it
                    if isOptimal
                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write cplex solution in fout") 
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Start a chronometer 
                    startingTime = time()
                    
                    # While the grid is not solved and less than 100 seconds are elapsed
                    while !isOptimal && resolutionTime < 100
                        
                        # TODO 
                        println("In file resolution.jl, in method solveDataSet(), TODO: fix heuristicSolve() arguments and returned values")
                        
                        # Solve it and get the results
                        isOptimal, resolutionTime = heuristicSolve()

                        # Stop the chronometer
                        resolutionTime = time() - startingTime
                        
                    end

                    # Write the solution (if any)
                    if isOptimal

                        # TODO
                        println("In file resolution.jl, in method solveDataSet(), TODO: write the heuristic solution in fout")
                        
                    end 
                end

                println(fout, "solveTime = ", resolutionTime) 
                println(fout, "isOptimal = ", isOptimal)
                
                # TODO
                println("In file resolution.jl, in method solveDataSet(), TODO: write the solution in fout") 
                close(fout)
            end


            # Display the results obtained with the method on the current instance
            include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end
