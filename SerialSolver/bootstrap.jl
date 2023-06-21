(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

using SerialSolver
const UserApp = SerialSolver
SerialSolver.main()
