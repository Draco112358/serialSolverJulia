include("FFT_solver_QS_S_type.jl")
include("create_volumes_mapping_v2.jl")
include("create_volume_centers.jl")
include("create_Grids_externals.jl")
include("compute_FFT_mutual_coupling_mats.jl")
include("mesher_FFT.jl")

using JSON, DelimitedFiles, JSON3, MAT, JLD2, ProfileView
using MLUtils: unsqueeze


function encode_complex(z)
    if z isa Complex
        return (z.re, z.im)
    else
        error("Object of type $typeof(z) is not JSON serializable")
    end
end
       
function dump_json_data(matrix_Z,matrix_S,matrix_Y, num_ports)

    z = [[[[0.1, 0.0]]]]
    pop!(z)

    matrix_Z = convert(Array{ComplexF64, 3}, matrix_Z)
    matrix_S = convert(Array{ComplexF64, 3}, matrix_S)
    matrix_Y = convert(Array{ComplexF64, 3}, matrix_Y)

    for i in range(1,num_ports)
        for j in range(1, num_ports)
            elements = map(v -> reinterpret(Float64, [v]), matrix_Z[i,j,:])
            push!(z, [elements])
        end
    end

    s = [[[[0.1, 0.0]]]]
    pop!(s)
    
    for i in range(1,num_ports)
        for j in range(1, num_ports)
            elements = map(v -> reinterpret(Float64, [v]), matrix_S[i,j,:])
            push!(s, [elements])
        end
    end

    y = [[[[0.1, 0.0]]]]
    pop!(y)
    
    for i in range(1,num_ports)
        for j in range(1, num_ports)
            elements = map(v -> reinterpret(Float64, [v]), matrix_Y[i,j,:])
            push!(y, [elements])
        end
    end

    solver_matrices_dict = Dict(
        "matrix_Z" => JSON.json(z),
        "matrix_S" => JSON.json(s),
        "matrix_Y" => JSON.json(y)
    )
    
    return solver_matrices_dict
end


Base.@kwdef struct signal
    dict_element
    value = complex(float(dict_element["Re"]),float(dict_element)["Im"])
end
        
Base.@kwdef struct geom_attributes
    dict_element
    radius = dict_element["radius"]
    segments = dict_element["segments"]
end
        
Base.@kwdef struct transf_params
    dict_element
    position = dict_element["position"]
    rotation = dict_element["rotation"]  
    scale = dict_element["scale"]   
end

Base.@kwdef struct element
    dict_element
    name = dict_element["name"]
    type = dict_element["type"]
    keyComponent = dict_element["keyComponent"]
    geometryAttributes = geom_attributes(dict_element["geometryAttributes"])
    transformationParams = transf_params(dict_element["transformationParams"])
end
        
            
        
Base.@kwdef struct port
    dict_element
    name = dict_element["name"]
    type = dict_element["type"]
    inputElement = element(dict_element["inputElement"])
    outputElement = element(dict_element["outputElement"])
    rlcParams = dict_element["rlcParams"]
    isSelected = dict_element["isSelected"]
end

        
Base.@kwdef struct lumped_element
    dict_element
    name = dict_element["name"]
    type = dict_element["type"]
    value = dict_element["value"]
    inputElement = element(dict_element["inputElement"])
    outputElement = element(dict_element["outputElement"])
    rlcParams = dict_element["rlcParams"]
    isSelected = dict_element["isSelected"]
end

Base.@kwdef mutable struct material
    dict_element
    name = dict_element["name"]
    color = dict_element["color"]
    permeability = dict_element["permeability"]
    tangent_delta_permeability = dict_element["tangent_delta_permeability"]
    custom_permeability = dict_element["custom_permeability"]
    permittivity = dict_element["permittivity"]
    tangent_delta_permittivity = dict_element["tangent_delta_permittivity"]
    custom_permittivity = dict_element["custom_permittivity"]
    conductivity = dict_element["conductivity"]
    tangent_delta_conductivity = dict_element["tangent_delta_conductivity"]
    custom_conductivity = dict_element["custom_conductivity"]
    sigmar = dict_element["conductivity"]
    tan_D = dict_element["tangent_delta_permittivity"]
    eps_re = dict_element["permittivity"]
    mur = dict_element["permeability"]
    epsr = 1
    Rx = nothing
    Ry = nothing
    Rz = nothing
    Cx = nothing
    Cy = nothing
    Cz = nothing
end
        
mutable struct port_def
    port_start
    port_end
    port_voxels
    port_nodes
    surf_s_port_nodes
    surf_e_port_nodes
end

mutable struct le_def
    value
    type
    le_start
    le_end
    le_voxels
    le_nodes
    surf_s_le_nodes
    surf_e_le_nodes
end
                

function read_ports(inputData, escal)
    #@assert inputData isa Dict
    ports = inputData["ports"]
    port_objects = [el for el in ports]
    input_positions = []
    output_positions = []
    N_PORTS = length(port_objects)
    
    for port_object in port_objects
        @assert length(port_object.inputElement.transformationParams.position)==3
        ipos = zeros((1,3))
        ipos[1, 1] = port_object.inputElement.transformationParams.position[1]*escal
        ipos[1, 2] = port_object.inputElement.transformationParams.position[2]*escal
        ipos[1, 3] = port_object.inputElement.transformationParams.position[3]*escal
        push!(input_positions, ipos)
        @assert length(port_object.outputElement.transformationParams.position)==3
        opos = zeros((1, 3))
        opos[1, 1] = port_object.outputElement.transformationParams.position[1]*escal
        opos[1, 2] = port_object.outputElement.transformationParams.position[2]*escal
        opos[1, 3] = port_object.outputElement.transformationParams.position[3]*escal
        push!(output_positions, opos)
    end
    @assert length(input_positions)==N_PORTS && length(output_positions)==N_PORTS
    inp_pos = []
    for i in input_positions
        push!(inp_pos, unsqueeze([i], dims=2))
    end
    out_pos = []
    for i in output_positions
        push!(out_pos, unsqueeze([i], dims=2))
    end
    ports_out = port_def(inp_pos, out_pos,zeros(Int64, (N_PORTS, 2)),zeros(Int64,(N_PORTS, 2)), Array{Any}(undef,0), Array{Any}(undef,0))
    
    return ports_out
end


function read_lumped_elements(inputData, escal)
    
    #@assert inputData isa Dict
    
    lumped_elements = inputData["lumped_elements"]
    
    lumped_elements_objects = [el for el in lumped_elements]
    input_positions = []
    output_positions = []
    values = []
    types = []
    N_LUMPED_ELEMENTS = length(lumped_elements_objects)
    if N_LUMPED_ELEMENTS == 0
        lumped_elements_out = le_def(zeros(0),zeros(Int64, 0),zeros((0, 3)),zeros((0, 3)),zeros(Int64, (0, 2)),zeros(Int64, (0, 2)), Array{Any}(undef,0), Array{Any}(undef,0))
        @assert length(input_positions)==N_LUMPED_ELEMENTS && length(output_positions)==N_LUMPED_ELEMENTS && length(values)==N_LUMPED_ELEMENTS && length(types)==N_LUMPED_ELEMENTS
    else
        for lumped_element_object in lumped_elements_objects
            @assert length(lumped_element_object.inputElement.transformationParams.position)==3
            ipos = zeros((1,3))
            ipos[1, 1] = lumped_element_object.inputElement.transformationParams.position[1]*escal
            ipos[1, 2] = lumped_element_object.inputElement.transformationParams.position[2]*escal
            ipos[1, 3] = lumped_element_object.inputElement.transformationParams.position[3]*escal
            push!(input_positions, ipos)
            @assert length(lumped_element_object.outputElement.transformationParams.position)==3
            opos = zeros((1, 3))
            opos[1, 1] = lumped_element_object.outputElement.transformationParams.position[1]*escal
            opos[1, 2] = lumped_element_object.outputElement.transformationParams.position[2]*escal
            opos[1, 3] = lumped_element_object.outputElement.transformationParams.position[3]*escal
            push!(output_positions, opos)            
            lvalue = zeros(1)
            lvalue[1] = lumped_element_object.value
            append!(values, lvalue)
            
            ltype = zeros(Int64, 1)
            ltype[1] = lumped_element_object.type
            push!(types, ltype)
        end
    end
            
        @assert length(input_positions)==N_LUMPED_ELEMENTS && length(output_positions)==N_LUMPED_ELEMENTS && length(values)==N_LUMPED_ELEMENTS && length(types)==N_LUMPED_ELEMENTS
    
        value = []
        for i in values
            push!(value, i[1])
        end

        type = []
        for i in types
            push!(type, i[1])
        end

        in_pos = []
        for i in input_positions
            push!(in_pos, unsqueeze([i], dims=2))
        end

        out_pos = []
        for i in output_positions
            push!(out_pos, unsqueeze([i], dims=2))
        end

        lumped_elements_out = le_def(value,type,in_pos, out_pos,zeros(Int64, (N_LUMPED_ELEMENTS, 2)), (Int64, (N_LUMPED_ELEMENTS, 2)), Array{Any}(undef,0), Array{Any}(undef,0))

    return lumped_elements_out
end

function read_materials(inputData)
    #@assert inputData isa Dict
    materials = inputData["materials"]
    # for el in materials
    #     el["sigmar"] = 9.4e6
    #     el["tan_D"] = 0
    #     el["eps_re"] = 1
    #     el["mur"] = 1
    # end
    materials_objects = [material(dict_element=el) for el in materials]
    return materials_objects
end

function read_signals(inputData)
    #@assert inputData isa Dict
    signals = inputData["signals"]
    signals_objects = [el for el in signals]
    return signals_objects
end

struct GMRES_set
    Inner_Iter
    Outer_Iter
    tol
end

function getEscalFrom(unit)
    escal = 1
    if (unit == "m")
        escal = 1
    end
    if (unit == "dm")
        escal = 1e-1
    end
    if (unit == "cm")
        escal = 1e-2
    end
    if (unit == "mm")
        escal = 1e-3
    end
    if (unit == "microm")
        escal = 1e-6
    end
    if (unit == "nanom")
        escal = 1e-9
    end
    return escal
end

    
function doSolving(mesherOutput, solverInput, solverAlgoParams)    

    mesherDict = mesherOutput
    inputDict = solverInput
    unit = solverInput["unit"]
    escal = getEscalFrom(unit)
    #JSON3.write("/Users/edgardovittoria/Downloads/mesherOutput.json", mesherDict)
    #JSON3.write("/Users/edgardovittoria/Downloads/solverInput.json", inputDict)
    #@save "/Users/edgardovittoria/Downloads/wptTest.jl" mesherDict inputDict unit

    
    sx, sy, sz = mesherDict["cell_size"]["cell_size_x"]*1000*escal,mesherDict["cell_size"]["cell_size_y"]*1000*escal,mesherDict["cell_size"]["cell_size_z"]*1000*escal
   

    origin_x = mesherDict["origin"]["origin_x"]
    origin_y = mesherDict["origin"]["origin_y"]
    origin_z = mesherDict["origin"]["origin_z"]

    origin = (origin_x,origin_y,origin_z)
    Nx = Int64(mesherDict["n_cells"]["n_cells_x"])
    Ny = Int64(mesherDict["n_cells"]["n_cells_y"])
    Nz = Int64(mesherDict["n_cells"]["n_cells_z"])

    testarray = []
    for (index, value) in mesherDict["mesher_matrices"]
        push!(testarray, copy(value))
    end


    grids = []

    for values in testarray
        if (length(testarray) == 1)
            grids = unsqueeze([values], dims=2)
        else
            push!(grids, unsqueeze(values, dims=2))
        end
    end


    frequencies = inputDict["frequencies"]
    freq = Array{Float64}(undef, 1, length(frequencies))
    for i in range(1, length(frequencies))
        freq[1,i] = frequencies[i]
    end
    #freq = convert(Array{Float64}, freq)

    n_freq = length(freq)
    
    PORTS = read_ports(inputDict, escal)

    L_ELEMENTS = read_lumped_elements(inputDict, escal)

    MATERIALS = read_materials(inputDict) 
    SIGNALS = read_signals(inputDict)
    
    
    
    
    # # START SETTINGS--------------------------------------------
    
    inner_Iter = solverAlgoParams["innerIteration"]
    outer_Iter = solverAlgoParams["outerIteration"]
    tol = solverAlgoParams["convergenceThreshold"]*ones((n_freq))
    # ind_low_freq= filter(i -> !iszero(freq[i]), findall(f -> f<1e5, frequencies))
    # tol[ind_low_freq] .= 1e-7
    

    GMRES_settings = GMRES_set(inner_Iter,outer_Iter,tol)
    QS_Rcc_FW=1; # 1 QS, 2 Rcc, 3 Taylor
    use_escalings=1;

    
   
    mapping_vols,num_centri=create_volumes_mapping_v2(grids)
    
    centri_vox,id_mat=create_volume_centers(grids,mapping_vols,num_centri,sx,sy,sz,origin);


    externals_grids=create_Grids_externals(grids);
    escalings,incidence_selection,circulant_centers,diagonals,expansions,ports,lumped_elements,li_mats,Zs_info=mesher_FFT(use_escalings,MATERIALS,sx,sy,sz,grids,centri_vox,externals_grids,mapping_vols,PORTS,L_ELEMENTS, origin);
    
    # matwrite("/Users/edgardovittoria/Downloads/matfile.mat", Dict(
    #     "grids" =>  grids,
    #     "num_full_vox" => num_centri,
    #     "sx" => sx,
    #     "sy" => sy,
    #     "sz" => sz
    # ))
   
    FFTCP,FFTCLp= @time compute_FFT_mutual_coupling_mats(circulant_centers,escalings,Nx,Ny,Nz,QS_Rcc_FW);

    #inputSolve = matread("/Users/edgardovittoria/Downloads/matfile.mat")
    println("time for solver")
    #Profile.clear()
    out = @time FFT_solver_QS_S_type(freq,escalings,incidence_selection,FFTCP,FFTCLp,diagonals,ports,lumped_elements,expansions,GMRES_settings,Zs_info,QS_Rcc_FW);
    #ProfileView.view()
    return dump_json_data(out["Z"],out["S"],out["Y"], length(inputDict["ports"]))
end


function doSolvingTest(inputDict, mesherDict)    

    unit = "mm"
    escal = getEscalFrom(unit)
    
    
    sx, sy, sz = mesherDict["cell_size"]["cell_size_x"]*1000*escal,mesherDict["cell_size"]["cell_size_y"]*1000*escal,mesherDict["cell_size"]["cell_size_z"]*1000*escal
   

    origin_x = mesherDict["origin"]["origin_x"]
    origin_y = mesherDict["origin"]["origin_y"]
    origin_z = mesherDict["origin"]["origin_z"]

    origin = (origin_x,origin_y,origin_z)
    Nx = Int64(mesherDict["n_cells"]["n_cells_x"])
    Ny = Int64(mesherDict["n_cells"]["n_cells_y"])
    Nz = Int64(mesherDict["n_cells"]["n_cells_z"])

    testarray = []
    for (index, value) in mesherDict["mesher_matrices"]
        push!(testarray, copy(value))
    end


    grids = []

    for values in testarray
        if (length(testarray) == 1)
            grids = unsqueeze([values], dims=2)
        else
            push!(grids, unsqueeze(values, dims=2))
        end
    end


    frequencies = inputDict["frequencies"]
    freq = []
    foreach(x->push!(freq, x), frequencies)
    freq = convert(Array{Float64}, freq)
    
    n_freq = length(freq)
    
    PORTS = read_ports(inputDict, escal)

    L_ELEMENTS = read_lumped_elements(inputDict, escal)

    MATERIALS = read_materials(inputDict) 
    SIGNALS = read_signals(inputDict)
    
    
    
    
    # # START SETTINGS--------------------------------------------
    
    inner_Iter = 100
    outer_Iter = 1
    tol = 0.0001*ones((n_freq))
    # ind_low_freq= filter(i -> !iszero(freq[i]), findall(freq -> freq<1e5, freq))
    # tol[ind_low_freq] .= 1e-7
    

    GMRES_settings = GMRES_set(inner_Iter,outer_Iter,tol)
    QS_Rcc_FW=1; # 1 QS, 2 Rcc, 3 Taylor
    use_escalings=1;

    mapping_vols,num_centri=create_volumes_mapping_v2(grids)
    #centri_vox,id_mat=create_volume_centers(grids,mapping_vols,num_full_vox,sx,sy,sz);
    centri_vox,id_mat=create_volume_centers(grids,mapping_vols,num_centri,sx,sy,sz,origin);


    externals_grids=create_Grids_externals(grids);
    #escalings,incidence_selection,circulant_centers,diagonals,expansions,ports,lumped_elements,li_mats,Zs_info=mesher_FFT(use_escalings,MATERIALS,sx,sy,sz,grids,centri_vox,externals_grids,mapping_vols,PORTS,L_ELEMENTS, origin);
    #ports,lumped_elements=mesher_FFT(use_escalings,MATERIALS,sx,sy,sz,grids,centri_vox,externals_grids,mapping_vols,PORTS,L_ELEMENTS, origin);
    escalings,incidence_selection,circulant_centers,diagonals,expansions,ports,lumped_elements,li_mats,Zs_info=mesher_FFT(use_escalings,MATERIALS,sx,sy,sz,grids,centri_vox,externals_grids,mapping_vols,PORTS,L_ELEMENTS, origin);
    FFTCP,FFTCLp=compute_FFT_mutual_coupling_mats(circulant_centers,escalings,Nx,Ny,Nz,QS_Rcc_FW);

    inputSolve = matread("/Users/edgardovittoria/Downloads/matfile.mat")
    #println("time for solver")
    #port_def(inp_pos, out_pos,zeros(Int64, (N_PORTS, 2)),zeros(Int64,(N_PORTS, 2)), Array{Any}(undef,0), Array{Any}(undef,0))
    #out = @time FFT_solver_QS_S_type(freq,inputSolve["escalings"],inputSolve["incidence_selection"],inputSolve["FFTCP"],inputSolve["FFTCLp"],inputSolve["diagonals"],ports,lumped_elements,inputSolve["expansions"],GMRES_settings,inputSolve["Zs_info"],QS_Rcc_FW);

    println("time for solver")
    out= FFT_solver_QS_S_type(freq,escalings,incidence_selection,FFTCP,FFTCLp,diagonals,ports,lumped_elements,expansions,GMRES_settings,Zs_info,QS_Rcc_FW);
    #println(out)
    ProfileView.view()
    return dump_json_data(out["Z"],out["S"],out["Y"], length(inputDict["ports"]))
end

function test() 
    json_string = read("/Users/edgardovittoria/Downloads/mesherOutput.json", String)
    mesherOutput = JSON3.read(json_string)
    json_string_2 = read("/Users/edgardovittoria/Downloads/solverInput.json", String)
    solverInput = JSON3.read(json_string_2)
    # # json_string_3 = read("/tmp/solverAlgoParams.json", String)
    # # solverAlgoParams = JSON3.read(json_string_3)
    # doSolvingTest(solverInput, mesherOutput)
    
    doSolvingTest(solverInput, mesherOutput)
    
end

#test()
