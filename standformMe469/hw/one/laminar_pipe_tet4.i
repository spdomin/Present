Simulations:
  - name: sim1
    time_integrator: ti_1
    optimizer: opt1
    error_estimator: errest_1

linear_solvers:

  - name: solve_scalar
    type: tpetra
    method: gmres
    preconditioner: sgs
    tolerance: 1e-5
    max_iterations: 50
    kspace: 50
    output_level: 0

  - name: solve_cont
    type: tpetra
    method: gmres 
    preconditioner: muelu
    tolerance: 1e-5
    max_iterations: 100
    kspace: 100
    output_level: 0
    recompute_preconditioner: no
    muelu_xml_file_name: milestone.xml

realms:

  - name: fluidRealm
    mesh: pipe_tet4.g 
    use_edges: no
    automatic_decomposition_type: rcb

    time_step_control:
     target_courant: 2.0
     time_step_change_factor: 1.2

    equation_systems:
      name: theEqSys
      max_iterations: 2
   
      solver_system_specification:
        velocity: solve_scalar
        pressure: solve_cont

      systems:

        - LowMachEOM:
            name: myLowMach
            max_iterations: 1
            convergence_tolerance: 1e-8

    initial_conditions:

      - constant: ic_1
        target_name: block_1
        value:
          pressure: 0.0
          velocity: [0.0,0.0,0.0]
       
    material_properties:

      target_name: block_1

      specifications:
 
        - name: density
          type: constant
          value: tbd

        - name: viscosity
          type: constant
          value: tbd

    boundary_conditions:

    - open_boundary_condition: bc_left
      target_name: surface_2
      open_user_data:
        pressure: tbd
        velocity: [0.0,0.0,0.0]

    - open_boundary_condition: bc_right
      target_name: surface_1
      open_user_data:
        pressure: 0.0
        velocity: [0.0,0.0,0.0]

    - wall_boundary_condition: bc_wall
      target_name: surface_3
      wall_user_data:
        velocity: [0.0,0.0,0.0]

    solution_options:
      name: myOptions
      turbulence_model: laminar
  
      options:
            
        - projected_nodal_gradient:
            velocity: element
            pressure: element

    output:
      output_data_base_name: output/laminar_pipe_tet4.e
      output_frequency: 100
      output_node_set: no 
      output_variables:
       - velocity
       - pressure
       - viscosity
       - density

Time_Integrators:
  - StandardTimeIntegrator:
      name: ti_1
      start_time: 0
      termination_step_count: 200
      time_step: 0.01
      time_stepping_type: adaptive
      time_step_count: 0
      second_order_accuracy: no

      realms:
        - fluidRealm
