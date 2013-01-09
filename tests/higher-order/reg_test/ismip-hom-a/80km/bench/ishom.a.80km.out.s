 Layout(EW,NS) =            41           41  total procs =             1
  
 Compute higher-order ice velocities, time =    0.000000000000000     
  
 Running Payne/Price higher-order dynamics with JFNK solver
NOXINIT CALLED    for nelem=43706
NOXInit: param list is: (delete this debug line)
Jacobian Operator = Matrix-Free   [unused]
Matrix-Free Perturbation = 0.0001   [unused]
CISM: Number of Time Steps To Use LOCA = 0   [unused]
Lean Matrix Free = 1   [unused]
NOX -> 
 Nonlinear Solver = Line Search Based   [unused]
 Direction -> 
  Method = Newton   [unused]
  Newton -> 
   Forcing Term Method = Type 2   [unused]
   Rescue Bad Newton Solve = 1   [unused]
   Stratimikos Linear Solver -> 
    NOX Stratimikos Options -> 
     [empty list]
    Stratimikos -> 
     Linear Solver Type = Belos   [unused]
     Preconditioner Type = None   [unused]
     Linear Solver Types -> 
      Belos -> 
       Solver Type = Block GMRES   [unused]
       Solver Types -> 
        Block GMRES -> 
         Convergence Tolerance = 0.0001   [unused]
         Output Frequency = 1   [unused]
         Output Style = 1   [unused]
         Verbosity = 33   [unused]
         Maximum Iterations = 100   [unused]
         Block Size = 1   [unused]
         Num Blocks = 100   [unused]
         Flexible Gmres = 1   [unused]
 Line Search -> 
  [empty list]
 Printing -> 
  Output Information -> 
   Error = 1   [unused]
   Warning = 1   [unused]
   Outer Iteration = 1   [unused]
   Parameters = 0   [unused]
   Details = 0   [unused]
   Stepper Iteration = 1   [unused]
   Linear Solver Details = 0   [unused]
 Status Tests -> 
  Test Type = Combo   [unused]
  Combo Type = OR   [unused]
  Number of Tests = 2   [unused]
  Test 0 -> 
   Test Type = NormF   [unused]
   Tolerance = 0.0001   [unused]
  Test 1 -> 
   Test Type = MaxIters   [unused]
   Maximum Iterations = 100   [unused]

NOXSolve called

************************************************************************
-- Nonlinear Solver Step 0 -- 
||F|| = 6.789e+02  step = 0.000e+00  dx = 0.000e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 26 iterations with total CPU time of 10.0719 sec

************************************************************************
-- Nonlinear Solver Step 1 -- 
||F|| = 8.283e+04  step = 1.000e+00  dx = 6.514e-03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88974 sec

************************************************************************
-- Nonlinear Solver Step 2 -- 
||F|| = 6.203e+04  step = 1.000e+00  dx = 2.065e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 16 iterations with total CPU time of 6.1285 sec

************************************************************************
-- Nonlinear Solver Step 3 -- 
||F|| = 4.808e+04  step = 1.000e+00  dx = 1.762e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.51051 sec

************************************************************************
-- Nonlinear Solver Step 4 -- 
||F|| = 3.737e+05  step = 1.000e+00  dx = 8.910e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.89575 sec

************************************************************************
-- Nonlinear Solver Step 5 -- 
||F|| = 2.172e+05  step = 1.000e+00  dx = 3.966e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88236 sec

************************************************************************
-- Nonlinear Solver Step 6 -- 
||F|| = 1.541e+05  step = 1.000e+00  dx = 2.925e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 11 iterations with total CPU time of 4.27423 sec

************************************************************************
-- Nonlinear Solver Step 7 -- 
||F|| = 4.695e+05  step = 1.000e+00  dx = 6.481e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.89308 sec

************************************************************************
-- Nonlinear Solver Step 8 -- 
||F|| = 3.670e+05  step = 1.000e+00  dx = 6.738e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 10 iterations with total CPU time of 3.79303 sec

************************************************************************
-- Nonlinear Solver Step 9 -- 
||F|| = 2.415e+05  step = 1.000e+00  dx = 2.471e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 10 iterations with total CPU time of 3.901 sec

************************************************************************
-- Nonlinear Solver Step 10 -- 
||F|| = 1.183e+06  step = 1.000e+00  dx = 9.246e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.99147 sec

************************************************************************
-- Nonlinear Solver Step 11 -- 
||F|| = 4.387e+06  step = 1.000e+00  dx = 2.813e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.756009 sec

************************************************************************
-- Nonlinear Solver Step 12 -- 
||F|| = 3.815e+06  step = 1.000e+00  dx = 7.146e-03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88951 sec

************************************************************************
-- Nonlinear Solver Step 13 -- 
||F|| = 1.356e+06  step = 1.000e+00  dx = 2.128e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 13 iterations with total CPU time of 5.02159 sec

************************************************************************
-- Nonlinear Solver Step 14 -- 
||F|| = 2.364e+06  step = 1.000e+00  dx = 4.668e-01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.13476 sec

************************************************************************
-- Nonlinear Solver Step 15 -- 
||F|| = 2.090e+06  step = 1.000e+00  dx = 1.837e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88114 sec

************************************************************************
-- Nonlinear Solver Step 16 -- 
||F|| = 1.427e+06  step = 1.000e+00  dx = 3.578e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 10 iterations with total CPU time of 3.81059 sec

************************************************************************
-- Nonlinear Solver Step 17 -- 
||F|| = 7.671e+05  step = 1.000e+00  dx = 1.250e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 13 iterations with total CPU time of 5.02741 sec

************************************************************************
-- Nonlinear Solver Step 18 -- 
||F|| = 1.085e+06  step = 1.000e+00  dx = 2.659e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 12 iterations with total CPU time of 4.57427 sec

************************************************************************
-- Nonlinear Solver Step 19 -- 
||F|| = 9.225e+05  step = 1.000e+00  dx = 1.000e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.989 sec

************************************************************************
-- Nonlinear Solver Step 20 -- 
||F|| = 8.465e+06  step = 1.000e+00  dx = 3.645e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88479 sec

************************************************************************
-- Nonlinear Solver Step 21 -- 
||F|| = 6.354e+05  step = 1.000e+00  dx = 3.227e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.73875 sec

************************************************************************
-- Nonlinear Solver Step 22 -- 
||F|| = 2.751e+07  step = 1.000e+00  dx = 5.538e-01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.51799 sec

************************************************************************
-- Nonlinear Solver Step 23 -- 
||F|| = 6.668e+06  step = 1.000e+00  dx = 2.097e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 12 iterations with total CPU time of 4.59712 sec

************************************************************************
-- Nonlinear Solver Step 24 -- 
||F|| = 5.282e+06  step = 1.000e+00  dx = 4.622e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.757239 sec

************************************************************************
-- Nonlinear Solver Step 25 -- 
||F|| = 2.813e+06  step = 1.000e+00  dx = 1.537e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 12 iterations with total CPU time of 4.64981 sec

************************************************************************
-- Nonlinear Solver Step 26 -- 
||F|| = 3.579e+06  step = 1.000e+00  dx = 6.356e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.37813 sec

************************************************************************
-- Nonlinear Solver Step 27 -- 
||F|| = 3.175e+06  step = 1.000e+00  dx = 2.210e-03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.753855 sec

************************************************************************
-- Nonlinear Solver Step 28 -- 
||F|| = 2.385e+06  step = 1.000e+00  dx = 1.719e-02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 10 iterations with total CPU time of 3.91458 sec

************************************************************************
-- Nonlinear Solver Step 29 -- 
||F|| = 3.204e+07  step = 1.000e+00  dx = 2.441e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.52784 sec

************************************************************************
-- Nonlinear Solver Step 30 -- 
||F|| = 1.765e+08  step = 1.000e+00  dx = 6.832e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.753525 sec

************************************************************************
-- Nonlinear Solver Step 31 -- 
||F|| = 1.208e+08  step = 1.000e+00  dx = 2.503e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.75109 sec

************************************************************************
-- Nonlinear Solver Step 32 -- 
||F|| = 1.057e+09  step = 1.000e+00  dx = 4.503e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.75417 sec

************************************************************************
-- Nonlinear Solver Step 33 -- 
||F|| = 1.625e+08  step = 1.000e+00  dx = 8.137e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.30569 sec

************************************************************************
-- Nonlinear Solver Step 34 -- 
||F|| = 1.189e+08  step = 1.000e+00  dx = 1.958e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.52115 sec

************************************************************************
-- Nonlinear Solver Step 35 -- 
||F|| = 2.748e+08  step = 1.000e+00  dx = 4.895e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.376566 sec

************************************************************************
-- Nonlinear Solver Step 36 -- 
||F|| = 1.403e+08  step = 1.000e+00  dx = 3.411e-01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.75208 sec

************************************************************************
-- Nonlinear Solver Step 37 -- 
||F|| = 9.106e+07  step = 1.000e+00  dx = 4.058e-01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.7321 sec

************************************************************************
-- Nonlinear Solver Step 38 -- 
||F|| = 3.277e+08  step = 1.000e+00  dx = 4.721e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.13602 sec

************************************************************************
-- Nonlinear Solver Step 39 -- 
||F|| = 1.534e+08  step = 1.000e+00  dx = 4.407e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.12743 sec

************************************************************************
-- Nonlinear Solver Step 40 -- 
||F|| = 9.584e+07  step = 1.000e+00  dx = 2.762e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 3.12252 sec

************************************************************************
-- Nonlinear Solver Step 41 -- 
||F|| = 4.276e+08  step = 1.000e+00  dx = 3.136e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.381912 sec

************************************************************************
-- Nonlinear Solver Step 42 -- 
||F|| = 2.996e+08  step = 1.000e+00  dx = 8.383e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.755549 sec

************************************************************************
-- Nonlinear Solver Step 43 -- 
||F|| = 1.516e+08  step = 1.000e+00  dx = 7.216e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88712 sec

************************************************************************
-- Nonlinear Solver Step 44 -- 
||F|| = 7.563e+07  step = 1.000e+00  dx = 2.464e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.76027 sec

************************************************************************
-- Nonlinear Solver Step 45 -- 
||F|| = 1.110e+09  step = 1.000e+00  dx = 5.409e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.38171 sec

************************************************************************
-- Nonlinear Solver Step 46 -- 
||F|| = 9.411e+08  step = 1.000e+00  dx = 1.061e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.752704 sec

************************************************************************
-- Nonlinear Solver Step 47 -- 
||F|| = 3.122e+08  step = 1.000e+00  dx = 2.544e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 3.12778 sec

************************************************************************
-- Nonlinear Solver Step 48 -- 
||F|| = 2.072e+09  step = 1.000e+00  dx = 5.711e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.377864 sec

************************************************************************
-- Nonlinear Solver Step 49 -- 
||F|| = 1.809e+09  step = 1.000e+00  dx = 1.488e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.759864 sec

************************************************************************
-- Nonlinear Solver Step 50 -- 
||F|| = 5.140e+08  step = 1.000e+00  dx = 5.436e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.35313 sec

************************************************************************
-- Nonlinear Solver Step 51 -- 
||F|| = 3.354e+09  step = 1.000e+00  dx = 8.571e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.50335 sec

************************************************************************
-- Nonlinear Solver Step 52 -- 
||F|| = 8.734e+08  step = 1.000e+00  dx = 3.321e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88737 sec

************************************************************************
-- Nonlinear Solver Step 53 -- 
||F|| = 4.884e+08  step = 1.000e+00  dx = 6.271e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.59543 sec

************************************************************************
-- Nonlinear Solver Step 54 -- 
||F|| = 1.565e+10  step = 1.000e+00  dx = 1.207e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.755862 sec

************************************************************************
-- Nonlinear Solver Step 55 -- 
||F|| = 1.224e+10  step = 1.000e+00  dx = 4.860e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.5035 sec

************************************************************************
-- Nonlinear Solver Step 56 -- 
||F|| = 1.080e+10  step = 1.000e+00  dx = 2.250e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.26802 sec

************************************************************************
-- Nonlinear Solver Step 57 -- 
||F|| = 5.339e+09  step = 1.000e+00  dx = 2.485e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.41712 sec

************************************************************************
-- Nonlinear Solver Step 58 -- 
||F|| = 3.604e+09  step = 1.000e+00  dx = 1.374e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.28379 sec

************************************************************************
-- Nonlinear Solver Step 59 -- 
||F|| = 4.172e+09  step = 1.000e+00  dx = 5.009e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.382611 sec

************************************************************************
-- Nonlinear Solver Step 60 -- 
||F|| = 9.763e+08  step = 1.000e+00  dx = 4.487e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.41097 sec

************************************************************************
-- Nonlinear Solver Step 61 -- 
||F|| = 8.061e+08  step = 1.000e+00  dx = 5.780e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 11 iterations with total CPU time of 4.29704 sec

************************************************************************
-- Nonlinear Solver Step 62 -- 
||F|| = 1.623e+09  step = 1.000e+00  dx = 2.332e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.13677 sec

************************************************************************
-- Nonlinear Solver Step 63 -- 
||F|| = 1.306e+09  step = 1.000e+00  dx = 6.032e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.28771 sec

************************************************************************
-- Nonlinear Solver Step 64 -- 
||F|| = 7.767e+08  step = 1.000e+00  dx = 6.318e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.41717 sec

************************************************************************
-- Nonlinear Solver Step 65 -- 
||F|| = 4.288e+08  step = 1.000e+00  dx = 1.291e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 12 iterations with total CPU time of 4.68541 sec

************************************************************************
-- Nonlinear Solver Step 66 -- 
||F|| = 1.518e+09  step = 1.000e+00  dx = 5.619e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.13376 sec

************************************************************************
-- Nonlinear Solver Step 67 -- 
||F|| = 1.303e+09  step = 1.000e+00  dx = 3.645e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.99133 sec

************************************************************************
-- Nonlinear Solver Step 68 -- 
||F|| = 2.241e+10  step = 1.000e+00  dx = 1.147e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.14787 sec

************************************************************************
-- Nonlinear Solver Step 69 -- 
||F|| = 6.156e+09  step = 1.000e+00  dx = 1.043e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.13781 sec

************************************************************************
-- Nonlinear Solver Step 70 -- 
||F|| = 4.376e+09  step = 1.000e+00  dx = 4.379e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.7307 sec

************************************************************************
-- Nonlinear Solver Step 71 -- 
||F|| = 5.227e+09  step = 1.000e+00  dx = 4.860e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.13844 sec

************************************************************************
-- Nonlinear Solver Step 72 -- 
||F|| = 2.782e+09  step = 1.000e+00  dx = 3.038e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.51599 sec

************************************************************************
-- Nonlinear Solver Step 73 -- 
||F|| = 2.101e+09  step = 1.000e+00  dx = 7.320e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.97534 sec

************************************************************************
-- Nonlinear Solver Step 74 -- 
||F|| = 1.188e+10  step = 1.000e+00  dx = 7.605e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.761217 sec

************************************************************************
-- Nonlinear Solver Step 75 -- 
||F|| = 5.904e+09  step = 1.000e+00  dx = 6.030e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.765621 sec

************************************************************************
-- Nonlinear Solver Step 76 -- 
||F|| = 3.058e+09  step = 1.000e+00  dx = 2.669e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.73987 sec

************************************************************************
-- Nonlinear Solver Step 77 -- 
||F|| = 7.732e+08  step = 1.000e+00  dx = 1.563e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.50738 sec

************************************************************************
-- Nonlinear Solver Step 78 -- 
||F|| = 2.386e+10  step = 1.000e+00  dx = 2.521e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.90076 sec

************************************************************************
-- Nonlinear Solver Step 79 -- 
||F|| = 2.097e+10  step = 1.000e+00  dx = 7.899e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 3.02642 sec

************************************************************************
-- Nonlinear Solver Step 80 -- 
||F|| = 1.378e+10  step = 1.000e+00  dx = 3.218e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 3.02535 sec

************************************************************************
-- Nonlinear Solver Step 81 -- 
||F|| = 2.963e+09  step = 1.000e+00  dx = 1.958e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.73752 sec

************************************************************************
-- Nonlinear Solver Step 82 -- 
||F|| = 4.694e+09  step = 1.000e+00  dx = 7.397e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.377355 sec

************************************************************************
-- Nonlinear Solver Step 83 -- 
||F|| = 4.107e+09  step = 1.000e+00  dx = 4.673e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.60999 sec

************************************************************************
-- Nonlinear Solver Step 84 -- 
||F|| = 3.482e+09  step = 1.000e+00  dx = 7.895e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.2787 sec

************************************************************************
-- Nonlinear Solver Step 85 -- 
||F|| = 2.959e+09  step = 1.000e+00  dx = 2.787e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.88737 sec

************************************************************************
-- Nonlinear Solver Step 86 -- 
||F|| = 2.096e+09  step = 1.000e+00  dx = 3.362e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.73443 sec

************************************************************************
-- Nonlinear Solver Step 87 -- 
||F|| = 1.702e+10  step = 1.000e+00  dx = 7.309e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.12895 sec

************************************************************************
-- Nonlinear Solver Step 88 -- 
||F|| = 1.509e+10  step = 1.000e+00  dx = 1.638e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 3.12166 sec

************************************************************************
-- Nonlinear Solver Step 89 -- 
||F|| = 1.558e+10  step = 1.000e+00  dx = 9.271e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.62098 sec

************************************************************************
-- Nonlinear Solver Step 90 -- 
||F|| = 5.895e+10  step = 1.000e+00  dx = 3.520e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.65251 sec

************************************************************************
-- Nonlinear Solver Step 91 -- 
||F|| = 3.476e+10  step = 1.000e+00  dx = 1.438e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.67579 sec

************************************************************************
-- Nonlinear Solver Step 92 -- 
||F|| = 2.545e+10  step = 1.000e+00  dx = 1.014e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.76417 sec

************************************************************************
-- Nonlinear Solver Step 93 -- 
||F|| = 1.671e+10  step = 1.000e+00  dx = 1.239e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.61138 sec

************************************************************************
-- Nonlinear Solver Step 94 -- 
||F|| = 1.251e+11  step = 1.000e+00  dx = 2.080e+08
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.50904 sec

************************************************************************
-- Nonlinear Solver Step 95 -- 
||F|| = 1.328e+11  step = 1.000e+00  dx = 5.513e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.90105 sec

************************************************************************
-- Nonlinear Solver Step 96 -- 
||F|| = 9.911e+10  step = 1.000e+00  dx = 1.061e+08
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.27363 sec

************************************************************************
-- Nonlinear Solver Step 97 -- 
||F|| = 6.106e+10  step = 1.000e+00  dx = 6.758e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.43256 sec

************************************************************************
-- Nonlinear Solver Step 98 -- 
||F|| = 3.233e+10  step = 1.000e+00  dx = 1.356e+08
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.37452 sec

************************************************************************
-- Nonlinear Solver Step 99 -- 
||F|| = 9.773e+10  step = 1.000e+00  dx = 1.841e+08
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.761926 sec

************************************************************************
-- Nonlinear Solver Step 100 -- 
||F|| = 8.181e+10  step = 1.000e+00  dx = 6.325e+07 (Failed!)
************************************************************************

************************************************************************
-- Final Status Test Results --
Failed.......OR Combination -> 
  **...........F-Norm = 8.181e+10 < 1.000e-04
               (Unscaled Two-Norm, Absolute Tolerance)
  Failed.......Number of Iterations = 100 < 100
************************************************************************
Nonlinear solver failed to converge!
Convergence Stats: for step  #1 : Newton, Krylov, Kr/Ne; LastKrylov, LastTol: 100  591  5.91  2 0.837
NOXFinish called
  
 Writing diagnostics to log file, time(yr) =    0.000000000000000     
 idiag_global, jdiag_global:            1            1
 local rank, i, j:            0            3            3
Application 14119425 resources: utime ~263s, stime ~2s Rss ~85436 inblocks ~72554 outblocks ~313633
