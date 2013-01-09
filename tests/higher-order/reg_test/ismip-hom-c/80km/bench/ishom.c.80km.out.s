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
||F|| = 1.358e+02  step = 0.000e+00  dx = 0.000e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 49 iterations with total CPU time of 16.7861 sec

************************************************************************
-- Nonlinear Solver Step 1 -- 
||F|| = 1.151e+08  step = 1.000e+00  dx = 4.689e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.969223 sec

************************************************************************
-- Nonlinear Solver Step 2 -- 
||F|| = 9.701e+07  step = 1.000e+00  dx = 1.172e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.28987 sec

************************************************************************
-- Nonlinear Solver Step 3 -- 
||F|| = 7.547e+07  step = 1.000e+00  dx = 2.180e-01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.58587 sec

************************************************************************
-- Nonlinear Solver Step 4 -- 
||F|| = 3.400e+07  step = 1.000e+00  dx = 7.943e-01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.25938 sec

************************************************************************
-- Nonlinear Solver Step 5 -- 
||F|| = 1.871e+07  step = 1.000e+00  dx = 1.766e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.29235 sec

************************************************************************
-- Nonlinear Solver Step 6 -- 
||F|| = 6.892e+06  step = 1.000e+00  dx = 4.845e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.29714 sec

************************************************************************
-- Nonlinear Solver Step 7 -- 
||F|| = 3.435e+06  step = 1.000e+00  dx = 9.053e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.30643 sec

************************************************************************
-- Nonlinear Solver Step 8 -- 
||F|| = 1.928e+06  step = 1.000e+00  dx = 4.001e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.970618 sec

************************************************************************
-- Nonlinear Solver Step 9 -- 
||F|| = 9.402e+05  step = 1.000e+00  dx = 2.653e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.265 sec

************************************************************************
-- Nonlinear Solver Step 10 -- 
||F|| = 6.819e+05  step = 1.000e+00  dx = 8.820e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.29785 sec

************************************************************************
-- Nonlinear Solver Step 11 -- 
||F|| = 6.221e+05  step = 1.000e+00  dx = 1.164e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.32557 sec

************************************************************************
-- Nonlinear Solver Step 12 -- 
||F|| = 3.428e+05  step = 1.000e+00  dx = 4.023e+00
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.26019 sec

************************************************************************
-- Nonlinear Solver Step 13 -- 
||F|| = 6.373e+05  step = 1.000e+00  dx = 9.415e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.327312 sec

************************************************************************
-- Nonlinear Solver Step 14 -- 
||F|| = 4.997e+05  step = 1.000e+00  dx = 1.042e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.643734 sec

************************************************************************
-- Nonlinear Solver Step 15 -- 
||F|| = 2.945e+05  step = 1.000e+00  dx = 1.147e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.66043 sec

************************************************************************
-- Nonlinear Solver Step 16 -- 
||F|| = 1.478e+06  step = 1.000e+00  dx = 1.197e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.644621 sec

************************************************************************
-- Nonlinear Solver Step 17 -- 
||F|| = 1.286e+06  step = 1.000e+00  dx = 1.460e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.28718 sec

************************************************************************
-- Nonlinear Solver Step 18 -- 
||F|| = 9.670e+05  step = 1.000e+00  dx = 4.780e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 1.94486 sec

************************************************************************
-- Nonlinear Solver Step 19 -- 
||F|| = 6.875e+05  step = 1.000e+00  dx = 1.153e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.25919 sec

************************************************************************
-- Nonlinear Solver Step 20 -- 
||F|| = 3.564e+05  step = 1.000e+00  dx = 3.861e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 2.91095 sec

************************************************************************
-- Nonlinear Solver Step 21 -- 
||F|| = 2.054e+05  step = 1.000e+00  dx = 3.043e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 2.91501 sec

************************************************************************
-- Nonlinear Solver Step 22 -- 
||F|| = 1.133e+05  step = 1.000e+00  dx = 8.962e+01
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.03435 sec

************************************************************************
-- Nonlinear Solver Step 23 -- 
||F|| = 2.234e+05  step = 1.000e+00  dx = 5.582e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.646092 sec

************************************************************************
-- Nonlinear Solver Step 24 -- 
||F|| = 1.074e+05  step = 1.000e+00  dx = 9.201e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.645392 sec

************************************************************************
-- Nonlinear Solver Step 25 -- 
||F|| = 5.685e+04  step = 1.000e+00  dx = 4.210e+02
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.37659 sec

************************************************************************
-- Nonlinear Solver Step 26 -- 
||F|| = 4.721e+04  step = 1.000e+00  dx = 2.036e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.965932 sec

************************************************************************
-- Nonlinear Solver Step 27 -- 
||F|| = 2.709e+04  step = 1.000e+00  dx = 7.069e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.69001 sec

************************************************************************
-- Nonlinear Solver Step 28 -- 
||F|| = 1.208e+04  step = 1.000e+00  dx = 9.085e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 10 iterations with total CPU time of 3.32115 sec

************************************************************************
-- Nonlinear Solver Step 29 -- 
||F|| = 1.679e+04  step = 1.000e+00  dx = 3.861e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.021 sec

************************************************************************
-- Nonlinear Solver Step 30 -- 
||F|| = 6.813e+04  step = 1.000e+00  dx = 3.503e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.28971 sec

************************************************************************
-- Nonlinear Solver Step 31 -- 
||F|| = 5.527e+04  step = 1.000e+00  dx = 9.677e+03
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.01282 sec

************************************************************************
-- Nonlinear Solver Step 32 -- 
||F|| = 4.959e+04  step = 1.000e+00  dx = 1.324e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 10 iterations with total CPU time of 3.32302 sec

************************************************************************
-- Nonlinear Solver Step 33 -- 
||F|| = 8.722e+04  step = 1.000e+00  dx = 8.446e+04
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.68055 sec

************************************************************************
-- Nonlinear Solver Step 34 -- 
||F|| = 2.365e+05  step = 1.000e+00  dx = 8.690e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.03548 sec

************************************************************************
-- Nonlinear Solver Step 35 -- 
||F|| = 1.095e+06  step = 1.000e+00  dx = 4.750e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.61739 sec

************************************************************************
-- Nonlinear Solver Step 36 -- 
||F|| = 6.647e+05  step = 1.000e+00  dx = 7.519e+05
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.35349 sec

************************************************************************
-- Nonlinear Solver Step 37 -- 
||F|| = 4.382e+05  step = 1.000e+00  dx = 6.111e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.67435 sec

************************************************************************
-- Nonlinear Solver Step 38 -- 
||F|| = 6.117e+05  step = 1.000e+00  dx = 1.977e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 1.05483 sec

************************************************************************
-- Nonlinear Solver Step 39 -- 
||F|| = 2.089e+06  step = 1.000e+00  dx = 1.509e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.28901 sec

************************************************************************
-- Nonlinear Solver Step 40 -- 
||F|| = 1.308e+06  step = 1.000e+00  dx = 4.634e+06
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.02707 sec

************************************************************************
-- Nonlinear Solver Step 41 -- 
||F|| = 5.639e+06  step = 1.000e+00  dx = 1.302e+08
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.70361 sec

************************************************************************
-- Nonlinear Solver Step 42 -- 
||F|| = 9.832e+06  step = 1.000e+00  dx = 1.057e+08
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.29973 sec

************************************************************************
-- Nonlinear Solver Step 43 -- 
||F|| = 3.108e+06  step = 1.000e+00  dx = 5.966e+07
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.3855 sec

************************************************************************
-- Nonlinear Solver Step 44 -- 
||F|| = 1.932e+08  step = 1.000e+00  dx = 1.710e+09
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.70077 sec

************************************************************************
-- Nonlinear Solver Step 45 -- 
||F|| = 2.248e+10  step = 1.000e+00  dx = 7.670e+10
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.02063 sec

************************************************************************
-- Nonlinear Solver Step 46 -- 
||F|| = 8.901e+10  step = 1.000e+00  dx = 7.045e+12
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.645972 sec

************************************************************************
-- Nonlinear Solver Step 47 -- 
||F|| = 7.493e+10  step = 1.000e+00  dx = 2.251e+13
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.966982 sec

************************************************************************
-- Nonlinear Solver Step 48 -- 
||F|| = 5.319e+10  step = 1.000e+00  dx = 2.099e+13
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.0069 sec

************************************************************************
-- Nonlinear Solver Step 49 -- 
||F|| = 2.707e+10  step = 1.000e+00  dx = 1.886e+13
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.03418 sec

************************************************************************
-- Nonlinear Solver Step 50 -- 
||F|| = 2.404e+11  step = 1.000e+00  dx = 3.482e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.970124 sec

************************************************************************
-- Nonlinear Solver Step 51 -- 
||F|| = 2.060e+11  step = 1.000e+00  dx = 2.651e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.67481 sec

************************************************************************
-- Nonlinear Solver Step 52 -- 
||F|| = 5.491e+12  step = 1.000e+00  dx = 5.015e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.28901 sec

************************************************************************
-- Nonlinear Solver Step 53 -- 
||F|| = 4.072e+12  step = 1.000e+00  dx = 5.469e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.61846 sec

************************************************************************
-- Nonlinear Solver Step 54 -- 
||F|| = 2.112e+12  step = 1.000e+00  dx = 5.528e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.27383 sec

************************************************************************
-- Nonlinear Solver Step 55 -- 
||F|| = 7.454e+11  step = 1.000e+00  dx = 1.193e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.04373 sec

************************************************************************
-- Nonlinear Solver Step 56 -- 
||F|| = 6.634e+11  step = 1.000e+00  dx = 1.632e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.647133 sec

************************************************************************
-- Nonlinear Solver Step 57 -- 
||F|| = 4.948e+11  step = 1.000e+00  dx = 5.306e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 1.95395 sec

************************************************************************
-- Nonlinear Solver Step 58 -- 
||F|| = 4.659e+11  step = 1.000e+00  dx = 9.325e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.970478 sec

************************************************************************
-- Nonlinear Solver Step 59 -- 
||F|| = 3.596e+11  step = 1.000e+00  dx = 5.505e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 1.9443 sec

************************************************************************
-- Nonlinear Solver Step 60 -- 
||F|| = 4.533e+11  step = 1.000e+00  dx = 9.152e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.29199 sec

************************************************************************
-- Nonlinear Solver Step 61 -- 
||F|| = 9.965e+11  step = 1.000e+00  dx = 6.625e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.651347 sec

************************************************************************
-- Nonlinear Solver Step 62 -- 
||F|| = 8.820e+11  step = 1.000e+00  dx = 1.034e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.652282 sec

************************************************************************
-- Nonlinear Solver Step 63 -- 
||F|| = 8.217e+11  step = 1.000e+00  dx = 1.126e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.292 sec

************************************************************************
-- Nonlinear Solver Step 64 -- 
||F|| = 8.231e+11  step = 1.000e+00  dx = 1.621e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.647323 sec

************************************************************************
-- Nonlinear Solver Step 65 -- 
||F|| = 9.688e+11  step = 1.000e+00  dx = 5.504e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.301 sec

************************************************************************
-- Nonlinear Solver Step 66 -- 
||F|| = 8.484e+11  step = 1.000e+00  dx = 6.067e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.70928 sec

************************************************************************
-- Nonlinear Solver Step 67 -- 
||F|| = 2.527e+11  step = 1.000e+00  dx = 4.958e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.37085 sec

************************************************************************
-- Nonlinear Solver Step 68 -- 
||F|| = 1.617e+11  step = 1.000e+00  dx = 5.481e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.37663 sec

************************************************************************
-- Nonlinear Solver Step 69 -- 
||F|| = 6.820e+11  step = 1.000e+00  dx = 7.854e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.982294 sec

************************************************************************
-- Nonlinear Solver Step 70 -- 
||F|| = 4.000e+11  step = 1.000e+00  dx = 6.369e+15
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.72585 sec

************************************************************************
-- Nonlinear Solver Step 71 -- 
||F|| = 1.415e+12  step = 1.000e+00  dx = 2.743e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.96521 sec

************************************************************************
-- Nonlinear Solver Step 72 -- 
||F|| = 1.233e+12  step = 1.000e+00  dx = 1.462e+16
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 9 iterations with total CPU time of 3.05475 sec

************************************************************************
-- Nonlinear Solver Step 73 -- 
||F|| = 1.739e+13  step = 1.000e+00  dx = 3.199e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.329026 sec

************************************************************************
-- Nonlinear Solver Step 74 -- 
||F|| = 1.262e+13  step = 1.000e+00  dx = 2.088e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.980734 sec

************************************************************************
-- Nonlinear Solver Step 75 -- 
||F|| = 8.715e+12  step = 1.000e+00  dx = 1.306e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.61568 sec

************************************************************************
-- Nonlinear Solver Step 76 -- 
||F|| = 2.921e+12  step = 1.000e+00  dx = 7.896e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.658103 sec

************************************************************************
-- Nonlinear Solver Step 77 -- 
||F|| = 9.552e+11  step = 1.000e+00  dx = 2.837e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.37983 sec

************************************************************************
-- Nonlinear Solver Step 78 -- 
||F|| = 3.935e+12  step = 1.000e+00  dx = 5.490e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.647108 sec

************************************************************************
-- Nonlinear Solver Step 79 -- 
||F|| = 2.690e+12  step = 1.000e+00  dx = 1.890e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.61924 sec

************************************************************************
-- Nonlinear Solver Step 80 -- 
||F|| = 1.757e+12  step = 1.000e+00  dx = 4.074e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.658409 sec

************************************************************************
-- Nonlinear Solver Step 81 -- 
||F|| = 7.793e+11  step = 1.000e+00  dx = 1.937e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 1.9297 sec

************************************************************************
-- Nonlinear Solver Step 82 -- 
||F|| = 3.105e+11  step = 1.000e+00  dx = 2.350e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.34504 sec

************************************************************************
-- Nonlinear Solver Step 83 -- 
||F|| = 1.188e+12  step = 1.000e+00  dx = 8.287e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 2 iterations with total CPU time of 0.652017 sec

************************************************************************
-- Nonlinear Solver Step 84 -- 
||F|| = 6.623e+11  step = 1.000e+00  dx = 7.767e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.971302 sec

************************************************************************
-- Nonlinear Solver Step 85 -- 
||F|| = 2.763e+11  step = 1.000e+00  dx = 1.840e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.69878 sec

************************************************************************
-- Nonlinear Solver Step 86 -- 
||F|| = 3.299e+11  step = 1.000e+00  dx = 1.766e+17
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 6 iterations with total CPU time of 2.03104 sec

************************************************************************
-- Nonlinear Solver Step 87 -- 
||F|| = 1.918e+12  step = 1.000e+00  dx = 2.444e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.29243 sec

************************************************************************
-- Nonlinear Solver Step 88 -- 
||F|| = 1.761e+12  step = 1.000e+00  dx = 1.016e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.61827 sec

************************************************************************
-- Nonlinear Solver Step 89 -- 
||F|| = 5.221e+11  step = 1.000e+00  dx = 2.442e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 8 iterations with total CPU time of 2.681 sec

************************************************************************
-- Nonlinear Solver Step 90 -- 
||F|| = 9.094e+12  step = 1.000e+00  dx = 9.636e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.96542 sec

************************************************************************
-- Nonlinear Solver Step 91 -- 
||F|| = 5.741e+12  step = 1.000e+00  dx = 9.185e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 4 iterations with total CPU time of 1.29135 sec

************************************************************************
-- Nonlinear Solver Step 92 -- 
||F|| = 3.505e+12  step = 1.000e+00  dx = 6.196e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.36637 sec

************************************************************************
-- Nonlinear Solver Step 93 -- 
||F|| = 1.208e+12  step = 1.000e+00  dx = 2.746e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.35617 sec

************************************************************************
-- Nonlinear Solver Step 94 -- 
||F|| = 8.497e+12  step = 1.000e+00  dx = 1.918e+19
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.97015 sec

************************************************************************
-- Nonlinear Solver Step 95 -- 
||F|| = 7.644e+12  step = 1.000e+00  dx = 7.808e+18
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 7 iterations with total CPU time of 2.3537 sec

************************************************************************
-- Nonlinear Solver Step 96 -- 
||F|| = 1.524e+14  step = 1.000e+00  dx = 8.894e+20
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 1 iterations with total CPU time of 0.328072 sec

************************************************************************
-- Nonlinear Solver Step 97 -- 
||F|| = 1.214e+14  step = 1.000e+00  dx = 3.730e+20
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 3 iterations with total CPU time of 0.969467 sec

************************************************************************
-- Nonlinear Solver Step 98 -- 
||F|| = 9.107e+13  step = 1.000e+00  dx = 3.204e+20
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.61823 sec

************************************************************************
-- Nonlinear Solver Step 99 -- 
||F|| = 4.567e+13  step = 1.000e+00  dx = 3.031e+20
************************************************************************

 
 The Belos solver of type "Belos::BlockGmresSolMgr<...,double>{Variant='Flexible', Ortho Type='DGKS', Block Size=1, Num Blocks=100, Max Restarts=20}" returned a solve status of "SOLVE_STATUS_CONVERGED" in 5 iterations with total CPU time of 1.61769 sec

************************************************************************
-- Nonlinear Solver Step 100 -- 
||F|| = 7.473e+12  step = 1.000e+00  dx = 2.616e+20 (Failed!)
************************************************************************

************************************************************************
-- Final Status Test Results --
Failed.......OR Combination -> 
  **...........F-Norm = 7.473e+12 < 1.000e-04
               (Unscaled Two-Norm, Absolute Tolerance)
  Failed.......Number of Iterations = 100 < 100
************************************************************************
Nonlinear solver failed to converge!
Convergence Stats: for step  #1 : Newton, Krylov, Kr/Ne; LastKrylov, LastTol: 100  533  5.33  5 0.165
NOXFinish called
  
 Writing diagnostics to log file, time(yr) =    0.000000000000000     
 idiag_global, jdiag_global:            1            1
 local rank, i, j:            0            3            3
Application 13924814 resources: utime ~213s, stime ~1s Rss ~100372 inblocks ~72586 outblocks ~313633
