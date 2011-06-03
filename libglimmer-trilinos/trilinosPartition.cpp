#include <iostream>
#include "Epetra_LocalMap.h"
#ifdef GLIMMER_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Teuchos_RCP.hpp"
#include "config.inc"

// Define global variables.
static Teuchos::RCP<const Epetra_Map> partitionMap;


extern "C" {

  // doPartition and getPartition use Epetra to partition the global
  // problem in parallel. 
  // This is needed for serial-glimmer and parallel-trilinos,
  // These function will not be needed with a distributed glimmer
  // since the PDE code will pick the partitioning.


  void FC_FUNC(dopartition,DOPARTITION)
    (int& matrixSize, int& mySize) {


#ifdef GLIMMER_MPI
    // On Linux, Jaguar, the MPI_Init in Fortran is recopgnized by C++
    // On Bill's Mac, it is not, so this extra MPI_Init is needed
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
      int    argc;
      char** argv;
      MPI_Init(&argc, &argv);
    }
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    partitionMap = Teuchos::rcp(new Epetra_Map(matrixSize,1,comm) );
    mySize = partitionMap->NumMyElements();

    cout << "Trilinos Partition: doPartition has mySize = " << mySize << endl;
  }


  void FC_FUNC(getpartition,GETPARTITION) (int& mySize, int* myIndicies) {


      // Copy indices into array to send back to glimmer
      partitionMap->MyGlobalElements(myIndicies);

  }
} // extern"C"
