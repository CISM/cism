#include <iostream>
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_CombineMode.h"
#include "matrixInterface.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Time.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

// Define variables that are global to this file.
Teuchos::RCP<TrilinosMatrix_Interface> interface;
Teuchos::RCP<Epetra_CrsMatrix> savedMatrix_A;
Teuchos::RCP<Epetra_CrsMatrix> savedMatrix_C;
Teuchos::RCP<Epetra_Vector> soln;
Teuchos::RCP<Teuchos::ParameterList> pl;
Teuchos::RCP<Teuchos::FancyOStream> out;
Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows;
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory;
Teuchos::RCP<const Thyra::LinearOpBase<double> > thyraOper;

extern "C" {
  //================================================================
  //================================================================
  // RN_20091215: This needs to be called only once in the beginning
  // to set up the problem.
  //================================================================

#ifdef _xlC_
  void inittrilinos(int& bandwidth, int& mySize, int* myIndicies) {
#else
  void inittrilinos_(int& bandwidth, int& mySize, int* myIndicies) {
#endif

#ifdef GLIMMER_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    Teuchos::RCP<const Epetra_Map> rowMap = 
      Teuchos::rcp(new Epetra_Map(-1,mySize,myIndicies,1,comm) );

    soln = Teuchos::rcp(new Epetra_Vector(*rowMap));
    pl = Teuchos::getParametersFromXmlFile("strat1.xml");

    out = Teuchos::VerboseObjectBase::getDefaultOStream();

    // Create an interface that holds a CrsMatrix instance and some useful methods.
    interface = Teuchos::rcp(new TrilinosMatrix_Interface(rowMap, bandwidth, comm));

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(pl);
    lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");
    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

    lows=Teuchos::null;
    thyraOper=Teuchos::null;
  }

  //============================================================
  // RN_20091118: This is to update the matrix with new entries.
  //============================================================
#ifdef _xlC_
  void putintotrilinosmatrix(int& rowInd, int& colInd, double& val) {
#else
  void putintotrilinosmatrix_(int& rowInd, int& colInd, double& val) {
#endif

    const Epetra_Map& map = interface->getRowMap();
    // If this row is not owned on this processor, then do nothing
    if (!map.MyGID(rowInd)) return;

    Epetra_CrsMatrix& matrix = *(interface->getOperator());

    if (!interface->isSparsitySet()) {
      // The matrix has not been "FillComplete()"ed. First fill of time step.
      int ierr = matrix.InsertGlobalValues(rowInd, 1, &val, &colInd);
      if (ierr<0) {cout << "Error Code for " << rowInd << "  " << colInd << "  = ("<< ierr <<")"<<endl; exit(1);}
      else if (ierr>0) cout << "Warning Code for " << rowInd << "  " << colInd << "  = ("<< ierr <<")"<<endl;
    }
    else {
      // Subsequent matrix fills of each time step.
      int ierr = matrix.ReplaceGlobalValues(rowInd, 1, &val, &colInd);
    
      if (ierr != 0) { // Sparsity pattern has changed. Create fresh matrix
	cout << "Warning: Trilinos matrix has detected a new entry (" 
             << rowInd << ", " << colInd << ", " << val 
             << ")\n\t that did not exist before. A new matrix will be formed!"
             << "\n\t This is expensive, and we should figure out why this is"
             << "\n\t happening and avoid it! -AGS" << endl;

	int matrixSize = interface->matrixOrder();
	int bandwidth = interface->bandwidth();
	
	Teuchos::RCP<Epetra_CrsMatrix> newMatrix =
	  Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, bandwidth) );
	
	int numEntries;
	double *values = new double[bandwidth];
	int *indices = new int[bandwidth];
	
	// Copy the old matrix to the new matrix.
	for (int j=0; j<matrixSize; ++j) {
	  if (map.MyGID(j) ) {
	    int aNumber = bandwidth;
	    ierr = matrix.ExtractGlobalRowCopy(j, aNumber, numEntries,
						values, indices);
	    assert(ierr >= 0);
	    ierr = newMatrix->InsertGlobalValues(j, numEntries, &(values[0]),
						 &(indices[0]) );
	    assert(ierr >= 0);
	  }
	}
	
	// Insert the new entry.
	if (map.MyGID(rowInd) ) {
	  ierr = newMatrix->InsertGlobalValues(rowInd, 1, &val, &colInd);
	}

	interface->updateOperator(newMatrix);
	
	delete[] values;
	delete[] indices;
      }
    }
  }

  //========================================================
  // RN_20091118: This is to make calls to Trilinos solvers.
  //========================================================
#ifdef _xlC_
  void solvewithtrilinos(double* rhs, double* answer, double& elapsedTime) {
#else
  void solvewithtrilinos_(double* rhs, double* answer, double& elapsedTime) {
#endif
    //Teuchos::Time linearTime("LinearTime"); linearTime.start();

    // Lock in sparsity pattern
    if (!interface->isSparsitySet()) interface->finalizeSparsity();

    Teuchos::RCP<Epetra_Vector> epetraSol = soln;
    Teuchos::RCP<Epetra_Vector> epetraRhs = interface->getPartitionedVec(rhs);

    thyraOper = Thyra::epetraLinearOp(interface->getOperator());
    Teuchos::RCP<Thyra::VectorBase<double> >
      thyraRhs = Thyra::create_Vector(epetraRhs, thyraOper->range() );
    Teuchos::RCP<Thyra::VectorBase<double> >
      thyraSol = Thyra::create_Vector(epetraSol, thyraOper->domain() );

    lows = Thyra::linearOpWithSolve(*lowsFactory, thyraOper);

    Thyra::SolveStatus<double>
      status = Thyra::solve(*lows, Thyra::NOTRANS, *thyraRhs, &*thyraSol);

    interface->spreadVector(*soln, answer);

    //elapsedTime = linearTime.stop();*out << "Total time elapsed for calling Solve(): " << elapsedTime << endl;
  }

#ifdef _xlC_
  void savetrilinosmatrix(int* i) {
#else
  void savetrilinosmatrix_(int* i) {
#endif
    if (!interface->isSparsitySet()) interface->finalizeSparsity();
    if (*i==0)
      savedMatrix_A = Teuchos::rcp(new Epetra_CrsMatrix(*(interface->getOperator())));
    else if (*i==1)
      savedMatrix_C = Teuchos::rcp(new Epetra_CrsMatrix(*(interface->getOperator())));
    else
      assert(false);
  }

#ifdef _xlC_
  void restoretrilinosmatrix(int* i) {
#else
  void restoretrilinosmatrix_(int* i) {
#endif
    if (*i==0)
      interface->updateOperator(savedMatrix_A); 
    else if (*i==1)
      interface->updateOperator(savedMatrix_C); 
    else
      assert(false);
  }

#ifdef _xlC_
  void matvecwithtrilinos(double* x, double* answer) {
#else
  void matvecwithtrilinos_(double* x, double* answer) {
#endif
    const Epetra_Map& map = interface->getRowMap(); 
    Teuchos::RCP<Epetra_Vector> epetra_x = interface->getPartitionedVec(x);
    Epetra_Vector y(map);
    interface->getOperator()->Multiply(false, *epetra_x, y);
    interface->spreadVector(y, answer);
  }

} // extern"C"
