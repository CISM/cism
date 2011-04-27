#include "Teuchos_ConfigDefs.hpp"

#ifdef GLIMMER_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_Version.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Trilinos_Util.h"

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

extern "C" {

#ifdef _xlC_
  void solvetriadmatrixwithtrilinos(int& nnz, int& order, int* row, 
              int* col, double* val, double* rhs, double* solution) {
#else
  void solvetriadmatrixwithtrilinos_(int& nnz, int& order, int* row, 
              int* col, double* val, double* rhs, double* solution) {
#endif

    try{
    
#ifdef GLIMMER_MPI
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    int i, j, ierr;
    int MyPID = Comm.MyPID();
    bool verbose = (MyPID == 0);
    Epetra_Map RowMap(order, 0, Comm);
    int NumMyElements = RowMap.NumMyElements();
    int *MyGlobalElements = new int[NumMyElements];
    RowMap.MyGlobalElements(&MyGlobalElements[0]);

#ifdef GLIMMER_MPI
    int nPEs;
    MPI_Comm_size(MPI_COMM_WORLD, &nPEs);
#endif

    int anEst = nnz / order + 1;
    Epetra_CrsMatrix A(Copy, RowMap, anEst);
    
    for (j=0; j<nnz; ++j) {
      if (RowMap.MyGID(row[j]) ) {
	ierr = A.InsertGlobalValues(row[j], 1, &(val[j]), &(col[j]) );
	assert(ierr >= 0);
      }
    }
    
    ierr = A.FillComplete();
    assert(ierr == 0);
    
    //-------------------------------------------------------------------------
    // RN_20091221: Taking care of the rhs
    //-------------------------------------------------------------------------
    Epetra_Vector b(RowMap);

    // Inserting values into the rhs
    double *MyGlobalValues = new double[NumMyElements];
    for (j=0; j<NumMyElements; ++j) {
      MyGlobalValues[j] = rhs[MyGlobalElements[j] ];
    }
    ierr = b.ReplaceGlobalValues(NumMyElements, &MyGlobalValues[0],
				 &MyGlobalElements[0]);

    //-------------------------------------------------------------------------
    // RN_20091221: Taking care of the solution
    //-------------------------------------------------------------------------
    Epetra_Vector x(RowMap);

    Teuchos::ParameterList paramList;

    Teuchos::RCP<Teuchos::ParameterList>
      paramList1 = Teuchos::rcp(&paramList, false);
    Teuchos::updateParametersFromXmlFile("./strat1.xml", paramList1.get() );

    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    Teuchos::RCP<Epetra_CrsMatrix> epetraOper = Teuchos::rcp(&A, false);
    Teuchos::RCP<Epetra_Vector> epetraRhs = Teuchos::rcp(&b, false);
    Teuchos::RCP<Epetra_Vector> epetraSol = Teuchos::rcp(&x, false);

    Teuchos::RCP<const Thyra::LinearOpBase<double> >
      thyraOper = Thyra::epetraLinearOp(epetraOper);
    Teuchos::RCP<Thyra::VectorBase<double> >
      thyraRhs = Thyra::create_Vector(epetraRhs, thyraOper->range() );
    Teuchos::RCP<Thyra::VectorBase<double> >
      thyraSol = Thyra::create_Vector(epetraSol, thyraOper->domain() );

    linearSolverBuilder.setParameterList(Teuchos::rcp(&paramList, false) );

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> >
      lows = Thyra::linearOpWithSolve(*lowsFactory, thyraOper);

    Thyra::SolveStatus<double>
      status = Thyra::solve(*lows, Thyra::NOTRANS, *thyraRhs, &*thyraSol);

    thyraSol = Teuchos::null;

    // For debugging =)
    //    cout << "A: " << A << endl;
    //    cout << "b: " << b << endl;
    //    cout << "x: " << x << endl;

    //    Epetra_Vector temp(RowMap);
    //    double residualNorm;
    //    A.Multiply(false, x, temp);
    //    temp.Update(-1, b, 1);
    //    temp.Norm2(&residualNorm);

    Epetra_LocalMap localMap(order, 0, Comm);
    Epetra_Vector xExtra(localMap); // local vector in each processor

    // RN_20091218: Create an import map and then import the data to the
    // local vector.
    Epetra_Import import(localMap, RowMap);

    xExtra.Import(x, import, Add);
    xExtra.ExtractCopy(solution);

    delete[] MyGlobalElements;
    delete[] MyGlobalValues;
    }
  
    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (string& s) {
      cout << s << endl;
    }
    catch (char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception!" << endl;
    }
  }

} // extern "C"
