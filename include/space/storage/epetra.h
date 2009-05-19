/* Trilinos */
#include <Epetra_Map.h>
#include <Epetra_FEVector.h>
/* /Trilinos */

namespace Epetra {

  class PointFunction {
    Epetra_Vector coefficients;
  }

  class FEFunction 
  {
    Epetra_FEVector coefficients;

    FEFunction(Epetra_FEVector& cs) : coefficients(cs) {}
    FEFunction(){};

    
  };

  class FEOperator {
    
  };
};
