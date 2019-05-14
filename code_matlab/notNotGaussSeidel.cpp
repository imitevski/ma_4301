/* MyMEXFunction
 * Adds second input to each
 * element of first input
 * a = MyMEXFunction(a,b);
 */

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <cmath>

using namespace matlab::data;
using namespace std;
using matlab::mex::ArgumentList;
using matlab::engine::convertUTF8StringToUTF16String;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        //checkArguments(outputs, inputs);
        TypedArray<double> u = move(inputs[0]);
        TypedArray<double> F = move(inputs[1]);
        const double h = inputs[2][0];
        const int m = inputs[3][0];
        const int n = inputs[4][0];
        
        double A_xy, A_vw;
        
        //mexPrintf("The reciprocal of h is %f\n",3.5);
        //mexEvalString("drawnow;");
        
        for (int i = 1; i < m - 1; i++) {
            
            for (int j = 1; j < n - 1; j++) {
                
                A_xy = (1/h)*(1/h)*(1/h)*(1/h)
                      *(u[i-1][j]+u[i+1][j]-2*u[i][j])
                      *(u[i][j-1]+u[i][j+1]-2*u[i][j]);
                
                A_vw = 0.25*(1/h)*(1/h)*(1/h)*(1/h)
                      *(u[i-1][j-1]+u[i+1][j+1]-2*u[i][j])
                      *(u[i+1][j-1]+u[i-1][j+1]-2*u[i][j]);
                
                if (A_xy <= A_vw) {
                    
                    u[i][j] = 0.25*(
                                    u[i+1][j]+u[i-1][j]+u[i][j-1]+u[i][j+1]
                                   )
                    -0.5*sqrt(
                              0.25*(
                                    (u[i+1][j]+u[i-1][j]-u[i][j-1]-u[i][j+1])
                                   *(u[i+1][j]+u[i-1][j]-u[i][j-1]-u[i][j+1])
                                   )
                    +h*h*h*h*F[i][j]);
                
                }   else   {
                    
                    u[i][j] = 0.25*(
                                    u[i-1][j+1]+u[i+1][j-1]
                                   +u[i+1][j+1]+u[i-1][j-1]
                                   )
                    -0.5*sqrt(
                              0.25*(
                                    (u[i-1][j+1]+u[i+1][j-1]-u[i+1][j+1]-u[i-1][j-1])
                                   *(u[i-1][j+1]+u[i+1][j-1]-u[i+1][j+1]-u[i-1][j-1])
                                   )
                    +4*h*h*h*h*F[i][j]);
                
                }
            
            }

        }
        
        outputs[0] = u;
        outputs[1][0] = 1/h;
        
    /*
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        
        // Get array factory
        ArrayFactory factory;
        
        // Check first input argument
        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getType() == ArrayType::COMPLEX_DOUBLE ||
            inputs[0].getNumberOfElements() != 1)
        {
            matlabPtr->feval(convertUTF8StringToUTF16String("error"),
                             0,
                             std::vector<Array>({ factory.createScalar("First input must scalar double") }));
        }
        
        // Check second input argument
        if (inputs[1].getType() != ArrayType::DOUBLE ||
            inputs[1].getType() == ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(convertUTF8StringToUTF16String("error"),
                             0,
                             std::vector<Array>({ factory.createScalar("Input must double array") }));
        }
        // Check number of outputs
        if (outputs.size() > 1) {
            matlabPtr->feval(convertUTF8StringToUTF16String("error"),
                             0,
                             std::vector<Array>({ factory.createScalar("Only one output is returned") }));
        }
    }
     */
    }
    
};
