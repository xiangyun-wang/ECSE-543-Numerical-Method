import Foundation
// print("Hello world!")
// var a = 1
// var b: Int
// b = a
// print(sqrt(4))

/*
 * Input: Lower Triangular Matrix, RHS vector
 * Output: Matrix solution
 */
print ("Hello Swift")

func forward(Lower : [[Double]], b : [Double]) -> [Double]? {
    let n = b.count                                     // matrix dimension
    var output = Array(repeating: 0.0, count: n)        // output vector initialization
    var j: Int                                          // column tracker
    for i in 0...n-1 {                                  // loop over all rows
        output[i] = b[i]                                // set solution to the corresponding RHS vector entry
        j = 0                                           // reset column tracker
        while(j < i) {                                  // only have to compute when column number < row number (because matrix is lower triangular)
            output[i] -= Lower[i][j] * output[j]        // within the row, subtract the multiplication of knows variables and corresponding coefficients
            j += 1                                      // update the column tracker
        }
        // for j in 0...i-1 {
        //     output[i] -= Lower[i][j] * output[j]
        // }
        if Lower[i][i] == 0{                            // if the diagonal element of current row is 0, return nil
            return nil
        }
        output[i] /= Lower[i][i]                        // othereise, divide the current answer entry with the corresponding coefficient
    }
    return output
}

/*
 * Input: Upper Triangular Matrix, RHS vector
 * Output: Matrix solution
 *
 * Comments: Extremely similar to forward substitution. The only difference is, 
 * instead of looping from the first row to the last row, backward substitution 
 * loops fro the last row to the first row since the matrix is upper triangular
 * (no extra comments will be added below)
 */
func backward(Upper : [[Double]], b : [Double]) -> [Double]? {
    let n = b.count
    var output = Array(repeating: 0.0, count: n)
    var j: Int
    for i in stride (from: n-1, through: 0, by: -1) {
        output[i] = b[i]
        j = n-1
        while(j > i) {
            output[i] -= Upper[i][j]*output[j]
            j -= 1
        }
        if Upper[i][i] == 0{
            return nil
        }
        output[i] /= Upper[i][i]
    }
    return output  
}

/*
 * Input: Matrix A to be decomposed
 * Output: Matrix L, where A = L*transpose(L)
 * 
 */
func cholesky_decomposition(A : [[Double]]) -> [[Double]]? {
    let n = A.count
    var extra_term: Double
    var output = Array(repeating: Array(repeating: 0.0, count: n), count: n)
    var k: Int
    for i in 0...n-1 {
        for j in 0...i {
            extra_term = 0
            k = 0
            while (k<j){
                extra_term += output[i][k] * output[j][k]
                k += 1;
            }
            // ----- implies that A is not a symmetric positive definite matrix -----
            if A[i][j] - extra_term <= 0 {
                return nil
            }
            // ----- end of check for SPD ------
            if i == j {
                output[i][j] = sqrt(A[i][i] - extra_term)
            } else {
                output[i][j] = 1.0 / output[j][j] * (A[i][j] - extra_term)
            }
            // if (output[i][j] == 0) {
            //     return nil
            // }
        }
    }
    return output
}

func find_half_band(A: [[Double]]) -> Int {
    let n = A.count
    var band = 0
    for i in 0...n-1 {
        for j in 0...i {
            
        }
    }
}

func chelesky_decomposition_optimized(A: [[Double]]) -> [Double]? {
   let n = A.count
    var extra_term: Double
    var output = Array(repeating: Array(repeating: 0.0, count: n), count: n)
    var k: Int
    for i in 0...n-1 {
        for j in 0...i {
            extra_term = 0
            k = 0
            while (k<j){
                extra_term += output[i][k] * output[j][k]
                k += 1;
            }
            // ----- implies that A is not a symmetric positive definite matrix -----
            if A[i][j] - extra_term <= 0 {
                return nil
            }
            // ----- end of check for SPD ------
            if i == j {
                output[i][j] = sqrt(A[i][i] - extra_term)
            } else {
                output[i][j] = 1.0 / output[j][j] * (A[i][j] - extra_term)
            }
            // if (output[i][j] == 0) {
            //     return nil
            // }
        }
    }
    return output 
}

func cheleskySolver(A: [[Double]], b: [Double]) -> [Double]? {
    return nil
}

// func PSD_Check(A: [[Double]]) -> Bool {
//     return false
// }

func dotproduct_vector(A: [[Double]], b:[Double]) -> [Double]? {
    return nil
}

func transpose(A: [[Double]]) -> [[Double]]? {
    return nil
}

func PSD_Generation(size: Int) -> [[Double]]? {
    return nil
}

func Read_Circuit_File(Dir: String) {

}