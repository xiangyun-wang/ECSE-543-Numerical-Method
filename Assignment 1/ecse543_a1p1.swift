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
    var band: Int = 0
    var N_Zero: Int = 0
    for i in 0...A.count-1 {
        N_Zero = 0
        for j in 0...i{
            if A[i][j] == 0{
                N_Zero += 1;
            }else{
                break
            }
        }
        if i+1-N_Zero > band {
            band = i+1-N_Zero        
        }
    }

    return band
    // initially set band = 0, 
    // loop over all rows
    // counter add until first non-zero element is met in a row (N_zero)
    // if i+1-N_zero > band, then band = i+1-N_zero
    // return band
}

func chelesky_decomposition_optimized_bandwidth_forward(A: [[Double]], RHS: [Double], band: Int) -> ([[Double]], [Double])? {
    let n = A.count
    var A_a = A
    var RHS_a = RHS
    var i_range: Int

    for j in 0...n-1 {
        if A_a[j][j] <= 0 {
            return nil
        }
        A_a[j][j] = sqrt(A_a[j][j])
        RHS_a[j] = RHS_a[j]/A_a[j][j]
        if j+band-1 > n-1 {
            i_range = n-1
        }else{
            i_range = j+band-1
        }
        for i in j+1...i_range {   // for optimized, it should be i in j+1...j+(band-1) or j+1...n-1 if j+(band-1) > n-1
            A_a[i][j] = A_a[i][j]/A_a[j][j]
            RHS_a[i] = RHS_a[i] - A_a[i][j]*RHS_a[j]
            for k in j+1...i {
                A_a[i][k] = A_a[i][k] - A_a[i][j]*A_a[k][j]
            }
        }
    }
    return (A_a, RHS_a)
}

func cheleskySolver(A: [[Double]], b: [Double]) -> [Double]? {
    return nil
}

// func PSD_Check(A: [[Double]]) -> Bool {
//     return false
// }

func dotproduct_vector(A: [[Double]], b:[Double]) -> [Double]? {
    var output = Array(repeating: 0.0, count: b.count)
    for i in 0...b.count-1 {
        for j in 0...b.count-1 {
            output[i] += A[i][j] * b[j]
        }
    }
    return output
}

func transpose(A: [[Double]]) -> [[Double]]{
    var output = A
    for i in 0...A.count-1 {
        for j in 0...A.count-1 {
            output[j][i] = A[i][j]
        }
    }
    return output
}

func PSD_Generation(size: Int) -> [[Double]]? {
    return nil
}

func Read_Circuit_File(Dir: String) -> [(String, String, Double, Double, Double)] {
    var output: [String]
    var tmp: [String]
    var new_output: [(String, String, Double, Double, Double)] = []
    do {
        let file_read = try String(contentsOfFile: Dir)
        output = (file_read.split(separator: "\n").map{String($0)})
        for element in output {
            tmp = element.split(separator: " ").map{String($0)}
            new_output.append((tmp[0], tmp[1], Double(tmp[2])!, Double(tmp[3])!, Double(tmp[4])!))
        }
    } catch{
        return []
    }
    return new_output
}

let read_test = Read_Circuit_File(Dir: "test.txt")
print(read_test[0])