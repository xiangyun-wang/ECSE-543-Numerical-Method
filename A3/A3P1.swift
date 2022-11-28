import Foundation

/*
 * Input: Lower Triangular Matrix, RHS vector
 * Output: Matrix solution
 */

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

// element-wise subtraction for two matrices
func subtract(A: [Double], B: [Double]) -> [Double]{
    var output = A
    for i in 0...A.count-1{
        output[i] = A[i] - B[i]
    }
    return output
}

func addition(A: [Double], B: [Double]) -> [Double]{
    var output = A
    for i in 0...A.count-1{
        output[i] = A[i] + B[i]
    }
    return output
}

func multiplication_vector(a: [Double], c: Double) -> [Double]{
    var output = a
    for i in 0...a.count-1{
        output[i] = c * a[i]
    }
    return output
}

// dot product of matrix and a vector
func dotproduct_vector(A: [[Double]], b:[Double]) -> [Double] {
    var output = Array(repeating: 0.0, count: A.count)
    for i in 0...A.count-1 {
        for j in 0...b.count-1 {
            output[i] += A[i][j] * b[j]
        }
    }
    return output
}

// dot product with two matrices
func dotproduct_matrix(A: [[Double]], B: [[Double]]) -> [[Double]]{
    var output:[[Double]] = []
    let B_transpose = transpose(A: B)
    for i in 0...B_transpose.count-1{
        output.append(dotproduct_vector(A:A, b:B_transpose[i]))
    }
    return transpose(A: output)
}

// transpose of a matrix
func transpose(A: [[Double]]) -> [[Double]]{
    var output = Array(repeating: Array(repeating: 0.0, count: A.count), count: A[0].count)
    for i in 0...A.count-1 {
        for j in 0...A[0].count-1 {
            output[j][i] = A[i][j]
        }
    }
    return output
}

func transpose_vector(a: [Double]) -> [[Double]]{
    var output = Array(repeating: Array(repeating: 0.0, count: 1), count: a.count)
    for i in 0...a.count-1 {
        output[i][1] = a[i]
    }
    return output
}

func find_half_band(A: [[Double]]) -> Int {
    // initially set band = 0, 
    // loop over all rows
    // counter add until first non-zero element is met in a row (N_zero)
    // if i+1-N_zero > band, then band = i+1-N_zero
    // return band
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
}

func cholesky_decomposition(A: [[Double]], RHS: [Double]) -> ([[Double]]?, [Double]?) {
    let n = A.count
    var A_a = A
    var RHS_a = RHS

    for j in 0...n-1 {
        if A_a[j][j] <= 0 {     // not symmetric positive definite or singular
            return (nil, nil)
        }
        A_a[j][j] = sqrt(A_a[j][j])     // update diagonal element first
        RHS_a[j] = RHS_a[j]/A_a[j][j]   // update corresponding right hand size
        if j == n-1 {   // if it's the last column, only the diagonal element need to be updated, break directly
            break
        }

        for i in j+1...n-1 {   // update elements below current diagonal element (column)
            A_a[i][j] = A_a[i][j]/A_a[j][j]
            RHS_a[i] = RHS_a[i] - A_a[i][j]*RHS_a[j]
            for k in j+1...i {
                A_a[i][k] = A_a[i][k] - A_a[i][j]*A_a[k][j]
            }
        }
    }
    return (A_a, RHS_a)
}


func choleskySolver(A: [[Double]], b: [Double]) -> [Double]? {
    // decompose matrix
    let (L, y) = cholesky_decomposition(A: A, RHS: b)
    // check if decomposition is successful (not symmetric positive definite or singular)
    if L == nil {
        return nil
    }
    // the returned "y" above has already been done forward substitution
    // only backward substitution is needed
    return backward(Upper: transpose(A: L!), b: y!) 
}

func cholesky_decomposition_optimized_bandwidth_forward(A: [[Double]], RHS: [Double], band: Int) -> ([[Double]]?, [Double]?) {
    let n = A.count
    var A_a = A
    var RHS_a = RHS
    var i_range: Int

    for j in 0...n-1 {
        if A_a[j][j] <= 0 {
            return (nil, nil)
        }
        A_a[j][j] = sqrt(A_a[j][j])
        RHS_a[j] = RHS_a[j]/A_a[j][j]
        if j == n-1 {
            break
        }
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

func choleskySolver_Optimized(A: [[Double]], b: [Double]) -> [Double]? {
    let band = find_half_band(A: A) 
    let (L, y) = cholesky_decomposition_optimized_bandwidth_forward(A: A, RHS: b, band: band)
    if L == nil {
        print("Error: Matrix not SPD!!");
        return nil
    }
    return backward(Upper: transpose(A: L!), b: y!) 
}

func two_norm(A: [Double]) -> Double{
    var norm = 0.0
    for element in A{
        norm = norm + element*element
    }
    return sqrt(norm)
}

func infinity_norm(A: [Double]) -> Double{
    var max: Double = 0.0
    for element in A {
        if max < element * element{
            max = element * element
        }
    }
    return sqrt(max)
}

func Lagrange_L(current_order: Int, x_vector: [Double], x_value: Double) -> Double{
    var numerator_L = 1.0;
    var denominator_L = 1.0;
    for i in 0...x_vector.count-1{
        if i != current_order {
            numerator_L *= x_value - x_vector[i]
            denominator_L *= x_vector[current_order] - x_vector[i]
        }
    }
    return numerator_L/denominator_L; 
}

func Lagrange_interpolate(x_vector: [Double], y_vector: [Double], x_value: Double) -> Double{
    var y_value = 0.0
    for i in 0...x_vector.count-1{
        y_value += y_vector[i] * Lagrange_L(current_order: i, x_vector: x_vector, x_value: x_value)
    }
    return y_value
}

func Cubic_Hermite(){
    
}


