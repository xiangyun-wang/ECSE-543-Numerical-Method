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

func CG_solve(A:[[Double]],b:[Double], threshold: Double) -> ([Double],[Double],[Double]){
    var x = Array(repeating: 0.0, count: b.count)
    var r = subtract(A: b, B: dotproduct_vector(A: A, b:x))
    var p = r
    var two_norm_array: [Double] = []
    var infinity_norm_array: [Double] = []
    var infinity_error = infinity_norm(A: r)
    var error = two_norm(A: r)
    two_norm_array.append(error)
    infinity_norm_array.append(infinity_error)
    var alpha = 0.0
    var beta = 0.0
    var pt = [p]
    while error > threshold {
        pt = [p]
        alpha = dotproduct_vector(A: pt, b: r)[0] / dotproduct_vector(A: dotproduct_matrix(A: pt, B: A), b: p)[0]
        x = addition(A: multiplication_vector(a: p, c: alpha), B: x)
        r = subtract(A: b, B: dotproduct_vector(A: A, b:x))
        beta = (-1) * dotproduct_vector(A: dotproduct_matrix(A: pt, B: A), b: r)[0] / 
            dotproduct_vector(A: dotproduct_matrix(A: pt, B: A), b: p)[0] 
        p = addition(A: multiplication_vector(a: p, c: beta), B: r)
        error = two_norm(A: r)
        infinity_error = infinity_norm(A: r)
        two_norm_array.append(error)
        infinity_norm_array.append(infinity_error)

    }
    return (x,two_norm_array,infinity_norm_array)
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

var grid = Array(repeating: Array(repeating: 0.0, count: 19), count: 19)
grid[0][0] = -4
grid[0][1] = 1
grid[0][2] = 2
grid[1][1] = -4
grid[1][0] = 1
grid[1][3] = 2
grid[2][2] = -4
grid[2][0] = 1
grid[2][4] = 1
grid[2][3] = 1
grid[3][3] = -4
grid[3][1] = 1
grid[3][2] = 1
grid[3][5] = 1
grid[4][4] = -4
grid[4][2] = 1
grid[4][5] = 1
grid[4][9] = 1
grid[5][5] = -4
grid[5][3] = 1
grid[5][4] = 1
grid[5][6] = 1
grid[5][10] = 1
grid[6][6] = -4
grid[6][5] = 1
grid[6][7] = 1
grid[6][11] = 1
grid[7][7] = -4
grid[7][6] = 1
grid[7][8] = 1
grid[7][12] = 1
grid[8][8] = -4
grid[8][7] = 2
grid[8][13] = 1
grid[9][9] = -4
grid[9][4] = 1
grid[9][10] = 1
grid[9][14] = 1
grid[10][10] = -4
grid[10][5] = 1
grid[10][9] = 1
grid[10][11] = 1
grid[10][15] = 1
grid[11][11] = -4
grid[11][6] = 1
grid[11][10] = 1
grid[11][12] = 1
grid[11][16] = 1
grid[12][12] = -4
grid[12][7] = 1
grid[12][11] = 1
grid[12][13] = 1
grid[12][17] = 1
grid[13][13] = -4
grid[13][8] = 1
grid[13][12] = 2
grid[13][18] = 1
grid[14][14] = -4
grid[14][9] = 1
grid[14][15] = 1
grid[15][15] = -4
grid[15][10] = 1
grid[15][14] = 1
grid[15][16] = 1
grid[16][16] = -4
grid[16][11] = 1
grid[16][15] = 1
grid[16][17] = 1
grid[17][17] = -4
grid[17][12] = 1
grid[17][16] = 1
grid[17][18] = 1
grid[18][18] = -4
grid[18][13] = 1
grid[18][17] = 2

var b = Array(repeating: 0.0, count: 19)

b[1] = -110
b[3] = -110
b[6] = -110
b[7] = -110
b[8] = -110 

var test1 = choleskySolver_Optimized(A: grid, b: b)
//print(test1)

// cholesky
var cholesky_ans = choleskySolver_Optimized(A: dotproduct_matrix(A: transpose(A: grid), B: grid), 
    b: dotproduct_vector(A: transpose(A: grid), b: b))

print("Cholesky Decomposition: ")
print(((cholesky_ans!)))
// Conjugate Gradient

var (CG_ans_x, CG_ans_two_norm, CG_ans_infinity_norm) = CG_solve(A:dotproduct_matrix(A: transpose(A: grid), B: grid),b:dotproduct_vector(A: transpose(A: grid), b: b), threshold: 1e-6)
print("Conjugate Gradient: ")
print(CG_ans_x)
print(CG_ans_two_norm)
print(CG_ans_infinity_norm)
print(grid)
print(b)