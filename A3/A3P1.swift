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

// ================== A3 Q1 ab =============================

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
// =============== A3 Q1 c ==============================
func Cubic_Hermite_interpolate(x_vector: [Double], y_vector: [Double], x_value: Double) -> Double{
    // find subdomain
    var sub_domain_index = 0;
    while(true){
        if x_value == x_vector[sub_domain_index]{
            return y_vector[sub_domain_index]
        }
        if x_value < x_vector[sub_domain_index] {
            break;
        }else{
            sub_domain_index += 1
            if sub_domain_index > 5{
                sub_domain_index = 5
                break
            }
        }
    }
    let L1 = (x_value-x_vector[sub_domain_index])/(x_vector[sub_domain_index-1] - x_vector[sub_domain_index])
    let L2 = (x_value-x_vector[sub_domain_index-1])/(x_vector[sub_domain_index] - x_vector[sub_domain_index-1])
    let L1_d = 1/(x_vector[sub_domain_index-1] - x_vector[sub_domain_index])
    let L2_d = 1/(x_vector[sub_domain_index] - x_vector[sub_domain_index-1])
    
    let U1 = (1-2*L1_d*(x_value-x_vector[sub_domain_index-1])) * (L1 * L1)
    let U2 = (1-2*L2_d*(x_value-x_vector[sub_domain_index])) * (L2 * L2)
    let V1 = (x_value-x_vector[sub_domain_index-1]) * (L1 * L1)
    let V2 = (x_value-x_vector[sub_domain_index]) * (L2 * L2)

    var b1 = 0.0
    if sub_domain_index != 1 {
        b1 = (y_vector[sub_domain_index-1]-y_vector[sub_domain_index-2]) / (x_vector[sub_domain_index-1]-x_vector[sub_domain_index-2])
    }
    let b2 = (y_vector[sub_domain_index]-y_vector[sub_domain_index-1]) / (x_vector[sub_domain_index]-x_vector[sub_domain_index-1])

    return y_vector[sub_domain_index-1] * U1 + b1 * V1 + y_vector[sub_domain_index] * U2 + b2 * V2
}

// =============== A3 Q1 d e f =====================

func Newton_Raphson(guess: Double, threshold: Double, X: [Double], Y: [Double]) -> (Double,Int){
    var output = guess
    let f_0 = Evaluate_Nonlinear_Function(input: 0.0, X: X, Y: Y)
    var iteration = 0
    while (abs(Evaluate_Nonlinear_Function(input: output, X: X, Y: Y)/f_0) > threshold){
        output -= Evaluate_Nonlinear_Function(input: output, X: X, Y: Y)/Nonlinear_Derivative(input: output, X: X, Y: Y)
        iteration += 1
    }
    return (output,iteration)
}

func Evaluate_Nonlinear_Function(input: Double, X: [Double], Y: [Double]) -> Double{
    return 3.97887 * 1e7 * input + 0.3 * Piecewise(input: input, X: X, Y: Y) - 8000
}

func Evaluate_Nonlinear_Function_Derivative(input: Double, X: [Double], Y: [Double]) -> Double{
    return 3.97887 * 1e7 + 0.3 * Nonlinear_Derivative(input: input, X: X, Y: Y)
}

func Nonlinear_Derivative(input: Double, X: [Double], Y: [Double]) -> Double{
    var index = X.count-1
    for i in 1...X.count-1 {
        if input < X[i]{
            index = i
            break
        }
    }
    return (Y[index]-Y[index-1])/(X[index]-X[index-1])
}

func Piecewise(input: Double, X: [Double], Y: [Double]) -> Double{
    var index = X.count-1
    for i in 1...X.count-1 {
        if input < X[i]{
            index = i-1
            break
        }
    }
    let slope = Nonlinear_Derivative(input: input, X: X, Y: Y)
    return Y[index] + (input-X[index]) * slope    
}

func Successive_Sub(guess: Double, threshold: Double, X: [Double], Y: [Double]) -> (Double,Int){
    let f_0 = Evaluate_Nonlinear_Function(input: 0.0, X: X, Y: Y)
    var output = guess
    var iteration = 0
    while(abs(Evaluate_Nonlinear_Function(input: output, X: X, Y: Y)/f_0) > threshold){
        output = output - Evaluate_Nonlinear_Function(input: output, X: X, Y: Y)/(3.97887*1e7) * 1e-3 // smaller step size
        iteration += 1
    }
    return (output,iteration)
}

func Q2_Jaccobian(X: [Double]) -> [[Double]]{
    var output = Array(repeating: Array(repeating: 0.0, count: 2), count: 2)
    output[0][0] = 1/512 + (1/(25e-3))*0.8e-6 * exp((X[0]-X[1])/(25e-3))
    output[0][1] = (-1)*(1/(25e-3))*0.8e-6 * exp((X[0]-X[1])/(25e-3))
    output[1][0] = (-1) * (1/(25e-3))*0.8e-6 * exp((X[0]-X[1])/(25e-3))
    output[1][1] = (1/(25e-3))*1.1e-6 * exp((X[1])/(25e-3)) + (1/(25e-3))*0.8e-6 * exp((X[0]-X[1])/(25e-3))
    return output
}

func Q2_F_Eva(X: [Double]) -> [Double]{
    var output = [0.0, 0.0]
    output[0] = (X[0]-0.2)/512 + 0.8e-6 * (exp((X[0]-X[1])/(25e-3))-1)
    output[1] = 1.1e-6 * (exp((X[1])/(25e-3))-1) - 0.8e-6 * (exp((X[0]-X[1])/(25e-3))-1)

    return output
}

func Q2_NR(guess: [Double], threshold: Double) -> ([Double],[Double],[[Double]],[[Double]]){
    var X = guess
    var f_value = [0.0,0.0]
    var f_d = [[0.0,0.0],[0.0,0.0]]
    var delta_x = [0.0,0.0]
    var f_set = [[Double]]()
    var x_set = [[Double]]()
    var delta_set = [Double]()
    while(true){
        f_value = Q2_F_Eva(X: X)
        f_set.append(f_value)
        f_d = Q2_Jaccobian(X: X)
        delta_x = choleskySolver(A: f_d,b: f_value)!
        X = subtract(A: X, B: delta_x)
        x_set.append(X)
        delta_set.append(two_norm(A: delta_x))
        if two_norm(A: delta_x) < threshold{
            return (X,delta_set, f_set, x_set)
        }
    }
}



// =============== A3 Q2 ===========================


// ============= A3 Q1 abc testing =================

var first_six_x = [0.0,0.2,0.4,0.6,0.8,1.0];
var first_six_y = [0.0,14.7,36.5,71.7,121.4,197.4]

var six_x = [0.0, 1.3, 1.4, 1.7, 1.8, 1.9]
var six_y = [0.0, 540.6, 1062.8, 8687.4, 13924.3, 22650.2]

var continuous_x = [0.00]
for i in 1...190 {
    continuous_x.append(0.01 * Double(i))
}

var y_Lagrange_first_six = [Double]()
var y_Lagrange_six = [Double]()
var y_Cubic_Hermite_six = [Double]()
var counter = 0
for element in continuous_x{
    y_Lagrange_first_six.append(Lagrange_interpolate(x_vector: first_six_x, y_vector: first_six_y, x_value: element))
    y_Lagrange_six.append(Lagrange_interpolate(x_vector: six_x, y_vector: six_y, x_value: element))
    y_Cubic_Hermite_six.append(Cubic_Hermite_interpolate(x_vector: six_x, y_vector: six_y, x_value: element))
}

print(y_Lagrange_first_six)
print(continuous_x)
print(y_Lagrange_six)
print(y_Cubic_Hermite_six)

// ============== A3 Q1 def testing ==============
let B = [0.0e-4,0.2e-4,0.4e-4,0.6e-4,0.8e-4,1.0e-4,1.1e-4,1.2e-4,1.3e-4,1.4e-4,1.5e-4,1.6e-4,1.7e-4,1.8e-4,1.9e-4]
let H = [0.0,14.7,36.5,71.7,121.4,197.4,256.2,348.7,540.6,1062.8,2318.0,4781.9,8687.4,13924.3,22650.2]
var (answer,iteration) = Newton_Raphson(guess: 0.0, threshold: 1e-6, X: B, Y: H)
print("Newton Raphson: ")
print(answer)
print("Iteration: ")
print(iteration)
print("Successive Substitution: ")
var (answer_2,iteration_2) = Successive_Sub(guess: 0.0, threshold: 1e-6, X: B, Y: H)
print(answer_2)



// ============== A3 Q2 testing ==================
var (answer,delta,f_values, x_values) = Q2_NR(guess: [0.0,0.0],threshold: 1e-6)
print("Final Converged Node Voltages: [V1, V2]")
print(answer)
print("Delta X: ")
print(delta)
print("f(Vn) Values: ([f1(Vn), f2(Vn)])")
print(f_values)
print("Vn Values: ([V1, V2])")
print(x_values)