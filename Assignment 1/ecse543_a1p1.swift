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
            if (i == j) && (A[i][j] - extra_term < 0) {
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

// -----------------------------------------------------------------
func cholesky_decomposition_optimized(A: [[Double]], RHS: [Double]) -> ([[Double]]?, [Double]?) {
    let n = A.count
    var A_a = A
    var RHS_a = RHS

    for j in 0...n-1 {
        if A_a[j][j] <= 0 {
            return (nil, nil)
        }
        A_a[j][j] = sqrt(A_a[j][j])
        RHS_a[j] = RHS_a[j]/A_a[j][j]
        if j == n-1 {
            break
        }

        for i in j+1...n-1 {   // for optimized, it should be i in j+1...j+(band-1) or j+1...n-1 if j+(band-1) > n-1
            A_a[i][j] = A_a[i][j]/A_a[j][j]
            RHS_a[i] = RHS_a[i] - A_a[i][j]*RHS_a[j]
            for k in j+1...i {
                A_a[i][k] = A_a[i][k] - A_a[i][j]*A_a[k][j]
            }
        }
    }
    return (A_a, RHS_a)
}
// -----------------------------------------------------------------


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

func choleskySolver(A: [[Double]], b: [Double]) -> [Double]? {
    // let L = cholesky_decomposition(A: A)
    // if L == nil {
    //     return nil
    // }
    // // forward
    // let y = forward(Lower: L!, b: b)  
    // let solution = backward(Upper: transpose(A: L!), b: y!)
    // return solution
    
    let (L, y) = cholesky_decomposition_optimized(A: A, RHS: b)
    if L == nil {
        return nil
    }

    return backward(Upper: transpose(A: L!), b: y!) 
}

func choleskySolver_Optimized(A: [[Double]], b: [Double]) -> [Double]? {
    let band = find_half_band(A: A) 
    print("band: " + String(band))
    let (L, y) = cholesky_decomposition_optimized_bandwidth_forward(A: A, RHS: b, band: band)
    if L == nil {
        return nil
    }
    return backward(Upper: transpose(A: L!), b: y!) 
}

// func PSD_Check(A: [[Double]]) -> Bool {
//     return false
// }

func dotproduct_vector(A: [[Double]], b:[Double]) -> [Double] {
    var output = Array(repeating: 0.0, count: A.count)
    for i in 0...A.count-1 {
        for j in 0...b.count-1 {
            output[i] += A[i][j] * b[j]
        }
    }
    return output
}

// need to check
func dotproduct_matrix(A: [[Double]], B: [[Double]]) -> [[Double]]{
    var output:[[Double]] = []
    let B_transpose = transpose(A: B)
    for i in 0...B_transpose.count-1{
        output.append(dotproduct_vector(A:A, b:B_transpose[i]))
    }
    return transpose(A: output)
}

func transpose(A: [[Double]]) -> [[Double]]{
    var output = Array(repeating: Array(repeating: 0.0, count: A.count), count: A[0].count)
    for i in 0...A.count-1 {
        for j in 0...A[0].count-1 {
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

func Read_Matrix_File(Dir: String) -> [[Double]] {
    var output: [String]
    var tmp: [Double]
    var new_output: [[Double]] = []
    do {
        let file_read = try String(contentsOfFile: Dir)
        output = (file_read.split(separator: "\n").map{String($0)})
        for element in output {
            tmp = element.split(separator: " ").map{Double($0)!}
            new_output.append(tmp)
        }
    } catch{
        return []
    }
    return new_output
}

func print_vector(vector: [Double]){
    for i in 0...vector.count-1{
        print(vector[i])
    }
}

func Extract_Circuit(circuit: [(String, String, Double, Double, Double)]) -> ([[Double]],[Double],[Double]){
    var y = Array(repeating: Array(repeating: 0.0, count: circuit.count), count: circuit.count)
    var E = Array(repeating: 0.0, count: circuit.count)
    var J = Array(repeating: 0.0, count: circuit.count)
    for i in 0...circuit.count-1 {
        let (_,_,voltage,resistance,current) = circuit[i]
        E[i] = voltage
        if resistance < 0{
            y[i][i] = 0
        }else{
            y[i][i] = 1.0/resistance
        }
        J[i] = current
    }
    return (y,E,J)
}

func subtract(A: [Double], B: [Double]) -> [Double]{
    var output = A
    for i in 0...A.count-1{
        output[i] = A[i] - B[i]
    }
    return output
}

func Solve_Simple_Circuit(Matrix_Dir: String, Circuit_Dir: String) -> [Double]{
    
    let circuit = Read_Circuit_File(Dir: Circuit_Dir)
    //print("Circuit Read --- Passed")
    let reduced_A = Read_Matrix_File(Dir: Matrix_Dir)
    //print("Matrix Read --- Passed")
    let (y, E, J) = Extract_Circuit(circuit: circuit)    
    //print("Circuit Extract --- Passed")
    // let tmp1 = dotproduct_matrix(A: reduced_A, B: y)
    // print("temp1 --- Passed" + String(tmp1.count) + " " + String(tmp1[0].count))
    // let tmp2 = transpose(A: reduced_A)
    // print("temp2 --- Passed" + String(tmp2.count) + " " + String(tmp2[0].count))
    let A = dotproduct_matrix(A: dotproduct_matrix(A: reduced_A, B: y), B: transpose(A: reduced_A))
    //print("A --- Passed")
    let b = dotproduct_vector(A: reduced_A,b: subtract(A: J,B: dotproduct_vector(A: y,b: E)))
    //print("b --- Passed")
    //let output = choleskySolver_Optimized(A: A, b: b)
    let output = choleskySolver(A: A, b: b)
    return output!
}


func Find_Eqv_Resistance(N: Int, R: Double) -> [[Double]] {
    let N_R = N*(N-1)*2
    var A = Array(repeating: Array(repeating: 0.0, count: (N_R+1)), count: (N*N))
    var y = Array(repeating: Array(repeating: 0.0, count: N_R+1), count: N_R+1)
    let J = Array(repeating: 0.0, count: N_R+1)
    var E = Array(repeating: 0.0, count: N_R+1)
    for i in 0...y.count-1 {
        y[i][i] = 1.0/R
    }
    E[N_R] = -100
    var column_to_write = 0
    for i in 1...N*N-1 {
        if (i < N-1 && i > 0) {
            // lower boarder
            A[i][column_to_write] = -1
            A[i][column_to_write+1] = 1
            A[i][column_to_write+N] = 1
            column_to_write += 1
        } else if ( i == N-1 || i == N*(N-1)){
            // corner
            A[i][column_to_write] = -1
            A[i][column_to_write+N] = 1
            column_to_write += 1
        }else if (i == N*N-1){
            // top right corner
            A[i][column_to_write] = -1
            A[i][column_to_write+N-1] = -1
            A[i][column_to_write+N] = 1
        }  else if (i%N==0){
            A[i][column_to_write] = -1
            A[i][column_to_write+N] = 1
            A[i][column_to_write+2*N-1] = 1
            // left boarder
            column_to_write += 1
        } else if (i%N==N-1){
            A[i][column_to_write] = -1
            A[i][column_to_write+N-1] = -1
            A[i][column_to_write+2*N-1] = 1
            // right boarder
            column_to_write += N
        } else if (i>N*(N-1) && i < N*N-1){
            // top boarder
            A[i][column_to_write] = -1
            A[i][column_to_write+N-1] = -1
            A[i][column_to_write+N] = 1
            column_to_write += 1
        }else{
            // middle
            A[i][column_to_write] = -1
            A[i][column_to_write+N-1] = -1
            A[i][column_to_write+N] = 1
            A[i][column_to_write+2*N-1] = 1
            column_to_write += 1
        }
    }
    A.remove(at: 0)

    let A_a = dotproduct_matrix(A: dotproduct_matrix(A: A, B: y), B: transpose(A: A))
    //print("A --- Passed")
    let b_b = dotproduct_vector(A: A,b: subtract(A: J,B: dotproduct_vector(A: y,b: E)))
    //print("b --- Passed")
    //let output = choleskySolver_Optimized(A: A, b: b)
    //let output = choleskySolver(A: A_a, b: b_b)
    let start = DispatchTime.now()
    //let output = choleskySolver(A: A_a, b: b_b)
    let output = choleskySolver_Optimized(A: A_a, b: b_b)
    let end = DispatchTime.now()
    let nanoTime = end.uptimeNanoseconds - start.uptimeNanoseconds
    print(Double(nanoTime)/1000000) // milisecon
    //choleskySolver_Optimized
    //let Req = R*output![N*N-2]/(1-output![N*N-2])
    let Req = output![output!.count-1]*R/(100-output![output!.count-1])
    return [[Req]]
    // A completed
}

// let test3 = dotproduct_matrix(A: dotproduct_matrix(A: test1, B: test2), B: transpose(A: test1))

// print(find_half_band(A: test3))

// var solution = Solve_Simple_Circuit(Matrix_Dir: "test_5_A_reduced.txt", Circuit_Dir: "test_5_circuit.txt")
// print(solution)
// .removeat()


let A = [[6.0, 0, 0],[0, 55, 225],[0, 225, 979]]

let b = [1.0,2.0,3.0]


// var solution = choleskySolver(A: A, b: b)
// print_vector(vector: solution!)
// var solution1 = choleskySolver_Optimized(A:A, b:b)
// print_vector(vector: solution1!)

// let read_test = Read_Circuit_File(Dir: "test.txt")
// print(read_test[0])
let z = Find_Eqv_Resistance(N: 15, R: 100.0)
print(z)
// for element in z {
//     print(element)
// }
//print(Find_Eqv_Resistance(N: 20, R: 100.0))