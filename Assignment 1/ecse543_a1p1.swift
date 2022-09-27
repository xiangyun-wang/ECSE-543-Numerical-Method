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

func Simple_Circuit_A_b(Matrix_Dir: String, Circuit_Dir: String) -> ([[Double]], [Double]){
    
    let circuit = Read_Circuit_File(Dir: Circuit_Dir)
    let reduced_A = Read_Matrix_File(Dir: Matrix_Dir)
    let (y, E, J) = Extract_Circuit(circuit: circuit)    
    let A = dotproduct_matrix(A: dotproduct_matrix(A: reduced_A, B: y), B: transpose(A: reduced_A))
    let b = dotproduct_vector(A: reduced_A,b: subtract(A: J,B: dotproduct_vector(A: y,b: E)))
    //let output = choleskySolver(A: A, b: b)
    return (A, b)
}

func Solve_Circuit(Circuit: ([[Double]], [Double]), Optimize: Bool) -> [Double]{
    let (A, b) = Circuit
    if Optimize {
        return choleskySolver_Optimized(A: A, b: b)!
    } else {
        return choleskySolver(A: A, b: b)!
    }
}


func Mesh_Circuit_Generation(N: Int, R: Double, V: Double) -> ([[Double]], [Double]) {
    let N_R = N*(N-1)*2
    var A = Array(repeating: Array(repeating: 0.0, count: (N_R+1)), count: (N*N))
    var y = Array(repeating: Array(repeating: 0.0, count: N_R+1), count: N_R+1)
    let J = Array(repeating: 0.0, count: N_R+1)
    var E = Array(repeating: 0.0, count: N_R+1)
    for i in 0...y.count-1 {
        y[i][i] = 1.0/R
    }
    E[N_R] = -V
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

    return (A_a, b_b)
}

func find_Req(N: Int, R: Double, V: Double, Optimize: Bool) -> (Double, Double){
    let mesh_circuit = Mesh_Circuit_Generation(N: N, R: R, V: V)
    let start = DispatchTime.now()
    let v_nodes = Solve_Circuit(Circuit: mesh_circuit, Optimize: Optimize)
    let end = DispatchTime.now()
    //print(Double(end.uptimeNanoseconds - start.uptimeNanoseconds)/1000000) // milisecon
    let Req = v_nodes[v_nodes.count-1]*R/(V-v_nodes[v_nodes.count-1])
    return (Req, Double(end.uptimeNanoseconds - start.uptimeNanoseconds)/1000000)
}

// let test3 = dotproduct_matrix(A: dotproduct_matrix(A: test1, B: test2), B: transpose(A: test1))

// print(find_half_band(A: test3))

// var solution = Solve_Circuit(Circuit: Simple_Circuit_A_b(Matrix_Dir: "test_5_A_reduced.txt", Circuit_Dir: "test_5_circuit.txt"), Optimize: false)
// print(solution)

// let (Req, time_consumed) = find_Req(N: 10, R: 100.0, V: 100, Optimize: false)
// print(Req)
// print(time_consumed)


//-------------------------- part 3 start ---------------------------
func find_residue(matrix : [[Double]]) -> Double {  // need to change
    var residue = 0.0;
    var possible_residue = 0.0;
    for i in 1...matrix.count-2 {
        for j in 1...matrix.count-2 {
            possible_residue = -4 * matrix[i][j] + matrix[i][j-1] + matrix[i][j+1] + matrix[i-1][j] + matrix[i+1][j]
            if possible_residue > residue {
                residue = possible_residue
            }
        }
    }
    return residue
}

func finite_difference_SOR(h: Double, w: Double, threshold: Double) -> [[Double]]{
    let inner_voltage = 110.0
    let outer_voltage = 0.0
    let inner_x = 0.06/h
    let inner_y = 0.08/h
    let free_top = 0.1
    let free_right = 0.1
    let outer_x = 0
    let outer_y = 0
    let x_range_discret = 0.1/h-1
    let y_range_discret = 0.1/h-1

    var grid_potential = Array(repeating: Array(repeating: 0.0, count: y_range_discret+1), count: x_range_discret+1)
    var iteration = 0
    var residue = Double.greatestFiniteMagnitude
    var tmp = 0.0

    for j in inner_y...grid_potential[0].count-1 {
        grid_potential[inner_x][j] = 110.0
    }

    for i in inner_x...grid_potential.count-1 {
        grid_potential[i][inner_y] = 110.0
    }

    while (residue > threshold) {
        for i in 1...x_range_discret{
            for j in 1...y_range_discret{
                if !(i >=  inner_x && j >= inner_y) {
                    if i == x_range_discret {       // right boundary
                        tmp = 2 * grid_potential[i-1][j] + grid_potential[i][j-1] + grid_potential[i][j+1] 
                    } else if j == y_range_discret {    // top boundary
                        tmp = 2 * grid_potential[i][j-1] + grid_potential[i-1][j] + grid_potential[i+1][j] 
                    } else {
                        tmp = grid_potential[i][j+1] + grid_potential[i][j-1] + grid_potential[i-1][j] + grid_potential[i+1][j] 
                    }
                    grid_potential[i][j] = (1-w) * grid_potential[i][j] + w/4.0*tmp
                }
            }
        }

        residue  = find_residue(grid_potential)
        iteration += iteration
    }
    return grid_potential
}
//-------------------------- part 3 end ---------------------------