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

var S = Array(repeating: Array(repeating: 0.0, count: 34), count: 34)
var potential = Array(repeating: Array(repeating: 0.0, count: 1), count: 34)

S[0][0] = 1.0000
S[1][0] = -0.5000
S[6][0] = -0.5000
S[0][1] = -0.5000
S[1][1] = 2.0000
S[2][1] = -0.5000
S[7][1] = -1.0000
S[1][2] = -0.5000
S[2][2] = 2.0000
S[3][2] = -0.5000
S[8][2] = -1.0000
S[2][3] = -0.5000
S[3][3] = 2.0000
S[4][3] = -0.5000
S[9][3] = -1.0000
S[3][4] = -0.5000
S[4][4] = 2.0000
S[5][4] = -0.5000
S[10][4] = -1.0000
S[4][5] = -0.5000
S[5][5] = 1.0000
S[11][5] = -0.5000
S[0][6] = -0.5000
S[6][6] = 2.0000
S[7][6] = -1.0000
S[12][6] = -0.5000
S[1][7] = -1.0000
S[6][7] = -1.0000
S[7][7] = 4.0000
S[8][7] = -1.0000
S[13][7] = -1.0000
S[2][8] = -1.0000
S[7][8] = -1.0000
S[8][8] = 4.0000
S[9][8] = -1.0000
S[14][8] = -1.0000
S[3][9] = -1.0000
S[8][9] = -1.0000
S[9][9] = 4.0000
S[10][9] = -1.0000
S[15][9] = -1.0000
S[4][10] = -1.0000
S[9][10] = -1.0000
S[10][10] = 4.0000
S[11][10] = -1.0000
S[16][10] = -1.0000
S[5][11] = -0.5000
S[10][11] = -1.0000
S[11][11] = 2.5000
S[17][11] = -1.0000
S[6][12] = -0.5000
S[12][12] = 2.0000
S[13][12] = -1.0000
S[18][12] = -0.5000
S[7][13]  = -1.0000
S[12][13] = -1.0000
S[13][13] = 4.0000
S[14][13] = -1.0000
S[19][13] = -1.0000
S[8][14] = -1.0000
S[13][14] = -1.0000
S[14][14] = 4.0000
S[15][14] = -1.0000
S[20][14] = -1.0000
S[9][15] = -1.0000
S[14][15] = -1.0000
S[15][15] = 4.0000
S[16][15] = -1.0000
S[21][15] = -1.0000
S[10][16] = -1.0000
S[15][16] = -1.0000
S[16][16] = 3.5000
S[17][16] = -1.0000
S[22][16] = -0.5000
S[11][17] = -1.0000
S[16][17] = -1.0000
S[17][17] = 2.5000
S[23][17] = -0.5000
S[12][18] = -0.5000
S[18][18] = 2.0000
S[19][18] = -1.0000
S[24][18] = -0.5000
S[13][19] = -1.0000
S[18][19] = -1.0000
S[19][19] = 4.0000
S[20][19] = -1.0000
S[25][19] = -1.0000
S[14][20] = -1.0000
S[19][20] = -1.0000
S[20][20] = 4.0000
S[21][20] = -1.0000
S[26][20] = -1.0000
S[15][21] = -1.0000
S[20][21] = -1.0000
S[21][21] = 4.0000
S[22][21] = -1.0000
S[27][21] = -1.0000
S[16][22] = -0.5000
S[21][22] = -1.0000
S[22][22] = 3.5000
S[23][22] = -1.0000
S[28][22] = -1.0000
S[17][23] = -0.5000
S[22][23] = -1.0000
S[23][23] = 2.0000
S[29][23] = -0.5000
S[18][24] = -0.5000
S[24][24] = 2.0000
S[25][24] = -1.0000
S[30][24] = -0.5000
S[19][25] = -1.0000
S[24][25] = -1.0000
S[25][25] = 4.0000
S[26][25] = -1.0000
S[31][25] = -1.0000
S[20][26] = -1.0000
S[25][26] = -1.0000
S[26][26] = 4.0000
S[27][26] = -1.0000
S[32][26] = -1.0000
S[21][27] = -1.0000
S[26][27] = -1.0000
S[27][27] = 3.0000
S[28][27] = -0.5000
S[33][27] = -0.5000
S[22][28] = -1.0000
S[27][28] = -0.5000
S[28][28] = 2.0000
S[29][28] = -0.5000
S[23][29] = -0.5000
S[28][29] = -0.5000
S[29][29] = 1.0000
S[24][30] = -0.5000
S[30][30] = 1.0000
S[31][30] = -0.5000
S[25][31] = -1.0000
S[30][31] = -0.5000
S[31][31] = 2.0000
S[32][31] = -0.5000
S[26][32] = -1.0000
S[31][32] = -0.5000
S[32][32] = 2.0000
S[33][32] = -0.5000
S[27][33] = -0.5000
S[32][33] = -0.5000
S[33][33] = 1.0000

potential[0][0] = 0
potential[1][0] = 0
potential[2][0] = 0
potential[3][0] = 0
potential[4][0] = 0
potential[5][0] = 0
potential[6][0] = 0
potential[7][0] = 6.8275
potential[8][0] = 13.1569
potential[9][0] = 18.0785
potential[10][0] = 20.7092
potential[11][0] = 24.8164
potential[12][0] = 0
potential[13][0] = 14.1532
potential[14][0] = 27.7217
potential[15][0] = 38.4477
potential[16][0] = 39.9421
potential[17][0] = 41.3318
potential[18][0] = 0
potential[19][0] = 22.0635
potential[20][0] = 45.1289
potential[21][0] = 68.0485
potential[22][0] = 78.6175
potential[23][0] = 77.1417
potential[24][0] = 0
potential[25][0] = 28.9720
potential[26][0] = 62.6819
potential[27][0] = 110.0000
potential[28][0] = 110.0000
potential[29][0] = 110.0000
potential[30][0] = 0
potential[31][0] = 31.1427
potential[32][0] = 66.6266
potential[33][0] = 110.0000
let epsilon = 8.854e-12
let total_energy = 0.5 * dotproduct_matrix(A: dotproduct_matrix(A: transpose(A: potential), B: S), B: potential)[0][0];
let V = 110.0
let C = 2.0*total_energy*epsilon/(V*V) * 4
print(C);