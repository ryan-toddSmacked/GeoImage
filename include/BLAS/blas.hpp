#pragma once

#include <array>
#include <cmath>

namespace blas {

    //--------------------------------------------------------------------------
    // Fixed-size matrix types for transformation calculations
    //--------------------------------------------------------------------------
    
    template<size_t Rows, size_t Cols>
    using Matrix = std::array<std::array<double, Cols>, Rows>;
    
    template<size_t N>
    using Vector = std::array<double, N>;

    //--------------------------------------------------------------------------
    // Basic matrix operations
    //--------------------------------------------------------------------------
    
    // Matrix-vector multiplication: y = A * x
    template<size_t Rows, size_t Cols>
    inline Vector<Rows> matVecMul(const Matrix<Rows, Cols>& A, const Vector<Cols>& x) {
        Vector<Rows> y{};
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                y[i] += A[i][j] * x[j];
            }
        }
        return y;
    }
    
    // Matrix transpose: B = A^T
    template<size_t Rows, size_t Cols>
    inline Matrix<Cols, Rows> transpose(const Matrix<Rows, Cols>& A) {
        Matrix<Cols, Rows> B{};
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                B[j][i] = A[i][j];
            }
        }
        return B;
    }
    
    // Matrix-matrix multiplication: C = A * B
    template<size_t M, size_t N, size_t P>
    inline Matrix<M, P> matMul(const Matrix<M, N>& A, const Matrix<N, P>& B) {
        Matrix<M, P> C{};
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < P; ++j) {
                for (size_t k = 0; k < N; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    }
    
    // A^T * b for matrix A and vector b
    template<size_t Rows, size_t Cols>
    inline Vector<Cols> matTransposeVecMul(const Matrix<Rows, Cols>& A, const Vector<Rows>& b) {
        Vector<Cols> y{};
        for (size_t j = 0; j < Cols; ++j) {
            for (size_t i = 0; i < Rows; ++i) {
                y[j] += A[i][j] * b[i];
            }
        }
        return y;
    }

    //--------------------------------------------------------------------------
    // 2x2 matrix inverse (for completeness)
    //--------------------------------------------------------------------------
    
    inline bool inverse2x2(const Matrix<2, 2>& A, Matrix<2, 2>& Ainv) {
        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if (std::abs(det) < 1e-12) return false;
        double invDet = 1.0 / det;
        Ainv[0][0] =  A[1][1] * invDet;
        Ainv[0][1] = -A[0][1] * invDet;
        Ainv[1][0] = -A[1][0] * invDet;
        Ainv[1][1] =  A[0][0] * invDet;
        return true;
    }

    //--------------------------------------------------------------------------
    // 3x3 matrix inverse
    //--------------------------------------------------------------------------
    
    inline bool inverse3x3(const Matrix<3, 3>& A, Matrix<3, 3>& Ainv) {
        // Calculate determinant
        double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                   - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                   + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
        
        if (std::abs(det) < 1e-12) return false;
        double invDet = 1.0 / det;
        
        // Adjugate matrix / determinant
        Ainv[0][0] =  (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * invDet;
        Ainv[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * invDet;
        Ainv[0][2] =  (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * invDet;
        Ainv[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) * invDet;
        Ainv[1][1] =  (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * invDet;
        Ainv[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) * invDet;
        Ainv[2][0] =  (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * invDet;
        Ainv[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) * invDet;
        Ainv[2][2] =  (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * invDet;
        
        return true;
    }

    //--------------------------------------------------------------------------
    // 6x6 matrix inverse (for exact 3-point affine solution)
    //--------------------------------------------------------------------------
    
    inline bool inverse6x6(const Matrix<6, 6>& A, Matrix<6, 6>& Ainv) {
        // Gauss-Jordan elimination with partial pivoting
        Matrix<6, 6> aug{};
        
        // Initialize augmented matrix [A | I]
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                aug[i][j] = A[i][j];
            }
        }
        
        // Identity portion tracked separately
        Matrix<6, 6> inv{};
        for (size_t i = 0; i < 6; ++i) {
            inv[i][i] = 1.0;
        }
        
        for (size_t col = 0; col < 6; ++col) {
            // Find pivot
            size_t maxRow = col;
            double maxVal = std::abs(aug[col][col]);
            for (size_t row = col + 1; row < 6; ++row) {
                if (std::abs(aug[row][col]) > maxVal) {
                    maxVal = std::abs(aug[row][col]);
                    maxRow = row;
                }
            }
            
            if (maxVal < 1e-12) return false;  // Singular matrix
            
            // Swap rows
            if (maxRow != col) {
                std::swap(aug[col], aug[maxRow]);
                std::swap(inv[col], inv[maxRow]);
            }
            
            // Scale pivot row
            double pivot = aug[col][col];
            for (size_t j = 0; j < 6; ++j) {
                aug[col][j] /= pivot;
                inv[col][j] /= pivot;
            }
            
            // Eliminate column
            for (size_t row = 0; row < 6; ++row) {
                if (row != col) {
                    double factor = aug[row][col];
                    for (size_t j = 0; j < 6; ++j) {
                        aug[row][j] -= factor * aug[col][j];
                        inv[row][j] -= factor * inv[col][j];
                    }
                }
            }
        }
        
        Ainv = inv;
        return true;
    }

    //--------------------------------------------------------------------------
    // Solve 6x6 linear system Ax = b
    //--------------------------------------------------------------------------
    
    inline bool solve6x6(const Matrix<6, 6>& A, const Vector<6>& b, Vector<6>& x) {
        Matrix<6, 6> Ainv;
        if (!inverse6x6(A, Ainv)) return false;
        x = matVecMul(Ainv, b);
        return true;
    }

    //--------------------------------------------------------------------------
    // Least squares solution for overdetermined systems
    // Solves A*x = b in least squares sense: x = (A^T * A)^-1 * A^T * b
    //--------------------------------------------------------------------------
    
    // 4 points -> 8 equations, 6 unknowns (A is 8x6)
    inline bool leastSquares8x6(const Matrix<8, 6>& A, const Vector<8>& b, Vector<6>& x) {
        // Compute A^T * A (6x6)
        Matrix<6, 6> AtA{};
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                for (size_t k = 0; k < 8; ++k) {
                    AtA[i][j] += A[k][i] * A[k][j];
                }
            }
        }
        
        // Compute A^T * b (6x1)
        Vector<6> Atb{};
        for (size_t i = 0; i < 6; ++i) {
            for (size_t k = 0; k < 8; ++k) {
                Atb[i] += A[k][i] * b[k];
            }
        }
        
        // Solve (A^T * A) * x = A^T * b
        return solve6x6(AtA, Atb, x);
    }
    
    // 5 points -> 10 equations, 6 unknowns (A is 10x6)
    inline bool leastSquares10x6(const Matrix<10, 6>& A, const Vector<10>& b, Vector<6>& x) {
        // Compute A^T * A (6x6)
        Matrix<6, 6> AtA{};
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                for (size_t k = 0; k < 10; ++k) {
                    AtA[i][j] += A[k][i] * A[k][j];
                }
            }
        }
        
        // Compute A^T * b (6x1)
        Vector<6> Atb{};
        for (size_t i = 0; i < 6; ++i) {
            for (size_t k = 0; k < 10; ++k) {
                Atb[i] += A[k][i] * b[k];
            }
        }
        
        // Solve (A^T * A) * x = A^T * b
        return solve6x6(AtA, Atb, x);
    }

    //--------------------------------------------------------------------------
    // Affine transformation from point correspondences
    // Given N point pairs (px, py) -> (gx, gy), compute transformation matrix
    //
    // Affine transform: gx = a*px + b*py + c
    //                   gy = d*px + e*py + f
    //
    // Returns 6 coefficients [a, b, c, d, e, f]
    //--------------------------------------------------------------------------
    
    // Exact solution from 3 points
    inline bool affineFrom3Points(
        const std::array<double, 3>& px, const std::array<double, 3>& py,  // pixel coords
        const std::array<double, 3>& gx, const std::array<double, 3>& gy,  // geo coords
        Vector<6>& coeffs)
    {
        // Build 6x6 system
        // [px0 py0 1  0   0   0] [a]   [gx0]
        // [px1 py1 1  0   0   0] [b]   [gx1]
        // [px2 py2 1  0   0   0] [c] = [gx2]
        // [0   0   0  px0 py0 1] [d]   [gy0]
        // [0   0   0  px1 py1 1] [e]   [gy1]
        // [0   0   0  px2 py2 1] [f]   [gy2]
        
        Matrix<6, 6> A{};
        Vector<6> b{};
        
        for (size_t i = 0; i < 3; ++i) {
            A[i][0] = px[i];
            A[i][1] = py[i];
            A[i][2] = 1.0;
            b[i] = gx[i];
            
            A[i + 3][3] = px[i];
            A[i + 3][4] = py[i];
            A[i + 3][5] = 1.0;
            b[i + 3] = gy[i];
        }
        
        return solve6x6(A, b, coeffs);
    }
    
    // Least squares solution from 4 points
    inline bool affineFrom4Points(
        const std::array<double, 4>& px, const std::array<double, 4>& py,
        const std::array<double, 4>& gx, const std::array<double, 4>& gy,
        Vector<6>& coeffs)
    {
        // Build 8x6 system (overdetermined)
        Matrix<8, 6> A{};
        Vector<8> b{};
        
        for (size_t i = 0; i < 4; ++i) {
            A[i][0] = px[i];
            A[i][1] = py[i];
            A[i][2] = 1.0;
            b[i] = gx[i];
            
            A[i + 4][3] = px[i];
            A[i + 4][4] = py[i];
            A[i + 4][5] = 1.0;
            b[i + 4] = gy[i];
        }
        
        return leastSquares8x6(A, b, coeffs);
    }
    
    // Least squares solution from 5 points
    inline bool affineFrom5Points(
        const std::array<double, 5>& px, const std::array<double, 5>& py,
        const std::array<double, 5>& gx, const std::array<double, 5>& gy,
        Vector<6>& coeffs)
    {
        // Build 10x6 system (overdetermined)
        Matrix<10, 6> A{};
        Vector<10> b{};
        
        for (size_t i = 0; i < 5; ++i) {
            A[i][0] = px[i];
            A[i][1] = py[i];
            A[i][2] = 1.0;
            b[i] = gx[i];
            
            A[i + 5][3] = px[i];
            A[i + 5][4] = py[i];
            A[i + 5][5] = 1.0;
            b[i + 5] = gy[i];
        }
        
        return leastSquares10x6(A, b, coeffs);
    }
    
    //--------------------------------------------------------------------------
    // Convert 6 affine coefficients to 4x4 transformation matrix (row-major)
    // Input: [a, b, c, d, e, f] where gx = a*px + b*py + c, gy = d*px + e*py + f
    // Output: 4x4 matrix for GeoTIFF ModelTransformation
    //--------------------------------------------------------------------------
    
    inline std::array<double, 16> affineCoeffsToMatrix4x4(const Vector<6>& coeffs) {
        return {
            coeffs[0], coeffs[1], 0.0, coeffs[2],  // row 0: a, b, 0, c
            coeffs[3], coeffs[4], 0.0, coeffs[5],  // row 1: d, e, 0, f
            0.0,       0.0,       0.0, 0.0,        // row 2: 0, 0, 0, 0
            0.0,       0.0,       0.0, 1.0         // row 3: 0, 0, 0, 1
        };
    }

}  // namespace blas
