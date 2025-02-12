#ifndef DYCKERHOFF_MATRIX_H
#define DYCKERHOFF_MATRIX_H

#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <memory>
#include <functional>
#include <stdexcept>

//#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#define new DEBUG_NEW


namespace dyMatrixClass {

	std::ostream& shortMat(std::ostream & os);
	std::ostream& longMat(std::ostream & os);

	enum class matrix_direction { mLeft = 1, mRight = 2, mLeftRight = 3};

	struct definiteness_error : std::runtime_error {
		definiteness_error(std::string s) : std::runtime_error{ s } {}
	};

	class cMatrix {
	protected:
		int m, n;
		std::unique_ptr<double[]> elems;
	public:
		cMatrix() = default;
		cMatrix(cMatrix&& m);															// Move constructor
		cMatrix& operator=(cMatrix&& m);												// Move assignment
		cMatrix(const cMatrix&);														// Copy constructor
		cMatrix& operator=(const cMatrix&);												// Copy assignment

		cMatrix(int m, int n)															// Constructor
			: m{ m }, n{ n }, elems{ new double[m * n]{} } {};
		cMatrix(std::initializer_list<std::initializer_list<double>> args);				// List initializer constructor
		cMatrix(std::initializer_list<double> args);									// List initializer constructor (column vector)

		void SetSize(int _m, int _n);
		void Resize(int _m, int _n, bool keepData = false);
		long index(int i, int j) const { return i * n + j; }
		double& operator()(int i, int j) const { return elems[i * n + j]; }
		double& operator()(int i) const { return elems[i]; }
		double* operator[](int i) const { return elems.get() + i * n; }

		double* data() { return elems.get(); }
		const double* data() const { return elems.get(); }
		double* begin() { return elems.get(); }
		const double* begin() const { return elems.get(); }
		double* end() { return elems.get() + m * n; }
		const double* end() const { return elems.get() + m * n; }

		int Rows() const { return m; }
		int Cols() const { return n; }

		cMatrix Row(int i) const;
		cMatrix Col(int j) const;
		cMatrix Diag() const;
		cMatrix SubMatrix(int i1, int i2, int j1, int j2) const;

		operator double() const;

		void XChangeRows(int i1, int i2) { std::swap_ranges(elems.get() + i1*n, elems.get() + (i1 + 1)*n, elems.get() + i2*n); }
		void XChangeCols(int j1, int j2) { for (int i = 0; i < m; i++) std::swap(elems[i*n + j1], elems[i*n + j2]); }

		cMatrix& apply(std::function<double(double)> f);

		cMatrix& operator*=(double c);
		cMatrix& operator/=(double c);

		cMatrix& operator+=(const cMatrix& B);
		cMatrix& operator-=(const cMatrix& B);
		cMatrix& operator*=(const cMatrix& B);

		cMatrix operator-();

		friend cMatrix operator+(const cMatrix& A, const cMatrix& B);
		friend cMatrix operator-(const cMatrix& A, const cMatrix& B);
		friend cMatrix operator*(const cMatrix& A, const cMatrix& B);

		friend bool operator==(const cMatrix& A, const cMatrix& B);
		friend bool operator!=(const cMatrix& A, const cMatrix& B);
		friend bool operator<=(const cMatrix& A, const cMatrix& B);
		friend bool operator>=(const cMatrix& A, const cMatrix& B);
		friend bool operator< (const cMatrix& A, const cMatrix& B);
		friend bool operator> (const cMatrix& A, const cMatrix& B);

		friend std::ostream& operator<< (std::ostream& os, const cMatrix& A);
		friend std::istream& operator>> (std::istream& is, const cMatrix& A);
		friend cMatrix Trans(const cMatrix& A);

		friend cMatrix operator*(const cMatrix& A, double c);
		friend cMatrix operator*(double c, const cMatrix& A);
		friend cMatrix operator/(const cMatrix& A, double c);

		friend double FrobeniusNorm(const cMatrix& A);
		friend double RowSumNorm(const cMatrix& A);
		friend double ColSumNorm(const cMatrix& A);

		friend double Sum(const cMatrix& A);

		cMatrix& SetEpsToZero(double eps);

		cMatrix& SetDiag(double D[]) { for (int i = 0; i < std::min(m, n); i++) elems[index(i, i)] = D[i]; return *this; };

		void Householder1(double& New, int first, int last, int pos, bool transpose = false);
		void Householder2(double& New, int first, int last, int pos, bool Symm = false);
		void ApplyHH(double w[], const cMatrix& Q, cMatrix& B, bool Transpose = false) const;
		void ApplyModifiedHouseholder(const double w[], const double y[], int first, int last, int pos, matrix_direction Direction);
		void ApplyGivens(double c, double s, int k1, int k2, matrix_direction Direction);
		void AccumulateHH(double D[], int Ofs, int nh) const;
		void AccumulateHH(double D[], int Ofs, int nh, cMatrix& Q, bool Transpose = false, int c1 = 0, int c2 = -1) const;

		void BiDiag (double HD[], double ND[], cMatrix* U = nullptr, cMatrix* V = nullptr) const;
		void TriDiag(double HD[], double ND[], cMatrix* U = nullptr) const;
		void Hessenberg(cMatrix& H, cMatrix* U = nullptr) const;
		void Hessenberg(cMatrix& H, cMatrix& U) const { return Hessenberg(H, &U); }

		cMatrix& Rref();

		void LU_Factorization(cMatrix& LU, int p[], double& det) const;
		void QR_Factorization(cMatrix& QR, double D[], int p[], bool Transpose = false, int* rank = nullptr) const;

		cMatrix Nullspace_QR();
		cMatrix OrthogonalComplement();

		void Cholesky(cMatrix& G, int* pRank = nullptr) const;
		void InversePosDef(cMatrix& InvA, int* pRank = nullptr) const;
		//friend cMatrix InversePosDef(const cMatrix& A, int* pRank = nullptr);
        friend cMatrix InversePosDef(const cMatrix& A, int* pRank);
		//friend cMatrix LGS_Solve_PosDef(const cMatrix& A, const cMatrix& B, int* pRank = nullptr);
        friend cMatrix LGS_Solve_PosDef(const cMatrix& A, const cMatrix& B, int* pRank);

		void LGS_Solve_LU(int p[], const cMatrix& B, cMatrix& X);
		friend cMatrix LGS_Solve_LU(const cMatrix& A, const cMatrix& B);
		void LGS_Solve_QR(double D[], int p[], int rank, const cMatrix& B, cMatrix& X, bool transpose = false);
		friend cMatrix LGS_Solve_QR(const cMatrix& A, const cMatrix& B);

		cMatrix& Invert();
		friend cMatrix Inverse(cMatrix A);

		int EigenSymm(double EW[], cMatrix* U = nullptr);
		int EigenSymm(double EW[], cMatrix& U) { return EigenSymm(EW, &U); }

		int EigenQR(std::complex<double> EW[], cMatrix* H = nullptr, cMatrix* U = nullptr);
		int EigenQR(std::complex<double> EW[], cMatrix& H)             { return EigenQR(EW, &H); }
		int EigenQR(std::complex<double> EW[], cMatrix& H, cMatrix& U) { return EigenQR(EW, &H, &U); };

		int SVD(double w[], cMatrix* U = nullptr, cMatrix* V = nullptr);
		int SVD(double w[], cMatrix& U)             { return SVD(w, &U);     }
		int SVD(double w[], cMatrix& U, cMatrix& V) { return SVD(w, &U, &V); }
	};

	class cVector : public cMatrix {
	public:
		cVector() = default;

		cVector& operator=(cMatrix&& m);												// Move assignment
		cVector& operator=(const cMatrix&);												// Copy assignment

		cVector(int n) : cMatrix(n, 1) {};												// Constructor
			

		cVector(std::initializer_list<double> args);									// List initializer constructor (column vector)

		double& operator[](int i) const { return elems[i]; }
		void SetSize(int _n) { cMatrix::SetSize(1, _n); }

		friend double InnerProduct(const cVector& A, const cVector& B);
	};

}

#endif
