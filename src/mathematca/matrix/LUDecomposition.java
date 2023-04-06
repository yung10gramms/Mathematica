package mathematca.matrix;


public class LUDecomposition {

	private double[][] A;
	private double[] P;
	
	/** Class is a Java implementation of the LU functions 
	 * cited in {@link https://en.wikipedia.org/wiki/LU_decomposition#Algorithms}
	 * 
	 * @param A array of pointers to rows of a square matrix having dimension N
	 * @param P array (must be initialized to length N+1)
	 *  */
	protected LUDecomposition(double[][] a, double[] p) {
		if(a == null) {
			System.out.println("a must not be null\n\n");
			throw new NullPointerException();
		}
		if(a.length == 0) {
			System.out.println("a must have dimension greater than 0\n\n");
			throw new NullPointerException();
		}
		if(a.length != a[0].length) {
			System.out.println("Class only works with square matrices\n\n");
			throw new NullPointerException();
		} 
		if(p == null || p.length != a.length + 1) {
			System.out.println("p array must be initialized to length N+1\n\n");
			throw new NullPointerException();
		}
		
		A = a;
		P = p;
	}

	/** 
	 * 
	 * 
	 * <p>
	 *  Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
	 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
	 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
	 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
	 * <p>
	 * see also {@link https://en.wikipedia.org/wiki/LU_decomposition#Algorithms}
	 * 
	 * @param Tol small tolerance number to detect failure when the matrix is near degenerate
	 * 
	 */
	int LUPDecompose(double Tol) {

		int N = A.length;
	    int i, j, k, imax; 
	    double maxA, absA;
	    double[] ptr;

	    for (i = 0; i <= N; i ++)
	    	P[i] = i;
	    
	    for (i = 0; i < N; i++) {
	        maxA = 0.0;
	        imax = i;

	        for (k = i; k < N; k++)
	            if ((absA = Math.abs(A[k][i])) > maxA) { 
	                maxA = absA;
	                imax = k;
	            }

	        if (maxA < Tol) return 0; //failure, matrix is degenerate

	        if (imax != i) {
	            //pivoting P
	            j = (int) P[i];
	            P[i] = P[imax];
	            P[imax] = j;

	            //pivoting rows of A
	            ptr = A[i];
	            A[i] = A[imax];
	            A[imax] = ptr;

	            //counting pivots starting from N (for determinant)
	            P[N]++;
	        }

	        for (j = i + 1; j < N; j++) {
	            A[j][i] /= A[i][i];

	            for (k = i + 1; k < N; k++)
	                A[j][k] -= A[j][i] * A[i][k];
	        }
	    }

	    return 1;  //decomposition done 
	}

	/*
	 * 
	 * inputs are A,P, filled in LUPDecompose; b is the RHS vector;
	 * OUTPUT: x - solution vector of A*x=b
	 */
	void LUPSolve(double[] b, double[] x) {
		
		
		int N = A.length;
		
	    for (int i = 0; i < N; i++) {
	        x[i] = b[(int) P[i]];

	        for (int k = 0; k < i; k++)
	            x[i] -= A[i][k] * x[k];
	    }

	    for (int i = N - 1; i >= 0; i--) {
	        for (int k = i + 1; k < N; k++)
	            x[i] -= A[i][k] * x[k];

	        x[i] /= A[i][i];
	    }
	}

	/* INPUT: A,P filled in LUPDecompose; N - dimension
	 * OUTPUT: IA is the inverse of the initial matrix
	 */
	void LUPInvert(double[][] IA) {
		int N = A.length;
		
	    for (int j = 0; j < N; j++) {
	        for (int i = 0; i < N; i++) {
	            IA[i][j] = P[i] == j ? 1.0 : 0.0;

	            for (int k = 0; k < i; k++)
	                IA[i][j] -= A[i][k] * IA[k][j];
	        }

	        for (int i = N - 1; i >= 0; i--) {
	            for (int k = i + 1; k < N; k++)
	                IA[i][j] -= A[i][k] * IA[k][j];

	            IA[i][j] /= A[i][i];
	        }
	    }
	}

	/* INPUT: A,P filled in LUPDecompose; N - dimension. 
	 * OUTPUT: Function returns the determinant of the initial matrix
	 */
	double LUPDeterminant() {
		int N = A.length;
	    double det = A[0][0];

	    for (int i = 1; i < N; i++)
	        det *= A[i][i];

	    return (P[N] - N) % 2 == 0 ? det : -det;
	}
}
