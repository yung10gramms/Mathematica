package mathematca.matrix;

import java.util.Random;

import mathematca.Complex;

public class LinearAlgebra {

	private static final int MINMODE = 1;
	private static final int MAXMODE = -1;
				
	private static int COUNTER = 0;
	
	protected static void initcounter() {
		COUNTER = 0;
	}
	
	protected static int getcounter() {
		return COUNTER;
	}
	
	protected static void incrcounter() {
		COUNTER++;
	}
	
	public LinearAlgebra() {
		// TODO Auto-generated constructor stub
	}

	protected static double[][] matmul2(double[][] A, double[][] B) {
		int lenAy = A.length;
		int lenAx = A[0].length;
		int lenBy = B.length;
		int lenBx = B[0].length;
		
		if(lenAx != lenBy) {
			System.out.println("matrix dimensions do not match.\n");
			return null;
		}
		
		double[][] out = new double[lenAy][lenBx];
		
		for(int col = 0; col < lenBx; col ++) {
			for(int row = 0; row < lenAy; row ++) {
				out[row][col] = dotProduct(getRow(A, row), getColumn(B, col));
				incrcounter();
			}
		}
		
		return out;
	}
	
	public static double[][] matmul(double[][]... A) {
		MatrixMultiplication matrixMul = new MatrixMultiplication();
		if(A.length == 1)
			return A[0];
		if(A.length == 2)
			return matmul2(A[0], A[1]);
		return matrixMul.multiply(A);
	}
	
	public static double[][] randn(int dimA, int dimB) {
		Random rand = new Random();
		double[][] output = new double[dimA][dimB];
		for(int i = 0; i < dimA; i ++)
			for(int j = 0; j < dimB; j ++)
				output[i][j] = rand.nextGaussian();
		return output;
	}
	
	public static double[][] rand(int dimA, int dimB) {
		Random rand = new Random();
		double[][] output = new double[dimA][dimB];
		for(int i = 0; i < dimA; i ++)
			for(int j = 0; j < dimA; j ++)
				output[i][j] = rand.nextDouble();
		return output;
	}
	
	public static double dotProduct(double[] u, double[] v) {
		if(u.length != v.length){
			System.out.println("vector dimensions do not match.\n");
			return 0;
		}
		double out = 0;
		for(int i = 0; i < v.length; i ++) {
			out += u[i]*v[i];
		}
		return out;
	}
	
	public static double[] getRow(double[][] mat, int row) {
		double[] out = new double[mat[0].length];
		for(int i = 0; i < out.length; i ++) {
			out[i] = mat[row][i];
		}
		return out;
	}
	
	public static double[] getColumn(double[][] mat, int column) {
		double[] out = new double[mat.length];
		for(int i = 0; i < out.length; i ++) {
			out[i] = mat[i][column];
		}
		return out;
	}
	
	public static void printMat(double[][] mat) {
		for(int i = 0; i < mat.length; i ++) {
			for(int j = 0; j < mat[i].length; j ++) {
				System.out.print(""+mat[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	private static int[] argminmax(double[][] mat, int func) {
		int[] argmin = new int[2];
		double min = (func == MINMODE) ? Double.MAX_VALUE : Double.MIN_VALUE;
		int mul = (func == MINMODE) ? 1 : -1;
		
		for(int i = 0; i < mat.length; i ++) {
			for(int j = 0; j < mat[i].length; j ++) {
				if(mat[i][j]*mul < min*mul) {
					min = mat[i][j];
					argmin[0] = i;
					argmin[1] = j;
				}
			}
		}
		return argmin;
	}
	
	public static int[] argmin(double[][] mat) {
		return argminmax( mat, MINMODE);
	}
	
	public static int[] argmax(double[][] mat) {
		return argminmax(mat, MAXMODE);
	}
	
	public static int argmax(double[] mat) {
		return argmax(new double[][]{mat})[0]; // is this efficient?
	}
	
	public static double min(double[][] mat) {		
		return minmax(mat, MINMODE);
	}
	
	public static double max(double[][] mat) {		
		return minmax(mat, MAXMODE);
	}
	
	private static double minmax(double[][] mat, int func) {
		double min = (func == MINMODE) ? Double.MAX_VALUE : Double.MIN_VALUE;
		int mul = (func == MINMODE) ? 1 : -1;
		
		for(int i = 0; i < mat.length; i ++) {
			for(int j = 0; j < mat[i].length; j ++) {
				if(mat[i][j]*mul < min*mul) {
					min = mat[i][j];		
				}
			}
		}
		return min;
	}
	
	/**
	 * Element-wise division. The output is a new array where out[i] = a[i]/b[i]
	 * 
	 * @param a nominator
	 * @param b denominator
	 * */
	public static Complex[] arrdiv(Complex[] a, Complex[] b) {
		int size = a.length;
		if(size != b.length)
			return null;
		
		Complex[] out = new Complex[size];
		for(int i = 0; i < size; i ++) {
			if(b[i].mod() == 0) {
				System.out.println("division by 0 not possible");
				return null;
			}
			out[i] = a[i].div(b[i]);
		}
		return out;
	}
	
	/**
	 * Element-wise multiplication. The output is a new array where out[i] = a[i]*b[i]
	 * 
	 * @param a first term
	 * @param b second term
	 * */
	public static Complex[] arrmul(Complex[] a, Complex[] b) {
		int size = a.length;
		if(size != b.length)
			return null;
		
		Complex[] out = new Complex[size];
		for(int i = 0; i < size; i ++) {
			out[i] = a[i].times(b[i]);
		}
		return out;
	}
	
	public static double[] linspace(double start, double end, int noSamples) {
		double step = (end-start)/noSamples;
		double[] out = new double[noSamples];
		for(int i = 0; i < noSamples; i ++) {
			out[i] = i*step + start;
		}
		return out;
	}
	
	public static double[] range(double start, double end, double step) {
		int noSamples = (int) ((end-start)/step);
		double[] out = new double[noSamples];
		for(int i = 0; i < noSamples; i ++) {
			out[i] = i*step + start;
		}
		return out;
	}
	
	public static double max(double... numbers) {
	    double maxel = Double.MIN_VALUE;
	    for (double num : numbers) {
	        maxel = maxel > num ? maxel : num;
	    }
	    return maxel;
	}	
	
	public static double norm(double[] vec, double p) {
		double out = 0;
		for(double c : vec) {
			out += Math.pow(c, p);
		}
		return Math.sqrt(out);
	}
	
	public static double norm2(double[] vec) {
		return norm(vec, 2);
	}
	
	/**
	 * fuction that returns a*vec
	 * */
	public static double[] scale(double[] vec, double a) {
		double[] out = new double[vec.length];
		for(int i = 0; i < vec.length; i ++) {
			out[i] = vec[i]*a;
		}
		return out;
	}
	
	/**
	 * fuction that returns x/||x||
	 * */
	public static double[] normalize(double[] vec) {
		return scale(vec, norm2(vec));
	}
	
	public static double[][] eye(int n) {
		double[][] out = zeros(n, n);
		for(int i = 0; i < n; i ++) {
			out[i][i] = 1;
		}
		return null;
	}
	
	/**
	 * n x m zero matrix
	 * 
	 * */
	public static double[][] zeros(int n, int m) {
		return new double[n][m];
	}
	
	/**
	 * Function that returns a submatrix of A, from columns <b>startx</b> to <b>endx</b> and
	 * from rows <b>starty</b> to <b>endy</b>
	 * 
	 * @param A input matrix
	 * @param starty row to start
	 * @param endy row to end
	 * @param startx column to start
	 * @param endx column to end
	 * 
	 * */
	public static double[][] subMatrix(double[][] A, int starty, int endy, int startx, int endx) {
		if(endx < startx)
			return null;
		if(endy < starty)
			return null;
		
		double[][] out = new double[endy-starty+1][endx-startx+1];
		for(int i = starty; i <= endy; i ++)
			for(int j = startx; j <= endx; j ++)
				out[i-starty][j-startx] = A[i][j];
		return out;
	}
	
	public static double[][] swapRows(double[][] A, int r1, int r2) {
		
		double[][] out = A.clone();
		for(int i = 0; i <= A.length; i ++)
		{	
			out[i][r1] = A[i][r2];
			out[i][r2] = A[i][r1];
		}	
		return out;
	}
	
	public static void swapRowsInPlace(double[][] A, int r1, int r2) {			
		double current;
		for(int i = 0; i < A.length; i ++)
		{	
			current = A[i][r1];
			A[i][r1] = A[i][r2];
			A[i][r2] = current;
		}			
	}
	
	public static void swapColsInPlace(double[][] A, int r1, int r2) {			
		double current;
		for(int i = 0; i < A[0].length; i ++)
		{	
			current = A[r1][i];
			A[r1][i] = A[r2][i];
			A[r2][i] = current;
		}			
	}
	
	public static double[][] cat(double[][] A, double[][] B, final String mode) {
		if(mode.equals("columns")) {
			if(A.length != B.length) {
				System.out.println("Concatenation not possible due to mismatch in dimensions (" +A.length+ "x" + A[0].length+", "+B.length+"x"+B[0].length+").");
				return null;
			}
			double[][] out = new double[A.length][A[0].length + B[0].length];
			for(int i = 0; i < A.length; i ++) {
				for(int j = 0; j < A[i].length; j ++) {
					out[i][j] = A[i][j];
				}
				for(int j = 0; j < B[i].length; j ++) {
					out[i][A[i].length + j] = A[i][j];
				}
			}
			return out;	
		} else if(mode.equals("rows")) {
			if(A[0].length != B[0].length) {

				System.out.println("Concatenation not possible due to mismatch in dimensions (" +A.length+ "x" + A[0].length+", "+B.length+"x"+B[0].length+").");
				return null;
			}
			double[][] out = new double[A.length + B.length][A[0].length];
			for(int i = 0; i < A.length; i ++) {
				for(int j = 0; j < A[i].length; j ++) {
					out[i][j] = A[i][j];
				}
			}
			for(int i = 0; i < B.length; i ++) {
				for(int j = 0; j < B[i].length; j ++) {
					out[i][j] = B[i+A.length][j];
				}
			}
			
			return out;
		}
		System.out.println("Error in mode.");
		return null;
	}
	
	/**
	 * Convert vector to matrix. You can copy it either to a row or to a column.
	 * 
	 * @param vector vector to convert to matrix
	 * @param mode preferred mode. It can be "row" or "column"
	 * */
	public static double[][] toMatrix(double[] vector, final String mode) {
		if(mode.equals("column")) {
			double[][] out = new double[vector.length][1];
			for (int i = 0; i < vector.length; i ++) {
				out[i][0] = vector[i];
			} 
			return out;
		} else if(mode.equals("row")) {
			double[][] out = new double[1][vector.length];
			for (int i = 0; i < vector.length; i ++) {
				out[0][i] = vector[i];
			} 
			return out;
		}
		System.out.println("Error in mode.");
		return null;
	}
	
	/**
	 * In place implementation of Gaussian Elimination algorithm. The pseudo code of 
	 * the algorithm can be found here 
	 * {@link https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode} 
	 * 
	 * @param A matrix to perform the algorithm on
	 *  */
	public static void gaussianElimination(double[][] A) {
		int h = 0, k = 0, m = A.length-1, n = A[0].length-1;
		
		while(h <= m && k <= n) {
			int i_max = argmax(subMatrix(A, h, m, k, k))[0] + h;
			if(A[i_max][k] == 0) /* No pivot in this column, pass to next column */
				k ++;
			else {
				swapRowsInPlace(A, h, i_max);
				for (int i = h + 1; i <= m; i ++) {
					double f = A[i][k] / A[h][k];
					A[i][k] = 0;
					for(int j = k+1; j <= n; j ++) {
						A[i][j] = A[i][j] - A[h][j] * f;
					}
				}
				h ++;
				k ++;
			}
		}				
	}
	
	public static double[] linsolve(double[][] A, double[] b) {
		double[][] aug = cat(A, toMatrix(b, "column"), "columns"); 
		gaussianElimination(aug);
		return getRow(aug, 0);
	}
	
	/**
	 * function that calls the LU decomposition function of the class
	 * LUDecomposition made in this project. As I cite there as well, the code for the 
	 * implementation of the LU algorithm can be found (in C, C# and MATLAB) in this wikipedia
	 * page:
	 * {@link https://en.wikipedia.org/wiki/LU_decomposition#Algorithms}
	 * 
	 * <p> 
	 * Note that you cannot instantiate, or call directly the functions of LUDecomposition
	 * class outside of this package, since I thought that the best design for this 
	 * project would be to have one class responsible for all matrix operations
	 * 
	 * @param A matrix to decompose
	 * @return L U matrices (in one structure)
	 * @throws NullPointerException if demensions do not match
	 * */
	public static double[][] LUDecompose(double[][] A) {
		int N = A.length;
		double[] P = range(0, N+1, 1);
		double[][] Aclone = A.clone();
		LUDecomposition lud = new LUDecomposition(Aclone, P);
		lud.LUPDecompose(0.01*N);
		return Aclone;
	}
	
	/**
	 * function that calls the LU-linear-equation-solver function of the class
	 * LUDecomposition made in this project. As I cite there as well, the code for the 
	 * implementation of the LU algorithm can be found (in C, C# and MATLAB) in this wikipedia
	 * page:
	 * {@link https://en.wikipedia.org/wiki/LU_decomposition#Algorithms}
	 * 
	 * <p>
	 * The system to solve looks like this: Ax = b, where A is a square matrix of dimension N, 
	 * b is a vector of dimension N as well, and x is the vector to determine 
	 * 
	 * <p> 
	 * Note that you cannot instantiate, or call directly the functions of LUDecomposition
	 * class outside of this package, since I thought that the best design for this 
	 * project would be to have one class responsible for all matrix operations
	 * 
	 * 
	 * @param A square A matrix 
	 * @param b right hand side vector
	 * @return x the solution 
	 * @throws NullPointerException if demensions do not match
	 * 
	 * */
	public static double[] LULinsolve(double[][] A, double[] b) {
		int N = A.length;
		
		if(b == null) 
		{
			System.out.println("\nVector b cannot be null\n");
			throw new NullPointerException();
		}
		
		if(b.length != N)
		{
			System.out.println("\nVector b must have the same dimensions as A\n");
			throw new NullPointerException();
		}
		
		double[] P = range(0, N+1, 1);
		double[][] Aclone = A.clone();
		
		LUDecomposition lud = new LUDecomposition(Aclone, P);
		lud.LUPDecompose(0.01*N);
		
		double[] x = new double[N];
		lud.LUPSolve(b, x);
		
		return x;
	}
	
	/**
	 * function that calls the LU-invert function of the class
	 * LUDecomposition made in this project. As I cite there as well, the code for the 
	 * implementation of the LU algorithm can be found (in C, C# and MATLAB) in this wikipedia
	 * page:
	 * {@link https://en.wikipedia.org/wiki/LU_decomposition#Algorithms}
	 * 
	 * <p>
	 * The function finds the matrix IA, where IA * A = I, where A is the matrix inserted 
	 * in the input of the function and I is the identity matrix of demension N
	 * 
	 * <p> 
	 * Note that you cannot instantiate, or call directly the functions of LUDecomposition
	 * class outside of this package, since I thought that the best design for this 
	 * project would be to have one class responsible for all matrix operations
	 * 
	 * 
	 * @param A square A matrix 
	 * @return IA the inverse of the matrix
	 * @throws NullPointerException if demensions do not match
	 * 
	 * */
	public static double[][] LUPInvert(double[][] A) {
		int N = A.length;
		
		double[] P = range(0, N+1, 1);
		double[][] Aclone = A.clone();
		
		LUDecomposition lud = new LUDecomposition(Aclone, P);
		lud.LUPDecompose(0.01*N);
		
		double[][] IA = new double[N][N];
		
		lud.LUPInvert(IA);
		return IA;
	}
}
