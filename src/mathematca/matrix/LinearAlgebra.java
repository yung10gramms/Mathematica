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
		
		double[][] out = new double[endy-starty][endx-startx];
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
		for(int i = 0; i <= A.length; i ++)
		{	
			current = A[i][r1];
			A[i][r1] = A[i][r2];
			A[i][r2] = current;
		}			
	}
	
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
}
