package mathematca.matrix;

import java.util.HashMap;

public class MatrixMultiplication {

	protected int[][] m;
    protected int[][] s;

	private double[][] a_matrix;
    
	protected MatrixMultiplication() {
		
	}
	
    protected void matrixChainOrder(int[] dims) {
        int n = dims.length - 1;
        m = new int[n][n];
        s = new int[n][n];

        for (int lenMinusOne = 1; lenMinusOne < n; lenMinusOne++) {
            for (int i = 0; i < n - lenMinusOne; i++) {
                int j = i + lenMinusOne;
                m[i][j] = Integer.MAX_VALUE;
                for (int k = i; k < j; k++) {
                    int cost = m[i][k] + m[k+1][j] + dims[i]*dims[k+1]*dims[j+1];
                    if (cost < m[i][j]) {
                        m[i][j] = cost;
                        s[i][j] = k;
                    }
                }
            }
        }
    }

    public void printOptimalParenthesizations() {
        boolean[] inAResult = new boolean[s.length];
        printOptimalParenthesizations(s, 0, s.length - 1, inAResult);
    }

    void printOptimalParenthesizations(int[][]s, int i, int j,  /* for pretty printing: */ boolean[] inAResult) {
        if (i != j) {
            printOptimalParenthesizations(s, i, s[i][j], inAResult);
            printOptimalParenthesizations(s, s[i][j] + 1, j, inAResult);
            String istr = inAResult[i] ? "_result " : " ";
            String jstr = inAResult[j] ? "_result " : " ";
            System.out.println(" A_" + i + istr + "* A_" + j + jstr);
            inAResult[i] = true;
            inAResult[j] = true;
        }
    }


    private void multiplyHelp(int[][]s, int i, int j,  HashMap<Integer, double[][]> matMap, boolean[] inAResult) {
        if (i != j) {
        	multiplyHelp(s, i, s[i][j], matMap, inAResult);
        	multiplyHelp(s, s[i][j] + 1, j, matMap, inAResult);
            
            double[][] firstMat = inAResult[i] ? a_matrix : matMap.get(i);
            double[][] secondMat = inAResult[j] ? a_matrix : matMap.get(j);
            
            a_matrix = LinearAlgebra.matmul2(firstMat, secondMat);
            inAResult[i] = true;
            inAResult[j] = true;
        }
    }
    
    protected double[][] multiply(double[][]... matrices) {
    	int numOfMats = matrices.length;
    	int[] dimensions = new int[numOfMats+1];
    	
    	HashMap<Integer, double[][]> map = new HashMap<>();
    	int p = 1;
    	int previousDimension = -1;
    	for(double[][] mat : matrices) {
    		map.put(p-1, mat);
    		if(p == 1) {
    			dimensions[0] = mat.length; 
    			previousDimension = dimensions[0];
    		}
    		if(mat.length != previousDimension)
    		{
    			System.out.println("matrix dimensions do not agree");
    			return null;
    		}
    		dimensions[p] = mat[0].length;
    		previousDimension = dimensions[p];
    		p ++;
    	}
    	matrixChainOrder(dimensions);
    	boolean[] inAResult = new boolean[s.length];
        multiplyHelp(s, 0, s.length - 1, map, inAResult);
        return a_matrix;
    }
    
    protected void test1() {
    	MatrixMultiplication mat = new MatrixMultiplication();
    	int[] dims = {100, 12, 2932, 2138, 1238, 128};
    	mat.matrixChainOrder(dims);
    	mat.printOptimalParenthesizations();
    	System.out.println();
    	
    	System.out.println();
    	for(int[] c1 : mat.s) {
    		for(int c2 : c1) 
    			System.out.print(c2+" "); 
    		System.out.println();
    	}
    }
    
    protected static void test2() {
    	MatrixMultiplication mat = new MatrixMultiplication();
    	
    	double[][] mat1 = {{1.2, 3.21, 4.32},{1.2, 3.11, 4.1},{1.2, 2.32, 4.3}};
    	double[][] mat2 = LinearAlgebra.randn(3, 15);
    	double[][] mat3 = LinearAlgebra.randn(15, 8);
    	double[][] mat4 = LinearAlgebra.randn(8, 15);
    	
    	double[][] res = mat.multiply(mat1, mat2, mat3, mat4);
    	
    	System.out.println();
    	printMatrixMatformat(mat1);
    	printMatrixMatformat(mat2);
    	printMatrixMatformat(mat3);
    	printMatrixMatformat(res);
    	
    	System.out.println("\n\nsmart way: "+LinearAlgebra.getcounter());
    	
    	LinearAlgebra.initcounter();
    	double[][] buf1 = LinearAlgebra.matmul2(mat1, mat2);
    	buf1 = LinearAlgebra.matmul2(buf1, mat3);
    	buf1 = LinearAlgebra.matmul2(buf1, mat4);
    	System.out.println("\n\ndumb way : "+LinearAlgebra.getcounter());
    	System.out.println();
    	
    }
    
    public static void printMatrix(double[][] mat) {
    	System.out.println();
    	for(double[] c1 : mat) {
    		for(double c2 : c1) 
    			System.out.print(c2+" "); 
    		System.out.println();
    	}
    }
    
    public static void printMatrixMatformat(double[][] mat) {
    	System.out.println();
    	System.out.print("[");
    	for(int i = 0; i < mat.length; i ++) {
    		for(int j = 0; j < mat[i].length; j ++) {
    			if(j != 0)
    				System.out.print(", ");
    			System.out.print(mat[i][j]);   			
    		}
    		if(i != mat.length - 1)
    			System.out.print("; ");
    	}
    	System.out.print("];");
    }
    
    public static void main(String[] args) {
    	
    	 test2();
    }
    
}
