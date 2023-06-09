package mathematca;

import mathematca.matrix.LinearAlgebra;

public class TestLinA {

	public TestLinA() {
		// TODO Auto-generated constructor stub
	}
	
	public static void main(String[] args) {
		double[][] a = new double[4][3];
		double[][] b = new double[3][30];
		
		for(double[] ar : a)
			for(int i = 0; i < ar.length; i ++)
				ar[i] = i;
//		LinearAlgebra.printMat(a);
		System.out.println();
		
		for(double[] ar : b)
			for(int i = 0; i < ar.length; i ++)
				ar[i] = i;
//		LinearAlgebra.printMat(b);
		System.out.println();
		
		double[][] matmat = LinearAlgebra.matmul(a, b);
		LinearAlgebra.printMat(matmat);
		System.out.println("\n\n");
		LinearAlgebra.gaussianElimination(matmat);
		LinearAlgebra.printMat(matmat);
		
		System.out.println("\n\n");
		double[][] matrix = {{1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}};
		double[] vector = {1, 2, 3, 4};
		double[] solution = LinearAlgebra.linsolve(matrix, vector);
		LinearAlgebra.printMat(LinearAlgebra.toMatrix(solution, "row"));
		
	}

}
