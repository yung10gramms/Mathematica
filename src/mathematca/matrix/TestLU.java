package mathematca.matrix;

public class TestLU {

	private TestLU() {
		// TODO Auto-generated constructor stub
	}

	
	private static void test1(double[][] A) {
		LinearAlgebra.printMat(A);
		double[] P = new double[5];
		LUDecomposition lud = new LUDecomposition(A,P);	
		lud.LUPDecompose(0.1);
		System.out.println();
		LinearAlgebra.printMat(A);
	}
	
	private static void test2(double[][] A) {
		System.out.println("\n\n");
		double[][] a = LinearAlgebra.LUDecompose(A);
		LinearAlgebra.printMat(a);
		System.out.println("\n\n");
	}
	
	public static void main(String[] args) {
		double[][] A = LinearAlgebra.rand(4, 4);
		test1(A);
		test2(A);
	}
}
