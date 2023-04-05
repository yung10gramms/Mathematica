package mathematca;

import mathematca.matrix.LinearAlgebra;

public class TestComplex {

	public TestComplex() {
		
	}
	
	public static void main(String[] args) {
		double[] array = LinearAlgebra.range(-1, 1, 0.01);
		Complex[] array_complex = 
				ComplexOperations.scale(ComplexOperations.toComplex(array), 
							ComplexOperations.complex(2*Math.PI*100));
		
		Complex[] x_array = ComplexOperations.sin(array_complex);
		
		Complex[] out = Signals.fft(1024, x_array);

		double[] Y = ComplexOperations.abs(out);
		System.out.print("\nY = [");
		for(double p : Y)
			System.out.print(", "+p+" ");
		System.out.println("];\n\n");
	}

}
