package mathematca;

public class ComplexOperations {

	private ComplexOperations() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * Function that creates a new Complex instance, in the form of z=|z|*exp(theta).
	 * 
	 * @param magnitude magnitude of z
	 * @param angle angle of z
	 * */
	public static Complex getComplexPolar(double magnitude, double angle) {
		double x = magnitude*Math.cos(angle);
		double y = magnitude*Math.sin(angle);
		return new Complex(x, y);
	}
	
	public static double[] abs(Complex[] in) {
		assert(in != null);
		
		if(in.length == 0)
			return null;
		assert(in[0] != null);
		
		double[] out = new double[in.length];
		
		for(int i = 0; i < in.length; i ++) {
			out[i] = in[i].mod();
		}
		return out;
	}
	
	public static Complex[] toComplex(double[] in) {
		Complex[] out = new Complex[in.length];
		for(int i = 0; i < in.length; i ++) {
			out[i] = new Complex(in[i], 0);
		}
		return out;
	}
	
	public static Complex complex(double in) {
		return new Complex(in, 0);
	}
	
	public static Complex[] conj(Complex[] in) {
		Complex[] out = new Complex[in.length];
		for(int i = 0; i < in.length; i ++) {
			out[i] = in[i].conj();
		}
		return out;
	}
	
	public static Complex[] zeros(int size) {
		Complex[] out = new Complex[size];
		for(int i = 0; i < size; i ++) {
			out[i] = new Complex(0, 0);
		}
		return out;
	}
	
	public static Complex[] sin(Complex[] in) {
		Complex[] out = new Complex[in.length];
		for(int i = 0; i < in.length; i ++) {
			out[i] = in[i].sin();
		}
		return out;
	}
	
	public static Complex[] cos(Complex[] in) {
		Complex[] out = new Complex[in.length];
		for(int i = 0; i < in.length; i ++) {
			out[i] = in[i].cos();
		}
		return out;
	}
	
	public static Complex[] scale(Complex[] in, Complex a) {
		Complex[] out = new Complex[in.length];
		for(int i = 0; i < in.length; i ++) {
			out[i] = in[i].times(a);
		}
		return out;
	}
}
