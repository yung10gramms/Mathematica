package mathematca;

import mathematca.matrix.LinearAlgebra;

public class Signals {

	private Signals() {}
	
	/**
	 * Implementation of the fft algorithm. 
	 * 
	 * @param n size of output array (must be pow of 2)
	 * @param arr the input array 
	 * */
	public static Complex[] fft(int n, Complex[] arr) {
		Complex[] paddedArr = ComplexOperations.zeros(n);
	    System.arraycopy(arr, 0, paddedArr, 0, arr.length);
	    arr = paddedArr;
	    
		if(n == 1)
			return new Complex[] {arr[0]};
		Complex[][] oddev = getEvenGetOdds(arr);

		Complex[] even = oddev[0];
		Complex[] odds = oddev[1];
		
		Complex[] d = fft(n/2, odds);
		Complex[] e = fft(n/2, even);
		
		
		Complex[] y = new Complex[n];
		for(int k = 0; k < n/2; k ++) {
			Complex omega = ComplexOperations.getComplexPolar(1, 2*Math.PI*k/n);
			
			y[k] = calculateCurrentTerm(omega, d[k], e[k]);
			y[k+n/2] = calculateOppositeTerm(omega, d[k], e[k]);	
		}
		
		return y;
	}
	
	/**
	 * Implementation of the ifft algorithm. 
	 * 
	 * @param n size of output array (must be pow of 2)
	 * @param arr the input array 
	 * */
	public static Complex[] ifft(int n, Complex[] arr) {
		return fft(n, ComplexOperations.conj(arr));
	}

	/**
	 * Function that splits an array into two sub-arrays, one with all elements of
	 * the array that have even indeces and one with the odd ones. <p>This function will 
	 * not overwrite/change the elements of the array. 
	 * @param array the input array
	 * */
	private static Complex[][] getEvenGetOdds(Complex[] array) {
		Complex[][] oddsEvens = new Complex[2][array.length/2];
		
		for(int i = 0; i < array.length; i ++) {
			if(i%2 == 0)
				oddsEvens[0][i/2] = array[i];
			else
				oddsEvens[1][(i-1)/2] = array[i];
		}
		return oddsEvens;
	}
	
	public static Complex[] conv(Complex[] x, Complex[] y) {
		int n = findClosestlog(x.length+y.length-1);
		Complex[] X = fft(n, x);
		Complex[] Y = fft(n, y);
		
		Complex[] Z = LinearAlgebra.arrmul(X, Y);
		
		return ifft(n, Z);
	}
	
	private static int findClosestlog(int x) {
		int xu = (int) Math.ceil(log2(x));
		return (int) Math.pow(2, xu);
	}
	
	private static double log2(double x) {
		return Math.log(x)/Math.log(2);
	}
	
	/*
	 * omega.times(d[k]).plus(e[k]);
	 * e[k] + omega*d[k]
	 * */
	private static final Complex calculateCurrentTerm(Complex w, Complex d, Complex e) {
		return w.times(d).plus(e);
	}
	
	private static final Complex calculateOppositeTerm(Complex w, Complex d, Complex e) {
		return e.minus(w.times(d));
	}
}
