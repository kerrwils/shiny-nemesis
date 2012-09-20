package lib;

import Jama.Matrix;

public class QRDecomposition {
	
	public static void main(String[] args) {
		/*
		 *   0.0608  0.0233 -0.0327
 			-0.0060  0.0056 -0.0098
 			-0.0095  0.0010 -0.0269
		 */
		
		long start, stop;
		
		Matrix m = new Matrix(new double[][]{{12,-51,4},{6,167,-68},{-4,24,-41}});
		
		start = System.currentTimeMillis();
		Matrix sol = invert(m);
		stop = System.currentTimeMillis();
		System.out.println("Time (ms): " + (stop - start));
		sol.print(6, 4);
		
	}
	
	/**
	 * Performs pseudo inversion from the results of the householder decomposition
	 * @param A
	 * @return
	 */
	public static Matrix invert(Matrix A) {
		
		Matrix R = householder(A);
		
		//A+ = ((RT*R)^-1)*AT
		return (((R.transpose()).times(R)).inverse()).times(A.transpose());
		
	}
	
	/**
	 * Preforms QR Decomposition with a householder pivot method on matrix m, stores the
	 * decomposed elements in Q and R.
	 * @param m Matrix to factor
	 */
	public static Matrix householder(Matrix m) {
		
		Matrix Q = null, R, q;
		int mm = m.getRowDimension(), mn = m.getColumnDimension();
		Matrix z = m;
		
		double[] e = new double[mm], x = new double[mm];
		double a;
		
		for (int k = 0; k < mn - 1; k++) {
			
			z = matrix_minor(z,k);
			
			x = z.getColumn(k);
			a = vnorm(x, mm);
				
			if(m.get(k, k) > 0)
				a = -a;

			e = new double[mm];
			e[k] = 1;
			
			vmadd(x,e,a,e,mm);
			vdiv(e, vnorm(e,mm), e, mm);
			q = vmul(e, mm);
			
			if(k==0){
				Q = q;
			}
			else {
				Q = q.times(Q);
			}
		}
		
		R = Q.times(m);
		Q = Q.transpose();
		
		return R;
	}
	
	public static Matrix matrix_minor(Matrix x, int d) {
		
		Matrix m = new Matrix(x.getRowDimension(), x.getColumnDimension());
		for (int i=0; i < d; i++)			
			m.set(i, i, 1);
		for (int i=d; i < x.getRowDimension(); i++)
			for (int j=d; j< x.getColumnDimension();j++)
				m.set(i, j, x.get(i, j));
		return m;
	}
	
	public static double vnorm(double[] x, int n) {
		double sum = 0;
		for(int i=0; i < n; i++)
			sum += x[i]*x[i];
		return Math.sqrt(sum);
	}

	public static void vmadd(double[] a, double[] b, double s, double[] c, int n) {
		for(int i=0;i<n;i++){
			c[i] = a[i] + b[i]*s;
		}
	}
	
	public static void vdiv(double[] x, double d, double[] y, int n) {
		for(int i=0;i<n;i++)
			y[i] = x[i] / d;
	}
	
	public static Matrix vmul(double v[], int n)
	{
		Matrix x = new Matrix(n, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				if (i == j)
					x.set(i,j,(-2*v[i]*v[j]) + 1);
				else
					x.set(i,j,-2*v[i]*v[j]);
			}
	 
		return x;
	}
	
}
