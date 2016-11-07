#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define PI 3.141592653589793
#define a 7.0710678118654755
#define x_min -1200.0
#define x_max 1200.0
#define dt 0.05
#define sigma 10
#define N 4096 //power of 2, good for FFT
#define x0 -100

//declaration of functions
double complex Gaussian_wp(double);
void normalize();
void initialize();
void write_file(FILE *);
void evolve();
double get_transprob();
double V(double);

//global variables
double complex V_exp[N];
double complex T_exp[N];
double complex wave_packet[N];   // x representation of wave packet
double complex p_wave_packet[N]; // p representation of wave packet
double p0;
double x_step;
fftw_plan fft, ifft;    // fft plan for fftw 

int main(int argc, char const **argv)
{
	x_step = (double)(x_max-x_min)/(double)N;
	double E0[argc-1];
	FILE *f;
	int i, k;
	double T, pl, pr;
	if(argv){
		for(i=0; i<argc-1; i++)
			E0[i]=atof(argv[i+1]);
	}
	else{
		printf("needs at least one argument\n");
		return 1;
	}

	for (i=0; i<argc-1; i++){
		printf("Current E = %lf\n", E0[i]);
		p0 = sqrt(2*E0[i]);
		/*the plans have to be made befor initialize(), 
		  because FFTW_MEASURE will override the arrays*/
		fft = fftw_plan_dft_1d(N, wave_packet, p_wave_packet, FFTW_FORWARD, FFTW_MEASURE);
		ifft = fftw_plan_dft_1d(N, p_wave_packet, wave_packet, FFTW_BACKWARD, FFTW_MEASURE);
		initialize();
		// f=fopen("result", "w");	
		// write_file(f);
		k=0;
		while(1){
			k++;
			evolve();
			pl=wave_packet[0]*conj(wave_packet[0]);
			pr=wave_packet[N-1]*conj(wave_packet[N-1]);
			//Make sure the wavepacket will not spread out
			if(pr >= 1e-7 || pl >= 1e-7) break;   
			// if(i%10==0)
			// 	write_file(f);
		}
		printf("Steps take = %d\n", k);
		// fclose(f);
		
		T = get_transprob();
		f = fopen("transprob", "a");
		fprintf(f, "%lf, %lf\n", E0[i], T);
		fclose(f);
	}

	fftw_destroy_plan(fft);
	fftw_destroy_plan(ifft);
	return 0;
}

double V(double x)
{
	if(0<=x && x<=a) 
		return 1;
		// return 4*(a-x)*x/(a*a);
	else return 0;	
}

void initialize()
{
	int i;
	double x, p;
	for(i=0; i<N; i++){
		x = x_min + x_step*(i+1);
		V_exp[i]=cexp(-I*V(x)*dt);
		wave_packet[i] = Gaussian_wp(x);
	}

	for(i=0; i<N; i++){
		// Considering Sampling Theorem && Cyclic displacement
		if(i<N/2) 
			p = i*2.0*PI/x_step/N;
		else      
			p = (i-N)*2.0*PI/x_step/N;
		T_exp[i] = cexp(-I*p*p*dt/2.0);
	}
	normalize();

}

void write_file(FILE *f)
{
	int i;
	double wr,wi;
	for(i=0; i<N; i++){
		wr = creal(wave_packet[i]);
		wi = cimag(wave_packet[i]);
		// if(wr!=0) printf("%lf + %lfi\n", wr, wi);
		// if (wr*wr+wi*wi!=0)printf("!!!\n");
		fprintf(f, "%lf ", wr*wr+wi*wi);
		if(i==N-1) fprintf(f, "\n");
	}
}

double complex Gaussian_wp(double x)
{
	return exp(-(x-x0)*(x-x0)/(2.0*sigma*sigma))*cexp(I*p0*(x-x0));
}

void normalize()
// This function nomalizes the wave packet numerically
{
	double ww[N];
	double prod=0;
	double A;
	int i;
	for(i=0; i<N; i++){
		ww[i] = wave_packet[i]*conj(wave_packet[i]);
	}
	//trapezoidal integral
	for(i=0; i<N-1; i++){
		prod += (ww[i]+ww[i+1])*x_step/2.0;
	}
	A=sqrt(prod);
	for(i=0; i<N; i++){
		wave_packet[i] = wave_packet[i]/A;
	}
}

void evolve()
// Makes one move
{
	int i;
	for(i=0; i<N; i++){
		wave_packet[i] *= V_exp[i];
	}
	fftw_execute(fft);
	for(i=0; i<N; i++){
		p_wave_packet[i] *= T_exp[i];
	}
	fftw_execute(ifft);
	normalize();
}

double get_transprob()
// Calculates the transmission probability by trapzoidal integration
{
	int i;
	double ww1, ww2, x, T=0;
	for(i=0; i<N-1; i++){
		x = x_min+x_step*(i+1);
		if(x<a) continue;
		ww1 = wave_packet[i]*conj(wave_packet[i]);
		ww2 = wave_packet[i+1]*conj(wave_packet[i+1]);
		T += (ww1+ww2)*x_step/2.0;
	}
	return T;
}









