#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

/*
 * Ali Abu-afash Nayef
 * Martín Rodríguez Arévalo
 * Alejandro Patiño Jaramillo
*/

double microsegundos() { /* obtiene la hora del sistema en microsegundos */
	struct timeval t;
	if (gettimeofday(&t, NULL) < 0 )
		return 0.0;
	return (t.tv_usec + t.tv_sec * 1000000.0);
}

void inicializar_semilla() {
	srand(time(NULL));
/* se establece la semilla de una nueva serie de enteros pseudo-aleatorios */
}

void aleatorio(int v[], int n) {
	int i, m=2*n+1;
	for (i=0; i < n; i++)
		v[i] = (rand() % m) - n;
	/* se generan números pseudoaleatorio entre -n y +n */
}

int hibbard(int v[], n){
	int m = 0;
	for (int k = 0; k < n; k++){
		int h = (int)pow(2, k + 1) - 1;
		if (h > n) break;
		v[m++] = h;
	}
	return m;
}

int knuth(int v[], int n){
	for (int k = 0; k < n; k++){
		v[k] = (int)(pow(3, k + 1) - 1) / 2;
	}
}

int sedgewick(int v[]) {
	v[0] = 1;
	for (int k = 1; k < n; k++){
		v[k] = (int)pow(4, k) + (3 * pow(2, k - 1)) + 1;
	}
}

int ciura(int v[], int n) {
	int ciura[] = {1, 4, 10, 23, 57, 132, 301, 701, 1750};
	int ciura_long = sizeof(ciura) / sizeof(ciura[0]);
	int m = 0;
	for (int k = 0; k < ciura_long && ciura[k] <= n; k++){
		v[m++] = ciura[k];
	}
	while(v[m - 1] < n){
		int sig = (int)(v[m - 1] * 2.25);
		if (sig > n) break;
		v[m++] = sig;
	}
	return m;
}

void ord_ins(int v[], int n) {
	for (int i = 1; i < n; i++){
		int x = v[i];
		int j = i - 1;
		while (j >= 0 && v[j] >= x){
			v[j + 1] = v[j];
			j--;
		}
		v[j + 1] = x;
	}
}

void ord_shell(int v[], int n, int incr[], int m){
	// incr es el vector de incrementos y el ultimo debe ser 1
	for (int k = 0; k < m; k++){
		int h = incr[k];
		for (int i = h + 1; i < n; i++) {
			int x = v[i];
			int j = i;
			while (j > h && v[j - h] > x){
				v[j] = v[j - h];
				j -= h;
			}
			v[j] = x;
		}
	}
}


/*void contarTiempoAlg1(int n, int k, int m, int exp, double confianza){
	printf("Algoritmo 1:\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n^1.8", "t (n) / n^2", "t (n) / n^2.2");
	double cte = 0;
	for (int i = 0; i < m; i++) {
		int conf = 0;
		int *v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		inicializar_vector(v, n);
		double ta = microsegundos();
		suma1(v, n);
		double tb = microsegundos();
		double t = tb - ta;

		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				inicializar_vector(v, n);
				suma1(v, n);
			}
			tb = microsegundos();
			double t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				inicializar_vector(v, n);
			}
			tb = microsegundos();
			double t2 = tb - ta;
			t = (t1 - t2) / k;
			conf = 1;
		}
		if (conf == 1){
			printf("(*)");
			printf("%7d", n);

		} else {
			printf("%10d", n);
		}
		printf("%18lf", t);
		printf("%18lf", t / pow((double)n, 1.8));
		printf("%18lf", t / pow((double)n, 2));
		cte += t / pow((double)n, 2);
		printf("%18.7lf", t / pow((double)n, 2.2));
		printf("\n");
		free(v);
		n *=exp;
	}
	cte /= (double)m;
	printf("Cte = %lf\n", cte);
}

void contarTiempoAlg2(int n, int k, int m, int exp, double confianza){
	printf("Algoritmo 2:\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n^0.8", "t (n) / n", "t (n) / n*log(n)");
	double cte = 0;
	for (int i = 0; i < m; i++) {
		int conf = 0;
		int *v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		inicializar_vector(v, n);
		double ta = microsegundos();
		suma2(v, n);
		double tb = microsegundos();
		double t = tb - ta;

		if (t < confianza){
			conf = 1;
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				inicializar_vector(v, n);
				suma2(v, n);
			}
			tb = microsegundos();
			double t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				inicializar_vector(v, n);
			}
			tb = microsegundos();
			double t2 = tb - ta;
			t = (t1 - t2) / k;
		}
		if (conf == 1){
			printf("(*)");
			printf("%7d", n);

		} else {
			printf("%10d", n);
		}
		printf("%18lf", t);
		printf("%18lf", t / pow((double)n, 0.8));
		printf("%18lf", t / (double)n);
		cte += t / (double)n;
		printf("%18.7lf", t / ((double)n * log(n)));
		printf("\n");
		free(v);
		n *=exp;
	}
	cte /= (double)m;
	printf("Cte = %lf\n", cte);
}
*/

void test_ins(){
	int n = 15;
	int v[n];
	aleatorio(v, n);
	ord_ins(v, n);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	printf("\n");
}
void test_shell(){
	int n = 15, m = 5;
	int v[n], incr[m];
	aleatorio(v, n);
	ciura(incr, m);
	ord_shell(v, n, incr, m);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);

	}
	printf("\n");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);
	
	}
	printf("\n");
}

int main(void){
	inicializar_semilla();
	test_ins();
	test_shell();
	return 0;
}