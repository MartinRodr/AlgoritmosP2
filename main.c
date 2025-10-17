#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <stdbool.h>

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

void inicializar_ascendente(int v[], int n) {
	for (int i = 0; i < n; i++) {
		v[i] = i;
	}
}

void inicializar_descendente(int v[], int n) {
	for (int i = 0; i < n; i++) {
		v[i] = n - i;
	}
}

bool ordenado(int v[], int n){
	for (int i = 1; i < n; i++){
		if (v[i - 1] > v[i]) { return false; }
	}
	return true;
}

int *hibbard(int n, int *m){
	int k = 1, cnt = 0, *incr;
	while ((int)pow(2, k) - 1 <= n){
		cnt++;
		k++;
	}
	if (cnt == 0) { *m = 0; return NULL; }
	incr = malloc(cnt * sizeof(int));
	if (!incr) { perror("malloc"); exit(1); }
	*m = cnt;
	for (int i = 0; i < cnt; i++){
		incr[i] = (int)pow(2, cnt - i) - 1;
	}
	return incr;
}

int *knuth(int n, int *m){
	int k = 1, cnt = 0, *incr;
	while (((int)pow(3, k) - 1) / 2 <= n) {
		cnt++;
		k++;
	}
	if (cnt == 0) { *m = 0; return NULL; }
	incr = malloc(cnt * sizeof(int));
	if (!incr) { perror("malloc"); exit(1); }
	*m = cnt;
	for (int i = 0; i < cnt; i++) {
		incr[i] = ((int)pow(3, cnt - i) - 1) / 2;
	}
	return incr;
}

int *sedgewick(int n, int *m) {
	int k = 1, cnt = 0, *incr;
	while ((int)(pow(4, k) + (3 * pow(2, k - 1)) + 1) <= n) {
		cnt++;
		k++;
	}
	if (cnt == 0) {
        *m = 0; return NULL;
    }
    incr = malloc(cnt * sizeof(int));
    if (!incr) { perror("malloc"); exit(1); }
    *m = cnt;
    for (int i = 0; i < cnt; i++) 
    	incr[i] = (int)(pow(4, cnt - i) + (3 * pow(2, cnt - i - 1)) + 1);
    if (incr[cnt - 1] != 1){
    	incr = realloc(incr, (cnt + 1) * sizeof(int));
    	if (!incr){
    		perror("realloc"); exit(1);
    	}
    	incr[cnt] = 1;
    	*m = cnt + 1;
    }
    return incr;
}

int *ciura(int n, int *m) {
	int ciura[] = {1, 4, 10, 23, 57, 132, 301, 701, 1750};
	int ciura_long = sizeof(ciura) / sizeof(ciura[0]); 
	int last = 1, cnt = 0, next = 0, *incr;
	// Contar base de ciura
	while ((cnt < ciura_long && ciura[cnt] <= n)) {
		last = ciura[cnt++]; }
	// Contar extendidos
	next = last;
	while (next < n){
		next = (int)round(last * 2.25);
		if (next >= n) break;
		last = next;
		cnt++;
	}
	if (cnt == 0) { *m = 0; return NULL; }
	incr = malloc(cnt * sizeof(int));
    if (!incr) { perror("malloc"); exit(1); }
    // Secuencia base de ciura
	for (int i = 0; i < ciura_long && ciura[i] <= n; i++){
		incr[cnt - i - 1] = ciura[i];
	}
	// Secuencia extendida
	last = incr[ciura_long - 1];
    for (int j = ciura_long; j < cnt; j++) {
        next = (int)round(last * 2.25);
        if (next >= n) break;
        incr[cnt - j - 1] = next;
        last = next;
    }
    *m = cnt;
	return incr;
}

/*int *ciura(int n, int *m) {
    int ciura[] = {1, 4, 10, 23, 57, 132, 301, 701, 1750};
    int ciura_len = sizeof(base) / sizeof(base[0]);
    int *incr = malloc(100 * sizeof(int)); // capacidad inicial
    if (!incr) { perror("malloc"); exit(1); }
    int cnt = 0;

    for (int i = 0; i < ciura_len && ciura[i] < n; i++)
        incr[cnt++] = ciura[i];

    // 2️⃣ Extender mientras el último sea menor que n
    int last = incr[cnt - 1];
    while (last < n) {
        int next = (int)round(last * 2.25);
        if (next >= n) break;
        incr[cnt++] = next;
        last = next;
    }

    // 3️⃣ Invertir el orden (Shell sort usa gaps de mayor a menor)
    for (int i = 0; i < cnt / 2; i++) {
        int tmp = incr[i];
        incr[i] = incr[cnt - i - 1];
        incr[cnt - i - 1] = tmp;
    }

    *m = cnt;
    return incr;
}*/

void ord_ins(int v[], int n) {
	int x, j;
	for (int i = 1; i < n; i++){
		x = v[i];
		j = i - 1;
		while (j >= 0 && v[j] >= x){
			v[j + 1] = v[j];
			j--;
		}
		v[j + 1] = x;
	}
}

void ord_shell(int v[], int n, int incr[], int m){
	// incr es el vector de incrementos y el ultimo debe ser 1
	int h, x, j;
	for (int k = 0; k < m; k++){
		h = incr[k];
		for (int i = h; i < n; i++) {
			x = v[i];
			j = i;
			while (j >= h && v[j - h] > x){
				v[j] = v[j - h];
				j -= h;
			}
			v[j] = x;
		}
	}
}

void contarTiempoInsAscend(int n, int k, int m, int exp, double confianza){
	int conf, *v;
	double ta, tb, t, t1, t2, cte;
	printf("Algoritmo Insercion: Vector Ascendente\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n^1.8", "t (n) / n^2", "t (n) / n^2.2");
	for (int i = 0; i < m; i++) {
		conf = 0;
		v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		inicializar_ascendente(v, n);
		ta = microsegundos();
		ord_ins(v, n);
		tb = microsegundos();
		t = tb - ta;
		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				inicializar_ascendente(v, n);
				ord_ins(v, n);
			}
			tb = microsegundos();
			t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				inicializar_ascendente(v, n);
			}
			tb = microsegundos();
			t2 = tb - ta;
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
	printf("Cte = %lf\n\n", cte);
}

void contarTiempoInsDescend(int n, int k, int m, int exp, double confianza){
	int conf, *v;
	double ta, tb, t, t1, t2, cte;
	printf("Algoritmo Insercion: Vector Descendente\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n^1.8", "t (n) / n^2", "t (n) / n^2.2");
	for (int i = 0; i < m; i++) {
		conf = 0;
		v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		inicializar_descendente(v, n);
		ta = microsegundos();
		ord_ins(v, n);
		tb = microsegundos();
		t = tb - ta;
		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				inicializar_descendente(v, n);
				ord_ins(v, n);
			}
			tb = microsegundos();
			t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				inicializar_descendente(v, n);
			}
			tb = microsegundos();
			t2 = tb - ta;
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
	printf("Cte = %lf\n\n", cte);
}

void contarTiempoInsDesord(int n, int k, int m, int exp, double confianza){
	int conf, *v;
	double ta, tb, t, t1, t2, cte;
	printf("Algoritmo Insercion: Vector Desordenado\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n^1.8", "t (n) / n^2", "t (n) / n^2.2");
	for (int i = 0; i < m; i++) {
		conf = 0;
		v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		aleatorio(v, n);
		ta = microsegundos();
		ord_ins(v, n);
		tb = microsegundos();
		t = tb - ta;
		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				aleatorio(v, n);
				ord_ins(v, n);
			}
			tb = microsegundos();
			t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				aleatorio(v, n);
			}
			tb = microsegundos();
			t2 = tb - ta;
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
	printf("Cte = %lf\n\n", cte);
}

void contarTiempoShellHibb(int n, int k, int m, int exp, double confianza){
	int conf, *v, *incr, mm;
	double ta, tb, t, t1, t2, cte;
	printf("Algoritmo Shell: Hibbard\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n*log(n)", "t (n) / n^1.18", "t (n) / n^1.3");
	for (int i = 0; i < m; i++) {
		conf = 0;
		v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		incr = hibbard(n, &mm);
		printf("\n");	
		aleatorio(v, n);
		ta = microsegundos();
		ord_shell(v, n, incr, mm);
		tb = microsegundos();
		t = tb - ta;
		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				aleatorio(v, n);
				ord_shell(v, n, incr, mm);
			}
			tb = microsegundos();
			t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				aleatorio(v, n);
			}
			tb = microsegundos();
			t2 = tb - ta;
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
		printf("%18lf", t / ((double)n * log(n)));
		printf("%18.7lf", t / pow((double)n, 1.18));
		cte += t / pow((double)n, 1.18);
		printf("%18lf", t / pow((double)n, 1.3));
		printf("\n");
		free(v);
		free(incr);
		n *=exp;
	}
	cte /= (double)m;
	printf("Cte = %lf\n\n", cte);
}

void contarTiempoShellKnuth(int n, int k, int m, int exp, double confianza){
	int conf, *v, *incr, mm;
	double ta, tb, t, t1, t2, cte;
	printf("Algoritmo Shell: Knuth\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n^1.1", "t (n) / n^1.163", "t (n) / n^1.3");
	for (int i = 0; i < m; i++) {
		conf = 0;
		v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		incr = knuth(n, &mm);
		printf("\n");
		aleatorio(v, n);
		ta = microsegundos();
		ord_shell(v, n, incr, mm);
		tb = microsegundos();
		t = tb - ta;
		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				aleatorio(v, n);
				ord_shell(v, n, incr, mm);
			}
			tb = microsegundos();
			t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				aleatorio(v, n);
			}
			tb = microsegundos();
			t2 = tb - ta;
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
		printf("%18.7lf", t / pow((double)n, 1.1));
		printf("%18lf", t / pow((double)n, 1.163));
		cte += t / pow((double)n, 1.1635);
		printf("%18lf", t / pow((double)n, 1.3));
		printf("\n");
		free(v);
		free(incr);
		n *=exp;
	}
	cte /= (double)m * 10;
	printf("Cte = %lf\n\n", cte);
}

void contarTiempoShellSedg(int n, int k, int m, int exp, double confianza){
	int conf, *v, *incr, mm;
	double ta, tb, t, t1, t2, cte;
	printf("Algoritmo Shell: Sedgewick\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / n^1.1", "t (n) / n^1.118", "t (n) / n^1.3");
	for (int i = 0; i < m; i++) {
		conf = 0;
		v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		incr = sedgewick(n, &mm);
		printf("\n");
		aleatorio(v, n);
		ta = microsegundos();
		ord_shell(v, n, incr, mm);
		tb = microsegundos();
		t = tb - ta;
		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				aleatorio(v, n);
				ord_shell(v, n, incr, mm);
			}
			tb = microsegundos();
			t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				aleatorio(v, n);
			}
			tb = microsegundos();
			t2 = tb - ta;
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
		printf("%18lf", t / (pow((double)n, 1)));
		printf("%18lf", t / pow((double)n, 1.118));
		cte += t / (pow((double)n, 1.118));
		printf("%18.7lf", t / pow((double)n, 1.3));
		printf("\n");
		free(v);
		free(incr);
		n *=exp;
	}
	cte /= (double)m;
	printf("Cte = %lf\n\n", cte);
}

void contarTiempoShellCiura(int n, int k, int m, int exp, double confianza){
	int conf, *v, *incr, mm;
	double ta, tb, t, t1, t2, cte;
	printf("Algoritmo Shell: Ciura\n");
	printf("%10s%18s%18s%18s%18s\n", 
		"n", "t (n)", "t (n) / log", "t (n) / n^1.21", "t (n) / n^4/3");
	for (int i = 0; i < m; i++) {
		conf = 0;
		v = malloc(n * sizeof(int));
		if (!v){
			perror("malloc");
			exit(EXIT_FAILURE);
		}
		incr = ciura(n, &mm);
		for (int j = 0; j < mm; j++) {
			printf("%d  ", incr[j]);
		}
		printf("\n");
		aleatorio(v, n);
		ta = microsegundos();
		ord_shell(v, n, incr, mm);
		tb = microsegundos();
		t = tb - ta;
		if (t < confianza){
			ta = microsegundos();
			for (int i = 0; i < k; i++){
				aleatorio(v, n);
				ord_shell(v, n, incr, mm);
			}
			tb = microsegundos();
			t1 = tb - ta;
			ta = microsegundos();
			for (int j = 0; j < k; j++){
				aleatorio(v, n);
			}
			tb = microsegundos();
			t2 = tb - ta;
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
		printf("%18lf", t / ((double)n * pow(log(n), 1.5)));
		printf("%18lf", t / pow((double)n, 1));
		//printf("%18lf", t / pow((double)n, 1.21));
		cte += t / pow((double)n, 1.2124);
		printf("%18.7lf", t / pow((double)n, 1.1));
		printf("\n");
		free(v);
		free(incr);
		n *=exp;
	}
	cte /= (double)m;
	printf("Cte = %lf\n\n", cte);
}

void test_ins(){
	int n = 15;
	int v[n];
	aleatorio(v, n);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	printf("\n");
	ord_ins(v, n);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	printf("\n\n");
}

void test_shell_hibbard(){
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	if (!ordenado(v, n)){
		incr = hibbard(n, &m);
		ord_shell(v, n, incr, m);
	}
	printf("\n");
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	printf("\n hibbard: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);
	}
	printf("\n\n");
	free(incr);
}

void test_shell_knuth(){
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	if (!ordenado(v, n)){
		incr = knuth(n, &m);
		ord_shell(v, n, incr, m);
	}
	printf("\n");
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);

	}
	printf("\n knuth: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);
	
	}
	printf("\n\n");
	free(incr);
}

void test_shell_sedgewick(){
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	if (!ordenado(v, n)){
		incr = sedgewick(n, &m);
		ord_shell(v, n, incr, m);
	}
	printf("\n");
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	printf("\n sedgewick: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);
	
	}
	printf("\n\n");
	free(incr);
}

void test_shell_ciura(){
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	if (!ordenado(v, n)){
		incr = ciura(n, &m);
		ord_shell(v, n, incr, m);
	}
	printf("\n");
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);
	}
	printf("\n ciura: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);		
	}
	printf("\n\n");
	free(incr);
}

int main(void){
	int k = 1000, n = 500, m = 11, exp = 2;
	double confianza = 500.00;
	inicializar_semilla();
	//contarTiempoInsAscend(n, k, m, exp, confianza);
	//contarTiempoInsDescend(n, k, m, exp, confianza);
	//contarTiempoInsDesord(n, k, m, exp, confianza);

	contarTiempoShellHibb(n, k, m, exp, confianza);
	contarTiempoShellKnuth(n, k, m, exp, confianza);
	contarTiempoShellSedg(n, k, m, exp, confianza);
	contarTiempoShellCiura(n, k, m, exp, confianza);

	/*test_ins();
	test_shell_hibbard();
	test_shell_knuth();
	test_shell_sedgewick();
	test_shell_ciura();*/
	return 0;
}