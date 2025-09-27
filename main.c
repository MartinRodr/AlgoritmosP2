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

bool ordenado(int v[], int n){
	for (int i = 1; i < n; i++){
		if (v[i - 1] > v[i]) { return false; }
	}
	return true;
}

int *hibbard(int n, int *m){
	int k = 1, cnt = 0;
	while ((int)pow(2, k) - 1 <= n){
		cnt++;
		k++;
	}
	if (cnt == 0) { *m = 0; return NULL; }
	int *incr = malloc(cnt * sizeof(int));
	if (!incr) { perror("malloc"); exit(1); }
	*m = cnt;
	for (int i = 0; i < cnt; i++){
		incr[i] = (int)pow(2, cnt - i) - 1;
	}
	return incr;
}

int *knuth(int n, int *m){
	int k = 1, cnt = 0;
	while (((int)pow(3, k) - 1) / 2 <= n) {
		cnt++;
		k++;
	}
	if (cnt == 0) { *m = 0; return NULL; }
	int *incr = malloc(cnt * sizeof(int));
	if (!incr) { perror("malloc"); exit(1); }
	*m = cnt;
	for (int i = 0; i < cnt; i++) {
		incr[i] = ((int)pow(3, cnt - i) - 1) / 2;
	}
	return incr;
}
 /*
 * REVISAR SEDGEWICK NO FUNCIONA
 *
 */
int *sedgewick(int n, int *m) {
	int k = 1, cnt = 0;
	while ((int)(pow(4, k) + (3 * pow(2, k - 1)) + 1) <= n) {
		cnt++;
		k++;
	}
	if (cnt == 0) {
        *m = 0; return NULL;
    }
    int *incr = malloc(cnt * sizeof(int));
    if (!incr) { perror("malloc"); exit(1); }
    *m = cnt;
    for (int i = 0; i < cnt; i++) {
    	if (i == cnt - 1) { incr[i] = 1; }
    	else{
    		incr[i] = (int)(pow(4, cnt - i) + (3 * pow(2, cnt - i - 1)) + 1);
    	}
	}
    return incr;
}

int *ciura(int n, int *m) {
	int ciura[] = {1, 4, 10, 23, 57, 132, 301, 701, 1750};
	int ciura_long = sizeof(ciura) / sizeof(ciura[0]);
	int last = 1, cnt = 0;
	while ((cnt < ciura_long && ciura[cnt] <= n)) {
		last = ciura[cnt++];
		//cnt++;
	}
	while (last < n){
		int next = (int)round(last * 2.25);
		if (next > n) break;
		last = next;
		cnt++;
	}
	if (cnt == 0) { *m = 0; return NULL; }
	int *incr = malloc(cnt * sizeof(int));
    if (!incr) { perror("malloc"); exit(1); }
    *m = cnt;
	for (int i = 0; i < ciura_long && ciura[i] <= n; i++){
		incr[cnt - i - 1] = ciura[i];
	}
	last = (cnt > 0) ? incr[cnt-1] : 1;
    for (int j = ciura_long; j < cnt; j++) {
        last = (int)round(last * 2.25);
        incr[cnt - j - 1] = last;
    }
	return incr;	
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
		for (int i = h; i < n; i++) {
			int x = v[i];
			int j = i;
			while (j >= h && v[j - h] > x){
				v[j] = v[j - h];
				j -= h;
			}
			v[j] = x;
		}
	}
}

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
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	if (!ordenado(v, n)){
		incr = ciura(n, &m);
		ord_shell(v, n, incr, m);
	}
	for (int i = 0; i < n; i++){
			printf("%d  ", v[i]);
	}
	printf("\n ciura: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);		
	}
	printf("\n");
	free(incr);
}

void p1(){
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	if (!ordenado(v, n)){
		incr = hibbard(n, &m);
		ord_shell(v, n, incr, m);
	}
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);

	}
	printf("\n hibbard: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);
	
	}
	printf("\n");
	free(incr);
}

void p2(){
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	if (!ordenado(v, n)){
		incr = sedgewick(n, &m);
		ord_shell(v, n, incr, m);
	}
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);

	}
	printf("\n hawd2: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);
	
	}
	printf("\n");
	free(incr);
}

void p3(){
	int n = 15, m;
	int v[n], *incr;
	aleatorio(v, n);
	if (!ordenado(v, n)){
		incr = knuth(n, &m);
		ord_shell(v, n, incr, m);
	}
	for (int i = 0; i < n; i++){
		printf("%d  ", v[i]);

	}
	printf("\n gege3: ");
	for (int i = 0; i < m; i++){
		printf("%d  ", incr[i]);
	
	}
	printf("\n");
	free(incr);
}

int main(void){
	inicializar_semilla();
	test_ins();
	test_shell();
	p1();
	p2();
	p3();
	return 0;
}