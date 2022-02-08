#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define N 20
int main(void) {
	int n = 0,i=0;
	int a[N] = { 0 };
	FILE* file;
	file = fopen("5.txt", "r");
	if (file == 0) {
		printf("file not found");
		exit(0);
	}
	for (i = 0; i < N; i++) {
		fscanf(file, "%d ", &a[i]);
		printf( "%d ", a[i]);
	}
	printf("\n");
	for (i = 0; i < N; i++) {
		if (a[i] > 0) {
			a[i] = a[i] + a[i] % 10;
		}
		if (a[i] < 0) {
			n = a[i]%10;
			n = abs(n);
			a[i] = a[i] - n;
		}
		printf("%d ", a[i]);
	}
	return 0;
}