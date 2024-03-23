/*
 * HW8.c
 * Fibonacci
 * Ashton Cole
 * COE 379L: Parallel Computing
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

long fib(int n) {
	if (n < 2) return n;
	else {
		long f1, f2;

#pragma omp taskgroup
		{
			// Each task gets its own n
#pragma omp task shared(f1) firstprivate(n) // depend(out:f1)
			f1 = fib(n - 1);

			// Each task gets its own n
#pragma omp task shared(f2) firstprivate(n) // depend(out:f2)
			f2 = fib(n - 2);
		}

#pragma omp taskwait
		return f1 + f2;
	}
}

int main(int argc, char **argv) {
	int n = 20;

	if (argc == 2) n = atoi(argv[1]);

#pragma omp parallel
#pragma omp single
	printf("Fib(%d)=%ld\n", n, fib(n));

	return 0;
}
