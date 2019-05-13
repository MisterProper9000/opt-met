#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <iostream>
#define alpha 0.38196601125 // (3 - sqrt(5)) / 2 ; (1 - alpha) = phi
static int counter = 0;

double function(double x) {
	counter++;
	return x * x - 2 * x + exp(-x);
}

double findMinGS(double(*f)(double), double a, double b, double* x_opt, double eps) {
	double min_f = 0;
	double lambda = 0, mu = 0, f_mu, f_lambda;
	int n = 0;
	if (b - a < eps)
		return (a + b) / 2;

	mu = b - (b - a) * alpha;
	lambda = a + (b - a) * alpha;
	printf("lambda_%i = %lf; mu_%i = %lf \n", n, lambda, n, mu);
	f_mu = f(mu);
	f_lambda = f(lambda);

	while (true)
	{
		if (f_lambda >= f_mu)
		{
			a = lambda;
			lambda = mu;
			f_lambda = f_mu;
			mu = b - alpha * (b - a);
			if (b - lambda < eps)
				break;
			f_mu = f(mu);
		}
		else
		{
			b = mu;
			mu = lambda;
			f_mu = f_lambda;
			lambda = a + alpha * (b - a);
			if (mu - a < eps)
				break;
			f_lambda = f(lambda);
		}
		n++;
		printf("lambda_%i = %lf; mu_%i = %lf \n", n, lambda, n, mu);
	}
	*x_opt = (lambda + mu) / 2;
	return f(*x_opt);
}

double findMinLE(double(*f)(double), double a, double b, double* x_opt, double eps, int parts) {
	double min_f = 0;
	double* x = new double[parts + 1];
	double* new_x = new double[parts + 1];
	double* y = new double[parts + 1];
	double* new_y = new double[parts + 1];
	for (int i = 0; i <= parts; i++) {
		x[i] = a + i * (b - a) / parts;
		y[i] = f(x[i]);
	}
	printf("[%lf, %lf]\n", x[0], x[parts]);
	while (x[parts] - x[0] > eps)
	{
		double y_min = INT32_MAX, x_left = 0.0;
		int min_i = 0;
		// находим экстремум в текущей выборке x[]
		for (int i = 0; i <= parts; i++) {
			if (y[i] < y_min) {
				y_min = y[i];
				min_i = i;
			}
		}

		// переопределяем интервал
		double h = 0.0; //величина шага в новом интервале
		if (min_i == 0) { // если минимум оказался на левой границе
			h = (x[min_i + 1] - x[min_i]) / parts;
			x_left = x[min_i]; //новая левая граница
		}
		else if (min_i == parts) { // если минимум оказался на правой границе
			h = (x[min_i] - x[min_i - 1]) / parts;
			x_left = x[min_i - 1];
		}
		else {  // если минимум оказался внутри
			h = (x[min_i + 1] - x[min_i - 1]) / parts;
			x_left = x[min_i - 1];
		}

		for (int i = 0; i <= parts; i++) {
			new_x[i] = x_left + i * h;
			// если такой x_i уже был в прошлой сетке, то зачем его вычислять? найдем его из старых данных
			bool found = false;
			for (int j = 0; j <= parts; j++) {
				if (new_x[i] == x[j]) {
					new_y[i] = y[j];
					found = true;
				}
			}

			// если нужное значение среди старых не обнаружено, придется его вычислить
			if (!found) {
				new_y[i] = f(new_x[i]);
			}
		}

		// обновим сетку и сеточную функцию
		for (int i = 0; i <= parts; i++) {
			x[i] = new_x[i];
			y[i] = new_y[i];
		}
		printf("[%lf, %lf]\n", x[0], x[parts]);
	}
	*x_opt = (x[parts] + x[0]) / 2;
	return f(*x_opt);
}

int main(void)
{
	double a = 1, b = 2, eps = 0.01;
	double(*pf)(double) = &function;
	// Testing GS method
	double x_opt = 0.0;
	double min = findMinGS(pf, a, b, &x_opt, eps);
	printf("GS: Practical count = %i \n", counter);
	printf("GS: Theoretical count = %i \n", (int)floor(log2(eps / (b - a)) / log2(1 - alpha)) + 2);
	printf("GS: Minimum = %lf at x = %lf \n\n", min, x_opt);
	counter = 0;
	// Testing localMin method
	x_opt = 0.0;
	min = findMinLE(pf, a, b, &x_opt, eps, 4);
	printf("LE: Practical count = %d \n", counter);
	printf("LE: Minimum = %lf at x = %lf \n", min, x_opt);
	system("pause");
	return 0;
}