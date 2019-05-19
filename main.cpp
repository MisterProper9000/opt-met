//////////////////////////////////////////////////////////////
// LOCAL INCLUDES
//////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <conio.h>
#include "MatrixWork.h"
#include <iostream>

#define EPSILON 1.e-3
static void print(double** in) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf("%lf ", in[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

static double* scalarMulVect(double scalar, double* vect) {
	double* res = new double[3];
	for (int i = 0; i < 3; i++) {
		res[i] = vect[i] * scalar;
	}
	return res;
}

static double* vMinus(double* v1, double* v2) {
	double* v = new double[3];
	for (int i = 0; i < 3; i++) {
		v[i] = v1[i] - v2[i];
	}
	return v;
}

static double dot(double *a, double *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static double Hdot(double *x, double** H, double *y) {
	double *t = new double[3];
	for (int i = 0; i < 3; i++)
	{
		t[i] = 0;
		for (int j = 0; j < 3; j++)
			t[i] += x[j] * H[j][i];
	}
	return dot(t, y);
}


static double* normalize(double *in, double** H) {
	double n = 1.0 / sqrt(Hdot(in, H, in));
	return scalarMulVect(n, in);
}

static double** classic_gram_schmidt(double** in, double** H)
{
	double **out = new double*[3];
	for (int i = 0; i < 3; i++)
		out[i] = new double[3];
	out[0] = normalize(in[0], H);
	out[1] = normalize(vMinus(in[1], scalarMulVect(Hdot(in[1], H, out[0]), out[0])), H);

	double* v1 = scalarMulVect(Hdot(in[2], H, out[0]), out[0]);
	double* v2 = scalarMulVect(Hdot(in[2], H, out[1]), out[1]);
	double* v = vMinus(in[2], v1);
	out[2] = normalize(vMinus(v, v2), H);
	/*std::cout << "****" << Hdot(out[2], H, out[1]) <<std::endl;
	std::cout << "****" << Hdot(out[0], H, out[1]) << std::endl;
	std::cout << "****" << Hdot(out[0], H, out[2]) << std::endl;*/

	return out;
}



double Function(double x1, double x2, double x3)
{
	return x1*x1 + x2*x2 + x3*x3 - x1*x2 - x2*x3 - 3 * x1 - 2 * x2 - x3;
}

void Gradient(double x[], double y[])
{
	y[0] = 2 * x[0] - x[1] - 3;
	y[1] = -x[0] + 2 * x[1] - x[2] - 2;
	y[2] = -x[1] + 2 * x[2] - 1;
}

double** Hessian(double x[] = nullptr)
{
	double **H = new double*[3];
	for (int i = 0; i < 3; i++)
		H[i] = new double[3];

	H[0][0] = 2;
	H[0][1] = H[1][0] = -1;
	H[0][2] = H[2][0] = 0;
	H[1][1] = 2;
	H[1][2] = H[2][1] = -1;
	H[2][2] = 2;
	return H;
}


double g(double x1, double x2, double x3)
{
	return (3 / 2.) * x1*x1 + (3 / 2.) * x2*x2 + (3 / 2.) * x3 * x3 - 2 * x1 * x2 - 2 * x2 * x3 - x1 - 2 * x2 - 3 * x3 + (x2 - 14)*(x2 - 14)*(x2 - 14)*(x2 - 14);
}

void gradG(double x[], double y[])
{
	y[0] = 3 * x[0] - 2 * x[1] - 1;
	y[1] = -2 * x[0] + 3 * x[1] - 2 * x[2] - 2 + 4 * (x[1] - 14)*(x[1] - 14)*(x[1] - 14);
	y[2] = -2 * x[1] + 3 * x[2] - 3;
}

double** hessianG(double x[])
{
	double **H = new double*[3];
	for (int i = 0; i < 3; i++)
		H[i] = new double[3];
	H[0][0] = 3;
	H[0][1] = H[1][0] = -2;
	H[0][2] = H[2][0] = 0;
	H[1][1] = 3 + 12 * (x[1] - 14) * (x[1] - 14);
	H[1][2] = H[2][1] = -2;
	H[2][2] = 3;
	return H;
}

void Gradient2(double x[], double y[], bool isModification, int recalculate)
{
	double **H, **invH, grad[3], det = 0, **tmp;
	H = new double*[3];
	invH = new double*[3];

	for (int i = 0; i < 3; i++)
	{
		H[i] = new double[3];
		invH[i] = new double[3];
	}

	if (recalculate)
	{
		if (!isModification) {
			tmp = Hessian();
		}
		else {
			tmp = hessianG(x);
		}
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
				H[i][j] = tmp[i][j];
		}
		inverse(H, invH, 3);
	}
	if (isModification == false) {
		Gradient(x, grad);
	}
	else {
		gradG(x, grad);
	}
	y[0] = (invH[0][0] * grad[0] + invH[0][1] * grad[1] + invH[0][2] * grad[2]);
	y[1] = (invH[1][0] * grad[0] + invH[1][1] * grad[1] + invH[1][2] * grad[2]);
	y[2] = (invH[2][0] * grad[0] + invH[2][1] * grad[1] + invH[2][2] * grad[2]);
}

double Step_Newton(double x[], double y[], double *f, double norm)
{
	double newtonStep = 1;
	double x_k_1[3];
	x_k_1[0] = x[0] - newtonStep * y[0];
	x_k_1[1] = x[1] - newtonStep * y[1];
	x_k_1[2] = x[2] - newtonStep * y[2];
	double fk1 = Function(x_k_1[0], x_k_1[1], x_k_1[2]);
	double fk2 = Function(x[0], x[1], x[2]);
	double mult;
	double grad[3];
	Gradient(x, grad);
	mult = -0.5 * (grad[0] * y[0] + grad[1] * y[1] + grad[2] * y[2]);
	while ((fk1 - fk2) > (mult * newtonStep * EPSILON))
	{
		newtonStep /= 2;
		x_k_1[0] = x[0] - newtonStep * y[0];
		x_k_1[1] = x[1] - newtonStep * y[1];
		x_k_1[2] = x[2] - newtonStep * y[2];
		fk1 = Function(x_k_1[0], x_k_1[1], x_k_1[2]);
		fk2 = Function(x[0], x[1], x[2]);
		Gradient(x, grad);
		mult = -0.5 * (grad[0] * y[0] + grad[1] * y[1] + grad[2] * y[2]);
	}

	x[0] -= newtonStep * y[0];
	x[1] -= newtonStep * y[1];
	x[2] -= newtonStep * y[2];

	*f = Function(x[0], x[1], x[2]);
	return newtonStep;
}
double Step_Newton_Modification(double x[], double y[], double *f, double norm)
{
	double newtonStep = 1;
	double x_k_1[3];
	x_k_1[0] = x[0] - newtonStep * y[0];
	x_k_1[1] = x[1] - newtonStep * y[1];
	x_k_1[2] = x[2] - newtonStep * y[2];
	double fk1 = g(x_k_1[0], x_k_1[1], x_k_1[2]);
	double fk2 = g(x[0], x[1], x[2]);
	double mult;
	double grad[3];
	gradG(x, grad);
	mult = -0.5 * (grad[0] * y[0] + grad[1] * y[1] + grad[2] * y[2]);
	while ((fk1 - fk2) > (mult * newtonStep * EPSILON))
	{
		newtonStep /= 2;
		x_k_1[0] = x[0] - newtonStep * y[0];
		x_k_1[1] = x[1] - newtonStep * y[1];
		x_k_1[2] = x[2] - newtonStep * y[2];
		fk1 = g(x_k_1[0], x_k_1[1], x_k_1[2]);
		fk2 = g(x[0], x[1], x[2]);
		gradG(x, grad);
		mult = -0.5 * (grad[0] * y[0] + grad[1] * y[1] + grad[2] * y[2]);
	}

	x[0] -= newtonStep * y[0];
	x[1] -= newtonStep * y[1];
	x[2] -= newtonStep * y[2];

	*f = g(x[0], x[1], x[2]);
	return newtonStep;
}

double* vector_dot_matrix(double* vec, double** Matrix)
{
	double result[3];
	for (int i = 0; i < 3; i++)
	{
		result[i] = 0;
		for (int j = 0; j < 3; j++)
			result[i] += vec[j] * Matrix[j][i];
	}
	return result;
}

int    fp_loc_iter = 0;
double* S_prev = new double[3];
double stepFR = 0;
static double **hess = nullptr;
static bool autoGenerate = false;
double Step_FletcherReeves(double x[], double y[], double *f, double norm)
{
	double* S_cur = new double[3];
	double **H = new double*[3];
	for (int i = 0; i < 3; i++)
		H[i] = new double[3];

	double **hessian = Hessian();

	if (!hess) {
		hess = Hessian();
	}
	double** A = Hessian();
	A[0][0] = 1;
	A[0][1] = A[1][0] = 4;
	A[0][2] = A[2][0] = 2;
	A[1][1] = 1;
	A[1][2] = A[2][1] = 6;
	A[2][2] = 1;

	double **ss = classic_gram_schmidt(A, hessian);

	if (fp_loc_iter == 0)
	{
		if (!autoGenerate) {
			S_cur[0] = S_prev[0] = -y[0];
			S_cur[1] = S_prev[1] = -y[1];
			S_cur[2] = S_prev[2] = -y[2];
		}
		else
			S_cur = S_prev = ss[0];

		double nominator = y[0] * S_cur[0] + y[1] * S_cur[1] + y[2] * S_cur[2];
		double *temp;
		temp = vector_dot_matrix(S_cur, hessian);
		double denominator = temp[0] * S_cur[0] + temp[1] * S_cur[1] + temp[2] * S_cur[2];
		stepFR = -nominator / denominator;

		fp_loc_iter++;
	}

	double nominator = 0;
	double denominator = 0;
	double *temp1;

	x[0] += stepFR * S_prev[0];
	x[1] += stepFR * S_prev[1];
	x[2] += stepFR * S_prev[2];

	/* Finding beta */
	double grad[3];
	Gradient(x, grad);
	temp1 = vector_dot_matrix(grad, hessian);
	nominator = temp1[0] * S_prev[0] + temp1[1] * S_prev[1] + temp1[2] * S_prev[2];
	temp1 = vector_dot_matrix(S_prev, hessian);
	denominator = temp1[0] * S_prev[0] + temp1[1] * S_prev[1] + temp1[2] * S_prev[2];
	double beta = nominator / denominator;

	S_cur[0] = -grad[0] + beta * S_prev[0];
	S_cur[1] = -grad[1] + beta * S_prev[1];
	S_cur[2] = -grad[2] + beta * S_prev[2];

	if (autoGenerate) {
		if (fp_loc_iter == 3)
			fp_loc_iter = 0;
		S_cur = ss[fp_loc_iter];
	}
	/* Finding step for next step*/
	double nominator1 = grad[0] * S_cur[0] + grad[1] * S_cur[1] + grad[2] * S_cur[2];
	double *temp;
	temp = vector_dot_matrix(S_cur, hessian);
	double denominator1 = temp[0] * S_cur[0] + temp[1] * S_cur[1] + temp[2] * S_cur[2];
	stepFR = -nominator1 / denominator1;

	/* Setting current as prev*/

	S_prev[0] = S_cur[0];
	S_prev[1] = S_cur[1];
	S_prev[2] = S_cur[2];

	*f = Function(x[0], x[1], x[2]);
	fp_loc_iter++;
	return stepFR;
}


double findAlpha(double x[], double y[], double *f, double S[])
{
	double alpha = (sqrt(5) - 1) / 2;
	double a = 0, b = 1;
	double l, m, fl, fm;

	l = a + (1 - alpha) * (b - a);
	m = a + b - l;
	fl = g(x[0] + S[0] * l, x[1] + S[1] * l, x[2] + S[2] * l);
	fm = g(x[0] + S[0] * m, x[1] + S[1] * m, x[2] + S[2] * m);

	while (fabs(m - l) > .00000001)
	{
		if (fl < fm)
		{
			b = m;
			m = l;
			fm = fl;
			l = a + (1 - alpha) * (b - a);
			fl = g(x[0] + S[0] * l, x[1] + S[1] * l, x[2] + S[2] * l);
		}
		else
		{
			a = l;
			l = m;
			fl = fm;
			m = a + alpha * (b - a);
			fm = g(x[0] + S[0] * m, x[1] + S[1] * m, x[2] + S[2] * m);
		}
	}
	return l;

}
int    mfp_loc_iter = 0;
double S_Mprev[3];
double stepMFR = 0;
double Step_ModificationFletcherReeves(double x[], double y[], double *f, double norm)
{
	double S_cur[3];
	double** hessian = hessianG(x);

	if (mfp_loc_iter == 0)
	{
		S_cur[0] = S_Mprev[0] = -y[0];
		S_cur[1] = S_Mprev[1] = -y[1];
		S_cur[2] = S_Mprev[2] = -y[2];
		stepMFR = findAlpha(x, y, f, S_cur);
		mfp_loc_iter++;
	}

	double nominator = 0;
	double denominator = 0;
	double *temp1 = NULL;

	x[0] += stepMFR * S_Mprev[0];
	x[1] += stepMFR * S_Mprev[1];
	x[2] += stepMFR * S_Mprev[2];

	/* Finding beta */
	double grad[3];
	gradG(x, grad);
	double delta[3];
	delta[0] = grad[0] - y[0];
	delta[1] = grad[1] - y[1];
	delta[2] = grad[2] - y[2];

	nominator = grad[0] * delta[0] + grad[1] * delta[1] + grad[2] * delta[2];
	denominator = y[0] * y[0] + y[1] * y[1] + y[2] * y[2];
	double beta = 0;
	if (mfp_loc_iter % 4 != 0)
	{
		beta = nominator / denominator;
		mfp_loc_iter++;
	}
	else mfp_loc_iter = 0;

	/* Finding S_cur*/
	S_cur[0] = -grad[0] + beta * S_Mprev[0];
	S_cur[1] = -grad[1] + beta * S_Mprev[1];
	S_cur[2] = -grad[2] + beta * S_Mprev[2];

	/* Finding step for next step*/
	stepMFR = findAlpha(x, y, f, S_cur);

	/* Setting current as prev*/
	S_Mprev[0] = S_cur[0];
	S_Mprev[1] = S_cur[1];
	S_Mprev[2] = S_cur[2];
	*f = g(x[0], x[1], x[2]);
	return stepMFR;
}

int Test(double(*Step)(double[], double[], double *, double), bool isModification, bool isNewton)
{
	int step = 0;
	double x[3] = { 1, 2, 3 }, x_prev[3], y[3], f = Function(x[0], x[1], x[2]), norm;
	int iters = 0;
	static double solution[3];
	if (isModification) {
		f = g(x[0], x[1], x[2]);
	}
	fp_loc_iter = 0;
	do
	{
		if (!isModification)
			Gradient(x, y);
		else
			gradG(x, y);

		norm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
		printf("\n x=(%06.4lf,%06.4lf,%06.4lf) \n y=(%06.4lf, %06.4lf,%06.4lf) \n grad_norm=%08.6lf,"
			"\n f=%07.5lf \n\n", x[0], x[1], x[2], y[0], y[1], y[2], norm, f);
		if (norm <= EPSILON)
			break;

		if (isNewton) {
			if (isModification)
				Gradient2(x, y, true, 1);
			else
				Gradient2(x, y, false, 1);
		}


		memcpy(x_prev, x, 2 * sizeof(double));
		printf(" step=%6.5lf\n", Step(x, y, &f, norm));

		step++;
	} while (norm > EPSILON);


	printf("\nDone in %d steps\n", step);

	return step;
}

int main(void)
{
	setlocale(LC_ALL, "Russian");

	std::cout << "1 - Метод Флетчера-Ривса" << std::endl;
	Test(Step_FletcherReeves, false, false);
	std::cout << std::endl;

	std::cout << "2 - Модифицированный метод Флетчера-Ривса" << std::endl;
	Test(Step_ModificationFletcherReeves, true, false);
	std::cout << std::endl;

	std::cout << "3 - Метод Ньютона" << std::endl;
	Test(Step_Newton, false, true);
	std::cout << std::endl;

	std::cout << "4 - Метод Ньютона для модификации" << std::endl;
	Test(Step_Newton_Modification, true, true);
	std::cout << std::endl;

	system("pause");
	return 0;
}