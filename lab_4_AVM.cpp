#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include "point.h"

using namespace std;

const int spacing = 7;
const int precision = 3;


const double dx = 0.05;
const double x = -2.4;
const Point initial(-1.3, 0.2);
const double a[] = { 1.645, -0.461, 0.423, 0.376, 0.012 };


double f(Point point);

vector<Point> method_euler(double dx);
vector<Point> mod_method_euler(double dx);
vector<Point> method_runge_kutta(double dx);

void graphics_1(vector<Point> euler, vector<Point> mod_euler, vector<Point> runge_kutta, vector<Point> reference);
void graphics_2(vector<Point> reference, vector<Point> euler_1, vector<Point> euler_2, vector<Point> euler_3, vector<Point> euler_4, vector<Point> euler_5);
void graphics_3(vector<double> method, vector<Point> error_points);

int main()
{
	cout << "variant #21" << endl;
	cout << "y' = 1,645 - 0.461x + 0.423x^2 + 0.376y + 0.012xy" << endl;
	cout << "Initial condition x = " << initial.x << " y = " << initial.y << endl;
	cout << "x = " << x << " dx = " << dx << endl << endl;

	cout << setprecision(precision) << fixed;

	cout << "Method Euler's" << endl;
	vector<Point> euler = method_euler(dx);

	cout << "\n\n" << "Modified method Euler's" << endl;
	vector<Point> mod_euler = mod_method_euler(dx);

	cout << "\n\n" << "Method Runge-Kutta (dx = 0.05)" << endl;
	vector<Point> runge_kutta = method_runge_kutta(dx);

	cout << "\n\n" << "Method Runge-Kutta (dx = 0.005)" << endl;
	vector<Point> reference = method_runge_kutta(dx / 10);

	cout << "\n\n" << "Method Euler's (dx = 0.025)" << endl;
	vector<Point> euler_1 = method_euler(dx / 2);

	cout << "\n\n" << "Method Euler's (dx = 0.0125)" << endl;
	vector<Point> euler_2 = method_euler(dx / 4);

	cout << "\n\n" << "Method Euler's (dx = 0.00833)" << endl;
	vector<Point> euler_3 = method_euler(dx / 6);

	cout << "\n\n" << "Method Euler's (dx = 0.00625)" << endl;
	vector<Point> euler_4 = method_euler(dx / 8);


	vector<double> method_errors;
	method_errors.push_back(abs(euler.back().y - reference.back().y));
	method_errors.push_back(abs(mod_euler.back().y - reference.back().y));
	method_errors.push_back(abs(runge_kutta.back().y - reference.back().y));
	cout << setw(16) << "Method" << setw(16) << "Absolute error" << endl;
	cout << setw(16) << "E" << setw(16) << method_errors[0] << endl;
	cout << setw(16) << "ME" << setw(16) << method_errors[1] << endl;
	cout << scientific << setw(16) << "RK" << setw(16) << method_errors[2] << endl;
	cout << setprecision(5) << fixed;

	vector<Point> points_errors;
	points_errors.push_back(Point(euler.size(), abs(euler.back().y - reference.back().y)));
	points_errors.push_back(Point(euler_1.size(), abs(euler_1.back().y - reference.back().y)));
	points_errors.push_back(Point(euler_2.size(), abs(euler_2.back().y - reference.back().y)));
	points_errors.push_back(Point(euler_3.size(), abs(euler_3.back().y - reference.back().y)));
	points_errors.push_back(Point(euler_4.size(), abs(euler_4.back().y - reference.back().y)));
	cout << setw(35) << "Number of points of integration" << setw(16) << "Absolute error" << endl;
	cout << setw(35) << euler.size() << setw(16) << points_errors[0].y << endl;
	cout << setw(35) << euler_1.size() << setw(16) << points_errors[1].y << endl;
	cout << setw(35) << euler_2.size() << setw(16) << points_errors[2].y << endl;
	cout << setw(35) << euler_3.size() << setw(16) << points_errors[3].y << endl;
	cout << setw(35) << euler_4.size() << setw(16) << points_errors[4].y << endl;


	graphics_1(euler, mod_euler, runge_kutta, reference);
	graphics_2(reference, euler, euler_1, euler_2, euler_3, euler_4);
	graphics_3(method_errors, points_errors);

	return 0;
}

double f(Point point)
{
	return a[0] + a[1] * point.x + a[2] * pow(point.x, 2) + a[3] * point.y + a[4] * point.x * point.y;
}

vector<Point> method_euler(double dx)
{

	cout << setw(spacing) << 'x' << setw(spacing) << 'y' << setw(spacing) << "K1" << endl;

	double l = max(initial.x, x) - min(initial.x, x);
	int n = (int)(l / dx);
	n += n * dx < l ? 1 : 0;
	dx *= initial.x < x ? 1 : -1;

	vector<Point> points;
	points.push_back(initial);

	Point prevPoint = initial;
	for (int i = 0; i < n; i++)
	{
		Point point;
		point.x = prevPoint.x + dx;
		double k1 = dx * f(prevPoint);
		point.y = prevPoint.y + k1;

		cout << setw(spacing) << prevPoint.x << setw(spacing) << prevPoint.y << setw(spacing) << k1 << endl;

		points.push_back(point);
		prevPoint = point;
	}

	cout << setw(spacing) << prevPoint.x << setw(spacing) << prevPoint.y << endl;

	return points;
}

vector<Point> mod_method_euler(double dx)
{

	cout << setw(spacing) << 'x' << setw(spacing) << 'y' << setw(spacing) << "K1" << setw(spacing) << "K2" << endl;

	double l = max(initial.x, x) - min(initial.x, x);
	int n = (int)(l / dx);
	n += n * dx < l ? 1 : 0;
	dx *= initial.x < x ? 1 : -1;

	vector<Point> points;
	points.push_back(initial);

	Point prevPoint = initial;
	for (int i = 0; i < n; i++)
	{
		Point point;
		point.x = prevPoint.x + dx;
		double k1 = (dx / 2) * f(prevPoint);
		double k2 = dx * f(Point(prevPoint.x + dx / 2, prevPoint.y + k1));
		point.y = prevPoint.y + k2;

		cout << setw(spacing) << prevPoint.x << setw(spacing) << prevPoint.y << setw(spacing) << k1 << setw(spacing) << k2 << endl;

		points.push_back(point);
		prevPoint = point;
	}

	cout << setw(spacing) << prevPoint.x << setw(spacing) << prevPoint.y << endl;

	return points;

}

vector<Point> method_runge_kutta(double dx)
{
	cout << setw(spacing) << 'x' << setw(spacing) << 'y' << setw(spacing) << "K1" << setw(spacing) << "K2" << setw(spacing) << "K3" << setw(spacing) << "K4" << endl;


	double l = max(initial.x, x) - min(initial.x, x);
	int n = (int)(l / dx);
	n += n * dx < l ? 1 : 0;
	dx *= initial.x < x ? 1 : -1;

	vector<Point> points;
	points.push_back(initial);

	Point prevPoint = initial;
	for (int i = 0; i < n; i++)
	{
		Point point;
		point.x = prevPoint.x + dx;
		double k1 = dx * f(prevPoint);
		double k2 = dx * f(Point(prevPoint.x + dx / 2, prevPoint.y + k1 / 2));
		double k3 = dx * f(Point(prevPoint.x + dx / 2, prevPoint.y + k2 / 2));
		double k4 = dx * f(Point(prevPoint.x + dx, prevPoint.y + k3));
		point.y = prevPoint.y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;

		cout << setw(spacing) << prevPoint.x << setw(spacing) << prevPoint.y << setw(spacing) << k1 << setw(spacing) << k2 << setw(spacing) << k3 << setw(spacing) << k4 << endl;

		points.push_back(point);
		prevPoint = point;
	}

	cout << setw(spacing) << prevPoint.x << setw(spacing) << prevPoint.y << endl;

	return points;
}

void graphics_1(vector<Point> euler, vector<Point> mod_euler, vector<Point> runge_kutta, vector<Point> reference)
{

	FILE* gnuplotPipe = _popen("gnuplot -p", "w");
	ofstream gnuplot(gnuplotPipe);
	gnuplot << "set terminal windows" << endl;
	gnuplot << "set title 'График зависимостей y = y(x)'" << endl;
	gnuplot << "set xlabel 'x'" << endl << "set ylabel 'y'" << endl;
	gnuplot << "plot '-' title 'Э' lw 2  with lines,"
		<< " '-' title 'МЭ' lw 2  with lines,"
		<< " '-' title 'РК (dx = 0.05)' lw 2  with lines,"
		<< " '-' title 'РК (dx = 0.005)' lw 2  with lines"
		<< endl;
	for (auto p = reference.cbegin(); p != reference.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = euler.cbegin(); p != euler.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = mod_euler.cbegin(); p != mod_euler.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = runge_kutta.cbegin(); p != runge_kutta.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	gnuplot.flush();

	_pclose(gnuplotPipe);

}

void graphics_2(vector<Point> reference, vector<Point> euler_1, vector<Point> euler_2, vector<Point> euler_3, vector<Point> euler_4, vector<Point> euler_5)
{

	FILE* gnuplotPipe = _popen("gnuplot -p", "w");
	ofstream gnuplot(gnuplotPipe);
	gnuplot << "set terminal windows" << endl;
	gnuplot << "set title 'График зависимостей y = y(x)'" << endl;
	gnuplot << "set xlabel 'x'" << endl << "set ylabel 'y'" << endl;
	gnuplot << "plot '-' title 'РК (dx = 0.005)' lw 2 with lines,"
		<< " '-' title 'Э (dx = 0.05)' lw 2 with lines,"
		<< " '-' title 'Э (dx = 0.025)' lw 2 with lines,"
		<< " '-' title 'Э (dx = 0.0125)' lw 2 with lines,"
		<< " '-' title 'Э (dx = 0.00833)' lw 2 with lines,"
		<< " '-' title 'Э (dx = 0.00625)' lw 2 with lines"
		<< endl;
	for (auto p = reference.cbegin(); p != reference.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = euler_1.cbegin(); p != euler_1.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = euler_2.cbegin(); p != euler_2.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = euler_3.cbegin(); p != euler_3.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = euler_4.cbegin(); p != euler_4.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	for (auto p = euler_5.cbegin(); p != euler_5.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	gnuplot.flush();

	_pclose(gnuplotPipe);

}

void graphics_3(vector<double> method, vector<Point> points)
{

	FILE* gnuplotPipe = _popen("gnuplot -p", "w");
	ofstream gnuplot(gnuplotPipe);
	gnuplot << "set terminal windows" << endl;
	gnuplot << "set multiplot layout 1, 2" << endl;
	gnuplot << "set title 'Зависимость абсолютной ошибки от метода'" << endl;
	gnuplot << "set xlabel 'Метод'" << endl << "set ylabel 'Абсолютная ошибка'" << endl;
	gnuplot << "set boxwidth 0.9 relative" << endl;
	gnuplot << "set style fill solid 1.0" << endl;
	gnuplot << "plot '-' using 2:xtic(1) title '' with boxes," << endl;
	gnuplot << "Э" << ' ' << method[0] << endl;
	gnuplot << "МЭ" << ' ' << method[1] << endl;
	gnuplot << "РК" << ' ' << method[2] << endl;
	gnuplot << "e" << endl;
	gnuplot.flush();

	gnuplot << "set title 'Зависимость абсолютной ошибки от кол-ва точек интегрирования'" << endl;
	gnuplot << "set xlabel 'Кол-во точек'" << endl << "set ylabel 'Абсолютная ошибка'" << endl;
	gnuplot << "plot '-' title '' with lines," << endl;
	for (auto p = points.cbegin(); p != points.cend(); ++p)
		gnuplot << p->x << ' ' << p->y << endl;
	gnuplot << "e" << endl;
	gnuplot.flush();

	_pclose(gnuplotPipe);

}

