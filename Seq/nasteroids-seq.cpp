#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <locale>
#include <vector>
#include <string>
#include <chrono>
#include <climits>
#include <cmath>
#include <random>
#include <thread>
#include <iomanip>
#include <ctime>
#include <fstream>



using namespace std;
using namespace std::chrono;
using clk = chrono::high_resolution_clock;

const double gravity = 6.674e-5;
const double intTime = 0.1;
const double dmin = 5.0;
const double width = 200.0;
const double height = 200.0;
const int m = 1000;
const int stdev = 50;

typedef struct planeta{
	double ejX, ejY, masa;
}planeta;

typedef struct asteroide{
	double ejX, ejY, masa, VelX, VelY, FuX, FuY;
}asteroide;

void calculoFuerzasPl(asteroide *a, planeta *b);
void calculoFuerzasAst(asteroide *a, asteroide *b);
void reboteAst(asteroide *a);
void velocidadAsteroide(asteroide *a);
bool reboteEntreAst(asteroide *a, asteroide *b);


int main (int argc, char *argv[]) {

	//Primero comprobamos que el numero de argumentos introducidos es correcto
	if(argc != 5)
	{
		cout << "nasteroids-seq: Wrong arguments"<< endl << "Correct use:"<< endl << "nasteroids-seq nAsteroides nIteraciones nPlanetas seed" << endl;
		return -1;
	}

	//Comprobamos la validez de los datos introducidos
	istringstream arg1(argv[1]);
	istringstream arg2(argv[2]);
	istringstream arg3(argv[3]);
	istringstream arg4(argv[4]);
	int a;

	if (!(arg1 >> a && arg1.eof()) || !(arg2 >> a && arg2.eof()) || !(arg3 >> a && arg3.eof()) || !(arg4 >> a && arg4.eof()))
	{
		cout << "nasteroids-seq: Wrong arguments"<< endl << "Correct use:"<< endl << "nasteroids-seq nAsteroides nIteraciones nPlanetas seed" << endl;
		return -1;
	}

	// Guardamos los argumentos introducidos
	int nAsteroides = atoi(argv[1]);
	int nIteraciones = atoi(argv[2]);
	int nPlanetas = atoi(argv[3]);
	int seed = atoi(argv[4]);
	int max_i = numeric_limits<int>::max();

	//Comprobamos que los numeros introducidos son validos
	if( (0 > nAsteroides || nAsteroides > 1000) || (0 > nIteraciones || nIteraciones > 5000) || (0 > nPlanetas || nPlanetas > max_i) || (0 > seed || seed > max_i) )
	{
		cout << "nasteroids-seq: Wrong arguments"<< endl << "Correct use:"<< endl << "nasteroids-seq nAsteroides nIteraciones nPlanetas seed" << endl;
		return -1;
	}

	auto start = clk::now();

	// random distribution
	default_random_engine re{seed};
	uniform_real_distribution<double> xdist{0.0, std::nextafter(width, std::numeric_limits<double>::max())};
	uniform_real_distribution<double> ydist{0.0, std::nextafter(height, std::numeric_limits<double>::max())};
	normal_distribution<double> mdist{m, stdev};

	// Fichero de salida init_conf
	ofstream configFile ("init_conf.txt");
	configFile << fixed << setprecision(3) << nAsteroides << " " << nIteraciones << " " << nPlanetas << " " << seed << endl;

	// Inicializamos los asteroides y sus parametros
	vector<asteroide> vec_Ast(nAsteroides);

	for (int unsigned i=0 ; i<vec_Ast.size() ; i++)
	{
		vec_Ast[i] = {xdist(re), ydist(re), mdist(re), 0.0, 0.0, 0.0, 0.0};
		configFile << fixed << setprecision(3) << vec_Ast[i].ejX << " " << vec_Ast[i].ejY << " " << vec_Ast[i].masa << endl;
	}

	// Inicializamos los planetas y sus parametros
	vector<planeta> vec_Pl(nPlanetas);

	for (int unsigned i=0 ; i<vec_Pl.size() ; i++)
	{
		if (i==0 || i%4 == 0)
		{
			vec_Pl[i] = {0.0, ydist(re), mdist(re)*10};
		}
		else if (i%4 == 1)
		{
			vec_Pl[i] = {xdist(re), 0.0, mdist(re)*10};
		}
		else if (i%4 == 2)
		{
			vec_Pl[i] = {width, ydist(re), mdist(re)*10};
		}
		else if (i%4 == 3)
		{
			vec_Pl[i] = {xdist(re), height, mdist(re)*10};
		}

		configFile << fixed << setprecision(3) << vec_Pl[i].ejX << " " << vec_Pl[i].ejY << " " << vec_Pl[i].masa << endl;
	}

		configFile.close();


	//Se calculan las fuerzas entre cuerpos
	for (int n=0; n < nIteraciones; n++)
	{
		for (int unsigned i=0; i<vec_Ast.size(); i++)
		{
			vec_Ast[i].FuX = 0.0;
			vec_Ast[i].FuY = 0.0;
		}
		for (int unsigned i=0; i<vec_Ast.size(); i++) 
		{
			for (int unsigned j=i+1; j<vec_Ast.size(); j++)
			{
				calculoFuerzasAst(&vec_Ast[i], &vec_Ast[j]);
			}
		}
		for (int unsigned i=0; i<vec_Ast.size(); i++) 
		{
			for (int unsigned j=0; j<vec_Pl.size(); j++)
			{
				calculoFuerzasPl(&vec_Ast[i], &vec_Pl[j]);
			}
		}
		for (int unsigned i=0; i<vec_Ast.size(); i++)
		{
			velocidadAsteroide(&vec_Ast[i]);
			reboteAst(&vec_Ast[i]);
		}

		// Calculo de choques entre asteroides
		double auxVelX = 0.0, auxVelY = 0.0;
		for (int unsigned i=0; i<vec_Ast.size(); i++)
		{
			for (int unsigned j=i+1; j<vec_Ast.size(); j++)
			{
				if(reboteEntreAst(&vec_Ast[i],&vec_Ast[j]))
				{
					auxVelX = vec_Ast[i].VelX;
					auxVelY = vec_Ast[i].VelY;
					vec_Ast[i].VelX = vec_Ast[j].VelX;
					vec_Ast[i].VelY = vec_Ast[j].VelY;
					vec_Ast[j].VelX = auxVelX;
					vec_Ast[j].VelY = auxVelY;
				}
			}
		}
	}

	// Fichero de salida out.txt
	ofstream outFile ("out.txt");

	for (int i = 0; i < nAsteroides; ++i)
	{
		outFile << fixed << setprecision(3) << vec_Ast[i].ejX << " " << vec_Ast[i].ejY << " " << vec_Ast[i].VelX << " "  << vec_Ast[i].VelY << " " << vec_Ast[i].masa << endl;
	}

	outFile.close();

	auto end = clk::now();
	auto dif = duration_cast<milliseconds>(end-start);
	cout << "Tiempo ejecucion: " << dif.count() << " ms" << endl;
}


void calculoFuerzasPl(asteroide *a, planeta *b) {

	double difX = a->ejX - b->ejX; double difY = a->ejY - b->ejY;
	double distAB = sqrt( (pow(difX, 2)) + (pow(difY, 2)) );
	double slope = difY / difX;

	if  (slope > 1)
	{
		slope = 1;
	}

	if (slope < -1)
	{
		slope = -1;
	}

	double ang = atan(slope);
	double F = (gravity * a->masa * b->masa) / (pow(distAB, 2)) ;

	if (F > 200.000)
	{
		F = 200.000;
	}

	double FuX = F * cos(ang);
	double FuY = F * sin(ang);

	a->FuX = a->FuX + FuX;
	a->FuY = a->FuY + FuY;
}

void velocidadAsteroide(asteroide *a) {

	double Ax = a->FuX / a->masa;
	double Ay = a->FuY / a->masa;

	a->VelX = a->VelX + (Ax * intTime);
	a->VelY = a->VelY + (Ay * intTime);
	a->ejX = a->ejX + (a->VelX * intTime);
	a->ejY = a->ejY + (a->VelY * intTime);
}

void calculoFuerzasAst(asteroide *a, asteroide *b) {

	double difX = a->ejX - b->ejX;
	double difY = a->ejY - b->ejY;
	double distAB = sqrt( (pow(difX, 2)) + (pow(difY, 2)) );

	if (distAB > dmin)
	{
		double slope = difY / difX;

		if  (slope > 1)
		{
			slope = 1;
		}

		if (slope < -1)
		{
			slope = -1;
		}

		double ang = atan(slope);
		double F = (gravity * a->masa * b->masa) / (pow(distAB, 2)) ;

		if (F > 200.000) 
		{
			F = 200.000;
		}

		double FuX = F * cos(ang);
		double FuY = F * sin(ang);

		a->FuX = a->FuX + FuX;
		a->FuY = a->FuY + FuY;

		b->FuX = b->FuX - FuX;
		b->FuY = b->FuY - FuY;
	}
}


void reboteAst(asteroide *a){

	if (a->ejX <= 0.0)
	{
		a->ejX = dmin;
		a->VelX = -1 * a->VelX;
	}
	if (a->ejY <= 0.0)
	{
		a->ejY = dmin;
		a->VelY = -1 * a->VelY;
	}
	if (a->ejX >= width)
	{
		a->ejX = width-dmin;
		a->VelX = -1 * a->VelX;
	}
	if (a->ejY >= height)
	{
		a->ejY = height-dmin;
		a->VelY = -1 * a->VelY;
	}
}

bool reboteEntreAst(asteroide *a, asteroide *b){

	double distBtwAB = sqrt((pow(a->ejX - b->ejX, 2)) + (pow(a->ejY - b->ejY, 2)));
	if (distBtwAB <= dmin)
	{
		return true;
	}
	return false;
}
