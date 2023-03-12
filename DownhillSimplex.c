/*****************************************************
 *                                                   *
 * Candidate number: 24563                           *
 *                                                   *
 * Downhill Simplex algorithm for finding the        *
 * minimum of the rosenbrock function                *
 *                                                   *
 * Starting vertices:                                *
 * p0 = (0,0), p1 = (2,0), p2 = (0,2)                *
 *                                                   *
 * This code can by amended for functions with more  *
 * variables by changing N to the number of          *
 * variables in the new function and the number of   *
 * initial vertices to N + 1                         *
 *                                                   *
 *****************************************************/

#include <stdio.h>
#include <math.h>

#define N 2     //number of dimensions

#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0
#define RHO 0.5     //coefficients for reflection, contraction, expansion and shrinking

#define TOL 1e-8
#define MAX_ITER 1000   //standard deviation tolerance and maximum iterations

typedef struct {
    double x[N];        //struct containing the coordinates of the vertices
} Point;

double func(Point P) {
    return 100 * pow(P.x[1] - pow(P.x[0], 2), 2) + pow(1 - P.x[0], 2);  //rosenbrock function
}

void sortPoints(double Y[], Point P[]) {
    double temp;
    Point pTemp;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N + 1; j++) {
            if (Y[i] > Y[j]) {
                temp = Y[i];
                Y[i] = Y[j];    //reordering outputs from smallest to largest
                Y[j] = temp;

                pTemp = P[i];
                P[i] = P[j];    //reordering the corresponding vertices
                P[j] = pTemp;
            }
        }
    }
}

void centroid(Point P[], Point* Pbar) {
    for (int i = 0; i < N; i++) {
        Pbar->x[i] = 0;
        for (int j = 0; j < N; j++) {
            Pbar->x[i] += P[j].x[i] / N;    //the centroid is the average of all vertices except Ph
        }
    }
}

void reflect(Point* Ps, Point* Pbar, Point P[]) {
    for (int i = 0; i < N; i++) {
        Ps->x[i] = (1 + ALPHA) * Pbar->x[i] - ALPHA * P[N].x[i];    //reflection of Ph, denoted by P*
    }
}

void insideContract(Point* Pss, Point P[], Point* Pbar) {
    for (int i = 0; i < N; i++) {
        Pss->x[i] = BETA * P[N].x[i] + (1 - BETA) * Pbar->x[i]; //this is used to reduce the size of the simplex
    }
}

void outsideContract(Point* Pss, Point* Ps, Point* Pbar) {
    for (int i = 0; i < N; i++) {
        Pss->x[i] = BETA * Ps->x[i] + (1 - BETA) * Pbar->x[i];  //this can be used to help to move the simplex out of a 
    }                                                           //narrow valley or a region with a high curvature
}

void expand(Point* Pss, Point* Ps, Point* Pbar) {
    for (int i = 0; i < N; i++) {
        Pss->x[i] = GAMMA * Ps->x[i] + (1 - GAMMA) * Pbar->x[i];    //expansion of P*, denoted by P**
    }
}

void shrink(Point P[]) {
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N; j++) {
            P[i].x[j] = P[0].x[j] + RHO * (P[i].x[j] - P[0].x[j]);  //shrinking simplex about the vertex Pl
        }
    }
}

void replacePoint(Point* new, Point* orig) {
    *new = *orig;                               //replacing Ph with P* or P**
}

_Bool MinCon(double Y[]) {
    double Ybar = 0, sum = 0;
    for (int i = 0; i < N + 1; i++) {
        Ybar += Y[i] / (N + 1);
    }
    for (int i = 0; i < N + 1; i++) {
        sum += pow(Y[i] - Ybar, 2) / N;
    }
    return (sqrt(sum) < TOL);   //condition for breaking the loop containing the simplex method
}

void simplex(Point P[]) {
    Point Pbar, Ps, Pss;
    double Ys, Yss, Y[N + 1];
    
    int a;  //iterations
    
    for (a = 0; a < MAX_ITER; a++) {

        for (int i = 0; i < N + 1; i++) {
            Y[i] = func(P[i]);              //initialising function outputs
        }

        sortPoints(Y, P);
        centroid(P, &Pbar);
        reflect(&Ps, &Pbar, P);

        Ys = func(Ps);

        if (Ys >= Y[0] && Ys < Y[N - 1]) {      //if y* >= yl and y* <= second worst vertex
            replacePoint(&P[N], &Ps);
        }

        else if (Ys < Y[0]) {           //if y* < yl
            expand(&Pss, &Ps, &Pbar);
            Yss = func(Pss);

            if (Yss < Ys) {                 //if y** < y*
                replacePoint(&P[N], &Pss);
            }

            else {
                replacePoint(&P[N], &Ps);
            }
        }

        else {          //if y* >= second worst vertex
            
            if (Ys < Y[N]) {                        //if y* < yh
                outsideContract(&Pss, &Ps, &Pbar);
                Yss = func(Pss);

                if (Yss < Ys) {                 //if y** < y*
                    replacePoint(&P[N], &Pss);
                }
                
                else {
                    shrink(P);
                }
            }

            else if (Ys >= Y[N]) {              //if y* >= yh
                insideContract(&Pss, P, &Pbar);
                Yss = func(Pss);

                if (Yss < Y[N]) {               //if y** < yh
                    replacePoint(&P[N], &Pss);
                }

                else {
                    shrink(P);
                }
            }
        }
      
        if (MinCon(Y)) {
            break;
        }
    }

    printf("Final vertices:\n");
    for (int i = 0; i < N + 1; i++) {
        printf("(%lf, %lf)\n", P[i].x[0], P[i].x[1]);
    }

    printf("\nEvaluations:\n");
    for (int i = 0; i < N + 1; i++) {
        printf("%lf\n", func(P[i]));
    }

    printf("\nIterations: %d\n", a + 1);
}

int main() {
    Point P[N + 1] = { {0, 0}, {2, 0}, {0, 2} };    //initial vertices (each vertex has N variables)
    simplex(P);
    return 0;
}
