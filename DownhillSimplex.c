#include <stdio.h>
#include <math.h>
#include <string.h>

#define N 2	//number of dimensions (best to match to number of variable inputs)
#define ALPHA 1.0
#define BETA 0.5	//reflection, expansion and contraction coefficients
#define GAMMA 2.0

typedef struct {
    double x[N];  //array of variable inputs (x0, x1, ...)
} Point;

double rosenbrock(Point P) {    //rosenbrock function
    return 100 * pow(P.x[1] - pow(P.x[0], 2), 2) + pow(1 - P.x[0], 2);
}

void file_data(char* textfile, double xlow, double xhigh, int values) { //writing data into textfile

    double x0 = xlow;
    double step = (xhigh - xlow) / (values - 1);

    FILE* fp = fopen(textfile, "w");

    for (int i = 0; i < 100; i++, x0 += step) {
        Point P = { x0, 1 };
        fprintf(fp, "%.3lf, %.3lf\n", x0, rosenbrock(P));
    }

    fclose(fp);
}

void sortPoints(double Y[], Point P[]) {	//reordering arrays from small to big
    double temp;
    Point pTemp;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N + 1; j++) {
            if (Y[i] > Y[j]) {
                temp = Y[i];
                Y[i] = Y[j];	//rearranging evaluations
                Y[j] = temp;

                pTemp = P[i];
                P[i] = P[j];	//rearranging corresponding vertices
                P[j] = pTemp;
            }
        }
    }
}

void centroid(Point P[], Point* Pbar) { //function to find centroid (average of all vertices except Ph)
    for (int i = 0; i < N; i++) {
        Pbar->x[i] = 0;
        for (int j = 0; j < N; j++) {
            Pbar->x[i] += P[j].x[i] / N;
        }
    }
}

void reflect(Point* Ps, Point* Pbar, Point P[]) {	//reflection of Ph, denoted by P*
    for (int i = 0; i < N; i++) {
        Ps->x[i] = (1 + ALPHA) * Pbar->x[i] - ALPHA * P[N].x[i];
    }
}

void expand(Point* Pss, Point* Ps, Point* Pbar) {	//expansion of P*, denoted by P**
    for (int i = 0; i < N; i++) {
        Pss->x[i] = GAMMA * Ps->x[i] + (1 - GAMMA) * Pbar->x[i];
    }
}

void contract(Point* Pss, Point P[], Point* Pbar) {	//contraction of the centroid in the direction of Ph, denoted by P**
    for (int i = 0; i < N; i++) {
        Pss->x[i] = BETA * P[N].x[i] + (1 - BETA) * Pbar->x[i];
    }
}

void replacePoint(Point* new, Point* orig) {	//replace Ph with P* or P**
    memcpy(new, orig, sizeof(Point));
}

void replacePi(Point P[]) {	//replace all vertices with the midpoint of that vertex and Pl
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N; j++) {
            P[i].x[j] = (P[i].x[j] + P[0].x[j]) / 2;
        }
    }
}

_Bool MinCon(double Y[]) {	//condition for breaking loop
    double Ybar = 0, sum = 0;
    for (int i = 0; i < N + 1; i++) {
        Ybar += Y[i] / (N + 1);
    }
    for (int i = 0; i < N + 1; i++) {
        sum += pow(Y[i] - Ybar, 2) / N;
    }
    return (sqrt(sum) < pow(10, -8));
}

void algorithm(Point P[]) {	//loop containing algorithm
    int a;
    for (a = 0; a < 1000; a++) {	//stops after 1000 iterations

        double Y[N + 1];
        for (int i = 0; i < N + 1; i++) {
            Y[i] = rosenbrock(P[i]);	//positions evaluated
        }

        sortPoints(Y, P);	//sort values and positions small to big

        Point Pbar, Ps, Pss;

        centroid(P, &Pbar);	//centroid calculation
        reflect(&Ps, &Pbar, P);	//P* calculation

        double Ys = rosenbrock(Ps);	//y* evaluation

        if (Ys < Y[0]) {

            expand(&Pss, &Ps, &Pbar);	//P** expansion
            double Yss = rosenbrock(Pss);	//y** evaluation

            if (Yss < Ys) {
                replacePoint(&P[2], &Pss);	//replace Ph with P**
            }

            else {
                replacePoint(&P[2], &Ps);	//replace Ph with P*
            }
        }

        else {
            for (int i = 1; i < N; i++) {
                if (Ys >= Y[i]) {

                    if (Ys < Y[2]) {
                        replacePoint(&P[2], &Ps);	//replace Ph with P*
                    }
                    contract(&Pss, P, &Pbar);	//P** contraction
                    double Yss = rosenbrock(Pss);	//y** evaluation

                    if (Yss < Y[2]) {
                        replacePoint(&P[2], &Pss);	//replace Ph with P**
                    }

                    else {
                        replacePi(P);	//replace Pi with (Pi+Pl)/2
                    }
                }

                else {
                    replacePoint(&P[2], &Ps);	//replace Ph with P*
                }
            }
        }

        if (MinCon(Y)) {	//condition for breaking loop
            break;
        }
    }
    printf("Final vertices:\n");
    for (int i = 0; i < N + 1; i++) {
        printf("(%lf, %lf)\n", P[i].x[0], P[i].x[1]);
    }

    printf("\nEvaluations:\n");
    for (int i = 0; i < N + 1; i++) {
        printf("%lf\n", rosenbrock(P[i]));
    }

    printf("\nIterations: %d\n", a);
}

int main() {
    file_data("data.txt", -2., 2., 100);

    Point P[N + 1] = { {0, 0}, {2, 0}, {0, 2} };	//initial vertices
    algorithm(P);
    return 0;
}
