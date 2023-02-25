#include <stdio.h>
#include <math.h>

typedef struct {	//struct containing position inputs
    double x0;
    double x1;
} Point;

double rosenbrock(Point vertex) {	//rosenbock parabolic valley function
    double x0 = vertex.x0, x1 = vertex.x1;
    return 100 * pow(x1 - pow(x0, 2), 2) + pow(1 - x0, 2);
}

void file_data(char *textfile, double xlow, double xhigh, int values) { //writing data into textfile

    double x0 = xlow;
    double step = (xhigh - xlow) / (values - 1); 
    
    FILE *fp = fopen(textfile, "w");
    
    for (int i = 0; i < 100; i++, x0 += step) {
        Point P = {x0, 1};
        fprintf(fp, "%.3lf, %.3lf\n", x0, rosenbrock(P));
    }
    
    fclose(fp);
}

void sortPoints(double Y[], Point P[]) {	//reordering arrays from small to big
    double temp;
    Point pTemp;
    for (int i = 0; i < 2; i++) {
        for (int j = i + 1; j < 3; j++) {
            if (Y[i] > Y[j]) {
                temp = Y[i];
                Y[i] = Y[j];	//rearranging evaluations
                Y[j] = temp;
                
                pTemp = P[i];
                P[i] = P[j];	//rearranging positions
                P[j] = pTemp;
            }
        }
    }
}

void getCentroid(Point P[], Point *Pbar) {	//centroid
    Pbar->x0 = (P[0].x0 + P[1].x0) / 2;
    Pbar->x1 = (P[0].x1 + P[1].x1) / 2;
}

void reflect(Point *Ps, Point *Pbar, Point P[]) {	//P*
    Ps->x0 = 2 * Pbar->x0 - P[2].x0;
    Ps->x1 = 2 * Pbar->x1 - P[2].x1;
}

void expand(Point *Pss, Point *Ps, Point *Pbar) {	//P** expansion
    Pss->x0 = 2 * Ps->x0 - Pbar->x0;
    Pss->x1 = 2 * Ps->x1 - Pbar->x1;
}

void contract(Point *Pss, Point P[], Point *Pbar) {	//P** contraction
    Pss->x0 = (P[2].x0 + Pbar->x0) / 2;
    Pss->x1 = (P[2].x1 + Pbar->x1) / 2;
}

void replacePoint(Point *dest, Point *orig) {	//replace Ph with P* or P**
    dest->x0 = orig->x0;
    dest->x1 = orig->x1;
}

void getPi(Point P[]) {	//replace every point with midpoint of that point and Pl
    for (int i = 0; i < 3; i++) {
        P[i].x0 = (P[i].x0 + P[0].x0) / 2;
        P[i].x1 = (P[i].x1 + P[0].x1) / 2;
    }
}

_Bool MinCon(double Y[]) {	//condition for breaking loop
    double Ybar = 0, sum = 0;
    for (int i = 0; i < 3; i++) {
        Ybar += Y[i] / 3;
    }
    for (int i = 0; i < 3; i++) {
        sum += pow(Y[i] - Ybar, 2) / 2;
    }
    return (sqrt(sum) < pow(10, -8));
}

void algorithm(Point P[]) {	//loop containing algorithm
    int a;
    for (a = 0; a < 1000; a++) {	//stops after 1000 iterations
        
        double Y[3];
        for (int i = 0; i < 3; i++) {
            Y[i] = rosenbrock(P[i]);	//positions evaluated
        }
       
        sortPoints(Y, P);	//sort values and positions small to big
        
        Point Pbar, Ps, Pss;
        
        getCentroid(P, &Pbar);	//centroid calculation
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
        
        else if (Ys >= Y[1]) {
		    
		    if (Ys < Y[2]) {
		    	replacePoint(&P[2], &Ps);	//replace Ph with P*
		    }
		    contract(&Pss, P, &Pbar);	//P** contraction
		    double Yss = rosenbrock(Pss);	//y** evaluation
		    
		    if (Yss < Y[2]) {
		    	replacePoint(&P[2], &Pss);	//replace Ph with P**
		    }
		    
		    else {
		    	getPi(P);	//replace Pi with (Pi+Pl)/2
		    }
        }
        
        else {
        	replacePoint(&P[2], &Ps);	//replace Ph with P*
        }
        
        if (MinCon(Y)) {	//condition for breaking loop
        	break;
        }
    }
        printf("Final vertices:\n");
        for (int i=0; i<3; i++) {
        	printf("(%lf, %lf)\n", P[i].x0, P[i].x1);
        }
        
        printf("\nEvaluations:\n");
        for (int i=0; i<3; i++) {
        	printf("%lf\n", rosenbrock(P[i]));
        }
        
        printf("\nIterations: %d\n", a+1);
}

int main() {
	file_data("data.txt", -2., 2., 100);
    
    Point P[3] = {{0, 0}, {2, 0}, {0, 2}};	//initial positions
	algorithm(P);
	return 0;
}
