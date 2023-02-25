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

int main() {

	file_data("data.txt", -2., 2., 100);
	
	return 0;
}