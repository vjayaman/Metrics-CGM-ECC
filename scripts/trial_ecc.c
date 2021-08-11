#include <stdio.h>
#include <math.h>

void transformData2(void *dm_unknown, int size, char* dtype, double mind, double maxd) {
	double (*dm)[size] = dm_unknown;
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			dm[i][j] = log10(dm[i][j] + 10);
		}
	}
	double min_dist = log10(mind + 10);
	double max_dist = log10(maxd + 10);

	if (dtype == "temp") {
		printf("Temporal");
	}else {
		printf("Geographical");
	}
}

int main() {
	printf("Testing conversion: \n");

	double dm[4][4] = {{1, 1, 1, 1}, 
		       {2, 2, 2, 2}, 
		       {3, 3, 3, 3}, 
	               {4, 4, 4, 4}};

	transformData2(dm, 4, "temp", 2.0, 10.0);
	return 0;
}
