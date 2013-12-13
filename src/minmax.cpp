#include "utils.h"

int main(int argc, char* argv[]) 
{
	string file;
	file.assign(argv[1]);

	int min, max;
	getMinMaxC(file.c_str(),&min,&max);

	cout << min << "\t" << max << endl;
}
