
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>


#include "mcrt.h"




int main()
{
	std::filesystem::current_path();
	std::filesystem::create_directory("slns");

	
	int interval = 3;
	size_t N = 10000;
	double R = 1.2;

	Mcrt code(R, N);
	code.execute_mcrt(interval);

	
	
	return 0;
}