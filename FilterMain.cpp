#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <stdint.h>
#include <fstream>
#include "Filter.h"
#include <omp.h>
using namespace std;
#include "rtdsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int main(int argc, char **argv) {

	if (argc < 2) {
		fprintf(stderr, "Usage: %s filter inputfile1 inputfile2 .... \n",
				argv[0]);
	}

	//
	// Convert to C++ strings to simplify manipulation
	//
	string filtername = argv[1];

	//
	// remove any ".filter" in the filtername
	//
	string filterOutputName = filtername;
	string::size_type loc = filterOutputName.find(".filter");
	if (loc != string::npos) {
		//
		// Remove the ".filter" name, which should occur on all the provided filters
		//
		filterOutputName = filtername.substr(0, loc);
	}

	Filter *filter = readFilter(filtername);

	double sum = 0.0;
	int samples = 0;

	for (int inNum = 2; inNum < argc; inNum++) {
		string inputFilename = argv[inNum];
		string outputFilename = "filtered-" + filterOutputName + "-"
				+ inputFilename;
		struct cs1300bmp *input = new struct cs1300bmp;
		struct cs1300bmp *output = new struct cs1300bmp;
		int ok = cs1300bmp_readfile((char *) inputFilename.c_str(), input);

		if (ok) {
			double sample = applyFilter(filter, input, output);
			sum += sample;
			samples++;
			cs1300bmp_writefile((char *) outputFilename.c_str(), output);
		}
		delete input;
		delete output;
	}
	fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

struct Filter *
readFilter(string filename) {
	ifstream input(filename.c_str());

	if (!input.bad()) {
		int size = 0;
		input >> size;
		Filter *filter = new Filter(size);
		int div;
		input >> div;
		filter->setDivisor(div);
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				int value;
				input >> value;
				filter->set(i, j, value);
			}
		}
		return filter;
	}
}

#if defined(__arm__)
static inline unsigned int get_cyclecount (void)
{
	unsigned int value;
	// Read CCNT Register
	asm volatile ("MRC p15, 0, %0, c9, c13, 0\t\n": "=r"(value));
	return value;
}

static inline void init_perfcounters (int32_t do_reset, int32_t enable_divider)
{
	// in general enable all counters (including cycle counter)
	int32_t value = 1;

	// peform reset:
	if (do_reset)
	{
		value |= 2;     // reset all counters to zero.
		value |= 4;// reset cycle counter to zero.
	}

	if (enable_divider)
		value |= 8;     // enable "by 64" divider for CCNT.

	value |= 16;

	// program the performance-counter control-register:
	asm volatile ("MCR p15, 0, %0, c9, c12, 0\t\n" :: "r"(value));

	// enable all counters:
	asm volatile ("MCR p15, 0, %0, c9, c12, 1\t\n" :: "r"(0x8000000f));

	// clear overflows:
	asm volatile ("MCR p15, 0, %0, c9, c12, 3\t\n" :: "r"(0x8000000f));
}

#endif

double applyFilter(struct Filter *filter, cs1300bmp *input, cs1300bmp *output) {

#if defined(__arm__)
	init_perfcounters (1, 1);
#endif

	long long cycStart, cycStop;
	double start, stop;
#if defined(__arm__)
	cycStart = get_cyclecount();
#else
	cycStart = rdtscll();
#endif

	output->width = input->width;
	output->height = input->height;

	// Get the height and width of the filter
	short int in_width = (input->width)-1;
	short int in_height = (input->height)-1;
	short int filter_array[3][3];

	// Get the filter data into local int array
	#pragma omp parallel for
	for (short int i = 0; i < 3; i++) {

		*((*(filter_array+i))+0) = filter->get(i, 0);
		*((*(filter_array+i))+1) = filter->get(i, 1);
		*((*(filter_array+i))+2) = filter->get(i, 2);

	}

	// get the size into local variable
	short int filterSize = filter->getSize();

	//get the divison
	int filter_div = filter->getDivisor();

	//Hi line filter
	if (filter_array[0][1] == -2) {
	#pragma omp parallel for
		for (short int plane = 0; plane < 3; plane++) {
			for (short int rows = 1; rows < in_height; rows++) {
				const short int new_row1 = rows - 1;
				const short int new_row3 = rows + 1;
				for (short int c = 1; c < in_width; c++) {
					int pxl1 = 0, pxl2 = 0, pxl3 = 0;
					const short int col1 = c - 1;
					const short int col3 = c + 1;

						pxl1 = (input->color[plane][new_row1][col1]  * (*((*(filter_array+0))+0)));
						pxl2 = (input->color[plane][new_row1][c] * (*((*(filter_array+0))+1)));
						pxl3 = (input->color[plane][new_row1][col3] * (*((*(filter_array+0))+2)));

						pxl1 += (input -> color[plane][rows][col1] * (*((*(filter_array+1))+0)) );
						pxl2 += (input -> color[plane][rows][c] * (*((*(filter_array+1))+1)) );
						pxl3 += (input -> color[plane][rows][col3] * (*((*(filter_array+1))+2)) );

						pxl1 +=	(input->color[plane][new_row3][col1]* (*((*(filter_array+2))+0)));
						pxl2 += (input->color[plane][new_row3][c] * (*((*(filter_array+2))+1)));
						pxl3 += (input->color[plane][new_row3][col3] * (*((*(filter_array+2))+2)));

					
					pxl1 = pxl1 + pxl2 + pxl3;

					if (pxl1 > 255) {
						output->color[plane][rows][c] = 255;
					} else if (pxl1 < 0) {
						output->color[plane][rows][c] = 0;
					} else
						output->color[plane][rows][c] = pxl1;
				}
			}
		}
	}
	// Gauss
	else if (filter_array[1][1] == 8) {

#pragma omp parallel for

		for (short int plane = 0; plane < 3; plane++) {
			for (short int rows = 1; rows < in_height; rows++) {
				const short int new_row1 = rows - 1;
				const short int new_row3 = rows + 1;

				for (short int col2 = 1; col2 < in_width; col2++) {

					short int pxl1 = 0, pxl2 = 0, pxl3 = 0;
					const short int col1 = col2 - 1;
					const short int col3 = col2 + 1;

					pxl1 = (input->color[plane][new_row1][col1]  * (*((*(filter_array+0))+0)));
					pxl2 = (input->color[plane][new_row1][col2] * (*((*(filter_array+0))+1)));
					pxl3 = (input->color[plane][new_row1][col3] * (*((*(filter_array+0))+2)));

					pxl1 += (input -> color[plane][rows][col1] * (*((*(filter_array+1))+0)) );
					pxl2 += (input -> color[plane][rows][col2] * (*((*(filter_array+1))+1)) );
					pxl3 += (input -> color[plane][rows][col3] * (*((*(filter_array+1))+2)) );

					pxl1 +=	(input->color[plane][new_row3][col1]* (*((*(filter_array+2))+0)));
					pxl2 += (input->color[plane][new_row3][col2] * (*((*(filter_array+2))+1)));
					pxl3 += (input->color[plane][new_row3][col3] * (*((*(filter_array+2))+2)));

					// divide by 24
					pxl1 = ((pxl1 + pxl2 + pxl3) >> 3) / 3;

					if (pxl1 > 255) {
						output->color[plane][rows][col2] = 255;

					} else if (pxl1 < 0) {
						output->color[plane][rows][col2] = 0;

					} else
						output->color[plane][rows][col2] = pxl1;
				}
			}
		}

	}
	// Emboss
	else if (filter_array[1][2] == -1) {
#pragma omp parallel for
		for (short int p = 0; p < 3; p++) {
			for (short int r = 1; r < in_height; r++) {
				const short int new_row1 = r - 1;
				const short int new_row3 = r + 1;
				for (short int c = 1; c < in_width; c++) {
					short int pxl1 = 0, pxl2 = 0, pxl3 = 0;
					const short int col3 = c + 1;
					const short int col1 = c - 1;

					pxl1 = (input->color[p][new_row1][col1]  * (*((*(filter_array+0))+0)));
					pxl2 = (input->color[p][new_row1][c] * (*((*(filter_array+0))+1)));
					pxl3 = (~(input->color[p][new_row1][col3])+1);

					pxl1 += (input -> color[p][r][col1] * (*((*(filter_array+1))+0)) );
					pxl2 += (input -> color[p][r][c] * (*((*(filter_array+1))+1)) );
					pxl3 += (~(input -> color[p][r][col3])+1);

					pxl1 +=	(input->color[p][new_row3][col1]* (*((*(filter_array+2))+0)));
					pxl2 += (~(input->color[p][new_row3][c])+1);
					pxl3 += (~(input->color[p][new_row3][col3])+1);

					pxl1 += pxl2 + pxl3;
					if (pxl1 > 255) {
						output->color[p][r][c] = 255;

					} else if (pxl1 < 0) {
						output->color[p][r][c] = 0;

					} else
						output->color[p][r][c] = pxl1;
				}
			}
		}
	}
	// Average
	else {
#pragma omp parallel for
		for (short int p = 0; p < 3; p++) {
			for (short int r = 1; r < in_height; r++) {
				const short int new_row1 = r - 1;
				const short int new_row3 = r + 1;
				for (short int c = 1; c < in_width; c++) {
					short int pxl1 = 0, pxl2 = 0, pxl3 = 0;
					const short int col1 = c - 1;
					const short int col3 = c + 1;

					pxl1 += (input->color[p][new_row1][col1]  * (*((*(filter_array+0))+0)));

					pxl2 += (input->color[p][new_row1][c] * (*((*(filter_array+0))+1)));
					pxl3 += (input->color[p][new_row1][col3] * (*((*(filter_array+0))+2)));

					pxl1 += (input -> color[p][r][col1] * (*((*(filter_array+1))+0)) );
					pxl2 += (input -> color[p][r][c] * (*((*(filter_array+1))+1)) );
					pxl3 += (input -> color[p][r][col3] * (*((*(filter_array+1))+2)) );

					pxl1 +=	(input->color[p][new_row3][col1]* (*((*(filter_array+2))+0)));
					pxl2 += (input->color[p][new_row3][c] * (*((*(filter_array+2))+1)));
					pxl3 += (input->color[p][new_row3][col3] * (*((*(filter_array+2))+2)));


					pxl1=(pxl1+pxl2+pxl3)/9;



					if (pxl1 < 0) {
						output->color[p][r][c] = 0;

					} else if (pxl1 > 255) {
						output->color[p][r][c] = 255;

					} else
						output->color[p][r][c] = pxl1;
				}
			}
		}
	}

#if defined(__arm__)
	cycStop = get_cyclecount();
#else
	cycStop = rdtscll();
#endif
	double diff = cycStop - cycStart;
	double diffPerPixel = diff / (output->width * output->height);
	fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n", diff,
			diff / (output->width * output->height));
	return diffPerPixel;
}
