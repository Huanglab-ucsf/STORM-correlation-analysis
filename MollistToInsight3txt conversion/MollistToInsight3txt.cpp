// MollistToInsight3txt.cpp 
// Convert coordinate list to Insight3 txt format
// Input: tab seperated txt file with six columns: Frame X Y Width Background Intensity
// By Yina Wang 2017.11.20

#include "stdafx.h"
#include "stdlib.h"
#include "string.h"

int main(int argc, char* argv[])
{
	if(argc<2)
	{
		printf("Input: tab seperated coordinate list, column order: Frame X Y Width Background Intensity.\n");
		return 0;
	}

	// Default parameters, which is listed in Insight3 txt file but not used by correlation analysis
	float height = 500;
	float area = 5000;
	float phi = 0;
	float Ax = 1;
	int length = 1;
	int link = -1;
	int valid = 1;
	float z = 0;
	float zc = 0;
	int cat = 1;

	// Open molecule list (in text format)
	FILE *input;
	char* input_name = argv[1];
	fopen_s(&input,argv[1],"rt");
	if(!input)
	{
		printf("Cannot open input file %s\n",argv[1]);
		return 0;
	}

	// Get the total number of molecules
	int nmol = 0;
	char* buffer=new char[20000];	//line buffer
	int noeol = 0;
	while(fgets(buffer, 20000, input)){
		noeol = 0;
		if(!strchr(buffer,'\n')){
			noeol = 1;
			continue;
		}
		nmol++;
	}
	if (noeol) nmol++;
	fseek(input, 0, SEEK_SET); 

	// read and convert coordinate list
	int filename_len = strlen(input_name);
	char* output_name = new char[filename_len+14];	// 14 for "_Insight3format"
	strncpy(output_name, input_name, filename_len-4);
	output_name[filename_len-4] = '\0';
	strcat(output_name, "_Insight3format");
	strcat(output_name, ".txt");

	printf("Writing coordinate list to %s...", output_name);

	FILE *output;
	fopen_s(&output, output_name, "wt");
	// Write header line
	fprintf(output, "Cas%d\tX\tY\tXc\tYc\tHeight\tArea\tWidth\tPhi\tAx\tBG\tI\tFrame\tLength\tLink\tValid\tZ\tZc\n", nmol);
	float xc,yc,w,b,I;
	int frame;
	fgets(buffer,20000,input);
	for(int i=0; i<nmol; i++)
	{
		fscanf(input,"%d\t%f\t%f\t%f\t%f\t%f\n",&frame,&xc,&yc,&w,&b,&I);
		fprintf(output,"%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n",
			cat,xc,yc,xc,yc,height,area,w*2,phi,Ax,b,I,frame,length,link,valid,z,zc);
	}
	fclose(input);
	fclose(output);

	printf(" ...done.\n");

	return 1;
}

