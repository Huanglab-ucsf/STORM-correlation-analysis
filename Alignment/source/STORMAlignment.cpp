// STORMAlignment.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ppl.h"

#define PI 3.141592653589793
#define TWOPI 6.283185307179586

struct molecule
{
	double x0;
	double y0;
	double x;
	double y;
	double w;
	double b;
	double i;
	int f;
	double lp;
	double weight;
};

// Coordinates of a correlation point
struct correlation_vector
{
	double xi;
	double eta;
	double var;
	double weight;
};


// Stores a roto-translational transform
struct transform
{
	double dx;
	double dy;
	double da;
};


double gaussian(double x, double y, double s2)
{
	return exp(-(x*x+y*y)/(2*s2))/(TWOPI*s2);
}

int max_finder_counter = 0;

void find_corr_max_brute(correlation_vector* c, int C, double step, double range, int n, double angle, double &c_max, double &da, double &dx, double &dy)
{
	
	// DEBUG START
	/*
	char logname[64];
	sprintf(logname, "corrs\\%d.txt", max_finder_counter);
	FILE *log = fopen(logname, "w");
	fprintf(log, "xi\teta\tC\n");
	max_finder_counter++;
	*/
	// DEBUG END

	double xi, eta, c_value;
	// Loop through correlation x grid
	for (int u=0; u<2*n+1; u++)
	{
		xi = u*step - range;
		// Loop through correlation y grid to get maximum
		for (int v=0; v<2*n; v++)
		{
			eta = v*step - range;
			c_value = 0;
			// Loop through all correlation vectors
			for (int k=0; k<C; k++)
			{
				// Evaluate Gaussian at the grid position
				c_value += c[k].weight*gaussian(c[k].xi-xi, c[k].eta-eta, c[k].var);
			}
			
			// DEBUG START
			//fprintf(log, "%f\t%f\t%f\n", xi, eta, c_value);
			// DEBUG END
			
			// Check if we have a larger correlation value
			if (c_value > c_max)
			{
				c_max = c_value;
				// Store transformation for this frame
				dx = xi;
				dy = eta;
				da = angle; 
			}
		}
	}
	
	// DEBUG START
	//fclose(log);
	// DEBUG END
}


void write_outfile(char* input_name, molecule* m, int N, int n_iterations, double pixel_size)
{
	int filename_len = strlen(input_name);
	char* output_name = new char[filename_len+8+5+1+32];	// 8 for "_aligned", 5 for the iteration number, 1 for null character, 32 to save the day
	strncpy(output_name, input_name, filename_len-4);
	output_name[filename_len-4] = '\0';
	strcat(output_name, "_aligned");
	char iteration_string[32];
	sprintf(iteration_string, "%05d", n_iterations);
	strcat(output_name, iteration_string);
	strcat(output_name, ".txt");

	printf("Writing aligned molecule list to %s...", output_name);

	FILE *input = fopen(input_name, "r");
	FILE *output;
	fopen_s(&output, output_name, "wt");
	// Write header line
	fprintf(output, "Cas%d\tX\tY\tXc\tYc\tHeight\tArea\tWidth\tPhi\tAx\tBG\tI\tFrame\tLength\tLink\tValid\tZ\tZc\n", N);
	float x,y,xc,yc,h,a,w,phi,ax,b,I,z,zc;
	int cat,valid,frame,length,link;
	char* buffer=new char[20000];
	fgets(buffer,20000,input);
	for(int i=0; i<N; i++)
	{
		fscanf(input,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%f\t%f\n",
			&cat,&x,&y,&xc,&yc,&h,&a,&w,&phi,&ax,&b,&I,&frame,&length,&link,&valid,&z,&zc);
		fprintf(output,"%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%d\t%d\t%d\t%.5f\t%.5f\n",
			cat,m[i].x0/pixel_size+10,m[i].y0/pixel_size+10,m[i].x/pixel_size+10,m[i].y/pixel_size+10,h,a,w,phi,ax,b,I,m[i].f,length,link,valid,z,zc);
	}
	fclose(input);
	fclose(output);

	printf(" ...done.\n");
}

void help()
{
		printf("Coordinate-based model-free STORM particle alignment by Bo Huang lab 2017.11.04\n"
			" STORMAlignment molecule list (.txt) [parameters]\n"
			" Parameters\n"
			"   -p=xxx   Image pixel size (nm) (default = 107)\n"
			"   -noEM=x  1 for EMCCD camera mode, 0 for not (default = 1)\n"
			"   -wo=x    1 for localization density weighting, 0 for not (default = 1)\n"
			"   -ts=xx   Translation step size (nm) (default = 5)\n"
			"   -tr=xxx  Translation range (nm) (default = 25)\n"
			"   -a=xxx   Rotation angle step size (degree) (default = 5)\n"
			"   -n=xxx   Maximum interation number (default = 10)\n"
			" The molecule list is the text list saved by Insight3.\n"
			" Output:\n"
			"  *_alignedxxxxx.txt: aligned molecule list in text format\n"
			"\n"
			);
}



int main(int argc, char* argv[])
{

	if(argc<2)
	{
		help();
		return 0;
	}

	// Default analysis parameters
	double pixel_size = 107;  //DNA-origami pixelsize 107nm
	int EM=1;
	double trans_step = 5;
	double trans_range = 25;
	double angle_step = TWOPI/(360/5);
	int weight_origin = 1;
	int max_interation = 10;

	// Read parameters
	for (int i=2; i<argc; i++)
	{
		if(!strncmp(argv[i],"-p=",3))
			pixel_size = atof(argv[i]+3);
		else if(!strncmp(argv[i],"-noEM=",6))
			EM = atof(argv[i]+6);
		else if(!strncmp(argv[i],"-wo=",4))
			weight_origin = atof(argv[i]+4);
		else if(!strncmp(argv[i],"-ts=",4))
			trans_step = atof(argv[i]+4);
		else if(!strncmp(argv[i],"-tr=",4))
			trans_range = atof(argv[i]+4);
		else if(!strncmp(argv[i],"-a=",3))
			angle_step = TWOPI*atof(argv[i]+3)/360;
		else if(!strncmp(argv[i],"-n=",3))
			max_interation = atof(argv[i]+3);
		else
		{
			printf("Invalid parameter %s\n", argv[i]);
			return 0;
		}
	}

	// Open molecule list (in text format)
	FILE *input;
	fopen_s(&input,argv[1],"rt");
	if(!input)
	{
		printf("Cannot open input file %s\n",argv[1]);
		return 1;
	}

	// Get the total number of molecules
	int N=0;
	fscanf_s(input,"Cas%d",&N);
	if(N==0)
	{
		printf("Invalid input file format\n");
		return 1;
	}

	// Allocating memory for molecules
	molecule *m = new molecule[N];
	if(m == NULL)
	{
		printf("Not enough memory for molecules.\n");
		return 1;
	}

	// Prepare to read the molecule list line by line
	char* buffer=new char[20000];	//line buffer
	fgets(buffer,20000,input);
	
	for(int i=0; i<N; i++)
	{
		fgets(buffer,20000,input);
		sscanf_s(buffer,"%*d%*lf%*lf%lf%lf%*lf%*lf%lf%*lf%*lf%lf%lf%d%*d%*d%*d%*lf%*lf",&(m[i].x),&(m[i].y),&(m[i].w),&(m[i].b),&(m[i].i),&(m[i].f));
		m[i].x *= pixel_size;
		m[i].y *= pixel_size;
	}
	fclose(input);

	// Get maximum frame number
	int F=0;
	for (int i=0; i<N; i++)
	{
		if (m[i].f > F)
		{
			F = m[i].f;
		}
	}

	// Calculating #molecules per frame.
	int *N_frame = new int[F+1];
	for (int i=0; i<=F; i++)
	{
		N_frame[i] = 0;
	}
	for (int i=0; i<N; i++)
	{
		N_frame[m[i].f]++;
	}
	
	printf("Aligning %d molecules in %d frames.\n", N, F);

	// Retrieving molecule indices where a new frame starts
	int *frame_start = new int[F+1];
	frame_start[0] = 0;
	int i=0;
	for (int f=1; f<=F; f++)
	{
		frame_start[f] = i;
		while (m[i].f==f)
		{
			i++;
		}
	}


	// Align by COM
	for (int f=1; f<=F; f++)
	{
		double com_x=0;
		double com_y=0;
		for (int i=frame_start[f]; i<frame_start[f]+N_frame[f]; i++)
		{
			com_x += m[i].x;
			com_y += m[i].y;
		}
		com_x /= N_frame[f];
		com_y /= N_frame[f];
		for (int i=frame_start[f]; i<frame_start[f]+N_frame[f]; i++)
		{
			m[i].x -= com_x;
			m[i].y -= com_y;
			m[i].x0 = m[i].x;
			m[i].y0 = m[i].y;
		}
	}

	// Calculate localization precisions
	double s2=0;
	for (int i=0; i<N; i++)
		s2 += m[i].w/2;
	s2 = s2*s2/(N*N); // (mean of molecule width/2) squared
	double pixelsize2 = pixel_size*pixel_size;
	double sa2 = s2 + pixelsize2/12;
	double var;
	double lp_mean = 0;
	for (int i=0; i<N; i++)
	{
		var = sa2 * ( 16/9. + 8*3.14159265359*sa2*m[i].b / ( m[i].i*pixelsize2 ) ) / m[i].i;
		if (EM==1)
			var *= 2;
		m[i].lp = sqrt(var);
		lp_mean += m[i].lp;
		if (m[i].lp != m[i].lp)
		{
			printf("Aborted: At least one localization precision is NaN.\n");
			return 0;
		}
	}
	lp_mean /= N;
	printf("Average localization precision: %f nm\n", lp_mean);

	// Calculate weighting factor based on local density
	for (int f=1; f<=F; f++)
	{
		for (int i=frame_start[f]; i<frame_start[f]+N_frame[f]; i++)
		{
			if (weight_origin == 1)
			{
				double density = 0.0;
				for (int j=frame_start[f]; j<frame_start[f]+N_frame[f]; j++)
				{
					double distance = sqrt((m[i].x-m[j].x)*(m[i].x-m[j].x)+(m[i].y-m[j].y)*(m[i].y-m[j].y));
					if (distance < 2*lp_mean)
					{
						density++;
					}
				}
				m[i].weight = 1.0/density;
			}
			else
			{
				m[i].weight = 1.0;
				//printf("No density weighting!\n");
			}
		}
	}

	double angle_range = TWOPI;
	int n_angles = int(ceil(angle_range/angle_step));
	angle_range = n_angles*angle_step;
	printf("Evaluating %d angles over range of %g degrees.\n", n_angles, 360*angle_range/TWOPI);

	int n_trans_steps = int(ceil(trans_range/trans_step));
	trans_range = n_trans_steps*trans_step;
	printf("Evaluating %dx%d translations over range of +/-%g nm.\n", 2*n_trans_steps+1, 2*n_trans_steps+1, trans_range);

	// Definitions for main loop		
	int n_iterations = 0;
	transform *t = new transform[F+1];
	int frame_counter;
	Concurrency::critical_section cs;
	double x_old, y_old;
	double t_max, da_max, translation;

	do
	{
		n_iterations++;
		printf("--- Iteration %d ---\n", n_iterations);
		
		frame_counter = 0;
		
		// Reset all transformations
		for (int f=1; f<=F; f++)
		{
			t[f].dx = 0;
			t[f].dy = 0;
			t[f].da = 0;
		}
		
		printf("    Processed 0 frames out of %d.\r", F);

		// Loop through each frame to get maximizing transformation
		Concurrency::parallel_for (1, F+1, [&] (int f)															
		//for (int f=1; f<=F; f++)
		{
			int C = N_frame[f]*(N-N_frame[f]);	// The number of correlation vectors
			correlation_vector *c = new correlation_vector[C];	// The cross correlation of the current frame with all other frames
			double c_max = 0;	// Current correlation function maximum value

			for (int a=0; a<n_angles; a++)
			{
				double angle = (a+1)*angle_step - angle_range/2;
				int k = 0;
				double x_rot, y_rot;
				for (int i=frame_start[f]; i<frame_start[f]+N_frame[f]; i++)
				{
					// Calculate rotated coordinates
					x_rot = m[i].x*cos(angle) - m[i].y*sin(angle);
					y_rot = m[i].x*sin(angle) + m[i].y*cos(angle);
					// Loop through molecules in all other frames and save correlation vector
					for (int j=0; j<N; j++)
					{
						if (m[j].f!=f)
						{
							c[k].xi = m[j].x - x_rot;
							c[k].eta = m[j].y - y_rot;
							c[k].var = 2 * (m[j].lp*m[j].lp + m[i].lp*m[i].lp);
							c[k].weight = m[i].weight;
							k++;
						}
					}
				}
				// Search for translational maximum for this angle
				find_corr_max_brute(c, C, trans_step, trans_range, n_trans_steps, angle, c_max, t[f].da, t[f].dx, t[f].dy);
			}
			delete [] c;
			cs.lock();
			frame_counter++;
			printf("    Processed %d frames out of %d.\r", frame_counter, F);
			cs.unlock();
		});
		//}
		printf("\n");

		// Now we have a transformation for each frame. Apply it to the respective molecules:
		for (int i=0; i<N; i++)
		{
			x_old = m[i].x;
			y_old = m[i].y;
			m[i].x = x_old*cos(t[m[i].f].da) - y_old*sin(t[m[i].f].da) + t[m[i].f].dx;
			m[i].y = x_old*sin(t[m[i].f].da) + y_old*cos(t[m[i].f].da) + t[m[i].f].dy;
		}

		// Calculate maximum translation and rotation
		t_max = 0;
		da_max = 0;
		for (int f=1; f<=F; f++)
		{
			translation = sqrt(t[f].dx*t[f].dx + t[f].dy*t[f].dy);
			if (translation>t_max)
				t_max = translation;
			if (abs(t[f].da) > da_max)
				da_max = abs(t[f].da);
		}
		printf("    Maximum translation: %g\n", t_max);
		printf("    Maximum rotation: %g\n", 360*da_max/TWOPI);

		write_outfile(argv[1], m, N, n_iterations, pixel_size);
	} while ((t_max >= trans_step || da_max >= angle_step) && n_iterations < max_interation);

	return 0;
}