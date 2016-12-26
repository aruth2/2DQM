#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "aruthio.h"
#include "aruthmath.h"

#include <gtk/gtk.h>

/* 2DQM.c created by Anthony Ruth on March 10th, 2016
 * Solves the shroedinger equation on a 2 dimensional square or hexagonal tiled surface 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * */


typedef struct { 
		int *compositevectors;
		int numcomp;
		double *latticevectors;
		double *vertices;
		int numvertices;
		} tiling;

	GtkWidget *TopLevel;

	GtkWidget *drawingarea;
	GtkWidget *rowsspin;
	GtkWidget *columnsspin;
	GtkWidget *fixspin;
	GtkWidget *potentialspin;
	GtkWidget *chargespin;
	GtkWidget *eigenvaluelist;
	GtkWidget *treeview;
	
	
	GtkWidget *squarewavefunction;
	GtkWidget *tiletype;
	GtkWidget *interspin;
	GtkWidget *broadeningspin;
	GtkWidget *upperlimitspin;

	
	
	GtkWidget *roundnessthicknessspin;
	GtkWidget *maximumradiusspin;
	GtkWidget *chargepotentialspin;

	
	double *potential;
	double *charge;
	int *fixedpoints;
	int *displayedwavefunctions;
	double *eigenvectors;
	double *eigenvalues;
	double *condensedeigenvalues,*condensedeigenvectors;
	int rows = 10;
	int columns = 10;
	int numpoints;
	tiling *tiles;
	double *xgrid;
	double *ygrid;
	
	double *chargedensity;
	double *chargepotential;
	double *pairdensity;
	
	double chargepotentialscalefactor = 0.01;//Determines how strong the potential is from the charges.
	
	FILE *wavefunctionpipe;
	FILE *fixedpipe;
	FILE *potentialpipe;
	FILE *chargepotentialpipe;
	FILE *dospipe;
	FILE *pointspecpipe;
	FILE *baderchargepipe;
	
	gboolean save()
{
	//getdata();
	
	printf("Saving data\n");

	
	GtkWidget *dialog;
	GtkFileChooser *chooser;
	gint res;

dialog = gtk_file_chooser_dialog_new ("Save Data To Folder",
                                      TopLevel,
                                     GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER,
                                      "Cancel",
                                      GTK_RESPONSE_CANCEL,
                                      "Save",
                                      GTK_RESPONSE_ACCEPT,
                                      NULL);
chooser = GTK_FILE_CHOOSER (dialog);

gtk_file_chooser_set_do_overwrite_confirmation (chooser, TRUE);


res = gtk_dialog_run (GTK_DIALOG (dialog));
	if (res == GTK_RESPONSE_ACCEPT)
  {
    char *foldername;
	char filename[1000];
    foldername = gtk_file_chooser_get_filename (chooser);
	
	sprintf(filename,"%s/settings",foldername);
	FILE *settingsfile = fopen(filename,"w");
	int type = gtk_combo_box_get_active (tiletype);
	
	fprintf(settingsfile,"Tiling = %d\n",type);
	fprintf(settingsfile,"Rows = %d\n",rows);
	fprintf(settingsfile,"Columns = %d\n",columns);
	fclose(settingsfile);
	
	sprintf(filename,"%s/mode",foldername);
	FILE *modefile = fopen(filename,"w");
	printintmatrixfile(modefile,fixedpoints,rows,columns);
	fclose(modefile);
	
	sprintf(filename,"%s/potential",foldername);
	FILE *potentialfile = fopen(filename,"w");
	printdoublematrixfile(potentialfile,potential,rows,columns);
	fclose(potentialfile);
	
	sprintf(filename,"%s/charge",foldername);
	FILE *chargefile = fopen(filename,"w");
	printdoublematrixfile(chargefile,charge,rows,columns);
	fclose(chargefile);	
  }

	gtk_widget_destroy (dialog);
	
	
	
	return FALSE;
}

	gboolean load()
{
	//getdata();
	
	printf("Loading data\n");

	
	GtkWidget *dialog;
	GtkFileChooser *chooser;
	gint res;

dialog = gtk_file_chooser_dialog_new ("Load Data From Folder",
                                      TopLevel,
                                     GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER,
                                      "Cancel",
                                      GTK_RESPONSE_CANCEL,
                                      "Load",
                                      GTK_RESPONSE_ACCEPT,
                                      NULL);
chooser = GTK_FILE_CHOOSER (dialog);

gtk_file_chooser_set_do_overwrite_confirmation (chooser, TRUE);


res = gtk_dialog_run (GTK_DIALOG (dialog));
	if (res == GTK_RESPONSE_ACCEPT)
  {
    char *foldername;
	char filename[1000];
    char buffer[1000];
    foldername = gtk_file_chooser_get_filename (chooser);
	
	printf("Loading Settings\n");
	sprintf(filename,"%s/settings",foldername);
	FILE *settingsfile = fopen(filename,"r");
	
	strip(settingsfile,"Tiling",buffer," =",1);
	gtk_combo_box_set_active (tiletype,atoi(buffer));
	printf("Tiles %d\n",atoi(buffer));
	
	strip(settingsfile,"Rows",buffer," =",1);
	gtk_spin_button_set_value(rowsspin,atoi(buffer));
	printf("Rows %d\n",atoi(buffer));


	strip(settingsfile,"Columns",buffer," =",1);
	gtk_spin_button_set_value(columnsspin,atoi(buffer));
	fclose(settingsfile);
	printf("Columns %d\n",atoi(buffer));

	
		printf("Making Grid\n");
	makegrid();
	
		printf("Loading Modes\n");
	sprintf(filename,"%s/mode",foldername);
	FILE *modefile = fopen(filename,"r");
	readintmatrixfile(modefile,fixedpoints,rows,columns);
	fclose(modefile);
	
		printf("Loading Potential\n");
	sprintf(filename,"%s/potential",foldername);
	FILE *potentialfile = fopen(filename,"r");
	readdoublematrixfile(potentialfile,potential,rows,columns);
	fclose(potentialfile);
	
		printf("Loading Charge\n");
	sprintf(filename,"%s/charge",foldername);
	FILE *chargefile = fopen(filename,"r");
	readdoublematrixfile(chargefile,charge,rows,columns);
	fclose(chargefile);	
  }

	gtk_widget_destroy (dialog);
	
	
	
	return FALSE;
}
	
	
	double squarelv[] = {1,0,0,1};
	double squarev[] = {0,0,1,0,1,1,0,1};
	tiling *square()
	{
		tiling *square = malloc(sizeof(tiling));
		square->numcomp = 0;
		square->latticevectors = squarelv;
		square->vertices = squarev;
		square->numvertices = 4;
		return square;
	}
	
	
	int hexcv[] = {1,-1};
	double hexlv[] = {sqrt(3)/2,1.0/2,0,1};
//	double hexv[] = {0,0.5,sqrt(3)/4,0,3*sqrt(3)/4,0,sqrt(3),0.5,3*sqrt(3)/4,1,sqrt(3)/4,1};
	double hexv[] = {0,0.5,1/(2*sqrt(3)),0,3/(2*sqrt(3)),0,4/(2*sqrt(3)),0.5,3/(2*sqrt(3)),1,1/(2*sqrt(3)),1};

	tiling *hexagonaltiling()
	{
		tiling *hex = malloc(sizeof(tiling));
		hex->numcomp = 1;
		hex->compositevectors = hexcv;
		hex->latticevectors = hexlv;
		hex->vertices = hexv;
		hex->numvertices = 6;
		return hex;
	}
	
	double roundness(int centerindex, double startradius, double stopradius)
	{
		//This constructs an annular ring around the grid point specified by centerindex
		//The relative deviation of the charge density inside that annular ring is returned
		//as a measure of "aromaticity"
		
		double psisq=0,psi=0;
		int pointsinside = 0;
		int i;
		double distance;
		for(i=0;i<numpoints;i++)
		{
			distance = pow(pow(*(xgrid+centerindex) - *(xgrid+i),2) + pow(*(ygrid+centerindex)-*(ygrid+i),2),0.5);
			if(distance >startradius && distance <stopradius)
			{
				psisq += pow(*(chargedensity+i),2);
				psi += *(chargedensity+i);
				pointsinside++;
			}
		}
		if(pointsinside == 0)
		return 1;
		
		psi /= pointsinside;
		psisq /= pointsinside;

//		return sqrt((psisq - pow(psi,2)))/psi;

		return sqrt((psisq - pow(psi,2)));
	}
	
	void angularroundness(int centerindex, double *r, double *theta, double *z)
	{
		//This generates the charge density as a function of r and theta
		
		
		int i;
		double distance;
		double angle;
		
		for(i=0;i<numpoints;i++)
		{
			if(i == centerindex)
			{
				r[i] = theta[i] = 0;
				z[i] = NAN;
				continue;
			}
			
			
			distance = pow(pow(*(xgrid+centerindex) - *(xgrid+i),2) + pow(*(ygrid+centerindex)-*(ygrid+i),2),0.5);
			angle = atan2(*(ygrid+centerindex)-*(ygrid+i),*(xgrid+centerindex) - *(xgrid+i));
			r[i] = distance;
			theta[i] = angle;
			if((*(fixedpoints+i) == 1))
			z[i] = 0;
			else
			z[i] = *(chargedensity+i);	
		}
		
	}
	
	
	void getgrid(tiling *tile, double *x, double *y)
	{
		int vec1,vec2;
		
		for(vec2 = 0;vec2<rows;vec2++)
		for(vec1 = 0;vec1<columns;vec1++)
		{
			*(x+vec2*columns+vec1) = vec1 * *(tile->latticevectors)+vec2 * *(tile->latticevectors+2);
			*(y+vec2*columns+vec1) = vec1 * *(tile->latticevectors+1)+vec2 * *(tile->latticevectors+3);
		}

	}
	
	void broadening(double *xin, double *yin, double *xout, double *yout, int numin, int numout, double broadening, double upperlimit)
	{
		/***********************************************/
		/* Convolves the input data y(xi) with a gaussian*/
		/* Assumes *xin and *(xin+numin-1) are the extremea
		 * 
		 *                                              */
		 
		 
		 double extreme1 = *(xin)-3*broadening;
		 double extreme2 = fmin(*(xin+numin-1),upperlimit);
		 
		 int in,out;
		 for(out = 0;out<numout;out++)
		 {
			*(yout + out) = 0;
			*(xout + out) = extreme1 + (extreme2-extreme1) * out /(numout-1);
			for(in = 0;in<numin;in++)
			{
				*(yout+out) += *(yin+in) * exp( - pow((*(xout+out) - *(xin+in))/broadening,2)/2);
			}
		 }
		 
		 
		 
	}
	
void normalize(double *wavefunction, double norm, int square, int numberofpoints, int *mask)
{
	int index;
	
	double sum = 0;
	for(index = 0;index<numberofpoints;index++)
	if(mask == NULL || !(*(mask+index) == 1))
	if(square)
	sum += *(wavefunction+index) * *(wavefunction+index);
	else
	sum += *(wavefunction+index);
	//printf("Sum is %g\n",sum);
	for(index = 0;index<numberofpoints;index++)
	if(mask == NULL || !(*(mask+index) == 1))
	if(square)
	*(wavefunction+index) *= sqrt(norm/sum);
	else
	*(wavefunction+index) *= norm/sum;
}

int gridbasedbc(int rows, int columns, int index, int shift, int boundaryconditions)
{
	/****************************************/
	/* Applies  boundary conditions         */
	/* To find the grid number which is     */
	/* shifted by shift from starting point */
	/* index based upon the type of boundary*/
	/* 0 open boundary, 2 periodic boundary */
	/* 3 terminated                         */
	
	int period;
	int sign = fabs(shift)/shift;
	
	if(fabs(shift) == columns)
	period = columns*rows;
	else
	period = columns;
	
	int finalindex;
	
	if((index + shift+period)/period == (index+period)/period)
	finalindex = index + shift;
	else
	switch(boundaryconditions)
	{
		case 0:
		finalindex = index - shift;//Bounce Back
		break;
		
		case 2:
		finalindex = index + shift  + -1 * sign * period;//Most important line
		break;
		
		case 3:
		finalindex = -1;
		break;
		}
	
	//printf("Index %d shift %d new index %d\n",index,shift,finalindex);
	return finalindex;
}
void nearestneighbors(int index, int *numnn, int *nn)
{
		*numnn = 4;
		*nn = gridbasedbc(rows,columns,index,1,*(fixedpoints+index));

 		*(nn+1) = gridbasedbc(rows,columns,index,-1,*(fixedpoints+index));

 		*(nn+2) = gridbasedbc(rows,columns,index,columns,*(fixedpoints+index));

 		*(nn+3) = gridbasedbc(rows,columns,index,-columns,*(fixedpoints+index));

		int compvectorindex;
		for(compvectorindex=0;compvectorindex<tiles->numcomp;compvectorindex++)
		{	
 		*(nn+*numnn) = gridbasedbc(rows,columns,index,*(tiles->compositevectors+2*compvectorindex),*(fixedpoints+index));
		if(*(nn+*numnn) != -1) //Not implemented yet, but for terminated boundary conditions.
		*(nn+*numnn) = gridbasedbc(rows,columns,*(nn+*numnn),columns * *(tiles->compositevectors+2*compvectorindex+1),*(fixedpoints+index));
		(*numnn)++;
		
		*(nn+*numnn) = gridbasedbc(rows,columns,index,-*(tiles->compositevectors+2*compvectorindex),*(fixedpoints+index));
		if(*(nn+*numnn) != -1) //Not implemented yet, but for terminated boundary conditions.
		*(nn+*numnn) = gridbasedbc(rows,columns, *(nn+*numnn),-columns * *(tiles->compositevectors+2*compvectorindex+1),*(fixedpoints+index));
		(*numnn)++;
		}
}


void fixedges(int *fixedpoints, int fixvalue)
{
	int row,column;

	for(row=0;row<rows;row++)
	{
		for(column=0;column<columns;column++)
		{
			if(row == 0 || row == rows -1 || column == 0 || column == columns -1)
			{
			*(fixedpoints+row * columns + column)=fixvalue; 
			}
		}
	}
	
}

void draw()
{
	
//	printf("Drawing\n");
	GtkAllocation *allocation = (GtkAllocation *)malloc(sizeof(GtkAllocation));		//Create the drawable pixmap
	gtk_widget_get_allocation(drawingarea, allocation);
	guint width = allocation->width;
	guint height = allocation->height;
	


	cairo_t *cr;
	cairo_surface_t *surface;
	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,width, height);
	
	cr = cairo_create (surface);
	cairo_set_source_rgb (cr, 0, 0, 0);
	cairo_rectangle(cr,0,0,width,height);
	cairo_fill(cr);
	
	
	
	int row, column;
	char *buffer = malloc(100*sizeof(char));
	int xmax = fmax(fmax(columns * *(tiles->latticevectors),rows * *(tiles->latticevectors+2)),columns * *(tiles->latticevectors)+rows * *(tiles->latticevectors+2));
	int ymax = fmax(fmax(columns * *(tiles->latticevectors+1),rows * *(tiles->latticevectors+3)),columns * *(tiles->latticevectors+1)+rows * *(tiles->latticevectors+3));
	int size = fmin(width/xmax,height/ymax);

//	int size = fmin(width/columns,height/rows);

	cairo_set_font_size(cr,size);
	
	int vertex;
	double x,y;
	for(row = 0;row<rows;row++)
	for(column = 0;column<columns;column++)
	{
		cairo_set_source_rgb(cr,1,1,1);
		
		switch(*(fixedpoints+row*columns+column))
		{
			case 1:
			cairo_set_source_rgb(cr,1,0,0);
			break;
			case 2:
			cairo_set_source_rgb(cr,0,1,0);
			break;
			case 3:
			cairo_set_source_rgb(cr,0,0,1);
			break;
			case 4:
			cairo_set_source_rgb(cr,1,1,0);
			break;
			case 5:
			cairo_set_source_rgb(cr,1,0,1);
			break;
			case 6:
			cairo_set_source_rgb(cr,0,1,1);
			break;
			case 7:
			cairo_set_source_rgb(cr,0,0,0);
			break;
			case 8:
			cairo_set_source_rgb(cr,0.5,0,0);
			break;
			case 9:
			cairo_set_source_rgb(cr,0,0.5,0);
			break;
		}
		/*
		if(*(fixedpoints + row*columns+column) == 1)
		cairo_set_source_rgb(cr,1,0,0);
		
		if(*(fixedpoints + row*columns+column) == 2)
		cairo_set_source_rgb(cr,0,1,0);
		
		if(*(fixedpoints + row*columns+column) == 3)
		cairo_set_source_rgb(cr,0,0,1);
		*/

		x = column* *(tiles->latticevectors) + row * *(tiles->latticevectors+2);
		y = column* *(tiles->latticevectors+1) + row * *(tiles->latticevectors+3);

		cairo_move_to(cr,x*size + (size-1) * *(tiles->vertices),y*size + (size-1) * *(tiles->vertices+1));
		for(vertex = 0;vertex<tiles->numvertices;vertex++)
		cairo_line_to(cr,x*size + (size-1) * *(tiles->vertices+2*vertex),y*size + (size-1) * *(tiles->vertices+1+2*vertex));
		
		cairo_close_path (cr);
		
		//cairo_rectangle(cr,column*size,row*size,size-1,size-1);
		
		cairo_fill(cr);
		
		if(*(potential+row*columns+column) != 0 && !(*(fixedpoints+row*columns+column) == 1))
		{
			sprintf(buffer,"%2g",*(potential+row*columns+column));//Make something to distinguish these two
			//cairo_move_to(cr, column*size,(row+1.0/2)*size);
			cairo_move_to(cr, (x-0.1)*size,(y+0.9)*size);
			cairo_set_source_rgb (cr, 0, 0, 0);
			cairo_show_text(cr,buffer);
		}
		else
		if(*(charge+row*columns+column) != 0 && !(*(fixedpoints+row*columns+column) == 1))
		{
			sprintf(buffer,"%2g",*(charge+row*columns+column));
			//cairo_move_to(cr, column*size,(row+1.0/2)*size);
			cairo_move_to(cr, (x-0.1)*size,(y+0.9)*size);
			cairo_set_source_rgb (cr, 0, 0, 0);
			cairo_show_text(cr,buffer);
		}
	}
	
	cairo_t *cr2;	
	cr2 = gdk_cairo_create (gtk_widget_get_window(drawingarea));
	cairo_set_source_surface(cr2, surface, 0, 0);
	cairo_paint(cr2);
	
	cairo_destroy(cr2);
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void printpm3dfile(char *filename,double *matrix1,double *matrix2, int numberOfRows, int numberOfColumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	FILE *outfile = fopen(filename,"w");
	int i;
	for (i = 0;i<numberOfRows*numberOfColumns ;i++ )
	{
		
			fprintf(outfile,"%.08g\t%.08g\n",*(matrix1+i),*(matrix2+i));
			if ((i+1) % numberOfColumns == 0)
			{
				fprintf(outfile,"\n");
			}
			
	}
	fclose(outfile);
}

void printsplotfile(char *filename,double *x, double *y, double *z, int numrows,int numcolumns)
{
	/****************************************************/
	/* This program outputs the contents of a matrix in */
	/* nice neat format.								*/
	/****************************************************/
	FILE *outfile = fopen(filename,"w");
	int row,column;
	int i;
	for (row = 0;row<numrows ;row++ )//Criss cross pattern
	{
	for(column = 0;column<numcolumns;column++)
	{
		i = row*numcolumns+column;
			fprintf(outfile,"%g %g %g\n",*(x+i),*(y+i),*(z+i));
	}
			fprintf(outfile,"\n");
	}
	
	for(column = numcolumns-1;column>=0;column--)
	{
	for (row = numrows-1;row>=0 ;row-- )
	{
		i = row*numcolumns+column;
			fprintf(outfile,"%g %g %g\n",*(x+i),*(y+i),*(z+i));
	}
			fprintf(outfile,"\n");
	}
	
	
	//This sections needs to be cleaned up heavily to be made more general right now it is specific to hex
	int compvec;
	//int v1dir,v2dir;//The sign of the vector 
	int start;
	for(compvec = 0;compvec<tiles->numcomp;compvec++)
	{
		
		//v1dir = *(tiles->compositevectors + 2*compvector) / fabs(*(tiles->compositevectors + 2*compvec));
		//v2dir = *(tiles->compositevectors + 2*compvector+1) / fabs(*(tiles->compositevectors + 2*compvec+1));

		for(start = 0;start < numrows;start++)
		{	
		for(row = start,column = 0;(row >= 0) && (column < numcolumns);column++,row--)
		{
			//printf("Sub comp\n");
			i = row*numcolumns+column;
			fprintf(outfile,"%g %g %g\n",*(x+i),*(y+i),*(z+i));
		}
			fprintf(outfile,"\n");
		}
		
		
		for(start = 0;start < numcolumns;start++)
		{	
		for(row = rows-1,column = start;(row >= 0) && (column < numcolumns);column++,row--)
		{
		//	printf("Sub comp\n");
			i = row*numcolumns+column;
			fprintf(outfile,"%g %g %g\n",*(x+i),*(y+i),*(z+i));
		}
			fprintf(outfile,"\n");
		}
		
	}
	
	//printf("THere\n");
	fclose(outfile);
}

/*
void makesquare(double voltage,double *potentials, int *fixedpoints, int columns, int rows, int squaresize)
{
	//Makes a square
	int row,column;
	
	int index;//The square will be traced out clockwise
	column = columns/2.0-squaresize/2.0;
	row = rows/2.0-squaresize/2.0;
	//printf("Starting at %d %d\n",column,row);
	for(index=0;index<squaresize-1;index++,column++)
	{
	*(potentials+column+row*columns) = voltage;
	*(fixedpoints+column+row*columns)=1;
	}
	for(index=0;index<squaresize-1;index++,row++)
	{
	*(potentials+column+row*columns) = voltage;
	*(fixedpoints+column+row*columns)=1;
	}
	for(index=0;index<squaresize-1;index++,column--)
	{
	*(potentials+column+row*columns) = voltage;
	*(fixedpoints+column+row*columns)=1;
	}
	for(index=0;index<squaresize-1;index++,row--)
	{
	*(potentials+column+row*columns) = voltage;
	*(fixedpoints+column+row*columns)=1;
	}
}

void makecircle(double voltage,double *potentials, int *fixedpoints, int columns, int rows, int radius, int xcenter, int ycenter)
{
	//Makes a circle
	int row,column;
	
	int index;//The circle will be traced out counterclockwise
	
	for(index=0;index<4*radius*M_PI;index++)//It will repeat points sometimes but should hit every point
	{
	column = xcenter + radius*cos(index/2.0/radius);
	row = ycenter + radius*sin(index/2.0/radius);
	*(potentials+column+row*columns) = voltage;
	if(fixedpoints != NULL)
	*(fixedpoints+column+row*columns)=1;
	}
}

void removefixedpointsincircle(int *fixedpoints, int rows, int columns, double radius, int xcenter, int ycenter)
{
		int x, y;
		
		for(x =0;x<2*radius;x++)
		for(y = 0;y<2*radius;y++)
		if((x-radius)*(x-radius) + (y-radius) * (y-radius) < radius * radius)
		{
		*(fixedpoints + (y+ycenter-(int)radius)*columns + x + xcenter-(int)radius) = 0;
		}
}

void makehexagon(double voltage, double *potentials, int *fixedpoints, int columns, int rows, int radius, int xcenter, int ycenter)
{
	int row,column;
	int x1,x2, y1, y2;
	double x,y;
	double theta;
	x1 = radius;
	x2 = radius * cos(M_PI/3);
	y1 = 0;
	y2 = radius * sin(M_PI/3);
	
	int symindex;
	int index;
	int maxpoints = radius*8;

	for(symindex=0;symindex<6;symindex++)	
	for(index = 0;index <maxpoints;index++)
	{
	x = x1 + (double)(x2-x1)/(maxpoints-1)*index;
	y = y1 + (double)(y2-y1)/(maxpoints-1)*index;
	theta = 2* M_PI * (double)symindex/6.0;
	column = xcenter + cos(theta)*x - sin(theta) * y ;
	row = ycenter + sin(theta) * x + cos(theta) * y;
	
	//printf("theta %g %g x %g  %d x %d\n",theta,x,y,column,row);
	//column = columns/2.0 + radius*cos(index/2.0/radius);
	//row = rows/2.0 + radius*sin(index/2.0/radius);
	*(potentials+column+row*columns) = voltage;
	*(fixedpoints+column+row*columns)=1;
	}
}
*/

void opengnuplotpipe(FILE *pipe,char *outputfilename,int delay)
{
	fprintf(pipe,"set term gif animate delay %d\nset output '%s'\nset key noautotitles\nset contour both\nset cntrparam levels auto 10\nset hidden3d\nset pm3d\n",delay,outputfilename);
}
void opengnuplotpipepicture(FILE *pipe, char *title)
{
	int interp;
	if(interspin != NULL)
	interp =  gtk_spin_button_get_value(GTK_SPIN_BUTTON(interspin));
	else
	interp = 5;
	//fprintf(pipe,"set key noautotitles\nset contour both\n set view 70,330 \nset cntrparam levels auto 10\nset hidden3d\nset pm3d\n");
	fprintf(pipe,"set key noautotitles\nset contour both\n set view 70,330 \nset cntrparam levels auto 10\nset view equal xy\nset hidden3d\nset pm3d ftriangles interpolate %d,%d\n",interp,interp);
	
	fprintf(pipe,"set key title '%s'\n",title);

}

void setinterpolate()
{
	int interp =  gtk_spin_button_get_value(GTK_SPIN_BUTTON(interspin));
	char buffer[1000];
	sprintf(buffer,"set pm3d ftriangles interpolate %d,%d\n",interp,interp);
	fprintf(potentialpipe,buffer);
	fprintf(fixedpipe,buffer);
	fprintf(wavefunctionpipe,buffer);
	fprintf(chargepotentialpipe,buffer);
	fprintf(baderchargepipe,buffer);

}

void closegnuplotpipe(FILE *pipe)
{
	//fprintf(pipe,"reread\n");
	fprintf(pipe,"exit");
	pclose(pipe);
}
void plot3d(FILE *pipe,char *filename)
{
	fprintf(pipe,"splot '%s' w pm3d\n",filename);
	//fprintf(pipe,"splot '%s'\n",filename);
	//printf("Plotting\n");
}

void splot3d(FILE *pipe,char *filename)
{
	fprintf(pipe,"splot '%s' u 1:2:3 w pm3d\n",filename);
	
//	fprintf(pipe,"splot '%s' u 1:2:3 w lines\n",filename);
	//fprintf(pipe,"splot '%s'\n",filename);
	//printf("Plotting\n");
}
/*
makepicture(char *prefix,double *potentials,double *fields,int *fixedpoints,int numdimensions,int *dimensions,int numpoints)
{
	int numiterations = 1000;
	printf("Running %s\n",prefix);
	
	char *datafilename = malloc(1000*sizeof(char));
	char *picturefilename = malloc(1000*sizeof(char));
	char *deletecommand = malloc(50*numiterations*sizeof(char));
	sprintf(picturefilename,"%s.png",prefix);
	sprintf(deletecommand,"rm ");
	sprintf(datafilename,"%s.txt",prefix);
	
	FILE *pipe = popen("gnuplot","w");
	opengnuplotpipepicture(pipe,picturefilename);
	int index;
	for(index=0;index<numiterations;index++)
	{
	printf("%d of %d iterations\r",index,numiterations);
	LaplacianStep(potentials,fixedpoints,numdimensions,dimensions);
	}
	
	tofield(potentials,fields,numdimensions,dimensions);
	printpm3dfile(datafilename,potentials,fields, *(dimensions+1), *dimensions);
	plot3d(pipe,datafilename,index);
	strcat(deletecommand,datafilename);
	strcat(deletecommand," ");
	
	closegnuplotpipe(pipe);
	char *openpicture = malloc(500*sizeof(char));
	sprintf(openpicture,"eog %s &\n",picturefilename);
	system(openpicture);
	system(deletecommand);
}*/

void findeigenvalues(double *Hamiltonian, int rows,int columns, double *eigenvalues, double *eigenvectors)
{
	/****************************************************************************/
	/* This function uses the dgeev subroutine from lapack to find the          */
	/* eigenvalues and eigenvectors of the Hamiltonian matrix. It calls the     */
	/* bubblesort function to sort the eigenvalues and eigenvectors, and        */
	/* finally prints out the resulting values and vectors.                     */
	/****************************************************************************/
	
	printf("Finding the eigenvalues of the Hamiltonian Matrix\n");
	/*char jobvl ='v';
	char jobvr = 'n';
//This symbolizes that the left and right eigenvectors will be found
	double *realEigenValues,*imaginaryEigenValues;
	realEigenValues = malloc(numberOfRows*sizeof(double));
	imaginaryEigenValues = malloc(numberOfRows*sizeof(double));
	
	//double *rightEigenVectors = malloc(numberOfRows*numberOfColumns*sizeof(double));	
	//Using the passed eigenvector structure instead
	//double *leftEigenVectors = malloc(numberOfRows*numberOfColumns*sizeof(double));
	
	//For a system with real eigenvalues the leftEigenVectors should match the right eigenvectors so only the right eigenvectors are printed
	int lwork = numberOfRows * 4;
	double *work = malloc(lwork*sizeof(double));
//This allocates workspace to be used by dgeev. The recomended space is 4 times the dimension
	int info = 0;//After dgeev info = 0 if the calculation went correctly
	
	
	double *fortranMatrix = malloc(numberOfRows*numberOfColumns*sizeof(double));
	ctofortran(Hamiltonian, fortranMatrix, numberOfRows, numberOfColumns);
	
	dgeev_(&jobvl,&jobvr,&numberOfRows,fortranMatrix,&numberOfRows,eigenvalues,imaginaryEigenValues,eigenvectors,&numberOfRows,NULL,&numberOfRows,work,&lwork,&info);
	*/
	
	
	//double *eigenvalues = malloc(rows*sizeof(double));

	char jobz = 'v';
	char uplo = 'u';
	int lwork = rows * rows * 4;
	double *work = malloc(lwork*sizeof(double));
	int liwork = 6 * rows + 3;
	double *iwork = malloc(liwork*sizeof(double)); 
	int info = 0;//After dsyevd info = 0 if the calculation went correctly	
	
	dsyevd_( &jobz, &uplo, &rows, Hamiltonian, &columns, eigenvalues, work, &lwork, iwork,  &liwork, &info );
	
	memcpy(eigenvectors,Hamiltonian,rows*columns*sizeof(double));
	//double *leftEigenVectors = matrix;
	free(work);
	free(iwork);
	if (info != 0)
	{
		printf("There was a problem calculating the eigenvalues\n");
		exit(0);
	}

	/*
	for(index=0;index<numdospoints;index++)
	{
		*(dos+index) = 0;
		lowerbound = mineigenvalue + (maxeigenvalue-mineigenvalue)*index/(numdospoints);
		upperbound = mineigenvalue + (maxeigenvalue-mineigenvalue)*(index+1)/(numdospoints);
		
		*(dosenergies + index) = (lowerbound+upperbound)/2;
		for(index2=0;index2<numberOfRows;index2++)
		{
			if(*(eigenvalues+index2) >= lowerbound && *(eigenvalues+index2) <= upperbound)
			(*(dos+index))++;
		}
	}*/	
}

void	finitedifferencematrix(double *FDM)
{
	int index;
	
	for(index = 0;index<numpoints*numpoints;index++)
	*(FDM+index) = 0;
	
	int compvectorindex;
	int row, column;
	int newindex;
	for(row = 0;row<rows;row++)
	for(column = 0;column<columns;column++)
	{
		index = row * columns+column;
		if(*(fixedpoints+index) == 1)
		continue;
		
 		newindex = gridbasedbc(rows,columns,index,1,*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		*(FDM + (index) * numpoints + newindex) += -1; 

 		newindex = gridbasedbc(rows,columns,index,-1,*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		*(FDM + (index) * numpoints + newindex) += -1; 

 		newindex = gridbasedbc(rows,columns,index,columns,*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		*(FDM + (index) * numpoints + newindex) += -1; 

 		newindex = gridbasedbc(rows,columns,index,-columns,*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		*(FDM + (index) * numpoints + newindex) += -1; 

		for(compvectorindex=0;compvectorindex<tiles->numcomp;compvectorindex++)
		{	
 		//printf("Here the number %d %d %d\n",compvectorindex,*(tiles->compositevectors+2*compvectorindex),*(tiles->compositevectors+2*compvectorindex+1));
 		newindex = gridbasedbc(rows,columns,index,*(tiles->compositevectors+2*compvectorindex),*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		newindex = gridbasedbc(rows,columns,newindex,columns * *(tiles->compositevectors+2*compvectorindex+1),*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		*(FDM + (index) * numpoints + newindex) += -1;  
		
		newindex = gridbasedbc(rows,columns,index,-*(tiles->compositevectors+2*compvectorindex),*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		newindex = gridbasedbc(rows,columns, newindex,-columns * *(tiles->compositevectors+2*compvectorindex+1),*(fixedpoints+index));
		if(newindex != -1) //Not implemented yet, but for terminated boundary conditions.
		*(FDM + (index) * numpoints + newindex) += -1;  
		}
		
		*(FDM + index* numpoints + index) = 4.0 + 2 * tiles->numcomp + *(potential+index) + *(chargepotential+index);


	}
	
}
FILE *openpipe()
    {
    
    FILE *pipe = popen("gnuplot","w");
    
    char *gnuplotscript = 
    "set xlabel 'X'\n"
    "set ylabel 'Y'\n";
    
    fprintf(pipe,gnuplotscript);
    
    return pipe;
    }

void plotline(char *filename,FILE *pipe,double fermienergy)
    {
    char command[100000];
    if(fermienergy != 0)
    sprintf(command,"set zeroaxis\nplot '%s' u ($1 - %g):($2) with lines\n",fermienergy,filename);
    else
    sprintf(command,"plot '%s' u 1:2 with lines\n",filename);

    fprintf(pipe,command);
    fflush(pipe);
    }
void plotyy(char *filename,FILE *pipe,double fermienergy)
    {
    char command[100000];
    if(fermienergy != 0)
    sprintf(command,"set zeroaxis\nset y2tics\nplot '%s' u ($1 - %g):($2) axes x1y1 with lines,'%s' u ($1 - %g):($3) axes x1y2 with lines\n",fermienergy,filename,fermienergy,filename);
    else
    sprintf(command,"set y2tics\nplot '%s' u 1:2 axes x1y1 with lines,'%s' u 1:3 axes x1y2 with lines\n",filename,filename);

    fprintf(pipe,command);
    fflush(pipe);
    }
    
    void plotlines(char *filename,FILE *pipe,double fermienergy, int numlines)
    {
    char command[100000];
    int index;
    char buffer1[1000],buffer2[1000];
    sprintf(buffer1,"");
    
    for(index = 0;index<numlines;index++)
    {
    if(fermienergy != 0)
    sprintf(buffer2,"%s '%s' u ($1 - %g):($%d) with lines,",buffer1,fermienergy,filename,index+2);
    else
    sprintf(buffer2,"%s '%s' u 1:%d with lines,",buffer1,filename,index+2);
	strcpy(buffer1,buffer2);
	
    }
    if(fermienergy != 0)
    sprintf(command,"set zeroaxis\nplot %s\n",fermienergy,buffer1);
    else
    sprintf(command,"plot %s\n",buffer1);
	
    fprintf(pipe,command);
    fflush(pipe);
    }

void openpipes()
{
	wavefunctionpipe = popen("gnuplot","w");
	opengnuplotpipepicture(wavefunctionpipe,"Wavefunctions");

	fixedpipe = popen("gnuplot","w");
	opengnuplotpipepicture(fixedpipe,"Boundary");
	
	potentialpipe = popen("gnuplot","w");
	opengnuplotpipepicture(potentialpipe,"Potential");
	
	chargepotentialpipe = popen("gnuplot","w");
	opengnuplotpipepicture(chargepotentialpipe,"Charge Potential");
	
	baderchargepipe = popen("gnuplot","w");
	opengnuplotpipepicture(baderchargepipe,"Bader Charge");
	fprintf(baderchargepipe,"unset pm3d\nunset hidden3d\nset view map\n");
	
	dospipe = openpipe();
	pointspecpipe = openpipe();

}

void runplot(char *filename, double *matrix, FILE *pipe)
{
		//printpm3dfile(filename,matrix, matrix, rows, columns);
		//plot3d(pipe,filename);
		printsplotfile(filename,xgrid,ygrid,matrix,rows,columns);
		splot3d(pipe,filename);
		
		fflush(pipe);
}

void runspheres(char *baderfilename, char *chargedensityfilename, double *matrix, FILE *pipe)
{
		//printpm3dfile(filename,matrix, matrix, rows, columns);
		//plot3d(pipe,filename);
		printsplotfile(baderfilename,xgrid,ygrid,matrix,rows,columns);
		//splot3d(pipe,filename);
		//fprintf(pipe,"splot '%s' u ($1):($2):(1):(($3)**(0.3333)*2) ps variable pt 7\n",filename);
		fprintf(pipe,"splot '%s' u 2:1:3 w pm3d, '%s' u ($2):($1):(0):(3*sqrt($3)) w points ps variable pt 7 lc rgb 'black'\n",chargedensityfilename,baderfilename);
		
		fflush(pipe);
}


static void
fixed_toggled (GtkCellRendererToggle *cell,
          gchar                 *path_str,
          gpointer               data)
{
  GtkTreeModel *model = (GtkTreeModel *)eigenvaluelist;
  GtkTreeIter  iter;
  GtkTreePath *path = gtk_tree_path_new_from_string (path_str);
  gboolean fixed;

  gtk_tree_model_get_iter (model, &iter, path);
  gtk_tree_model_get (model, &iter, 1, &fixed, -1);

	int *indices = gtk_tree_path_get_indices(path);
  fixed ^= 1;

 *(displayedwavefunctions + *indices) = fixed;
 
  gtk_list_store_set (GTK_LIST_STORE (model), &iter, 1, fixed, -1);

  gtk_tree_path_free (path);
  
	plotwavefunction();
}

/*void toggledisplay(int index)
{
	
	 printf("Toggling %d\n",index);
}*/

void calculate()
{
	

	int numiterations = 1;
//	printf("Running %s\n",prefix);
	int index;
	
	double *FDM = malloc(numpoints*numpoints*sizeof(double));
	finitedifferencematrix(FDM);
	
	//printmatrix(FDM,numpoints,numpoints);
	
	
	int numfixedpoints;
	for(index=0,numfixedpoints=0;index<numpoints;index++)
	if(*(fixedpoints+index) == 1)
	numfixedpoints++;
	
	double *condensedmatrix = malloc((numpoints-numfixedpoints) *(numpoints-numfixedpoints) * sizeof(double));
	maskcompress(FDM,condensedmatrix,2,numpoints,numfixedpoints,fixedpoints);
	//printf("Here is the compressed form \n");
	//printmatrix(condensedmatrix,numpoints-numfixedpoints,numpoints-numfixedpoints);
	
	//maskdecompress(FDM,condensedmatrix,2,numpoints,numfixedpoints,fixedpoints);
	//printf("Here is the decompressed form \n");
	//printmatrix(FDM,numpoints,numpoints);
	
	if(condensedeigenvalues != NULL)
	free(condensedeigenvalues);
	if(condensedeigenvectors != NULL)
	free(condensedeigenvectors); 
	condensedeigenvalues = malloc((numpoints-numfixedpoints)*sizeof(double));	
	condensedeigenvectors = malloc( (numpoints-numfixedpoints)*(numpoints-numfixedpoints)*sizeof(double));
	findeigenvalues(condensedmatrix,numpoints-numfixedpoints,numpoints-numfixedpoints,condensedeigenvalues,condensedeigenvectors);
	
	//printf("Decompressing eigenvalues\n");
	maskdecompress(eigenvalues,condensedeigenvalues,1,numpoints,numfixedpoints,fixedpoints);
	//printf("Decompressing eigenvectors\n");
	maskdecompress(eigenvectors,condensedeigenvectors,2,numpoints,numfixedpoints,fixedpoints);
	
	//printf("Compressed eigenvectors\n");
	//printmatrix(condensedeigenvectors,numpoints-numfixedpoints,numpoints-numfixedpoints);
	//printf("Decompressed eigenvectors\n");
	//printmatrix(eigenvectors,numpoints,numpoints);
	//findeigenvalues(FDM,numpoints,numpoints,eigenvalues,eigenvectors);
	
	//free(condensedeigenvalues);
//	free(condensedeigenvectors);
	free(condensedmatrix);
	
	
	//bubblesort(eigenvalues,eigenvectors,numberOfRows,numberOfColumns*sizeof(double));
	bubblesort(eigenvalues,eigenvectors,numpoints,numpoints*sizeof(double));


	dos();
	
	//printf("here\n");

	
	//vectorviewer(eigenvalues,eigenvectors,fixedpoints);

	free(FDM);
//	free(eigenvalues);
//	free(eigenvectors);

	for(index = 0;index<numpoints;index++)
	{
		normalize(eigenvectors + index*numpoints,1,1,numpoints,fixedpoints);
	}
	
	//printf("End of calculate function");
}

double localizationindex(int *usedgrid, int numgridpoints)
{
	double charge=0, numpairs=0;
	int i,j;
	
	for(i=0;i<numgridpoints;i++)
	{
		charge += *(chargedensity+*(usedgrid+i));
		for(j=i;j<numgridpoints;j++)
		numpairs += *(pairdensity + *(usedgrid+i) * numpoints + *(usedgrid+j));//The diagonal term
	}
	printf("The charge is %g\n",charge);
	printf("The number of pairs is %g\n",numpairs);
	double localization = charge-numpairs/charge;
	printf("localization is %g\n",localization);
	
	return localization;
}

void pd()
{
	/* Calculates the pair density of the system
	 * Also generates the global localization index
	 * 
	 * 
	 * */
	 int point1, point2,vector1,vector2,numused1,numused2;
	 for(point1 = 0;point1<numpoints;point1++)
	 for(point2 = 0;point2<numpoints;point2++)
	 *(pairdensity+point1*numpoints+point2) = 0;
	 
	 
	 for(point1 = 0;point1<numpoints;point1++)
	 {
	 if((*(fixedpoints + point1) == 1))
	 continue;
	 for(point2 = point1;point2<numpoints;point2++)
	 {
	  if((*(fixedpoints + point2) == 1))
	 continue;
	 \
	 numused1 = 0;
	 for(vector1 = 0;vector1<numpoints;vector1++)
	{
		if(*(eigenvalues + vector1) == 0)
		continue;
	numused2 = 0;
	for(vector2 = vector1;vector2<numpoints;vector2++)
	{	
		if(*(eigenvalues + vector2) == 0)
		continue;
		
		if(*(displayedwavefunctions + numused1) && *(displayedwavefunctions + numused2))
		{
		//printf("Adding to pair %d,%d, wavefunction %d, %g and wavefunction %d, %g\n",point1,point2,vector1,*(eigenvectors+numpoints*vector1 + point1),vector2,*(eigenvectors+numpoints*vector2 + point2));
		//*(pairdensity+point1*numpoints+point2) += *(eigenvectors+numpoints*vector1 + point1) * *(eigenvectors+numpoints*vector2 + point2);
		*(pairdensity+point1*numpoints+point2) += pow(*(eigenvectors+numpoints*vector1 + point1),2) * pow(*(eigenvectors+numpoints*vector2 + point2),2)
		+ pow(*(eigenvectors+numpoints*vector1 + point2),2) * pow(*(eigenvectors+numpoints*vector2 + point1),2)
		- 2 * *(eigenvectors+numpoints*vector1+point1) * *(eigenvectors+numpoints*vector1+point2) * *(eigenvectors+numpoints*vector2+point1) * *(eigenvectors+numpoints*vector2+point2);
		}
		numused2++;
	}
		numused1++;
	}
	*(pairdensity+point2*numpoints+point1) = *(pairdensity+point1*numpoints+point2);
}}
	
	
	int numfixedpoints,index;
	for(index=0,numfixedpoints=0;index<numpoints;index++)
	if(*(fixedpoints+index) == 1)
	numfixedpoints++;
	
	double compressedpairdensity[numpoints*numpoints];
	maskcompress(pairdensity,compressedpairdensity,2,numpoints,numfixedpoints,fixedpoints);
	printf("Here is the pair density matrix\n");
	printmatrix(compressedpairdensity,numpoints-numfixedpoints,numpoints-numfixedpoints);
	
	int fullrange[numpoints];
	int used = 0;
	for(point1=0;point1<numpoints;point1++)
	{
	if(*(fixedpoints+point1))
	continue;
	fullrange[used] = point1;
	used++;
	}
	localizationindex(fullrange,numpoints-numfixedpoints);
	
}

void dos()
{
		int numdospoints = 2000;
	double *dos = malloc(numdospoints*sizeof(double));
	double *integrateddos = malloc(numdospoints*sizeof(double));
	double *dosenergies = malloc(numdospoints*sizeof(double));
	//double charge = 0;
//	double mineigenvalue = *(eigenvalues);
//	double maxeigenvalue = *(eigenvalues + numpoints-1);
//	double lowerbound, upperbound;
	
	double one[numpoints];
	int index;
	for(index=0;index<numpoints;index++)
	one[index] = 1;	
	
	int numfixedpoints;
	for(index=0,numfixedpoints=0;index<numpoints;index++)
	if(*(fixedpoints+index) == 1)
	numfixedpoints++;
	
	double broad = gtk_spin_button_get_value(broadeningspin);
	double upperlimit = gtk_spin_button_get_value(upperlimitspin);
//	broadening(eigenvalues+numfixedpoints,one+numfixedpoints,dosenergies,dos,numpoints-numfixedpoints,numdospoints,broad);	
	broadening(condensedeigenvalues,one,dosenergies,dos,numpoints-numfixedpoints,numdospoints,broad,upperlimit);	

	//for(index = 0;index<numpoints;index++)
	//if(*(displayedwavefunctions + index))
	//charge += 2;
	
	normalize(dos,numpoints-numfixedpoints,0,numdospoints,NULL);
	
	*integrateddos = 0;
	for(index=1;index<numdospoints;index++)
	integrateddos[index] = dos[index] + integrateddos[index-1];


	double pointspec[numdospoints];
	double pointspecenergies[numdospoints];
	double pointldos[numpoints-numfixedpoints];
	
	for(index=0;index<numpoints;index++)
	pointldos[index] = 0;
	
	int pointspecindex = -1;
	int j;
	for(index=0;index<numpoints;index++)
	if(*(fixedpoints+index) == 4)//pointspec
	{
		pointspecindex = index;
		for(j=0;j<numpoints-numfixedpoints;j++)
		pointldos[j] += pow(*(condensedeigenvectors+j*(numpoints-numfixedpoints) + index),2);
	}
	if(pointspecindex != -1)
	broadening(condensedeigenvalues,pointldos,pointspecenergies,pointspec,numpoints-numfixedpoints,numdospoints,broad,upperlimit);	

	
	
	
	
	
	GtkTreeIter iter;
	gtk_list_store_clear (eigenvaluelist);
	//GtkWidget *togglebutton;
//	GtkAdjustment *occupancyadjustment = gtk_adjustment_new (0, 0, 2, 1, 1, 1);
	int numused = 0;
	for(index=0;index<numpoints;index++)
	{
		if(*(eigenvalues+index) != 0)
		{
		gtk_list_store_append (eigenvaluelist,&iter);
		//gtk_list_store_set (eigenvaluelist, &iter, 0,  *(eigenvalues+index),-1);
		//togglebutton = gtk_toggle_button_new();
		//g_signal_connect_swapped(togglebutton,"toggled",toggledisplay,index);
		gtk_list_store_set (eigenvaluelist, &iter, 0, *(eigenvalues+index),1,*(displayedwavefunctions+numused),-1);
//		gtk_list_store_set (eigenvaluelist, &iter, 0, *(eigenvalues+index),1,(double)1,-1);
//		gtk_list_store_set (eigenvaluelist, &iter, 0, *(eigenvalues+index),1,"1.0",-1);

//		gtk_list_store_set (eigenvaluelist, &iter, 0, *(eigenvalues+index),-1);
		numused++;
		}
	}
	//gtk_tree_view_columns_autosize(treeview);
	
	//for(index = 0;index<numpoints;index++)
	//printf("Eigenvalue is %g\n",*(eigenvalues+index));
	
	FILE *outfile = fopen("dos.dat","w");
	for(index=1;index<numdospoints;index++)
	{
	//fprintf(outfile,"%g %g\n",*(dosenergies+index),*(dos+index));
	fprintf(outfile,"%g %g %g\n",*(dosenergies+index),*(dos+index),*(integrateddos+index));
	}
	fclose(outfile);

	//plotline("dos.dat",openpipe(),0);
	//plotline("dos.dat",dospipe,0);
	plotyy("dos.dat",dospipe,0);
	
	//plotlines("dos.dat",dospipe,0,2);

	if(pointspecindex != -1)
	{
	outfile = fopen("pointspec.dat","w");
	for(index=1;index<numdospoints;index++)
	{
	fprintf(outfile,"%g %g\n",*(pointspecenergies+index),*(pointspec+index));
	}
	fclose(outfile);

	//plotline("dos.dat",openpipe(),0);
	plotline("pointspec.dat",pointspecpipe,0);
}


	free(dos);
	free(dosenergies);

	}
	
void bader(double totalcharge)
	{
		//Calculates the bader charges of the charge density
		int i,j;
		double epsilon = 1e-4;
		int maxbasins = 100;
		double *chargedensitycopy = malloc(numpoints*sizeof(double));
		double *weights = malloc(maxbasins*numpoints*sizeof(double));
		double *charges = malloc(maxbasins*sizeof(double));
		int *centers = malloc(maxbasins*sizeof(int));
		int *indices = malloc(numpoints*sizeof(int));
		int **neighbors = malloc(numpoints*sizeof(int *));
		
		
		int type = gtk_combo_box_get_active(tiletype);
		double a,l;
		int numneighbors;
		switch (type){
			case 0:
			a = l = 1;
			numneighbors = 4;
			break;
			
			case 1:
			a = 1/sqrt(3);
			l = 1;
			numneighbors = 6;
			break;
			
			default:
			a = l = 1;
			numneighbors = 4;
			break;
			};
		double sum;	
		double *fluxes = malloc(numneighbors*numpoints*sizeof(double));
		double *unnormalizedfluxes = malloc(numneighbors*sizeof(double));
			//printf("calculating fluxes\n");
		for(i=0;i<numpoints;i++)
		{
			*(neighbors+i) = malloc(numneighbors*sizeof(int));
			
			if(*(fixedpoints+i) == 1)
			{
				for(j=0;j<numneighbors;j++)
				{
					*(fluxes+i*numneighbors+j) == 0;
					
				}
				continue;
			}
			nearestneighbors(i,&numneighbors,*(neighbors+i));
			sum = 0;
			for(j=0;j<numneighbors;j++)
			{
				*(unnormalizedfluxes+j) = a/l * (*(chargedensity+*(*(neighbors+i)+j)) - *(chargedensity+i));
				if(*(unnormalizedfluxes+j)<0)
				*(unnormalizedfluxes+j) = 0;
				sum += *(unnormalizedfluxes+j);
			}
			for(j=0;j<numneighbors;j++)
			{
				if(sum > 0)
				*(fluxes + i*numneighbors+j) = *(unnormalizedfluxes+j)/sum;
				else
				*(fluxes + i*numneighbors+j) = 0;
				//printf("Flux from %d to %d is %g\n",i,*(*(neighbors+i)+j),*(fluxes + i*numneighbors+j));
			}
		}
		//printf("Allocating and freeing\n");
		
		
		int numbasins = 0;
		
		for(i=0;i<numpoints;i++)
		for(j=0;j<maxbasins;j++)
		*(weights+i*maxbasins+j) = 0;
		for(i=0;i<numpoints;i++)
		{
		*(chargedensitycopy+i) = *(chargedensity+i);
		*(indices+i) = i;
		}
		//printf("Sorting Charge density\n");
		bubblesortthreedatas(chargedensitycopy,fluxes,neighbors,indices,numpoints,numneighbors*sizeof(double),sizeof(int *),sizeof(int));
		
		int center;
		int k;
		int neighboriscenter;
		//printf("Calculating weights\n");
		for(i=0;i<numpoints;i++)
		{
			if(*(fixedpoints+*(indices+i)) == 1)
			continue;
				center = 1;
				//printf("Charge density is %g for index %d\n",*(chargedensitycopy+i),*(indices+i));
				for(j=0;j<numneighbors;j++)
				{
					//if(numbasins >0)
					//printf("Center %d neighbor %d\n",*(centers+numbasins-1),*(*(neighbors+i)+j));
					
					neighboriscenter = -1;
					for(k=0;k<numbasins;k++)
					if((*(centers+k) == *(*(neighbors+i)+j)))
					{
					neighboriscenter = k;
					//printf("neighbor %d is center\n",*(centers+k));
					}
					
					//This is a clusterfuck
					if((*(*(neighbors+i)+j) >= 0) && (*(*(neighbors+i)+j) < numpoints))
					if(neighboriscenter != -1 || (*(chargedensity+*(*(neighbors+i)+j)) - *(chargedensitycopy+i) > epsilon) || (*(chargedensitycopy+i) == 0) )
					{
						//printf("neighbor greater %d, thispoint zero %d\n",(*(chargedensity+*(*(neighbors+i)+j)) - *(chargedensitycopy+i) > epsilon),(*(chargedensitycopy+i) == 0) );
						center = 0;
						for(k=0;k<numbasins;k++)
						{
						if(neighboriscenter != -1 && !(*(chargedensity+*(*(neighbors+i)+j)) - *(chargedensitycopy+i) > epsilon))
						{
						if(k == neighboriscenter)
						{
						*(weights+k*numpoints+*(indices+i)) = 1;
						break;
						}
						else
						*(weights+k*numpoints+*(indices+i)) = 0;
						}
						else
						*(weights+k*numpoints + *(indices+i)) += *(fluxes+i*numneighbors+j) * *(weights+k*numpoints+*(*(neighbors+i)+j));
						//printf("Weight for basin %d is %g from neighbor weight %g and neighbor %d\n",k,*(weights+k*numpoints + *(indices+i)),*(weights+k*numpoints+*(*(neighbors+i)+j)),*(*(neighbors+i)+j));
						}
					}
					
				}
				if(center)
					{
						//printf("%d is a center \n",*(indices+i));
						*(centers+numbasins) = *(indices+i);
						*(weights+numbasins*numpoints+*(indices+i)) = 1;
						numbasins++;
					}
		}
		printf("Found %d basins\n",numbasins);
	//	printf("Calculating charges\n");
		double chargesum = 0;
		for(i=0;i<numbasins;i++)
		{
		*(charges+i) = 0;
		for(j=0;j<numpoints;j++)
		*(charges+i) += 2 * *(weights+i*numpoints+j) * *(chargedensity+j) * 1;//This 1 could be replaced by a Voronoi volume/area
		//A charge of 2 is used because at present we expect any occupied state to be doubly occupied
		printf("Bader charge at point %d location %g x %g is %g\n",*(centers+i),*(xgrid+*(centers+i)),*(ygrid+*(centers+i)),*(charges+i));
		chargesum += *(charges+i);
		}
		if(fabs(chargesum-totalcharge) > 0.01)
		printf("There is charge missing, bader charges total %g sum of charge density is %g\n",chargesum,totalcharge);
		
		double *badercharge = malloc(numpoints*sizeof(double));
		for(i=0;i<numpoints;i++)
		*(badercharge+i) = 0;
		for(i=0;i<numbasins;i++)
		*(badercharge + *(centers+i)) = *(charges+i);
		
		runspheres("badercharge.dat","wavefunction.dat",badercharge,baderchargepipe);
		
		free(fluxes);
		free(unnormalizedfluxes);
		free(charges);
		free(chargedensitycopy);
		free(weights);
		free(centers);
		free(indices);
		for(i=0;i<numpoints;i++)
		free(*(neighbors+i));
		free(neighbors);
		//printf("Done badering\n");
	}

void calculatenocharge()
{
	int index;
	for(index = 0;index<numpoints;index++)
	{
	*(chargedensity + index) = 0;
	*(chargepotential+index) = 0;
	}
	calculate();
}

double iterate()
{
	chargepotentialscalefactor = gtk_spin_button_get_value(GTK_SPIN_BUTTON(chargepotentialspin));
	//printf("Charge potential scale factor is %g\n",chargepotentialscalefactor);
	int pointindex;
	for(pointindex = 0;pointindex<numpoints;pointindex++)
	chargedensity[pointindex] = 0;
	
	int vectorindex;
	double oldenergy = 0;
	int numused = 0;
	double totalcharge = 0;
	for(vectorindex = 0;vectorindex<numpoints;vectorindex++)// Calculate charge density
	{
		if(*(eigenvalues + vectorindex) == 0)
		continue;
		
		if(*(displayedwavefunctions + numused))
		{
		oldenergy += *(eigenvalues+vectorindex);
		for(pointindex = 0;pointindex<numpoints;pointindex++)
		{
			if(! (*(fixedpoints+pointindex) == 1))
			{
			printf("Adding to location %d, wavefunction %d, %g\n",pointindex,vectorindex,*(eigenvectors + numpoints*vectorindex+pointindex));
			*(chargedensity+pointindex) += pow(*(eigenvectors + numpoints*vectorindex+pointindex),2);
			totalcharge += pow(*(eigenvectors + numpoints*vectorindex+pointindex),2);
			}
		}
		}
		numused++;
	}
	printf("Total charge is %g\n",totalcharge);
	bader(2*totalcharge);
	pd();
	
	int pointindex2;
	double distance;
	
	for(pointindex = 0;pointindex<numpoints;pointindex++)
	{
	*(chargepotential+pointindex) = 0;
	if(! (*(fixedpoints+pointindex) == 1))
	for(pointindex2 = 0;pointindex2 < numpoints;pointindex2++)
	//if(pointindex2 != pointindex)
	{
		if(*(fixedpoints+pointindex2) == 1)
		{
		*(chargepotential+pointindex) += -(2 * *(charge+pointindex2)) * chargepotentialscalefactor/(distance+1);//Factor of 2 here for double occupancy
		continue;
		}
		distance = pow(*(xgrid+pointindex)-*(xgrid+pointindex2),2) + pow(*(ygrid+pointindex)-*(ygrid+pointindex2),2);
		distance = sqrt(distance);
		*(chargepotential+pointindex) += (*(chargedensity+pointindex2)-2 * *(charge+pointindex2)) * chargepotentialscalefactor/(distance+1);//Factor of 2 here for double occupancy
	}
	}
	plotchargepotential();
	calculate();
	
	
	double newenergy = 0;
	
	for(vectorindex = 0,numused = 0;vectorindex<numpoints;vectorindex++)// Calculate charge density
	{
		if(*(eigenvalues + vectorindex) == 0)
		continue;
		
		if(*(displayedwavefunctions + numused))
		{
		printf("Adding eigenvalue %g\n",*(eigenvalues+vectorindex));
		newenergy += *(eigenvalues+vectorindex);
		}
		numused++;
	}
	printf("Old energy %g new energy %g convergence %g\n",oldenergy,newenergy,fabs(newenergy-oldenergy));
	
	return newenergy-oldenergy;
}

void converge()
{
	double convergencecriteria = 0.1;
	

	while(fabs(iterate()) > convergencecriteria || fabs(iterate()) > convergencecriteria);

//	while(fabs(iterate()) > convergencecriteria);

	printf("Done converging\n");
	//Make aromatic graph
	
}

void calculateroundness()
{
	
	int i,j;
	int *numnn = malloc(sizeof(int));
	int *nn = malloc(6*sizeof(int));
	int aromaticcenter;
	
	//double maxradius = 5;
	double radiusstep = 0.01;
	//double thickness = 1.5;
	
	double maxradius = gtk_spin_button_get_value(maximumradiusspin);
	double thickness = gtk_spin_button_get_value(roundnessthicknessspin);
	int numroundnesspoints = maxradius/radiusstep + 1;
	double *round = malloc(numroundnesspoints * sizeof(double));
	double *radii = malloc(numroundnesspoints * sizeof(double));
	
	FILE *outfile = fopen("roundness.dat","w");
	FILE *aromaticpipe = openpipe();
	
	FILE *roundness3Dpipe;
	roundness3Dpipe	= openpipe();
		fprintf(roundness3Dpipe,"set key noautotitles\nset contour both\n set view 70,330 \nset cntrparam levels auto 10\nset hidden3d\nset pm3d ftriangles interpolate 5,5\n");
	
		fprintf(roundness3Dpipe,"set xlabel 'R'\nset ylabel '{/Symbol f}'\n");
//		fprintf(roundness3Dpipe,"set xrange [:%g]\n",rmax);
		
	char *roundness3Dfilename = "roundness3D";
	
	int maxaromatics = 100;
	int numaromatics = 0;
	double *r = malloc(maxaromatics*numpoints*sizeof(double));
	double *theta = malloc(maxaromatics*numpoints*sizeof(double));
	double *z = malloc(maxaromatics*numpoints*sizeof(double));
	double rmax=0;
	double rdisplacement = 0;
//	char *setyrange = "set yrange [0:1]\n";
//	fprintf(aromaticpipe,setyrange);
	for(i=0;i<numpoints;i++)//Run aromatic calculations on a point which has fixedpoints = 1, but no fixed neighbors
	{
		if(*(fixedpoints+i) == 1)
		{
		nearestneighbors(i,numnn,nn);
		aromaticcenter = 1;
		for(j=0;j<*numnn;j++)
		if(*(fixedpoints+*(nn+j)) == 1)
		aromaticcenter = 0;
		
		if(aromaticcenter)
		{
			printf("Found that index %d is an aromatic center\n",i);
		for(j = 0;j<numroundnesspoints;j++)
		{
			*(radii+j) = j*radiusstep + thickness/2;
			*(round+j) = roundness(i,j*radiusstep,j*radiusstep+thickness);
			fprintf(outfile,"%g %g\n",*(radii+j),*(round+j));
		}
		//printf("About to plot\n");
		//fflush(outfile);
		//printf("A\n");
		fprintf(outfile,"\n");
		//printf("B\n");
		//outfile = fopen("roundness.dat","w");
		//printf("Done Plotting\n");
		
		
		
		//Aromatic Center Type 2
		//printf("a\n");
//		roundness3Dpipe = openpipe();
		//		printf("b\n");
		angularroundness(i,r+numaromatics*numpoints,theta+numaromatics*numpoints,z+numaromatics*numpoints);
		//printf("c\n");
		rmax = 50;
		for(j=0;j<numpoints;j++)
		if((r[j+numaromatics*numpoints] < rmax) && (*(fixedpoints+j) == 1) && (j!=i))
		rmax = r[j+numaromatics*numpoints];
		
		for(j=0;j<numpoints;j++)
		{
		if(r[j+numaromatics*numpoints] > rmax)
		z[j+numaromatics*numpoints] = NAN;
		r[j+numaromatics*numpoints] += rdisplacement;
		}
		rdisplacement += rmax + 5;
		
		printf("rmax is %g\n",rmax);
		//opengnuplotpipepicture(roundness3Dpipe,"Aromaticity");
		//printf("d\n");
		
		//printf("e\n");
		//usleep(1000);
		//fflush(roundness3Dpipe);
		//printf("f\n");
		numaromatics++;
		}	
		}
	}
	fclose(outfile);
	plotline("roundness.dat",aromaticpipe,0);


	return;
	//printsplotfile(roundness3Dfilename,r,theta,z,numaromatics,rows*columns);
	outfile = fopen(roundness3Dfilename,"w");
	
	//double *rnew = malloc(maxaromatics*numpoints*sizeof(double));
	//double *thetanew = malloc(maxaromatics*numpoints*sizeof(double));
	//double *znew = malloc(maxaromatics*numpoints*sizeof(double));
	
	bubblesorttwodatas(theta,r,z,numaromatics*numpoints,sizeof(double),sizeof(double));
	
	bubblesorttwodatas(r,theta,z,numaromatics*numpoints,sizeof(double),sizeof(double));
	
	int ilast = 0;
	for(i=0;i<numpoints*numaromatics;i++)
	{
	//if((i % numpoints) == 0 && i != 0)
	//fprintf(outfile,"\n\n");
	
	double epsilon = 0.0001;
	if(!isnan(z[i]))
	{
	if(fabs(theta[ilast] -theta[i]) > M_PI || fabs(r[ilast] - r[i]) > 2)
	fprintf(outfile,"\n\n");
	
	fprintf(outfile,"%g %g %g\n",r[i],theta[i],z[i]);
	ilast = i;
	}
	}
	bubblesorttwodatas(theta,r,z,numaromatics*numpoints,sizeof(double),sizeof(double));
	for(i=0;i<numpoints*numaromatics;i++)
	{
	//if((i % numpoints) == 0 && i != 0)
	//fprintf(outfile,"\n");
	
	double epsilon = 0.0001;
	if(!isnan(z[i]))
	{
	if(fabs(theta[ilast] -theta[i]) > M_PI || fabs(r[ilast] - r[i]) > 2)
	fprintf(outfile,"\n\n");
	
	fprintf(outfile,"%g %g %g\n",r[i],theta[i],z[i]);
	ilast = i;
	}
	}
	
	fclose(outfile);
	fprintf(roundness3Dpipe, "splot '%s' w lines\n",roundness3Dfilename);
	//splot3d(roundness3Dpipe,roundness3Dfilename);
	fflush(roundness3Dpipe);
	
	free(radii);
	free(round);
	free(r);
	free(theta);
	free(z);
}

void clear()
{
	int index;
	int fixset = 1;
	if(fixspin != NULL)
	fixset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(fixspin));

	
	for(index = 0;index<numpoints;index++)
	{
		*(potential+index) = 0;
		*(fixedpoints+index) = 0;
		*(charge+index) = 0;
		}
		fixedges(fixedpoints,fixset);
	
}

void corral()
{
	
	int index;
	int fixset = 1;
	if(fixspin != NULL)
	fixset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(fixspin));

	
	for(index = 0;index<numpoints;index++)
	{
		*(potential+index) = 0;
		*(fixedpoints+index) = 0;
		*(charge+index) = 0;
		}
		
	double xmajor,ymajor, xcent,ycent;
	
	xcent = *(xgrid+numpoints/2);
	ycent = *(ygrid+numpoints/2);		
	double x,y;
	if(gtk_combo_box_get_active (tiletype) == 0)
	{
	xmajor = columns/2;
	ymajor = rows/2;
	//xcent = columns/2;
	//ycent = rows/2;
	}
	if(gtk_combo_box_get_active (tiletype) == 1)
	{
		xmajor = columns/2-rows/12;// - rows/sqrt(3)/2;
		ymajor = fmin(rows*sqrt(3)/4-1,columns*sqrt(3)/4-1);
		
		x = xcent;
		y = ycent;
		xcent = x*sqrt(3)/2 + y/2;
		ycent = y*sqrt(3)/2 - x/2;
	}
	printf("center %g x %g, axis %g x %g\n",xcent,ycent,xmajor,ymajor);
	for(index = 0;index<numpoints;index++)
	{
		if(gtk_combo_box_get_active (tiletype) == 0)
		{
		x = *(xgrid+index);
		y = *(ygrid+index);
		}
		if(gtk_combo_box_get_active (tiletype) == 1)
		{
		x = *(xgrid+index)*sqrt(3)/2 + *(ygrid+index)/2;
		y = *(ygrid+index)*sqrt(3)/2- *(xgrid+index)/2;
		}
		printf("grid %g x %g primes %g x %g\n",*(xgrid+index),*(ygrid+index),x,y);
		if(!((pow((x-xcent)/xmajor,2) + pow((y-ycent)/ymajor,2)) < 1))
		*(fixedpoints+index) = fixset;
	}
	
	double e = sqrt(pow(xmajor,2) - pow(ymajor,2));
	
	int xfocus = columns/2-e,yfocus=rows/2;
	*(fixedpoints+xfocus+yfocus*columns) = fixset;
	xfocus = columns/2+e;
	*(fixedpoints+xfocus+yfocus*columns) = fixset;
	printf("focus %d x %d eccentricity %g\n",xfocus,yfocus,e);
	
}


void makegrid()
{
	
	int newrows = rows,newcolumns = columns;
	
	if(rowsspin != NULL)
	newrows = gtk_spin_button_get_value(GTK_SPIN_BUTTON(rowsspin));
	if(columnsspin != NULL)
	newcolumns = gtk_spin_button_get_value(GTK_SPIN_BUTTON(columnsspin));
	
	
	int type;
	if(tiletype != NULL)
	type = gtk_combo_box_get_active (tiletype);
	else
	type = 0;
	if(type == 0)
	tiles = square();
	else
	tiles = hexagonaltiling();

	numpoints = newrows * newcolumns;

	double *newpotential = malloc(numpoints*sizeof(double));
	double *newcharge = malloc(numpoints*sizeof(double));
	int *newfixedpoints = malloc(numpoints*sizeof(int));
	
	int index;

	if(potential != NULL)
	{
		for(index = 0;index<numpoints;index++)
		{
			*(newpotential+index) = 0;
			*(newcharge+index) = 0;
			*(newfixedpoints+index) = 0;			
		}
		
		int row,column;
		for(row = 0;(row<newrows) && (row<rows);row++)
		for(column=0;(column<newcolumns) && (column<columns);column++)
		{
//			printf("Row %d column %d\n",row,column);
			*(newpotential+row*newcolumns+column) = *(potential+row*columns+column);
			*(newcharge+row*newcolumns+column) = *(charge+row*columns+column);
			*(newfixedpoints+row*newcolumns+column) = *(fixedpoints+row*columns+column);			
		}
		
		free(potential);
		free(charge);
		free(fixedpoints);
		free(displayedwavefunctions);
		free(eigenvectors);
		free(eigenvalues);
		free(xgrid);
		free(ygrid);
		free(chargedensity);
		free(pairdensity);
		free(chargepotential);
		
		potential = newpotential;
		fixedpoints = newfixedpoints;
		charge = newcharge;
	}
	else
	{
		potential = newpotential;
		fixedpoints = newfixedpoints;
		charge = newcharge;
		clear();	
	}
	displayedwavefunctions = malloc(numpoints*sizeof(int));
	eigenvectors = malloc(numpoints*numpoints*sizeof(double));
	eigenvalues = malloc(numpoints*sizeof(double));
	xgrid = malloc(numpoints*sizeof(double));
	ygrid = malloc(numpoints*sizeof(double));
	rows = newrows;
	columns = newcolumns;
	chargedensity = malloc(numpoints*sizeof(double));
	pairdensity = malloc(numpoints*numpoints*sizeof(double));
	chargepotential = malloc(numpoints*sizeof(double));
	
	getgrid(tiles,xgrid,ygrid);

	
	for(index = 0;index<numpoints;index++)
	{
	*(displayedwavefunctions+index) = 0;
	*(chargedensity + index) = 0;
	*(chargepotential+index) = 0;
	}
		
	
}

void move(int displacement)
{
	
	int fixset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(fixspin));
	double potentialset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(potentialspin));
	double chargeset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(chargespin));
	
	/*
	if(fabs(potentialset) < 0.01)
	potentialset = 0;
	if(fabs(chargeset) < 0.01)
	chargeset = 0;
	*/
	double *newpotential = malloc(numpoints*sizeof(double));
	double *newcharge = malloc(numpoints*sizeof(double));
	int *newfixedpoints = malloc(numpoints*sizeof(int));
	
		
		int row,column, oldindex, newindex,edge;
		for(row = 0;row<rows;row++)
		for(column=0;column<columns;column++)
		{
//			printf("Row %d column %d\n",row,column);
			int newindex = row*columns+column;
			
			int oldindex = newindex-displacement;
			
			
			if(fabs(displacement) == 1)
			edge = (column -displacement <0) || (column-displacement >= columns);
			else
			edge = (oldindex < 0) || (oldindex >= rows*columns);
			
			
			if(edge)
			{
			*(newpotential+newindex) = potentialset;
			*(newcharge+newindex) = chargeset;
			*(newfixedpoints+newindex) = fixset;
			}
			else
			{
			*(newpotential+newindex) = *(potential+oldindex);
			*(newcharge+newindex) = *(charge+oldindex);
			*(newfixedpoints+newindex) = *(fixedpoints+oldindex);			
			}
		}
	free(potential);
	free(charge);
	free(fixedpoints);
	potential = newpotential;
	charge = newcharge;
	fixedpoints = newfixedpoints;
}

gboolean movecallback(int choice)
{
	int displacement = 0;
	switch(choice)
	{
		case 0://left
		displacement = -1;
		break;
		case 1://right
		displacement = 1;
		break;
		case 2://up
		displacement = -columns;
		break;
		case 3://down
		displacement = columns;
		break;
		default:
		printf("tried to perform an invalid move, choice is %d\n",choice);
		break;
	}
	move(displacement);
	draw();
}


void plotpotentialandfixed()
{
	char *potentialfilename = "potential.dat";
	char *fixedfilename = "fixed.dat";
	
	int index;
	double fixeddouble[numpoints];
	int nonzeropotential = 0;
	for(index = 0;index<numpoints;index++)
	{
	*(fixeddouble+index) = *(fixedpoints+index);
	if(*(potential+index) != 0)
	nonzeropotential = 1;
	}
	if(nonzeropotential)
	runplot(potentialfilename,potential,potentialpipe);
	runplot(fixedfilename,fixeddouble,fixedpipe);
}

void plotchargepotential()
{
	char *chargepotentialfilename = "chargepotential.dat";
	
	runplot(chargepotentialfilename,chargepotential,chargepotentialpipe);
}

void plotwavefunction()
{
	int square = gtk_toggle_button_get_active(squarewavefunction);
	double wavefunction[numpoints];
	int pointindex;
	for(pointindex = 0;pointindex<numpoints;pointindex++)
	wavefunction[pointindex] = 0;
	
	int vectorindex;
	int numused;
	for(vectorindex = 0,numused = 0;vectorindex<numpoints;vectorindex++)
	{
		if(*(eigenvalues+vectorindex) == 0)
		continue;
		
		if(*(displayedwavefunctions + numused))
		{
		//printf("Including wavefunction %d\n",vectorindex);
		for(pointindex = 0;pointindex<numpoints;pointindex++)
		{
			if(! (*(fixedpoints+pointindex) == 1))
			if(square)
			*(wavefunction+pointindex) += pow(*(eigenvectors + numpoints*vectorindex+pointindex),2);
			else
			*(wavefunction+pointindex) += *(eigenvectors + numpoints*vectorindex+pointindex);
		}
		}
		numused++;
	}
	runplot("wavefunction.dat",wavefunction,wavefunctionpipe);
}

void click(GtkWidget *widget, GdkEventMotion *event )
{
//	if(event-type == GDK_BUTTON_RELEASE_MASK)
//	return;
	GtkAllocation *allocation = (GtkAllocation *)malloc(sizeof(GtkAllocation));		//Create the drawable pixmap
	gtk_widget_get_allocation(drawingarea, allocation);
	guint width = allocation->width;
	guint height = allocation->height;
	
  int x, y;
  GdkModifierType state;



  if (event->is_hint)
    gdk_window_get_pointer (event->window, &x, &y, &state);
  else
    {
      x = event->x;
      y = event->y;
      state = event->state;

	}
	int left;
	if(state & GDK_BUTTON1_MASK)
	left = 1;
	else
	left = 0;
		
//		printf("click with state %d\n",state);

	int xmax = fmax(fmax(columns * *(tiles->latticevectors),rows * *(tiles->latticevectors+2)),columns * *(tiles->latticevectors)+rows * *(tiles->latticevectors+2));
	int ymax = fmax(fmax(columns * *(tiles->latticevectors+1),rows * *(tiles->latticevectors+3)),columns * *(tiles->latticevectors+1)+rows * *(tiles->latticevectors+3));
	int size = fmin(width/xmax,height/ymax);

	int min=-1;
	double mindist=10000,dist;
	int index;
	for(index = 0;index<numpoints;index++)
	{
		dist = pow(x-size * *(xgrid+index)-size/2.0,2)+pow(y-size * *(ygrid+index)-size/2.0,2);
		if(dist<mindist)
		{
			min = index;
			mindist = dist;
		}
	}
	//int row, column;
	//int size = fmin(width/columns,height/rows);
	//row = y/size;
	//column = x/size;

	int fixset = 0;
	double potentialset = 0;
	double chargeset =0;

	if(left)
	{
	fixset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(fixspin));
	potentialset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(potentialspin));
	chargeset = gtk_spin_button_get_value(GTK_SPIN_BUTTON(chargespin));
	}
//	*(potential + row*columns+column) = potentialset;
//	*(fixedpoints + row*columns+column) = fixset;
	*(potential + min) = potentialset;
	*(charge + min) = chargeset;
	*(fixedpoints + min) = fixset;
	
	gtk_widget_queue_draw(drawingarea);
	
	if(fixset == 4 && condensedeigenvalues != NULL)
	dos();
//	plotpotentialandfixed();
}
static gboolean
motion_notify_event( GtkWidget *widget, GdkEventMotion *event )
{

  int x, y;
  GdkModifierType state;

  if (event->is_hint)
    gdk_window_get_pointer (event->window, &x, &y, &state);
  else
    {
      x = event->x;
      y = event->y;
      state = event->state;
    }
	
	//printf("state is %d button1 %d button2 %d\n",state,GDK_BUTTON1_MASK,GDK_BUTTON2_MASK);
	if(state & GDK_BUTTON1_MASK || state & GDK_BUTTON2_MASK || state & GDK_BUTTON3_MASK)
	//printf("Click and drag\n");
	//if(state == 272)
	click(widget,event);
	
  return FALSE;
}

void main()
{
	 gtk_init(NULL,NULL);
	 makegrid();
	 openpipes();
	
//	printf("Here the number %g %g %g %g\n",*(tiles->latticevectors),*(tiles->latticevectors+1),*(tiles->latticevectors+2),*(tiles->latticevectors+3));

	GtkWidget *loadbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(loadbutton),gtk_image_new_from_stock(GTK_STOCK_OPEN,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(loadbutton, "clicked", G_CALLBACK (load),NULL);

	GtkWidget *savebutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(savebutton),gtk_image_new_from_stock(GTK_STOCK_SAVE,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(savebutton, "clicked", G_CALLBACK (save),NULL);

	 

	drawingarea = gtk_drawing_area_new ();
	gtk_widget_set_size_request(drawingarea,400, 400);
	gtk_widget_add_events (drawingarea, GDK_POINTER_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK |GDK_BUTTON_PRESS_MASK);
	g_signal_connect (drawingarea, "motion_notify_event",motion_notify_event, NULL);
	g_signal_connect (drawingarea, "button_press_event",click, NULL);
	g_signal_connect(drawingarea, "draw", G_CALLBACK(draw), NULL);

	tiletype = gtk_combo_box_text_new ();
	gtk_combo_box_text_append (tiletype,NULL,"Square");
	gtk_combo_box_text_append (tiletype,NULL,"Hexagonal");
	gtk_combo_box_set_active(tiletype,0);
	g_signal_connect_swapped(tiletype, "changed", G_CALLBACK (makegrid),NULL);
	g_signal_connect_swapped(tiletype, "changed", G_CALLBACK (draw),NULL);
	gtk_widget_set_tooltip_text(tiletype,"Select the tiling algorithm");


	GtkWidget *rowsbox = gtk_hbox_new(FALSE,0);
	GtkWidget *rowstext = gtk_label_new("Rows"); 
	rowsspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (10, 1, 100, 1, 10, 0.0)),1,0);
    gtk_box_pack_start (GTK_BOX(rowsbox),rowstext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(rowsbox),rowsspin, TRUE, TRUE, 10);                                  
	g_signal_connect_swapped(rowsspin, "value-changed", G_CALLBACK (makegrid),NULL);
	g_signal_connect_swapped(rowsspin, "value-changed", G_CALLBACK (draw),NULL);
	gtk_widget_set_tooltip_text(rowsspin,"Number of rows to tile");


	GtkWidget *columnsbox = gtk_hbox_new(FALSE,0);
	GtkWidget *columnstext = gtk_label_new("Columns"); 
	columnsspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (10, 1, 100, 1, 10, 0.0)),1,0);
    gtk_box_pack_start (GTK_BOX(columnsbox),columnstext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(columnsbox),columnsspin, TRUE, TRUE, 10);                                  
	g_signal_connect_swapped(columnsspin, "value-changed", G_CALLBACK (makegrid),NULL);
	g_signal_connect_swapped(columnsspin, "value-changed", G_CALLBACK (draw),NULL);
	gtk_widget_set_tooltip_text(columnsspin,"Number of columns to tile");

	GtkWidget *leftbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(leftbutton),gtk_image_new_from_stock(GTK_STOCK_GO_BACK,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(leftbutton, "clicked", G_CALLBACK (movecallback),0);
	
	GtkWidget *rightbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(rightbutton),gtk_image_new_from_stock(GTK_STOCK_GO_FORWARD,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(rightbutton, "clicked", G_CALLBACK (movecallback),1);
	
	GtkWidget *upbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(upbutton),gtk_image_new_from_stock(GTK_STOCK_GO_UP,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(upbutton, "clicked", G_CALLBACK (movecallback),2);
	
	GtkWidget *downbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(downbutton),gtk_image_new_from_stock(GTK_STOCK_GO_DOWN,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(downbutton, "clicked", G_CALLBACK (movecallback),3);
	
	GtkWidget *arrowsbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(arrowsbox),leftbutton, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(arrowsbox),rightbutton, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(arrowsbox),upbutton, FALSE, TRUE, 10);
	gtk_box_pack_start (GTK_BOX(arrowsbox),downbutton, FALSE, TRUE, 10);

	
	
	GtkWidget *fixbox = gtk_hbox_new(FALSE,0);
	GtkWidget *fixtext = gtk_label_new("Set Mode:"); 
	fixspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (1, 0, 9, 1, 1, 0.0)),1,0);
    gtk_box_pack_start (GTK_BOX(fixbox),fixtext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(fixbox),fixspin, TRUE, TRUE, 10);                                  
	gtk_widget_set_tooltip_text(fixspin,"Boundary Conditions for tiles\n 0 open (white), 1 blocked (red), 2 periodic (green), 3 closed (blue)"
                              ", 4 pointspec, 5-9 open (colored)\n 2 and 3 only affect edges");


	GtkWidget *potentialbox = gtk_hbox_new(FALSE,0);
	GtkWidget *potentialtext = gtk_label_new("Set Potential:"); 
	potentialspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (0, -100, 100, 0.01, 0.01, 0.01)),0.01,2);
    gtk_box_pack_start (GTK_BOX(potentialbox),potentialtext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(potentialbox),potentialspin, TRUE, TRUE, 10);
    gtk_widget_set_tooltip_text(potentialspin,"Applies potential to electrons within this tile.\n Only used in Fock calculations");

                                      

	GtkWidget *chargebox = gtk_hbox_new(FALSE,0);
	GtkWidget *chargetext = gtk_label_new("Set Charge:"); 
	chargespin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (0, -100, 100, 0.01, 0.01, 0.01)),0.01,2);
    gtk_box_pack_start (GTK_BOX(chargebox),chargetext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(chargebox),chargespin, TRUE, TRUE, 10);                                  
    gtk_widget_set_tooltip_text(chargespin,"Generates a potential based upon charge in this tile.\n Only used in Fock calculations");



	GtkWidget *clearbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(clearbutton),gtk_image_new_from_stock(GTK_STOCK_REFRESH,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(clearbutton, "clicked", G_CALLBACK (clear),NULL);
	g_signal_connect_swapped(clearbutton, "clicked", G_CALLBACK (draw),NULL);
    gtk_widget_set_tooltip_text(clearbutton,"Reset the grid");
    
    
    GtkWidget *corralbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(corralbutton),gtk_image_new_from_stock(GTK_STOCK_EXECUTE,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(corralbutton, "clicked", G_CALLBACK (corral),NULL);
	g_signal_connect_swapped(corralbutton, "clicked", G_CALLBACK (draw),NULL);
    gtk_widget_set_tooltip_text(corralbutton,"Create the best estimation of a quantum corral(ellipse)\n Use an odd number orf rows and columns");


	eigenvaluelist = gtk_list_store_new(2,G_TYPE_DOUBLE,G_TYPE_BOOLEAN);
	treeview = gtk_tree_view_new_with_model(eigenvaluelist);
    gtk_widget_set_tooltip_text(treeview,"Displays the eigenvalues of the system.\n The occupancy can be toggled to select the wavefunction for plotting and SCF");

	GtkCellRenderer *renderer = gtk_cell_renderer_text_new ();

	gtk_tree_view_insert_column_with_attributes
                               (treeview,
                                -1,
                                "Eigenvalues",
                                renderer,"text",0/*,"mode", GTK_CELL_RENDERER_MODE_EDITABLE*/,NULL);
	renderer = gtk_cell_renderer_toggle_new ();
//	renderer = gtk_cell_renderer_spin_new ();
//	renderer = gtk_cell_renderer_text_new ();
//	gtk_cell_renderer_toggle_set_radio(renderer,TRUE);


	GtkAdjustment *occupancyadjustment = gtk_adjustment_new (0, 0, 2, 1, 1, 1);
	gtk_tree_view_insert_column_with_attributes
                               (treeview,
                                -1,
                                "Occupancy",
                                renderer,"active",1/*,"editable",TRUE,"text",1*//*,"adjustment",occupancyadjustment*/,NULL);
	
//    g_signal_connect_swapped(renderer,"toggled",toggledisplay,10);
	
	g_signal_connect (renderer, "toggled",G_CALLBACK (fixed_toggled),gtk_tree_view_get_model(treeview));
//	g_signal_connect (renderer, "edited",G_CALLBACK (fixed_toggled),gtk_tree_view_get_model(treeview));

	GtkWidget *roundnessthicknessbox = gtk_hbox_new(FALSE,0);
	GtkWidget *roundnessthicknesstext = gtk_label_new("Thickness:"); 
	roundnessthicknessspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (1, 0.01, 10, 0.001, 0.001, 0.001)),0.001,3);
    gtk_box_pack_start (GTK_BOX(roundnessthicknessbox),roundnessthicknesstext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(roundnessthicknessbox),roundnessthicknessspin, TRUE, TRUE, 10); 
    gtk_widget_set_tooltip_text(roundnessthicknessspin,"The thickness of the annular region which is integrated over to calculated roundness");

	GtkWidget *maximumradiusbox = gtk_hbox_new(FALSE,0);
	GtkWidget *maximumradiustext = gtk_label_new("Maximum Radius:"); 
	maximumradiusspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (5, 1, 100, 0.001, 0.001, 0.001)),0.001,3);
    gtk_box_pack_start (GTK_BOX(maximumradiusbox),maximumradiustext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(maximumradiusbox),maximumradiusspin, TRUE, TRUE, 10); 
    gtk_widget_set_tooltip_text(maximumradiusspin,"The large-radius cutoff of roundness calculations.");

	
	GtkWidget *chargepotentialbox = gtk_hbox_new(FALSE,0);
	GtkWidget *chargepotentialtext = gtk_label_new("Charge Potential:"); 
	chargepotentialspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (0.01, 0.001, 10, 0.001, 0.001, 0.001)),0.001,3);
    gtk_box_pack_start (GTK_BOX(chargepotentialbox),chargepotentialtext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(chargepotentialbox),chargepotentialspin, TRUE, TRUE, 10); 
    gtk_widget_set_tooltip_text(chargepotentialspin,"The coefficient of coulomb potential in Fock calculations");

	GtkWidget *playbutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(playbutton),gtk_image_new_from_stock(GTK_STOCK_MEDIA_PLAY,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(playbutton, "clicked", G_CALLBACK (calculatenocharge),NULL);
	gtk_widget_set_tooltip_text(playbutton,"Hartree level calculation");


	GtkWidget *iteratebutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(iteratebutton),gtk_image_new_from_stock(GTK_STOCK_MEDIA_NEXT,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(iteratebutton, "clicked", G_CALLBACK (iterate),NULL);
    gtk_widget_set_tooltip_text(iteratebutton,"A single iteration of Fock calculation.");


	GtkWidget *convergebutton = gtk_button_new();
	gtk_button_set_image(GTK_BUTTON(convergebutton),gtk_image_new_from_stock(GTK_STOCK_REFRESH,GTK_ICON_SIZE_DND));
	g_signal_connect_swapped(convergebutton, "clicked", G_CALLBACK (converge),NULL);
	//g_signal_connect_swapped(convergebutton, "clicked", G_CALLBACK (calculateroundness),NULL);
    gtk_widget_set_tooltip_text(convergebutton,"Self-consistent iterative Fock");


	GtkWidget *playbox = gtk_hbox_new(FALSE,0);

	gtk_box_pack_start(GTK_BOX(playbox),playbutton, FALSE, TRUE, 10);
	gtk_box_pack_start(GTK_BOX(playbox),iteratebutton, FALSE, TRUE, 10);
	gtk_box_pack_start(GTK_BOX(playbox),convergebutton, FALSE, TRUE, 10);

	
	GtkWidget *scrolledwindow=  gtk_scrolled_window_new(NULL,NULL);
	gtk_scrolled_window_set_policy(scrolledwindow,GTK_POLICY_NEVER, GTK_POLICY_AUTOMATIC);
	gtk_container_add (GTK_CONTAINER (scrolledwindow),treeview);
	
	squarewavefunction = gtk_check_button_new_with_label("Square Wavefunction");
	
	GtkWidget *interbox = gtk_hbox_new(FALSE,0);
	GtkWidget *intertext = gtk_label_new("Interpolate"); 
	interspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (8, 1, 100, 1, 1, 1)),1,0);
    gtk_box_pack_start (GTK_BOX(interbox),intertext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(interbox),interspin, TRUE, TRUE, 10);                                  
	g_signal_connect_swapped(interspin, "value-changed", G_CALLBACK (setinterpolate),NULL);
	
	GtkWidget *broadeningbox = gtk_hbox_new(FALSE,0);
	GtkWidget *broadeningtext = gtk_label_new("Broadening"); 
	broadeningspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (0.1, 0.0001, 100, 0.001, 0.01, 0.01)),3,3);
    gtk_box_pack_start (GTK_BOX(broadeningbox),broadeningtext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(broadeningbox),broadeningspin, TRUE, TRUE, 10);                                  
	g_signal_connect_swapped(broadeningspin, "value-changed", G_CALLBACK (dos),NULL);
	
	
	GtkWidget *upperlimitbox = gtk_hbox_new(FALSE,0);
	GtkWidget *upperlimittext = gtk_label_new("UpperLimit"); 
	upperlimitspin = gtk_spin_button_new(GTK_ADJUSTMENT(gtk_adjustment_new (100, -100, 100, 0.01, 0.01, 0.1)),2,2);
    gtk_box_pack_start (GTK_BOX(upperlimitbox),upperlimittext, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(upperlimitbox),upperlimitspin, TRUE, TRUE, 10);                                  
	g_signal_connect_swapped(upperlimitspin, "value-changed", G_CALLBACK (dos),NULL);
	
	

	GtkWidget *filebox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(filebox),loadbutton, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(filebox),savebutton, FALSE, TRUE, 10);

	
	GtkWidget *topbox = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(topbox),tiletype, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(topbox),rowsbox, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(topbox),columnsbox, FALSE, TRUE, 10);
	
	GtkWidget *bottombox = gtk_hbox_new(FALSE,0);
    gtk_box_pack_start (GTK_BOX(bottombox),fixbox, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(bottombox),potentialbox, FALSE, TRUE, 10);
    
    GtkWidget *bottom2box = gtk_hbox_new(FALSE,0);
    gtk_box_pack_start (GTK_BOX(bottom2box),chargebox, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(bottom2box),clearbutton, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(bottom2box),corralbutton, FALSE, TRUE, 10);

	
	GtkWidget *leftbox = gtk_vbox_new(FALSE,0);
    gtk_box_pack_start (GTK_BOX(leftbox),filebox, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(leftbox),topbox, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(leftbox),arrowsbox, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(leftbox),drawingarea, TRUE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(leftbox),bottombox, FALSE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(leftbox),bottom2box, FALSE, TRUE, 10);

	GtkWidget *rightbox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start (GTK_BOX(rightbox),squarewavefunction, FALSE, TRUE, 0);
	gtk_box_pack_start (GTK_BOX(rightbox),interbox, FALSE, TRUE, 0);
	gtk_box_pack_start (GTK_BOX(rightbox),broadeningbox, FALSE, TRUE, 0);
	gtk_box_pack_start (GTK_BOX(rightbox),upperlimitbox, FALSE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(rightbox),scrolledwindow, TRUE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(rightbox),roundnessthicknessbox, FALSE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(rightbox),maximumradiusbox, FALSE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(rightbox),chargepotentialbox, FALSE, TRUE, 0);
    gtk_box_pack_start (GTK_BOX(rightbox),playbox, FALSE, TRUE, 0);

	GtkWidget *mainbox = gtk_hbox_new(FALSE,0);
    gtk_box_pack_start (GTK_BOX(mainbox),leftbox, TRUE, TRUE, 10);
    gtk_box_pack_start (GTK_BOX(mainbox),rightbox, FALSE, FALSE, 10);
	

	TopLevel = gtk_window_new (GTK_WINDOW_TOPLEVEL); 
	gtk_container_add (GTK_CONTAINER (TopLevel), mainbox);
	gtk_widget_show_all(TopLevel);
	
	 //		printf("Here the number %d %d\n",*(tiles->compositevectors),*(tiles->compositevectors+1));

	gtk_main();
}	
