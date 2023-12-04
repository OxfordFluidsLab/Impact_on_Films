/* 3phasedroponfilmonpool.c

	Ben Fudge, 04/12/2023
	benjamin.fudge@maths.ox.ac.uk

	3 Phase drop on film on top of pool pool impact

	Compile with 

		export OMP_NUM_THREADS=X

	where X is the number of threads you want to use

		qcc -O2 -w -fopenmp 3phasedroponfilmonpool.c -o 3phasedroponfilmonpool -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm

	And executed with

		./3phasedroponfilmonpool Re We Fr Real_Pool_Density Real_Drop_Viscosity Real_Drop_Air_ST Real_Pool_Air_ST Real_Drop_Pool_ST Film_Thickness MaxLevel
	
	where
		-Real_Pool_Density is the pool density in kg/m^3
		-Real_Pool_Viscosity is the pool kinematic viscosity in cSt
		-Real_Drop_Air_ST is the surface tension between the droplet and air in mNm (i.e. water is 72) and the same for the other two pairs
        -Film_Thickness is the thickness of the liquid film normalised to a droplet diameter of 1
		-MaxLevel is the maximum resolution level   

	e.g. 
	
		./3phasedroponfilmonpool 1000 100 10 934 2 15 20 5 0.5 10

	The parameters have to be in that order, if an incorrect number of parameters are provided them the simulation will abort.
*/

#define SLIP

#define harmu

#ifdef harmu
	#define mu(f1, f2, f3) (1./(clamp(f1,0,1)*1./mu1 + clamp(f2,0,1)*1./mu2 + clamp(f3,0,1)*1./mu3))
#endif

#include "axi.h"
#include "navier-stokes/centered.h"

#include "three-phaseCHIZARI.h"
#include "tension.h"

#include "view.h"
#include "tag.h"

#include "interface_area_axi.h"

// Needed to sort out directories for results
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// Grid refinement
#define LEVEL 6
int MAXLEVEL;

// Non dimensional droplet size, velocity
// and density as well as initial distance
// from the surface
#define D0 1.0
double V0;
#define RHO1 1.0
double S0;

#define SIZE 5.0

// Fixed Air Properties
#define REAL_AIR_DENSITY 1.2754 // Air Density in kg/m^3 (fixed)
#define REAL_AIR_VISC 18.1e-6 // Air Dynamic Viscosity in Pa.s (fixed)

// Fixed Droplet Properties for FC-770
#define REAL_DROP_DENSITY 1793 // Real Droplet Density in kg/m^3 (fixed)
#define REAL_DROP_VISCOSITY 1.4e-3 // Real Droplet Dynamic Viscosity in Pa.s (fixed)

// Fixed Droplet Properties for Water
//#define REAL_DROP_DENSITY 1000 // Real Droplet Density in kg/m^3 (fixed)
//#define REAL_DROP_VISCOSITY 1.0e-3 // Real Droplet Dynamic Viscosity in Pa.s (fixed)

// Varying Pool Properties
double REAL_POOL_DENSITY; // Real Pool Density in kg/m^3 (given)
double REAL_POOL_VISCOSITY; // Real Pool Kinematic Viscosity in cSt (convenient for SO) (given)

// Density Ratios
double RHO_AIR2DROP; // Air to Droplet Density Ratio (calculated)
double RHO_POOL2DROP; // Pool to Droplet Density Ratio (calculated)

// Viscosity Ratios
double MU_AIR2DROP; // Air to Droplet Dynamic Viscosity Ratio (calculated) 
double MU_POOL2DROP; // Pool to Droplet Dynamic Viscosity Ratio (calculated)

// Surface Tensions
double REAL_DROP_AIR_ST; // Droplet-Air Surface Tension in mN/m (given)
double REAL_POOL_AIR_ST; // Pool-Air Surface Tension in mN/m (given)
double REAL_DROP_POOL_ST; // Droplet-Pool Surface Tension in mN/m (given)

// Dimensionless Groups
double RE;
double WE;
double FR;

// For Surface Tensions
double surf_drop_air;
double sigmaSF_sigma_FA;
double surf_drop_pool;
double sigmaSA_sigma_FA;
double surf_pool_air;
double SIGMA1;
double SIGMA2;
double SIGMA3;

// For Adaptivity
double f1tol;
double f2tol;
double f3tol;
double ff1tol;
double ff2tol;
double ff3tol;
double vorttol;

// For Filtering
int filtersize;
double filtertol;

int filterflag=0;

// Film Thickness
double filmthickness;

// For saving stuff
char resultsfolder[500];

char restring[100];
char westring[100];
char frstring[100];
char viscratiostring[100];
char filmthicknessstring[100];
char maxlevelstring[100];

char interfacesfolder[500];
char snapshotsfolder[500];
char velsfolder[600];

char videoname[500];
char runstats[500];
char volumeslog[500];
char energylog[500];
char paramslog[500];
char surfacearealog[500];


scalar fupper[];
scalar flower[];

// Solid BCs for when using Gravity

// Bottom of Pool = LEFT of domain
// No slip, no permeability
u.n[left]=dirichlet(0); // Impermeable
u.t[left]=dirichlet(0); // No Slip
p[left]=neumann(0);
pf[left]=neumann(0);

// Side of Pool = TOP of domain
// Currently no slip and no permeability
// May change to slip condition
u.n[top]=dirichlet(0); // Impermeable
#ifdef SLIP
u.t[top]=neumann(0); // Slip 
#else
u.t[top]=dirichlet(0); // No Slip
#endif
p[top]=neumann(0);
pf[top]=neumann(0);

// Above the Pool = RIGHT of domain
// Outflow
u.n[right]=neumann(0);
p[right]=dirichlet(0);
pf[right]=dirichlet(0);

// Default for right edge is symmetry

int main(int argc, char * argv[])
{
	assert(argc==11);

	RE=atof(argv[1]);
	WE=atof(argv[2]);
	FR=atof(argv[3]);
	REAL_POOL_DENSITY=atof(argv[4]);
	REAL_POOL_VISCOSITY=atof(argv[5]);
	REAL_DROP_AIR_ST=atof(argv[6]);
	REAL_POOL_AIR_ST=atof(argv[7]);
	REAL_DROP_POOL_ST=atof(argv[8]);
	filmthickness=atof(argv[9]);
	MAXLEVEL=atoi(argv[10]);

	// Calculate Ratios
	RHO_AIR2DROP=REAL_AIR_DENSITY/REAL_DROP_DENSITY;
	RHO_POOL2DROP=REAL_POOL_DENSITY/REAL_DROP_DENSITY;

	MU_AIR2DROP=REAL_AIR_VISC/REAL_DROP_VISCOSITY;

	double REAL_POOL_VISCOSITY_Pas=REAL_POOL_VISCOSITY*REAL_POOL_DENSITY/1e6;
	MU_POOL2DROP=REAL_POOL_VISCOSITY_Pas/REAL_DROP_VISCOSITY;

	// For saving 
	sprintf(resultsfolder, "./Results_Pool_on_Film/");

	sprintf(restring, "Re_%1.03f_", RE);
	sprintf(westring, "We_%1.03f_", WE);
	sprintf(frstring, "Fr_%2.1f_", FR);
	sprintf(viscratiostring, "Visc_Ratio_%1.03f_", MU_POOL2DROP);
	sprintf(filmthicknessstring, "FILM_THICKNESS_%1.02f_", filmthickness);
	sprintf(maxlevelstring,"MAX_LEVEL_%i", MAXLEVEL);

	strcat(resultsfolder, restring);
	strcat(resultsfolder, westring);
	strcat(resultsfolder, frstring);
	strcat(resultsfolder, viscratiostring);
	strcat(resultsfolder, filmthicknessstring);
	strcat(resultsfolder, maxlevelstring);

	strcpy(interfacesfolder, resultsfolder);
	strcpy(snapshotsfolder, resultsfolder);
	strcpy(velsfolder, resultsfolder);
	strcpy(videoname, resultsfolder);
	strcpy(runstats, resultsfolder);
	strcpy(volumeslog, resultsfolder);
	strcpy(energylog, resultsfolder);
	strcpy(paramslog, resultsfolder);
	strcpy(surfacearealog, resultsfolder);

	strcat(interfacesfolder, "/Interfaces/");
	strcat(snapshotsfolder, "/Snapshots/");
	strcat(velsfolder, "/Velocities/");

	strcat(runstats, "/Simulation_Stats_");
	strcat(runstats, restring);
	strcat(runstats, westring);
	strcat(runstats, frstring);
	strcat(runstats, viscratiostring);
	strcat(runstats, filmthicknessstring);
	strcat(runstats, maxlevelstring);

	strcat(videoname, "/Video_");
	strcat(videoname, restring);
	strcat(videoname, westring);
	strcat(videoname, frstring);
	strcat(videoname, viscratiostring);
	strcat(videoname, filmthicknessstring);
	strcat(videoname, maxlevelstring);
	strcat(videoname, ".mp4");

	strcat(volumeslog, "/Volumes_Log_");
	strcat(volumeslog, restring);
	strcat(volumeslog, westring);
	strcat(volumeslog, frstring);
	strcat(volumeslog, viscratiostring);
	strcat(volumeslog, filmthicknessstring);
	strcat(volumeslog, maxlevelstring);

	strcat(energylog, "/Energy_Log_");
	strcat(energylog, restring);
	strcat(energylog, westring);
	strcat(energylog, frstring);
	strcat(energylog, viscratiostring);
	strcat(energylog, filmthicknessstring);
	strcat(energylog, maxlevelstring);

	strcat(surfacearealog, "/Surface_Area_Log_");
	strcat(surfacearealog, restring);
	strcat(surfacearealog, westring);
	strcat(surfacearealog, frstring);
	strcat(surfacearealog, viscratiostring);
	strcat(surfacearealog, filmthicknessstring);
	strcat(surfacearealog, maxlevelstring);			
		
	strcat(paramslog, "/Simulation_Parameters");

	// Folders for results
	struct stat st1 = {0};

		if (stat("./Results_Pool_on_Film", &st1) == -1) {
			mkdir("./Results_Pool_on_Film", 0700);
		}

	struct stat st2 = {0};
	int answer; 

		if (stat(resultsfolder, &st2) == -1) {
			mkdir(resultsfolder, 0700);
		}

		else {
			printf("Results already exist for this simulation: %s\n", resultsfolder);
			printf("Do you wish to delete the old results and replace with new ones? [1] Yes, [2] No ");
			scanf("%i", &answer);

			if (answer == 1) {
				printf("Deleting old results and running new simulation\n");
			}

			else {
				printf("Terminating Simulation\n");
				return 0;
			}
		}

	struct stat st3 = {0};

		if (stat(interfacesfolder, &st3) == -1) {
			mkdir(interfacesfolder, 0700);
		}

	struct stat st4 = {0};

		if (stat(snapshotsfolder, &st4) == -1) {
			mkdir(snapshotsfolder, 0700);
		}


	struct stat st5 = {0};

		if (stat(velsfolder, &st5) == -1) {
			mkdir(velsfolder, 0700);
		}


	// Delete text files if they already exist
	if ( access( videoname, F_OK ) != -1 ) {
		remove(videoname);
	}

	if ( access( runstats, F_OK ) != -1 ) {
		remove(runstats);
	}

	if ( access( volumeslog, F_OK) != -1) {
		remove(volumeslog);
	}

	if ( access( energylog, F_OK) != -1) {
		remove(energylog);
	}

	if ( access( surfacearealog, F_OK) != -1) {
		remove(surfacearealog);
	}

	// Simulation size and resolution
	size(SIZE*D0);
	origin(-0.5*SIZE,0);
	init_grid(1<<LEVEL);

	// Densities 
	rho1 = RHO1;
	rho2 = RHO1 * RHO_AIR2DROP;
	rho3 = RHO1 * RHO_POOL2DROP;

	// Defines Viscosity to get correct Re
	mu1 = 1.0/(double)RE;
	mu2 = mu1 * MU_AIR2DROP;
	mu3 = mu1 * MU_POOL2DROP;

	// Defines surface tensions to get correct We
	surf_drop_air = 1.0/WE;

	sigmaSF_sigma_FA = REAL_DROP_POOL_ST/REAL_DROP_AIR_ST;
	surf_drop_pool = surf_drop_air*sigmaSF_sigma_FA;

	sigmaSA_sigma_FA = REAL_POOL_AIR_ST/REAL_DROP_AIR_ST;
	surf_pool_air = surf_drop_air*sigmaSA_sigma_FA;

	SIGMA1 = 0.5*(surf_drop_air + surf_drop_pool - surf_pool_air);
	SIGMA2 = 0.5*(surf_drop_air + surf_pool_air - surf_drop_pool);
	SIGMA3 = 0.5*(surf_drop_pool + surf_pool_air - surf_drop_air);

	f1.sigma = SIGMA1;
	f2.sigma = SIGMA2;
	f3.sigma = SIGMA3;

	fupper.sigma=SIGMA3;
	flower.sigma=SIGMA3;

	printf("\nStarting Simulation Re %4.1f, We %3.1f, Fr %2.1f, Visc Ratio %1.02f & Max Level %i\n\n", 
	RE, WE, FR, MU_POOL2DROP, MAXLEVEL);

	run();
}

event defaults (i=0)
{
	interfaces = list_add (NULL, f1);
	interfaces = list_add (interfaces, f2);
	interfaces = list_add (interfaces, f3);
	interfaces = list_add (interfaces, fupper);
	interfaces = list_add (interfaces, flower);
}

double inittotalf1vol;
double inittotalf2vol;
double inittotalf3vol;

double levelwidth; 

event init (t=0)
{
	// Tolerances
	if (MAXLEVEL >= 10)
	{
		vorttol = 5e-1;
		levelwidth=0.005;
	}
	else
	{
		vorttol = 1e-1;
		levelwidth=0.015;
	}

	double offset=SIZE*pow(2,-MAXLEVEL)/3;

	S0=0.6;
	V0=1.0;

	// Have initial grid highly refined around droplet
	// as we know where it is
	refine(fabs(x) < 32.0*levelwidth && level < MAXLEVEL-5);
	refine(fabs(x) < 16.0*levelwidth && level < MAXLEVEL-4);
	refine(fabs(x) < 8.0*levelwidth && level < MAXLEVEL-3);
	refine(fabs(x) < 6.0*levelwidth && level < MAXLEVEL-2);
	refine(fabs(x) < 4.0*levelwidth && level < MAXLEVEL-1);
	refine(fabs(x) < 1.0*levelwidth && level < MAXLEVEL);

	refine(fabs(x-filmthickness) < 32.0*levelwidth && level < MAXLEVEL-5);
	refine(fabs(x-filmthickness) < 16.0*levelwidth && level < MAXLEVEL-4);
	refine(fabs(x-filmthickness) < 8.0*levelwidth && level < MAXLEVEL-3);
	refine(fabs(x-filmthickness) < 6.0*levelwidth && level < MAXLEVEL-2);
	refine(fabs(x-filmthickness) < 4.0*levelwidth && level < MAXLEVEL-1);
	refine(fabs(x-filmthickness) < 1.0*levelwidth && level < MAXLEVEL);

	refine(sq(x-S0+offset-filmthickness) + sq(y) < sq(0.52*D0) && sq(x-S0+offset-filmthickness) + sq(y) > sq(0.48*D0) && level < MAXLEVEL);

	// Initialises droplet with centre S0 from midline
	// Pool for negative x

	fraction(f1, union(-sq(x-S0+offset-filmthickness)-sq(y)+sq(0.5*D0), -x+offset)); // Drop
	fraction(f2, difference(x-offset-filmthickness, -sq(x-S0+offset-filmthickness)-sq(y)+sq(0.5*D0))); // Air
	fraction(f3, difference(-x+offset+filmthickness, -x+offset)); // Pool

	fraction(fupper, -x+offset+filmthickness);
	fraction(flower, -x+offset);

	// Droplet has initial velocity towards pool
	foreach()
	{
		u.x[] = -V0*f1[];
	}

	inittotalf1vol = statsf(f1).sum*2.0*pi;
	inittotalf2vol = statsf(f2).sum*2.0*pi;
	inittotalf3vol = statsf(f3).sum*2.0*pi;

	FILE * vollog = fopen(volumeslog, "a");
	// Iteration, Time, f_Total, f2_Total, f3_Total
	fprintf(vollog, "%i %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", i, t, inittotalf1vol, inittotalf1vol/inittotalf1vol,
	inittotalf2vol, inittotalf2vol/inittotalf2vol, inittotalf3vol, inittotalf3vol/inittotalf3vol);
	fclose(vollog);

	char snapshotname[400];
	char snapshottimestamp[100];
	sprintf(snapshottimestamp, "Initial.gfs");

	strcpy(snapshotname, snapshotsfolder);
	strcat(snapshotname, snapshottimestamp);

	output_gfs(file = snapshotname, translate = true, t=t);

	char interfacename[500];

	char interfacetimestamp[100];
	sprintf(interfacetimestamp, "Interface_Initial");

	strcpy(interfacename, interfacesfolder);
	strcat(interfacename, interfacetimestamp);

	char interface2name[500];
	strcpy(interface2name, interfacename);
	strcat(interface2name,"_2");

	char interface3name[500];
	strcpy(interface3name, interfacename);
	strcat(interface3name,"_3");

	char interfaceuppername[500];
	strcpy(interfaceuppername, interfacename);
	strcat(interfaceuppername,"_upper");

	char interfacelowername[500];
	strcpy(interfacelowername, interfacename);
	strcat(interfacelowername,"_lower");

	strcat(interfacename, "_1");

	FILE * fpp = fopen(interfacename, "w");
	output_facets(f1,fpp);
	fclose(fpp);

	FILE * fpp2 = fopen(interface2name, "w");
	output_facets(f2,fpp2);
	fclose(fpp2);

	FILE * fpp3 = fopen(interface3name, "w");
	output_facets(f3,fpp3);
	fclose(fpp3);

	FILE * fppupper = fopen(interfaceuppername, "w");
	output_facets(fupper,fppupper);
	fclose(fppupper);

	FILE * fpplower = fopen(interfacelowername, "w");
	output_facets(flower,fpplower);
	fclose(fpplower);

	f1tol=2e-2;
	f2tol=2e-2;
	f3tol=2e-2;
	ff1tol=2e-2;
	ff2tol=2e-2;
	ff3tol=2e-2;

	filtersize=4;
	filtertol=1e-3;

	// Prints all simulation parameters to a file
	FILE * parlog = fopen(paramslog, "w");
	fprintf(parlog,"Level %i\nMax Level %i\nDomain Size %1.1f\nDroplet Dimensionless Diameter %1.1f\nDroplet Initial Dimensionless Velocity %1.1f\n",
	LEVEL,MAXLEVEL,SIZE,D0,V0); 
	fprintf(parlog,"Droplet Dimensionless Density %1.1f\nDroplet Dimensionless Distance Above Pool %1.1f\nAir to Droplet Density Ratio %1.6f\n",
	RHO1,S0,RHO_AIR2DROP);
	fprintf(parlog,"Pool to Droplet Density Ratio %1.6f\nAir to Droplet Viscosity Ratio %1.6f\nPool to Droplet Viscosity Ratio %1.6f\n",
	RHO_POOL2DROP,MU_AIR2DROP,MU_POOL2DROP);
	fprintf(parlog,"Real Droplet-Air Surface Tension (mN/m) %1.2f\nReal Droplet-Pool Surface Tension (mN/m) %1.2f\n",REAL_DROP_AIR_ST,REAL_DROP_POOL_ST);
	fprintf(parlog,"Real Pool-Air Surface Tension (mN/m) %1.2f\nRe %4.1f\nWe %3.1f\nFr %1.1f\nf1 (Droplet) Density %1.6f\nf2 (Air) Density %1.6f\n",
	REAL_POOL_AIR_ST,RE,WE,FR,rho1,rho2);
	fprintf(parlog,"f3 (Pool) Density %1.6f\nf1 (Droplet) Viscosity %1.6f\nf2 (Air) Viscosity %1.6f\nf3 (Pool) Viscosity %1.6f\n",rho3,mu1,mu2,mu3);
	fprintf(parlog,"Drop-Air Dimensionless Surface Tension %1.6f\nDrop-Pool/Drop Air Surface Tension Ratio %1.6f\n",surf_drop_air,sigmaSF_sigma_FA);
	fprintf(parlog,"Drop-Pool Dimensionless Surface Tension %1.6f\nAir-Pool/Drop Air Surface Tension Ratio %1.6f\n",surf_drop_pool,sigmaSA_sigma_FA);
	fprintf(parlog,"Air-Pool Dimensionless Surface Tension %1.6f\nf1 (Droplet) Surface Tension %1.6f\nf2 (Air) Surface Tension %1.6f\n",
	surf_pool_air,SIGMA1,SIGMA2);
	fprintf(parlog,"f3 (Pool) Surface Tension %1.6f\n",SIGMA3);
	fprintf(parlog,"Adaptivity Level Widths %1.4f\nf1 (Droplet) Adaptivity Threshold %1.4f\nf2 (Air) Adaptivity Threshold %1.4f\n",
	levelwidth,f1tol,f2tol);
	fprintf(parlog,"f3 (Pool) Adaptivity Threshold %1.4f\nf1 (Droplet) Boundary Adaptivity Threshold %1.4f\nf2 (Air) Boundary Adaptivity Threshold %1.4f\n",
	f3tol,ff1tol,ff2tol);
	fprintf(parlog,"f3 (Pool) Boundary Adaptivity Threshold %1.4f\nVorticity Adaptivity Threshold %1.4f\nFilter Size %i\nFilter Threshold %1.4f\n",
	ff3tol,vorttol,filtersize,filtertol);
	#ifdef FILTERED
	fprintf(parlog, "Filtered Interface\n");
	#else
	fprintf(parlog, "Unfiltered Interface\n");
	#endif
	#ifdef mu
	fprintf(parlog, "Harmonic Mean for Viscosity\n");
	#else
	fprintf(parlog, "Arithmetic Mean for Viscosity\n");
	#endif
	#ifdef SLIP
	fprintf(parlog, "Slip Boundary Condition for Side Edge\n");
	#else
	fprintf(parlog, "No-Slip Boundary Condition for Side Edge\n");
	#endif
	fclose(parlog);

}

event acceleration (i++)
{
	// Defines gravity to get correct Fr
	// Also in this case downwards is in the 
	// negative x direction
	face vector av=a;
	foreach_face(x)
	{
		av.x[] -= 1.0/(FR*FR);
	}
}

scalar energydissipationrate[];
scalar energydissipationratedroplet[];
scalar energydissipationrateair[];
scalar energydissipationratepool[];
scalar energydissipationrate_dV[];
scalar dudx[];
scalar dvdy[];
scalar dvdx[];
scalar dudy[];

scalar kineticenergypervol[];
scalar kineticenergypervoldroplet[];
scalar kineticenergypervolair[];
scalar kineticenergypervolpool[];

scalar gravenergypervol[];
scalar gravenergypervoldroplet[];
scalar gravenergypervolair[];
scalar gravenergypervolpool[];

scalar viewingfield[];
scalar mymu[];
scalar myrho[];

scalar vel2[];
scalar myvort[];
scalar mylevel[];

double totaldissipatedenergy = 0;
double totaldissipateddropletenergy = 0;
double totaldissipatedpoolenergy = 0;
double totaldissipatedairenergy = 0;

double currentdissipatedenergy = 0;
double currentdissipateddropletenergy = 0;
double currentdissipatedpoolenergy = 0;
double currentdissipatedairenergy = 0;

double totalkineticenergy = 0;
double totalkineticdropletenergy = 0;
double totalkineticpoolenergy = 0;
double totalkineticairenergy = 0;

double totalgravenergy = 0;
double totalgravdropletenergy = 0;
double totalgravpoolenergy = 0;
double totalgravairenergy = 0;

double f1SurfaceArea;
double f2SurfaceArea;
double f3SurfaceArea;
double totalSurfaceArea;

double f1SurfaceEnergy;
double f2SurfaceEnergy;
double f3SurfaceEnergy;
double totalSurfaceEnergy;

double totalenergy;

event energiesandsurfaceareas (i++)
{
	foreach()
	{
		dudx[]=((u.x[1,0]-u.x[-1,0])/(2.0*Delta));
		dvdy[]=((u.y[0,1]-u.y[0,-1])/(2.0*Delta));
		dvdx[]=((u.y[1,0]-u.y[-1,0])/(2.0*Delta));
		dudy[]=((u.x[0,1]-u.x[0,-1])/(2.0*Delta));

		viewingfield[]=0.000*f1[]+0.5*f2[]+1.000*f3[];
		mymu[]=mu(f1[], f2[], f3[]);
		myrho[]=rho(f1[], f2[], f3[]);

		vel2[]=u.x[]*u.x[]+u.y[]*u.y[];
		mylevel[]=level;

		energydissipationrate[]=2.0*mu(f1[], f2[], f3[])*(dudx[]*dudx[]+dvdy[]*dvdy[]+(u.y[]*u.y[]/(cm[]*cm[])))
        +mu(f1[], f2[], f3[])*(dvdx[]+dudy[])*(dvdx[]+dudy[]);

		energydissipationrate_dV[]=energydissipationrate[]*Delta*Delta*2*M_PI*y;

		energydissipationratedroplet[]=energydissipationrate[]*f1[];
		energydissipationrateair[]=energydissipationrate[]*f2[];
		energydissipationratepool[]=energydissipationrate[]*f3[];

		kineticenergypervol[]=0.5*rho(f1[], f2[], f3[])*(u.x[]*u.x[] + u.y[]*u.y[]);

		kineticenergypervoldroplet[]=kineticenergypervol[]*f1[];
		kineticenergypervolair[]=kineticenergypervol[]*f2[];
		kineticenergypervolpool[]=kineticenergypervol[]*f3[];

		gravenergypervol[]=rho(f1[], f2[], f3[])*(x+0.5*SIZE*D0)/(FR*FR);

		gravenergypervoldroplet[]=gravenergypervol[]*f1[];
		gravenergypervolair[]=gravenergypervol[]*f2[];
		gravenergypervolpool[]=gravenergypervol[]*f3[];	
	}

	vorticity(u,myvort);

    currentdissipatedenergy = statsf(energydissipationrate).sum*2.0*pi;
	currentdissipateddropletenergy = statsf(energydissipationratedroplet).sum*2.0*pi;
	currentdissipatedpoolenergy = statsf(energydissipationratepool).sum*2.0*pi;
	currentdissipatedairenergy = statsf(energydissipationrateair).sum*2.0*pi;

    totaldissipatedenergy += currentdissipatedenergy*dt;
	totaldissipateddropletenergy += currentdissipateddropletenergy*dt;
	totaldissipatedpoolenergy += currentdissipatedpoolenergy*dt;
	totaldissipatedairenergy += currentdissipatedairenergy*dt;

	totalkineticenergy = statsf(kineticenergypervol).sum*2.0*pi;
	totalkineticdropletenergy = statsf(kineticenergypervoldroplet).sum*2.0*pi;
	totalkineticpoolenergy = statsf(kineticenergypervolpool).sum*2.0*pi;
	totalkineticairenergy = statsf(kineticenergypervolair).sum*2.0*pi;

	totalgravenergy = statsf(gravenergypervol).sum*2.0*pi;
	totalgravdropletenergy = statsf(gravenergypervoldroplet).sum*2.0*pi;
	totalgravpoolenergy = statsf(gravenergypervolpool).sum*2.0*pi;
	totalgravairenergy = statsf(gravenergypervolair).sum*2.0*pi;

	f1SurfaceArea=interface_area_axi(f1);
    f2SurfaceArea=interface_area_axi(f2);
    f3SurfaceArea=interface_area_axi(f3);
	totalSurfaceArea=f1SurfaceArea+f2SurfaceArea+f3SurfaceArea;

	f1SurfaceEnergy=f1SurfaceArea*SIGMA1;
	f2SurfaceEnergy=f2SurfaceArea*SIGMA2;
	f3SurfaceEnergy=f3SurfaceArea*SIGMA3;
	totalSurfaceEnergy=f1SurfaceEnergy+f2SurfaceEnergy+f3SurfaceEnergy;

	totalenergy=totaldissipatedenergy+totalkineticenergy+totalgravenergy+totalSurfaceEnergy;
}


event energyandsurfacearealogging (t += 0.001)
{
 	FILE * enlog = fopen(energylog, "a");
	fprintf(enlog, "%i %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f %1.5f\n"
	,i, t, currentdissipatedenergy, totaldissipatedenergy, currentdissipatedpoolenergy, totaldissipatedpoolenergy, 
	currentdissipateddropletenergy, totaldissipateddropletenergy, currentdissipatedairenergy, totaldissipatedairenergy, 
	totalkineticenergy, totalkineticdropletenergy, totalkineticpoolenergy, totalkineticairenergy,
	totalgravenergy, totalgravdropletenergy, totalgravpoolenergy, totalgravairenergy,
	totalSurfaceEnergy, f1SurfaceEnergy, f2SurfaceEnergy, f3SurfaceEnergy, totalenergy);
	fclose(enlog);

	FILE * arealog = fopen(surfacearealog, "a");
	fprintf(arealog, "%i %1.5f %1.5f %1.5f %1.5f %1.5f\n", i, t, f1SurfaceArea, f2SurfaceArea, f3SurfaceArea, totalSurfaceArea);
	fclose(arealog);   
}

// Adaptation is based on the errors on the tracer fields and 
// the vorticity
event controlledadapt (i++)
{
	scalar omega[];
	vorticity(u,omega);

	scalar ff1[], ff2[], ff3[];
		
	foreach()
	{
		ff1[]=f1[];
		ff2[]=f2[];
		ff3[]=f3[];
	}
	boundary({ff1, ff2, ff3});

	adapt_wavelet( {f1, f2, f3, ff1, ff2, ff3, u.x, u.y }, (double[])
	{f1tol, f2tol, f3tol, ff1tol, ff2tol, ff3tol, 3e-2, 3e-2}, MAXLEVEL, LEVEL-1);
}


event volumetracking (t += 0.01)
{
	if(i>0)
	{
		double totalf1vol = statsf(f1).sum*2.0*pi;
		double totalf2vol = statsf(f2).sum*2.0*pi;
		double totalf3vol = statsf(f3).sum*2.0*pi;

		FILE * vollog = fopen(volumeslog, "a");
		// Iteration, Time, f_Total, f2_Total, f3_Total
		fprintf(vollog, "%i %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", i, t, totalf1vol, totalf1vol/inittotalf1vol,
		totalf2vol, totalf2vol/inittotalf2vol, totalf3vol, totalf3vol/inittotalf3vol);
		fclose(vollog);
	}
}

// Remove small droplets (Noise)
event FilterDroplets (i++)
{
 		remove_droplets(f1, filtersize, filtertol);
 		remove_droplets(f1, filtersize, filtertol, true);
 		remove_droplets(f2, filtersize, filtertol);
 		remove_droplets(f2, filtersize, filtertol, true);
 		remove_droplets(f3, filtersize, filtertol);
 		remove_droplets(f3, filtersize, filtertol, true);
}

// Save video of simulation
event viewing (t += 0.005)
{
	view(width=1900, height=1050, fov=22.5, ty= 0, quat = { 0, 0, -0.707, 0.707 });
	
    clear();
	draw_vof("f1", lw=2);
	draw_vof("f2", lw=2);
	draw_vof("f3", lw=2);
	squares("viewingfield", map = cool_warm, min = -0.2, max = 1.2);
	mirror({0,1}) {
		draw_vof("f1", lw=2);	
		draw_vof("f2", lw=2);
		draw_vof("f3", lw=2);
		cells(lw=0.25);
		squares("mylevel", map = cool_warm, min = 5, max = 11);
	} 

	char timestring[100];
	sprintf(timestring, "t=%2.03f",t);
	draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);

	save(videoname);
}


// Useful information about simulation
event SimStats (i += 10) 
{
	timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
	FILE * runlog = fopen (runstats, "a");
	fprintf (runlog, "%i %g %g %ld %g %g\n", i, t, dt, grid->n, s.real, s.cpu);
	fclose(runlog);
	// i, time, no of cells, real time elapsed, cpu time
}


// Saves Gerris view snapshots which are quite useful
event snapshot (t += 0.1)
{
	char snapshotname[500];
	char snapshottimestamp[100];
	sprintf(snapshottimestamp, "Snapshot-%1.02f.gfs", t);

	strcpy(snapshotname, snapshotsfolder);
	strcat(snapshotname, snapshottimestamp);

	output_gfs(file = snapshotname, translate = true, t=t);
}

// Saves the interfaces 
event interfacelogging (t += 0.01)
{
	char interfacename[500];

	char interfacetimestamp[100];
	sprintf(interfacetimestamp, "Interface_%1.03f",t);

	strcpy(interfacename, interfacesfolder);
	strcat(interfacename, interfacetimestamp);

	char interface2name[500];
	strcpy(interface2name, interfacename);
	strcat(interface2name,"_2");

	char interface3name[500];
	strcpy(interface3name, interfacename);
	strcat(interface3name,"_3");

	char interfaceuppername[500];
	strcpy(interfaceuppername, interfacename);
	strcat(interfaceuppername,"_upper");

	char interfacelowername[500];
	strcpy(interfacelowername, interfacename);
	strcat(interfacelowername,"_lower");

	strcat(interfacename, "_1");

	FILE * fpp = fopen(interfacename, "w");
	output_facets(f1,fpp);
	fclose(fpp);

	FILE * fpp2 = fopen(interface2name, "w");
	output_facets(f2,fpp2);
	fclose(fpp2);

	FILE * fpp3 = fopen(interface3name, "w");
	output_facets(f3,fpp3);
	fclose(fpp3);

	FILE * fppupper = fopen(interfaceuppername, "w");
	output_facets(fupper,fppupper);
	fclose(fppupper);

	FILE * fpplower = fopen(interfacelowername, "w");
	output_facets(flower,fpplower);
	fclose(fpplower);
}

event linevels (t += 0.01)
{
	char velname00[200];
	char veltimestamp00[100];
	sprintf(veltimestamp00, "x001_Vels_%03.0f",100*t);
	strcpy(velname00, velsfolder);
	strcat(velname00, veltimestamp00);

	char velname05[200];
	char veltimestamp05[100];
	sprintf(veltimestamp05, "x005_Vels_%03.0f",100*t);
	strcpy(velname05, velsfolder);
	strcat(velname05, veltimestamp05);

	char velname10[200];
	char veltimestamp10[100];
	sprintf(veltimestamp10, "x010_Vels_%03.0f",100*t);
	strcpy(velname10, velsfolder);
	strcat(velname10, veltimestamp10);

	char velname15[200];
	char veltimestamp15[100];
	sprintf(veltimestamp15, "x015_Vels_%03.0f",100*t);
	strcpy(velname15, velsfolder);
	strcat(velname15, veltimestamp15);

	char velname20[200];
	char veltimestamp20[100];
	sprintf(veltimestamp20, "x020_Vels_%03.0f",100*t);
	strcpy(velname20, velsfolder);
	strcat(velname20, veltimestamp20);

	char velname25[200];
	char veltimestamp25[100];
	sprintf(veltimestamp25, "x025_Vels_%03.0f",100*t);
	strcpy(velname25, velsfolder);
	strcat(velname25, veltimestamp25);

	char velname30[200];
	char veltimestamp30[100];
	sprintf(veltimestamp30, "x030_Vels_%03.0f",100*t);
	strcpy(velname30, velsfolder);
	strcat(velname30, veltimestamp30);

	char velname35[200];
	char veltimestamp35[100];
	sprintf(veltimestamp35, "x035_Vels_%03.0f",100*t);
	strcpy(velname35, velsfolder);
	strcat(velname35, veltimestamp35);

	char velname40[200];
	char veltimestamp40[100];
	sprintf(veltimestamp40, "x040_Vels_%03.0f",100*t);
	strcpy(velname40, velsfolder);
	strcat(velname40, veltimestamp40);

	char velname45[200];
	char veltimestamp45[100];
	sprintf(veltimestamp45, "x045_Vels_%03.0f",100*t);
	strcpy(velname45, velsfolder);
	strcat(velname45, veltimestamp45);

	char velname50[200];
	char veltimestamp50[100];
	sprintf(veltimestamp50, "x050_Vels_%03.0f",100*t);
	strcpy(velname50, velsfolder);
	strcat(velname50, veltimestamp50);

	FILE * velslog00 = fopen(velname00, "w");
	FILE * velslog05 = fopen(velname05, "w");
	FILE * velslog10 = fopen(velname10, "w");
	FILE * velslog15 = fopen(velname15, "w");
	FILE * velslog20 = fopen(velname20, "w");
	FILE * velslog25 = fopen(velname25, "w");
	FILE * velslog30 = fopen(velname30, "w");
	FILE * velslog35 = fopen(velname35, "w");
	FILE * velslog40 = fopen(velname40, "w");
	FILE * velslog45 = fopen(velname45, "w");
	FILE * velslog50 = fopen(velname50, "w");

		for ( double xx = S0+filmthickness; xx > -1.5; xx += -0.001)
			{
				// x u.x u.y p 
				fprintf(velslog00, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.001), interpolate(u.y, xx, 0.001), interpolate(p, xx, 0.001));
				fprintf(velslog05, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.05), interpolate(u.y, xx, 0.05), interpolate(p, xx, 0.05));
				fprintf(velslog10, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.1), interpolate(u.y, xx, 0.1), interpolate(p, xx, 0.1));
				fprintf(velslog15, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.15), interpolate(u.y, xx, 0.15), interpolate(p, xx, 0.15));
				fprintf(velslog20, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.2), interpolate(u.y, xx, 0.2), interpolate(p, xx, 0.2));
				fprintf(velslog25, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.25), interpolate(u.y, xx, 0.25), interpolate(p, xx, 0.25));
				fprintf(velslog30, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.3), interpolate(u.y, xx, 0.3), interpolate(p, xx, 0.3));
				fprintf(velslog35, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.35), interpolate(u.y, xx, 0.35), interpolate(p, xx, 0.35));
				fprintf(velslog40, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.4), interpolate(u.y, xx, 0.4), interpolate(p, xx, 0.4));
				fprintf(velslog45, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.45), interpolate(u.y, xx, 0.45), interpolate(p, xx, 0.45));
				fprintf(velslog50, "%g %g %g %g\n", xx, interpolate(u.x, xx, 0.5), interpolate(u.y, xx, 0.5), interpolate(p, xx, 0.5));
			}

	fclose(velslog00);
	fclose(velslog05);
	fclose(velslog10);
	fclose(velslog15);
	fclose(velslog20);
	fclose(velslog25);
	fclose(velslog30);
	fclose(velslog35);
	fclose(velslog40);
	fclose(velslog45);
	fclose(velslog50);
}


// Defines end time and tasks to do at the end
event end (t=2.0)
{
	// End Timing
	timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
	FILE * runlog = fopen (runstats, "a");
	fprintf (runlog, "%i %g %g %ld %g %g\n", i, t, dt, grid->n, s.real, s.cpu);
	fclose(runlog);
	// i, time, no of cells, real time elapsed, cpu time

	// End Volume
	double totalf1vol = statsf(f1).sum*2*pi;
	double totalf2vol = statsf(f2).sum*2*pi;
	double totalf3vol = statsf(f3).sum*2*pi;

	FILE * vollog = fopen(volumeslog, "a");
	// Iteration, Time, f_Total, f2_Total, f3_Total
	fprintf(vollog, "%i %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n", i, t, totalf1vol, totalf1vol/inittotalf1vol,
	totalf2vol, totalf2vol/inittotalf2vol, totalf3vol, totalf3vol/inittotalf3vol);
	fclose(vollog);

	printf("\nFinished Simulation Re %4.1f, We %3.1f, Fr %2.1f, Visc Ratio %1.02f & Max Level %i\n\n", 
	RE, WE, FR, MU_POOL2DROP, MAXLEVEL);
}
