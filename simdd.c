/*******************************************************
 Drift-Diffusion Simulator (1-D)
 (Decoupled Method)
 Adapted from a Japanese Master thesis by Qikai Li
 Zhigang Shuai Group, ICCAS, Beijing 100190
 Email: qkli@iccas.ac.cn

 Last modified: April 12, 2009
 *******************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define true    1
#define false   0
#define bool    int

#define ECHARGE 1.602189e-19 	    // charge of an electron [c]
#define EMASS   9.109534e-31 	    // mass of an electron [kg] 
#define BOLTZ   1.380662e-23 	    // Boltzmann const. [J/K]
#define PI      3.141593e+00 	    // PI
#define EPS0    8.854188e-12 	    // dielectric constant in vacuum [F/m]
#define NM      10000 		        // Max number of nodes

/* eqType: 0 - Poisson equation, 
           1 - electron continuity equation
           1 - electron continuity equation */
enum {EQ_POISSON, EQ_ELECTRON_CONT, EQ_HOLE_CONT};

int    DDS_kernel(double vmax, double vdlt);
void   setup_tridiag_matrix(int eqType, double *at, double *bt, double *ct, double *ft, int ll2);
void   tridag(double *a, double *b, double *c, double *r, double *u, int n);
int    Poisson(double *at,double *bt,double *ct,double *ft,double *ut);
int    electron_continuity(double *at,double *bt,double *ct,double *ft,double *ut);
int    hole_continuity(double *at,double *bt,double *ct,double *ft,double *ut);
void   clear_vector(double *a, int n);
void   setup_env(double *acc,double *don,double *volt,double *cgm,double *elec,double *hole);
void   setup_recombi_array (int eqType, int ll, double *rcmb);
void   setup_mobility_array(int eqType, int ll, double *doud);
double Bernoulli(double x);
void   output(double *, double *, double *);

int    lnd, lpd, ltot;
int    grid_n, grid_n2;					// grid lattice count, & the count excluding two end points
int    ibmax;							// maximum bias voltage step number
double dltex;							// dx 
double rmass, eps;						// effective mass, dielectric constant
double cdin;							// intrinsic carrier density
double tau_p, tau_n;					// lift time of hole and electron, in [s]
double VT, Temp;						// VT = q/(kB*T), inversion of the normal thermal voltage definition
double vmax, vdlt, vold;				// voltage max, delta
double volt[NM], elec[NM], hole[NM];	// voltage, electron, & hole arrays
double don[NM], acc[NM], cgm[NM];       // donor, acceptor, & charge arrays
double ecurrt[NM], hcurrt[NM];			// electron current, & hole current arrays
double dl1[NM], dl2[NM], dudn[NM];		// electron density (n) arrays
double dp1[NM], dp2[NM], dudp[NM];		// hole density (p) arrays
double gam1[NM], gam2[NM], gam3[NM];	// discretization coefficients of 2nd derivative
double h_sp[NM], h_prm[NM];				// discretization distance h's
int    iter_max;						// maximum number of iteration
double Tol[3];							// converge limits for 3 equations


int main()
{
    // Parameter of Si
    rmass = 0.328*EMASS;    // effective mass
    eps   = 11.9*EPS0;      // absolute dielectric constant
    cdin  = 1.105e10*1e6;   // intrinsic carrier density [m^3]
    tau_p = 10.0e-9;		// in [s]
    tau_n = 10.0e-9;		// in [s]

    // setting up grid, impurity, temperature & initials
    setup_env(acc, don, cgm, volt, elec, hole);

    // max iteration
    iter_max = 1000;

	// max tolerance
    Tol[0] = 1.0e-5;		// Poisson equation
    Tol[1] = 1.0e-5;		// electron continuity equation
    Tol[2] = 1.0e-5;		// hole continuity equation

    // Applying voltage
    vmax = 0.5;				// max voltage
    vdlt = 0.1;				// voltage increment

    // Drift-Diffusion Simulator kernel
    DDS_kernel(vmax, vdlt);

    output(volt,elec,hole);

    return 0;
}


/*
  setting up grid, impurity, temperature & initial values
*/
void setup_env(double *acc,double *don,double *cgm,double *volt,double *elec,double *hole)
{
    double p_dope, n_dope;
    int i;

    printf("Device structure pn-diode\n");

    // grid
    lpd  = 100;				// number of p domain
    lnd  = 100;				// number of n domain
    ltot = lpd + lnd;		// total number of domains
    grid_n   = ltot + 1;	// total grid number
    grid_n2  = grid_n - 2;	// total grid number excluding two boundary points

    // dx [m]
    dltex = 10.0e-9;

    // temperature [K]
    Temp = 300.0;

    // same doping concentration
    p_dope = 1.0e16;		// [cm^-3]
    n_dope = 1.0e16;		// [cm^-3]

    // setup doping profile
    p_dope = p_dope*1.0e6;
    n_dope = n_dope*1.0e6;

    clear_vector(don,NM); 
	clear_vector(acc,NM);

    for (i = 0; i < grid_n; i++)
    {
        if (lpd <= i && i < grid_n)
        {
            // n-region
            don[i] = n_dope;
            acc[i] = 0.0;
        } else
        {
            // p-region
            don[i] = 0.0;
            acc[i] = p_dope;
        }

        cgm[i] = don[i] - acc[i];
    }

    VT = ECHARGE/(BOLTZ*Temp);		// inversion of thermal voltage, q/(kB*T)
	
	// grid space distance
    for (i = 0; i < grid_n; i++)
    {
        h_sp[i]  = dltex;
        h_prm[i] = (dltex+dltex)/2.0;
    }

    clear_vector(gam1,NM);
    clear_vector(gam2,NM);
    clear_vector(gam3,NM);

	// g1(n)V(n-1) + g2(n)V(n) + g3(n)V(n+1) = - q * [Nd(n) - Na(n) + p(n) - n(n)]/epsi0
	// coeffients of discretized 2nd derivative: g1(n) = g3(n) = 1/h^2, g2(n) = - 2 / h^2
    for (i = 1; i < grid_n-1; i++)
    {
        gam1[i] = 1.0 / h_sp[i-1] / h_prm[i];
        gam3[i] = 1.0 / h_sp[i] / h_prm[i];
        gam2[i] = -gam1[i] - gam3[i];
    }

    // initial values
    for (i = 0; i < grid_n; i++)
    {
        // voltage
        volt[i] = 0.0;

        if (cgm[i] < 0.0)
            volt[i] = log(-cdin/cgm[i])/VT;
        if (cgm[i] > 0.0)
            volt[i] = log(cgm[i]/cdin)/VT;

        vold = volt[0];

        // electron
        elec[i] = don[i];

        // hole
        hole[i] = acc[i];
    }

    hole[0] = -cgm[0] * (1.0+sqrt(1.0+(2.0*cdin/cgm[0])*(2.0*cdin/cgm[0]))) / 2.0;
    elec[0] = cdin*cdin/hole[0];
    elec[grid_n-1] = cgm[grid_n-1] * (1.0+sqrt(1.0+(2.0*cdin/cgm[grid_n-1])*(2.0*cdin/cgm[grid_n-1]))) / 2.0;
    hole[grid_n-1] = cdin*cdin/elec[grid_n-1];
}


/***************************************************
  Drift-Diffusion Simulator (DDS) kernel
  A(n) U(n-1) + B(n) U(n) + C(n) U(n+1) = F(n)
 ***************************************************/
int DDS_kernel(double vmax, double vdlt)
{
    double at[NM], bt[NM], ct[NM], ft[NM], ut[NM];     // TriDiagonal Matrix, solution in ut[]
    double built, vnew, sumcur, ecr, hcr;
    int ib, it, lj;
    int iter_ret[3];

    // check the voltage increment
    if (fabs(vdlt) > fabs(vmax))
    {
        if (vmax != 0.0)
        {
            printf("ERROR: V-step (vdlt) too large!\n");
            exit(1);
        }
    }

    // reverse bias voltage?
    if (vmax < 0.0)
        vdlt = -vdlt;

    ibmax = (int)(fabs(vmax/vdlt+1.0e-5)) + 1;
    printf("Maximum bias voltage step number: %d\n", ibmax);
    
    vnew = vold+0.0;

    for (ib = 0; ib < ibmax; ib++)
    {
        // loop for bias
        volt[0] = vnew;
        for (it = 0; it < iter_max; it++)
        {
            iter_ret[0] = Poisson            (at, bt, ct, ft, ut);
            iter_ret[1] = electron_continuity(at, bt, ct, ft, ut);
            iter_ret[2] = hole_continuity    (at, bt, ct, ft, ut);

            if (iter_ret[0] == 1 && iter_ret[1] == 1 && iter_ret[2] == 1)
                goto CONVERGED;
        }

        fprintf(stderr, "Failed to CONVERGE!\n");
        exit(1);

CONVERGED:
		// Update generation rate

        lj = grid_n/2;
        ecurrt[ib] = ecr = ECHARGE / h_prm[lj] * (dl1[lj]*elec[lj] + dl2[lj]*elec[lj+1]);
        hcurrt[ib] = hcr = ECHARGE / h_prm[lj] * (dp1[lj]*hole[lj] + dp2[lj]*hole[lj+1]);
        sumcur = ecr + hcr;
        built = volt[grid_n-1] - volt[0];

        /* printf("built-in: %9.5f \n", volt[grid_n-1]-volt[0]); */
        printf("%5d %9.3f %10.5f %15.5E %15.5E %15.5E\n", 
			   ib+1, ib*vdlt, built, ecr/1.0e4, hcr/1.0e4, (ecr+hcr)/1.0e4);

        vnew = vnew + vdlt;
    }
    
    return 0;
}


/*********************************************************************
 Poisson solver
 A(n) U(n-1) + B(n) U(n) + C(n) U(n+1) = F(n)
 *********************************************************************/
int Poisson(double *at, double *bt, double *ct, double *ft, double *ut)
{
    double error, max_error;
    int iter, lj, eqType;

    eqType = EQ_POISSON; 	// Type of the Poisson equation

    for (iter = 1; iter < iter_max; iter++)
    {
        setup_tridiag_matrix(eqType, at, bt, ct, ft, grid_n2);
        
        tridag(at, bt, ct, ft, ut, grid_n2);
        
        max_error = -1.0;
        for (lj = 1; lj < grid_n-1; lj++)
        {
			// Gummel's scheme:
            // update potential, electron density & hole density
            volt[lj] = volt[lj] + ut[lj-1];				// V(i+1) = V(i) + dV(i)
            elec[lj] = elec[lj] * exp(VT*ut[lj-1]);		// n(i+1) = n(i) * exp(VT*dV(i))
            hole[lj] = hole[lj] * exp(-VT*ut[lj-1]);	// p(i+1) = p(i) * exp(-VT*dV(i))

            // search max relative error
            if (fabs(volt[lj]) == 0.0)
            {
                fprintf(stderr, "ERROR: volt[] at lj = %d is zero!\n", lj);
                exit(1);
            } else
            {
                error = fabs(ut[lj-1]/volt[lj]);
                
                if (error > max_error)
                    max_error = error;
            }
        }

        if (max_error < Tol[eqType])
        {
            return 1;				// normal exit here if converged!
        }
    }

	// if reach here, oops!
    fprintf (stderr,"Failed to converge in Poisson()!\n");
    exit(1);
    
    return 0;
}


/***************************************************
 * Electron continuity solver
 ***************************************************/
int electron_continuity(double *at,double *bt,double *ct,double *ft,double *ut)
{
    double error, max_error;
    int iter, lj, eqType;

    eqType = EQ_ELECTRON_CONT;
    
    for (iter = 1; iter < iter_max; iter++)
    {
        setup_tridiag_matrix(eqType, at, bt, ct, ft, grid_n2);
        
        tridag(at, bt, ct, ft, ut, grid_n2);
        
        max_error = -1.0;
        for (lj = 1; lj < grid_n-1; lj++)
        {
            /* renew electron density */
            elec[lj] = elec[lj] + ut[lj-1];

            /* search max relative error */
            if (fabs(elec[lj]) == 0.0)
            {
				fprintf(stderr, "ERROR: elec[] at lj = %d is zero!\n", lj);
                exit(1);
            } else
            {
                error = fabs(ut[lj-1]/elec[lj]);
                if (error > max_error)
                    max_error = error;
            }
        }

        if (max_error < Tol[eqType])
        {
            return 1;				// normal exit here if converged!
        }
    }

	// if reach here, oops!
    fprintf (stderr, "Failed to converge in electron_continuity()!\n");
    exit(1);
    return 0;
}


/*******************************************************
 * Hole continuity solver
 *******************************************************/
int hole_continuity(double *at, double *bt, double *ct, double *ft, double *ut)
{
    double error, max_error;
    int iter, lj, eqType;

    eqType = EQ_HOLE_CONT;
    
    for (iter = 1; iter < iter_max; iter++)
    {
        setup_tridiag_matrix(eqType, at, bt, ct, ft, grid_n2);
        
        tridag(at, bt, ct, ft, ut, grid_n2);
        
        max_error = -1.0;
        for (lj = 1; lj < grid_n-1; lj++)
        {
            /* new potential, electron density, hole density... */
            hole[lj] = hole[lj] + ut[lj-1];

            /* search max relative error */
            if (fabs(hole[lj]) == 0.0)
            {
				fprintf(stderr, "ERROR: hole[] at lj = %d is zero!\n", lj);
                exit(1);
            } else
            {
                error = fabs(ut[lj-1]/hole[lj]);
                if (error > max_error)
                    max_error = error;
            }
        }

        if (max_error < Tol[eqType])
        {
            return 1;			// normal exit here if converged!
        }
    }

	// if reach here, oops!
    fprintf (stderr, "Failed to converge in hole_continuity()!\n");
    exit(1);

    return 0;
}


/**************************************************
 * Output results
 **************************************************/
void output(double *volt,double *elec,double *hole)
{
    FILE *fp0, *fp1, *fp2, *fp3;
    double cjn, cjp, dvl_el, dvl_hl, field[NM];
    double qro, sarea, cr1, cr2, ecr, hcr;
    int ib, lj;

    fp0 = fopen("cc17.txt", "w");
    fp1 = fopen("cc18.txt", "w");
    fp2 = fopen("cc19.txt", "w");
    fp3 = fopen("cc20.txt", "w");

    for (lj = 0; lj < grid_n; lj++)
    {
        qro = don[lj]-acc[lj]+hole[lj]-elec[lj];
        fprintf(fp0,"%8.1f %14.5e %14.5e %14.5e %14.5e\n", lj*dltex*1.0e9,
                qro/1.0e6, hole[lj]/1.0e6, elec[lj]/1.0e6, cgm[lj]/1.0e6);
    }

    clear_vector(field, NM);

    for (lj = 0; lj < grid_n; lj++)
    {
        cjn = cjp = 0.0;

        if (lj == 0)
            field[lj] = 0.0;

        if (lj != 0 && lj < grid_n-1)
        {
            field[lj] = (volt[lj-1] - volt[lj+1])/h_sp[lj]/2.0;
            cjn = ECHARGE/h_prm[lj]*(dl1[lj]*elec[lj]+dl2[lj]*elec[lj+1]);
            cjp = ECHARGE/h_prm[lj]*(dp1[lj]*hole[lj]+dp2[lj]*hole[lj+1]);
        }

        if (lj == grid_n-1)
            field[lj] = field[grid_n-1];

        fprintf(fp1, "%8.1f %9.5f ", lj*dltex*1.0e9, volt[lj]);
        fprintf(fp1, "%13.5E %13.5E ", field[lj]/1.0e5, cjn/1.0e4);
        fprintf(fp1, "%13.5E %13.5E\n",cjp/1.0e4, (cjn+cjp)/1.0e4);

        dvl_el = 0.0;
        if (elec[lj] != 0.0)
            dvl_el = cjn/elec[lj]/ECHARGE;

        dvl_hl = 0.0;
        if (hole[lj] != 0.0)
            dvl_hl = cjp/hole[lj]/ECHARGE;

        fprintf(fp2, "%8.1f %10.5f ", lj*dltex*1e9, volt[lj]);
        fprintf(fp2, "%13.5E %13.5E\n", elec[lj]/1e6,hole[lj]/1e6);
        // fprintf(fp2, "%13.5E %13.5E£¤n", dvl_el*1e2,dvl_hl*1e2);
    }

    /* at each step */
    sarea = 1.0;
    cr1 = cr2 = ecr = hcr = 0.0;
    for (ib = 0; ib < ibmax; ib++)
    {
        if (ib >= 0)
        {
            cr1 = ecurrt[ib];
            cr2 = hcurrt[ib];
            ecr = sarea*cr1/1e4;
            hcr = sarea*cr2/1e4;
        }
        fprintf(fp3,"%5d %9.3f %15.5E %15.5E %15.5E\n",ib+1, ib*vdlt, ecr, hcr, ecr+hcr);
    }

    fclose(fp0); 
    fclose(fp1); 
    fclose(fp2); 
    fclose(fp3);
}


/*
  Set all vector components to 0's
 */
void clear_vector(double *a, int n)
{
    int i;
    
    for (i = 0; i < n; i++)
    {
        a[i] = 0.0;
    }
}


/************************************************************
 * Setup tridiagonal matrix, i.e., at[], bt[], ct[], & ft[]
 * eqType == 0: Poisson equation 
 * eqType == 1: electron continuity equation 
 * eqType == 2: hole continuity equation 
 ************************************************************/
void setup_tridiag_matrix(int eqType, double *at, double *bt, double *ct, double *ft, int grid_n2)
{
    double beta[NM];	// handy array for Bernoulli function
	double doud[NM];	// mobility array
	double rcmb[NM];	// recombination array
    double dnmnt;
    int mt, nt;

    // for Poisson's equation
    if (eqType != EQ_POISSON)
    {
        setup_mobility_array(eqType, grid_n2+2, doud);		// mobility

        setup_recombi_array (eqType, grid_n2+2, rcmb);		// recombination

        clear_vector(beta,NM);							// holds [V(n)-V(n+1)]/VT

        for (mt = 0; mt < grid_n2+1; mt++)
        {
            beta[mt] = VT * (volt[mt] - volt[mt+1]);
        }
    }

    // for electron continuity equation
    if (eqType == EQ_ELECTRON_CONT)
    {
        clear_vector(dl1,NM);
        clear_vector(dl2,NM);
        clear_vector(dudn,NM);

        for (mt = 0; mt < grid_n2+1; mt++)
        {
            dnmnt    = tau_n*(hole[mt] + cdin) + tau_p*(elec[mt] + cdin);
            dudn[mt] = (elec[mt] - tau_p*rcmb[mt]) / dnmnt;

            dl1[mt]  = doud[mt] * (-Bernoulli(beta[mt])) / VT;
            dl2[mt]  = doud[mt] * Bernoulli(-beta[mt]) / VT;
        }
    }

    // for hole continuity equation
    if (eqType == EQ_HOLE_CONT)
    {
        clear_vector(dp1,NM);
        clear_vector(dp2,NM);
        clear_vector(dudp,NM);

        for (mt = 0; mt < grid_n2+1; mt++)
        {
            dnmnt    = tau_n*(hole[mt] + cdin) + tau_p*(elec[mt] + cdin);
            dudp[mt] = (elec[mt] - tau_n*rcmb[mt]) / dnmnt;

            dp1[mt]  = doud[mt] * Bernoulli(-beta[mt]) / VT;
            dp2[mt]  = doud[mt] * (-Bernoulli(beta[mt])) / VT;
        }
    }

    for (nt = 1; nt < grid_n2+1; nt++)
    {
        if (eqType == EQ_POISSON)		// make matrix for Poisson's equation
        {
            at[nt-1] = gam1[nt];
            bt[nt-1] = gam2[nt] - ECHARGE*VT/eps*(elec[nt]+hole[nt]);
            ct[nt-1] = gam3[nt];

            ft[nt-1] = -(cgm[nt]-elec[nt]+hole[nt])*ECHARGE/eps
                       -gam1[nt]*volt[nt-1] - gam2[nt]*volt[nt] - gam3[nt]*volt[nt+1];
        } 
		else if (eqType == EQ_ELECTRON_CONT)
        {
            at[nt-1] = -dl1[nt-1]*gam1[nt];
            bt[nt-1] = gam3[nt]*dl1[nt] - gam1[nt]*dl2[nt-1] - dudn[nt];
            ct[nt-1] = dl2[nt]*gam3[nt];

            ft[nt-1] = -(dl1[nt]*elec[nt]+dl2[nt]*elec[nt+1])*gam3[nt]
                       + (dl1[nt-1]*elec[nt-1] + dl2[nt-1]*elec[nt])*gam1[nt] + rcmb[nt];
        } 
		else if (eqType == EQ_HOLE_CONT)
        {
            at[nt-1] = -dp1[nt-1]*gam1[nt];
            bt[nt-1] = gam3[nt]*dp1[nt] - gam1[nt]*dp2[nt-1] + dudp[nt];
            ct[nt-1] = dp2[nt]*gam3[nt];

            ft[nt-1] = -(dp1[nt]*hole[nt]+dp2[nt]*hole[nt+1])*gam3[nt]
                       +(dp1[nt-1]*hole[nt-1] + dp2[nt-1]*hole[nt])*gam1[nt] - rcmb[nt];
        }
    }
}


/*
  Solve N x N tridiagonal matrix equation (M U = W -> U)
  Let a, b, and c be the vectors of the left, center and right diagonal elements of the matrix, respectively.
  Note that a1 and c_N are undefined, and can be conveniently set to zero.

  See Numerical Recipes, ch. 2.4
 */
void tridag(double *a, double *b, double *c, double *f, double *u, int n)
{
    double gam[NM];
	double bet;
    int j;

    if (b[0] == 0.0)
    {
        fprintf(stderr, "Error 1 in tridag!\n");
        exit(1);
    }

	// If this happens then you should rewrite your equation as a set of order N-1,
	// then u1 trivially eliminated.
    bet  = b[0];
    u[0] = f[0] / bet;

    for (j = 1; j < n; j++)
    {
        gam[j] = c[j-1] / bet;
        bet = b[j] - a[j]*gam[j];

        if (bet == 0.0)
        {
            fprintf(stderr, "Error 2 in tridag!\n");
            exit(1);
        }

        u[j] = (f[j] - a[j]*u[j-1]) / bet;
    }

    for (j = n-2; j >= 0; j--)
    {
        u[j] = u[j] - gam[j+1]*u[j+1];		// Back-substitution
    }
}


/*
  Bernoulli function, segmented to improve robustness.
 */
#define X1 -38.0
#define X2 -2.0e-6
#define X3 2.0e-6
#define X4 38.0
#define X5 40.0
double Bernoulli(double x)
{
    if (x < X1)
        return -x;
    else if (X1< x && x <= X2)
        return x/(exp(x)-1.0);
    else if (X2 < x && x <= X3)
        return 1.0-x/2.0;
    else if (X3 < x && x <= X4)
        return x*exp(-x)/(1.0 - exp(-x));
    else if (X4 < x && x <= X5)
        return x*exp(-x);
    else
        return 0.0;
}


/*
  Recombination model (Shockley-Read-Hall model) 

  U(x) = Up(x) = Un(x) = (p(x)*n(x)-ni^2) / (tau_p*[n(x) + nt(x)] + tau_n*[p(x) + pt(x)])

  See Semiconductor Devices: Physics and Technology (2nd, by Sze, p. 64£©
 */
void setup_recombi_array(int eqType, int ll, double *rcmb)
{
    double a1, b1;
    int mt;

    if (eqType == EQ_POISSON)
    {
        fprintf(stderr, "ERROR: eqType == 0 in setup_recombi_array()\n");
        exit(1);
    }

    if (ll > NM)
    {
        fprintf(stderr, "ERROR: ll > NM in setup_recombi_array()\n");
        exit(1);
    }

    if (ll < 3)
    {
        fprintf(stderr, "ERROR: ll < 3 in setup_recombi_array()\n");
        exit(1);
    }

    clear_vector(rcmb,ll);

    for (mt = 1; mt <= ll-1; mt++)
    {
        a1 = elec[mt]*hole[mt] - cdin*cdin;
        b1 = tau_p*(elec[mt]+cdin) + tau_n*(hole[mt]+cdin);

        rcmb[mt] = a1 / b1;
    }
}


/*
  Mobility model (Simplified model)
 */
void setup_mobility_array(int eqType, int ll, double *doud)
{
    int mt;

    if (eqType == EQ_POISSON)
    {
        printf("setup_mobility_array() eqType = 0\n");
        exit(1);
    }

    if (ll > NM)
    {
        printf("ll > NM in setup_mobility_array()\n");
        exit(1);
    }

    if (ll < 1)
    {
        printf("ll < 1 in setup_mobility_array()\n");
        exit(1);
    }

    clear_vector(doud,ll);

    for (mt = 0; mt <= ll-1; mt++)
    {
        if (eqType == EQ_ELECTRON_CONT)
			doud[mt] = 4.0e5 * pow(Temp, -2.6);
        else if (eqType == EQ_HOLE_CONT)
            doud[mt] = 2.5e4 * pow(Temp, -2.3);
    }

    doud[ll] = doud[ll-1];
}
