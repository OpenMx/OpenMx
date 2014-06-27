// From http://cran.r-project.org/web/packages/expm/index.html
// License: GPL2 / GPL3

/*  ===== File part of R package expm =====
 *
 *  Function to compute the matrix logarithm
 *
 *     log(M) = L such that
 *
 *     M  = exp(L) where
 *
 *     exp(L) = sum(n = 0:Inf; L^n / n!),
 *
 *  where M and L are an (n x n) matrix.
 *
 *  The functions therein use LAPACK and BLAS routines. Nicely
 *  formatted man pages for these can be found at
 *
 *    <http://www.mathkeisan.com/UsersGuide/E/>
 *
 *  AUTHORS: Christophe Dutang, based on code eigen,
 *
 *  i.e., function 'modLa_rg' and 'modLa_dgesv' in R's
 *  <Rsource>/src/modules/lapack/lapack.c, used in eigen()
 */

#include "glue.h"
#define _(s) s   // don't worry about localization

void logm_eigen(double *x, int n, double *z, double tol)
{
    if (n == 1)
        z[0] = log(x[0]);		/* scalar logarithm */
    else
    {
        const char *transa = "N";
        const int nsqr = n * n;

        const Rcomplex cone = {1., 0.}, czero = {0., 0.};
        int i, j;
        int info, lwork, is_conjug, is_diag;
        double onenorm, rcond, tmp;
        char jobVL[1], jobVR[1];

        /* Arrays */
        int *ipiv = (int *) R_alloc(n, sizeof(int)); /* permutation vector */
        double *left, *right, *workdiag; /* left and right eigenvectors and workspace for diagonalisation */
        double *wR = (double *) R_alloc(n, sizeof(double)); /* real part of eigenvalues */
        double *wI = (double *) R_alloc(n, sizeof(double)); /* imaginary part of eigenvalues */
        double *rworksing = (double *) R_alloc(2*n, sizeof(double)); /* working vector to test the singularity */
        Rcomplex *eigvect = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* (right) eigenvectors matrix */
        Rcomplex *eigvectinv = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* its inverse */
        Rcomplex *logeigval; /* complex matrix diag(log(eigenvalues)) */
        Rcomplex *ctmp = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex)); /* temp working variable */
        Rcomplex *worksing = (Rcomplex *) R_alloc(2*n, sizeof(Rcomplex)); /* workspace to test the singularity */

        Memcpy(z, x, nsqr);

        /* Test if x is diagonalisable by computing its eigenvalues and (right) eigenvectors */
        /* code based on modLa_rg in lapack.c, used in eigen.R */
        jobVL[0] = 'N';
        left = (double *) 0;
        jobVR[0] = 'V';
        right = (double *) R_alloc(nsqr, sizeof(double));

        /* 1 - ask for optimal size of work array */
        lwork = -1;
        F77_CALL(dgeev)(jobVL, jobVR, &n, z, &n, wR, wI,
			left, &n, right, &n, &tmp, &lwork, &info);
        if (info != 0)
            Rf_error(_("error code %d from Lapack routine dgeev"), info);
        lwork = (int) tmp;
        workdiag = (double *) R_alloc(lwork, sizeof(double));

        /* 2 - compute eigenvalues and (right) eigenvectors */
        F77_CALL(dgeev)(jobVL, jobVR, &n, z, &n, wR, wI,
			left, &n, right, &n, workdiag, &lwork, &info);
        if (info != 0)
            Rf_error(_("error code %d from Lapack routine dgeev"), info);

        /* try to invert the eigenvectors matrix */
        /* 1 - build the Rcomplex matrix with eigenvectors */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                is_conjug = 0;
                if(i < n-1)
                {	/* conjugate eigenvalues */
                    if(wR[i] == wR[i+1] && wI[i] == -wI[i+1] && wI[i] != 0.0)
                    {
                        is_conjug = 1;
                        eigvect[i * n + j].r = right[i * n + j];
                        eigvect[i * n + j].i = right[(i+1) * n + j];
                    }
                }
                if(i > 0)
                {	/* conjugate eigenvalues */
                    if(wR[i] == wR[i-1] && wI[i] == -wI[i-1] && wI[i] != 0.0)
                    {
                        is_conjug = 1;
                        eigvect[i * n + j].r = right[(i-1) * n + j];
                        eigvect[i * n + j].i = -right[i * n + j];
                    }
                }
                /* real eigenvalues */
                if(!is_conjug)
                {
                    eigvect[i * n + j].r = right[i * n + j];
                    eigvect[i * n + j].i = 0.0;
                }
                /* eigvectinv initialise with the identity matrix */
                eigvectinv[i * n +j].r = (i == j) ? 1.0 : 0.0;
                eigvectinv[i * n +j].i = 0.0;
            }
        }

        /* 2 - store the matrix eigvect (because function zgesv will change it) */
        Memcpy(ctmp, eigvect, nsqr);

        /* 3 - solve a linear complex equation system with eigvectinv equals
         * to matrix identity. hence, on exit eigvectinv contains the
         * inverse of complex matrix eigvect. code base on solve.R */
        F77_CALL(zgesv)(&n, &n, eigvect, &n, ipiv, eigvectinv, &n, &info);

        if (info > 0)
            is_diag = 0; //matrix eigvect is exactly singular.
        if (info < 0)
            Rf_error(_("argument %d of Lapack routine dgesv had invalid value"), -info);
        if (info == 0)
            is_diag = 1;

        /* check if matrix eigvectinv is numerically singular */
        if (is_diag)
        {
            /* compute the reciprocal condition number of eigvectinv. */

            /* 1 - compute the one norm of the matrix eigvectinv */
		char one = '1';
            onenorm = F77_CALL(zlange)(&one, &n, &n, eigvectinv, &n, (double*) NULL);

            /* 2 - estimates the reciprocal of the condition number
             * when the one norm is used. */
            F77_CALL(zgecon)(&one, &n, eigvectinv, &n, &onenorm, &rcond, worksing, rworksing, &info);

            if (rcond < tol)
                is_diag=0;
        }


        if (is_diag)
        {

            /* x is diagonalisable so
             * compute complex matrix operations :
             * eigvect %*% diag(log(eigenvalues)) %*% eigvectinv */

            /* 1 - logeigval is the complex matrix diag(log(eigenvalues)) */
            /* code based on z_log in complex.c */
            logeigval = (Rcomplex *) R_alloc(nsqr, sizeof(Rcomplex));
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if(i == j)
                    {
                        logeigval[i * n +j].r = log( sqrt(wR[i] * wR[i] + wI[i] * wI[i]) ) ;
                        logeigval[i * n +j].i = atan2( wI[i], wR[i] ) ;
                    }
                    else
                    {
                        logeigval[i * n +j].r = 0.0;
                        logeigval[i * n +j].i = 0.0;
                    }
                }
            }

            /* 2 - restore the matrix eigvect */
            Memcpy(eigvect, ctmp, nsqr);

            /* 3 - compute (complex) matrix product: ctmp <- eigvect * logeigval */
            F77_CALL(zgemm)(transa, transa, &n, &n, &n, &cone, eigvect, &n,
			    logeigval, &n, &czero, ctmp, &n);

            /* 4 - compute (complex) matrix product: logeigval <- ctmp * eigvectinv */
            F77_CALL(zgemm)(transa, transa, &n, &n, &n, &cone, ctmp, &n,
			    eigvectinv, &n, &czero, logeigval, &n);

            //TOCHECK
            /* store the real part in z */
            /* the matrix logarithm is always real,
             * even if x has complex eigen values. */
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
            {   z[i * n + j] = logeigval[i * n + j].r;
            }
            }


        }
	else
	    Rf_error("non diagonalisable matrix");
    }
}

/* Main function, the only one used by .Call(). */
SEXP do_logm_eigen(SEXP x, SEXP tolin)
{
    SEXP dims, z;
    int n, m;
    double *rx = REAL(x), *rz;
    double tol = Rf_asReal(tolin);

    if (!Rf_isNumeric(x) || !Rf_isMatrix(x))
        Rf_error(_("invalid argument"));

    dims = Rf_getAttrib(x, R_DimSymbol);
    n = INTEGER(dims)[0];
    m = INTEGER(dims)[0];
    if (n != m)
        Rf_error(_("non-square matrix"));
    if (n == 0)
        return(Rf_allocVector(REALSXP, 0));

    PROTECT(z = Rf_allocMatrix(REALSXP, n, n));
    rz = REAL(z);

    logm_eigen(rx, n, rz, tol);

    Rf_setAttrib(z, R_DimNamesSymbol, Rf_getAttrib(x, R_DimNamesSymbol));

    UNPROTECT(1);

    return z;
}
