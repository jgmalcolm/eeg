#include <mex.h>
#include <math.h>


// squared distance
static double dist2(const double *x, const double *y, int nch)
{
    double d = 0;
    for (int i = 0; i < nch; ++i) {
        double diff = x[i] - y[i];
        d += diff*diff;
    }
    return d;
}

static double min(double x, double y)
{
    return (x < y) ? x : y;
}

const double INF = 999999; // something relatively big

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // input already checked in m-file
    const mxArray *m_X = prhs[0], *m_Y = prhs[1];
    const mwSize nch = mxGetM(m_X), nx = mxGetN(m_X), ny = mxGetN(m_Y);
    const double *X = mxGetPr(m_X), *Y = mxGetPr(m_Y);

    mwSize pts = 0;
    double acc = 0;

    for (mwSize i = 0; i < nx; i++) {
        const double *x = X + nch*i;

        // determine closest element from same set (not self, unique)
        double dxx2 = INF; // something big
        for (mwSize j = 0; j < nx; j++) {
            if (i == j) continue; // skip self
            dxx2 = min(dxx2, dist2(x, X + nch*j, nch));
        }
        if (dxx2 == 0) continue; // two of equal distance

        // determine closest element from other set (unique)
        double dxy2 = INF; // something big
        for (mwSize j = 0; j < ny; j++) {
            dxy2 = min(dxy2, dist2(x, Y + nch*j, nch));
        }
        if (dxy2 == 0) continue; // two of equal distance

        // accumulate
        acc += log(dxy2 / dxx2);
        pts += 1;
    }

    plhs[0] = mxCreateDoubleScalar(acc / pts);
}
