#include "mex.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

double compute_border_median(const float* img, mwSize H, mwSize W, mwIndex x0, mwIndex y0, mwIndex R) {
    const mwSize win = 2 * R + 1;
    const mwSize nBorder = 4 * win - 4;
    std::vector<double> border;
    border.reserve(nBorder);

    const mwIndex xMin = x0 - R;
    const mwIndex xMax = x0 + R;
    const mwIndex yMin = y0 - R;
    const mwIndex yMax = y0 + R;

    for (mwIndex x = xMin; x <= xMax; ++x) {
        border.push_back(static_cast<double>(img[yMin + x * H]));
    }
    for (mwIndex x = xMin; x <= xMax; ++x) {
        border.push_back(static_cast<double>(img[yMax + x * H]));
    }
    for (mwIndex y = yMin + 1; y < yMax; ++y) {
        border.push_back(static_cast<double>(img[y + xMin * H]));
    }
    for (mwIndex y = yMin + 1; y < yMax; ++y) {
        border.push_back(static_cast<double>(img[y + xMax * H]));
    }

    const mwSize n = border.size();
    const mwSize mid = n / 2;
    std::nth_element(border.begin(), border.begin() + mid, border.end());
    const double upper = border[mid];
    if ((n & 1U) != 0U) {
        return upper;
    }
    std::nth_element(border.begin(), border.begin() + mid - 1, border.begin() + mid);
    const double lower = border[mid - 1];
    return 0.5 * (lower + upper);
}

double compute_border_median(const double* img, mwSize H, mwSize W, mwIndex x0, mwIndex y0, mwIndex R) {
    const mwSize win = 2 * R + 1;
    const mwSize nBorder = 4 * win - 4;
    std::vector<double> border;
    border.reserve(nBorder);

    const mwIndex xMin = x0 - R;
    const mwIndex xMax = x0 + R;
    const mwIndex yMin = y0 - R;
    const mwIndex yMax = y0 + R;

    for (mwIndex x = xMin; x <= xMax; ++x) {
        border.push_back(img[yMin + x * H]);
    }
    for (mwIndex x = xMin; x <= xMax; ++x) {
        border.push_back(img[yMax + x * H]);
    }
    for (mwIndex y = yMin + 1; y < yMax; ++y) {
        border.push_back(img[y + xMin * H]);
    }
    for (mwIndex y = yMin + 1; y < yMax; ++y) {
        border.push_back(img[y + xMax * H]);
    }

    const mwSize n = border.size();
    const mwSize mid = n / 2;
    std::nth_element(border.begin(), border.begin() + mid, border.end());
    const double upper = border[mid];
    if ((n & 1U) != 0U) {
        return upper;
    }
    std::nth_element(border.begin(), border.begin() + mid - 1, border.begin() + mid);
    const double lower = border[mid - 1];
    return 0.5 * (lower + upper);
}

template <typename T>
void localize_kernel(const T* img,
                     mwSize H,
                     mwSize W,
                     const double* candXY,
                     mwSize nCand,
                     mwIndex R,
                     T* xyOut,
                     double& nValidOut) {
    mwIndex nValid = 0;

    for (mwIndex i = 0; i < nCand; ++i) {
        const double xCand = candXY[i];
        const double yCand = candXY[i + nCand];
        const mwIndex x0 = static_cast<mwIndex>(xCand) - 1;
        const mwIndex y0 = static_cast<mwIndex>(yCand) - 1;

        if (x0 < R || x0 + R >= W || y0 < R || y0 + R >= H) {
            continue;
        }

        const double bg = compute_border_median(img, H, W, x0, y0, R);

        double sumW = 0.0;
        double sumDx = 0.0;
        double sumDy = 0.0;

        for (mwIndex dx = 0; dx <= 2 * R; ++dx) {
            const mwIndex x = x0 + dx - R;
            const double relX = static_cast<double>(dx) - static_cast<double>(R);
            for (mwIndex dy = 0; dy <= 2 * R; ++dy) {
                const mwIndex y = y0 + dy - R;
                const double relY = static_cast<double>(dy) - static_cast<double>(R);
                const double w = std::max(static_cast<double>(img[y + x * H]) - bg, 0.0);
                sumW += w;
                sumDx += relX * w;
                sumDy += relY * w;
            }
        }

        if (sumW <= 0.0) {
            continue;
        }

        const double xhat = xCand + sumDx / sumW;
        const double yhat = yCand + sumDy / sumW;
        xyOut[nValid] = static_cast<T>(xhat);
        xyOut[nValid + nCand] = static_cast<T>(yhat);
        ++nValid;
    }

    nValidOut = static_cast<double>(nValid);
}

} // namespace

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("localize_weighted_centroid_mex:nrhs",
                          "Expected 3 inputs: I, candXY, R.");
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt("localize_weighted_centroid_mex:nlhs",
                          "Expected at most 2 outputs: xyAll, nValid.");
    }

    const mxArray* Iin = prhs[0];
    const mxArray* candIn = prhs[1];
    const mxArray* Rin = prhs[2];

    if (mxGetNumberOfDimensions(Iin) != 2) {
        mexErrMsgIdAndTxt("localize_weighted_centroid_mex:image",
                          "Input image I must be 2-D.");
    }
    if (!mxIsSingle(Iin) && !mxIsDouble(Iin)) {
        mexErrMsgIdAndTxt("localize_weighted_centroid_mex:imageType",
                          "Input image I must be single or double.");
    }
    if (!mxIsDouble(candIn) || mxGetNumberOfDimensions(candIn) != 2 || mxGetN(candIn) != 2) {
        mexErrMsgIdAndTxt("localize_weighted_centroid_mex:candXY",
                          "candXY must be an n-by-2 double matrix.");
    }
    if (!mxIsDouble(Rin) || mxGetNumberOfElements(Rin) != 1) {
        mexErrMsgIdAndTxt("localize_weighted_centroid_mex:R",
                          "R must be a double scalar.");
    }

    const mwSize H = mxGetM(Iin);
    const mwSize W = mxGetN(Iin);
    const mwSize nCand = mxGetM(candIn);
    const mwIndex R = static_cast<mwIndex>(mxGetScalar(Rin));
    const double* candXY = mxGetPr(candIn);

    mxClassID classId = mxGetClassID(Iin);
    plhs[0] = mxCreateNumericMatrix(nCand, 2, classId, mxREAL);
    plhs[1] = mxCreateDoubleScalar(0.0);

    double nValid = 0.0;
    if (classId == mxSINGLE_CLASS) {
        localize_kernel(static_cast<const float*>(mxGetData(Iin)),
                        H, W, candXY, nCand, R,
                        static_cast<float*>(mxGetData(plhs[0])),
                        nValid);
    } else {
        localize_kernel(static_cast<const double*>(mxGetData(Iin)),
                        H, W, candXY, nCand, R,
                        static_cast<double*>(mxGetData(plhs[0])),
                        nValid);
    }

    *mxGetPr(plhs[1]) = nValid;
}
