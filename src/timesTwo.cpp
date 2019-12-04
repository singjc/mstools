#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

#include <cmath>
#include <climits>
#include <algorithm>
#include <vector>
#include "MSNumpress.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


static void encodeInt(
    const unsigned int x,
    unsigned char* res,
    size_t *res_length	
) {
  // get the bit pattern of a signed int x_inp
  unsigned int m;
  unsigned char i, l; // numbers between 0 and 9
  
  unsigned int mask = 0xf0000000;
  unsigned int init = x & mask;
  
  if (init == 0) {
    l = 8;
    for (i=0; i<8; i++) {
      m = mask >> (4*i);
      if ((x & m) != 0) {
        l = i;
        break;
      }
    }
    res[0] = l;
    for (i=l; i<8; i++) {
      res[1+i-l] = static_cast<unsigned char>( x >> (4*(i-l)) );
    }
    *res_length += 1+8-l;
    
  } else if (init == mask) {
    l = 7;
    for (i=0; i<8; i++) {
      m = mask >> (4*i);
      if ((x & m) != m) {
        l = i;
        break;
      }
    }
    res[0] = l + 8;
    for (i=l; i<8; i++) {
      res[1+i-l] = static_cast<unsigned char>( x >> (4*(i-l)) );
    }
    *res_length += 1+8-l;
    
  } else {
    res[0] = 0;
    for (i=0; i<8; i++) {
      res[1+i] = static_cast<unsigned char>( x >> (4*i) );
    }
    *res_length += 9;
    
  }
}

static void decodeInt(
    const unsigned char *data,
    size_t *di,
    size_t max_di,
    size_t *half,
    unsigned int *res
) {
  size_t n, i;
  unsigned int mask, m;
  unsigned char head;
  unsigned char hb;
  
  // Extract the first half byte, specifying the number of leading zero half
  // bytes of the final integer.
  // If half is zero, we look at the first half byte, otherwise we look at
  // the second (lower) half byte and advance the counter to the next char.
  if (*half == 0) {
    head = data[*di] >> 4;
  } else {
    head = data[*di] & 0xf;
    (*di)++;
  }
  
  *half = 1-(*half); // switch to other half byte
  *res = 0;
  
  if (head <= 8) {
    n = head;
  } else { // we have n leading ones, fill n half bytes in res with 0xf
    n = head - 8;
    mask = 0xf0000000;
    for (i=0; i<n; i++) {
      m = mask >> (4*i);
      *res = *res | m;
    }
  }
  
  if (n == 8) {
    return;
  }
  
  if (*di + ((8 - n) - (1 - *half)) / 2 >= max_di) {
    throw "[MSNumpress::decodeInt] Corrupt input data! ";
  }
  
  for (i=n; i<8; i++) {
    if (*half == 0) {
      hb = data[*di] >> 4;
    } else {
      hb = data[*di] & 0xf;
      (*di)++;
    }
    *res = *res | ( static_cast<unsigned int>(hb) << ((i-n)*4));
    *half = 1 - (*half);
  }
}

// This is only valid on systems were ints use more bytes than chars...

const int ONE = 1;
static bool is_little_endian() {
  return *((char*)&(ONE)) == 1;
}
bool IS_LITTLE_ENDIAN = is_little_endian();



/////////////////////////////////////////////////////////////

static void encodeFixedPoint(
    double fixedPoint, 
    unsigned char *result
) {
  int i;
  unsigned char *fp = (unsigned char*)&fixedPoint;
  for (i=0; i<8; i++) {
    result[i] = fp[IS_LITTLE_ENDIAN ? (7-i) : i];
  }
}



static double decodeFixedPoint(
    const unsigned char *data
) {
  int i;
  double fixedPoint;
  unsigned char *fp = (unsigned char*)&fixedPoint;
  
  for (i=0; i<8; i++) {
    fp[i] = data[IS_LITTLE_ENDIAN ? (7-i) : i];
  }
  
  return fixedPoint;
}


size_t encodePic(
    const double *data, 
    size_t dataSize, 
    unsigned char *result
) {
  size_t i, ri;
  unsigned int x;
  unsigned char halfBytes[10];
  size_t halfByteCount;
  size_t hbi;
  
  //printf("Encoding %d doubles\n", (int)dataSize);
  
  halfByteCount = 0;
  ri = 0;
  
  for (i=0; i<dataSize; i++) {
    
    if (THROW_ON_OVERFLOW && 
        (data[i] + 0.5 > INT_MAX || data[i] < -0.5)		){
      throw "[MSNumpress::encodePic] Cannot use Pic to encode a number larger than INT_MAX or smaller than 0.";
    }
    x = static_cast<unsigned int>(data[i] + 0.5);
    //printf("%d %d %d,   extrapol: %d    diff: %d \n", ints[0], ints[1], ints[2], extrapol, diff);
    encodeInt(x, &halfBytes[halfByteCount], &halfByteCount);
    
    for (hbi=1; hbi < halfByteCount; hbi+=2) {
      result[ri] = static_cast<unsigned char>(
        (halfBytes[hbi-1] << 4) | (halfBytes[hbi] & 0xf)
      );
      //printf("%x \n", result[ri]);
      ri++;
    }
    if (halfByteCount % 2 != 0) {
      halfBytes[0] = halfBytes[halfByteCount-1];
      halfByteCount = 1;
    } else {
      halfByteCount = 0;
    }
  }
  if (halfByteCount == 1) {
    result[ri] = static_cast<unsigned char>(halfBytes[0] << 4);
    ri++;
  }
  return ri;
}

//' @export
// [[Rcpp::export]]
void encodePic( const std::vector<double> &data,  
                std::vector<unsigned char> &result ) { 
  Rprintf( "encodePic" );
  for(int i=0; i<data.size(); ++i)
    std::cout << data[i] << ' ';
  
  size_t dataSize = data.size();
  result.resize(dataSize * 5);
  size_t encodedLength = encodePic(&data[0], dataSize, &result[0]);
  result.resize(encodedLength);
}

/////////////////////////////////////////////////////////////

double optimalLinearFixedPoint(
    const double *data, 
    size_t dataSize
) {
  /*
   * safer impl - apparently not needed though
   *
   if (dataSize == 0) return 0;
   
   double maxDouble = 0;
   double x;
   
   for (size_t i=0; i<dataSize; i++) {
   x = data[i];
   maxDouble = max(maxDouble, x);
   }
   
   return floor(0xFFFFFFFF / maxDouble);
   */
  if (dataSize == 0) return 0;
  if (dataSize == 1) return floor(0x7FFFFFFFl / data[0]);
  double maxDouble = max(data[0], data[1]);
  double extrapol;
  double diff;
  
  for (size_t i=2; i<dataSize; i++) {
    extrapol = data[i-1] + (data[i-1] - data[i-2]);
    diff = data[i] - extrapol;
    maxDouble = max(maxDouble, ceil(abs(diff)+1));
  }
  
  return floor(0x7FFFFFFFl / maxDouble);
}

double optimalLinearFixedPointMass(
    const double *data, 
    size_t dataSize,
    double mass_acc
) {
  if (dataSize < 3) return 0; // we just encode the first two points as floats
  
  // We calculate the maximal fixedPoint we need to achieve a specific mass
  // accuracy. Note that the maximal error we will make by encoding as int is
  // 0.5 due to rounding errors.
  double maxFP = 0.5 / mass_acc;
  
  // There is a maximal value for the FP given by the int length (32bit)
  // which means we cannot choose a value higher than that. In case we cannot
  // achieve the desired accuracy, return failure (-1).
  double maxFP_overflow = optimalLinearFixedPoint(data, dataSize);
  if (maxFP > maxFP_overflow) return -1;
  
  return maxFP;
}



size_t decodeLinear(
    const unsigned char *data,
    const size_t dataSize,
    double *result
) {
  size_t i;
  size_t ri = 0;
  unsigned int init, buff;
  int diff;
  long long ints[3];
  //double d;
  size_t di;
  size_t half;
  long long extrapol;
  long long y;
  double fixedPoint;
  
  //printf("Decoding %d bytes with fixed point %f\n", (int)dataSize, fixedPoint);
  
  if (dataSize == 8) return 0;
  
  if (dataSize < 8) 
    throw "[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read fixed point! ";
  
  fixedPoint = decodeFixedPoint(data);
  
  
  if (dataSize < 12) 
    throw "[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read first value! ";
  
  ints[1] = 0;
  for (i=0; i<4; i++) {
    ints[1] = ints[1] | ((0xff & (init = data[8+i])) << (i*8));
  }
  result[0] = ints[1] / fixedPoint;
  
  if (dataSize == 12) return 1;
  if (dataSize < 16) 
    throw "[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read second value! ";
  
  ints[2] = 0;
  for (i=0; i<4; i++) {
    ints[2] = ints[2] | ((0xff & (init = data[12+i])) << (i*8));
  }
  result[1] = ints[2] / fixedPoint;
  
  half = 0;
  ri = 2;
  di = 16;
  
  //printf("   di     ri      half    int[0]    int[1]    extrapol   diff\n");
  
  while (di < dataSize) {
    if (di == (dataSize - 1) && half == 1) {
      if ((data[di] & 0xf) == 0x0) {
        break;
      }
    }
    //printf("%7d %7d %7d %lu %lu %ld", di, ri, half, ints[0], ints[1], extrapol);
    
    ints[0] = ints[1];
    ints[1] = ints[2];
    decodeInt(data, &di, dataSize, &half, &buff);
    diff = static_cast<int>(buff);
    
    extrapol = ints[1] + (ints[1] - ints[0]);
    y = extrapol + diff;
    //printf(" %d \n", diff);
    result[ri++] 	= y / fixedPoint;
    ints[2] 		= y;
  }
  
  return ri;
}


//' @export
// [[Rcpp::export]]
void decodeLinear(
    const std::vector<unsigned char> &data,
    std::vector<double> &result
) {
  std::cout << "This is data: " << &data << '\n';
  std::cout << typeid(&data).name() << '\n';
  std::cout << "Data size: " << data.size() << '\n';
  std::cout << "This is result: " << &result << '\n';
  std::cout << typeid(&result).name() << '\n';
  std::cout << "Result size: " << result.size() << '\n';
  size_t dataSize = data.size();
  std::cout << "Line 407" << '\n';
  result.resize((dataSize - 8) * 2);
  std::cout << "Line 409" << '\n';
  std::cout << "Result size2: " << result.size() << '\n';
  std::cout << "Line 411" << '\n';
  size_t decodedLength = decodeLinear(&data[0], dataSize, &result[0]);
  std::cout << "Line 413" << '\n';
  result.resize(decodedLength);
  std::cout << "Result size3: " << result.size() << '\n';
  std::cout << "This is Result: " << &result << '\n';
}

/////////////////////////////////////////////////////////////


