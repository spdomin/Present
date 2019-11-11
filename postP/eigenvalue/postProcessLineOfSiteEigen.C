#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

class EigenValue {
public: 
  EigenValue() { w_ = 12;}
  ~EigenValue() {}

  void diagonalize( 
    const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3]);
  void output_matrix(
    const double (&A)[3][3], std::string message);
  void check_matrix_vector(
    const double (&D)[3][3], const double (&Q)[3][3]);
  void matrix_matrix_multiply(
    const double (&A)[3][3], const double (&B)[3][3], double (&C)[3][3]);
  void sort(const double (&D)[3][3], int (&row)[3]);

  int w_;
};

void
EigenValue::output_matrix(
  const double (&A)[3][3], std::string message)
{
  std::cout << message << std::setw(w_) << std::endl;
  for ( int i = 0; i < 3; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      std::cout << A[i][j] << " " << std::setw(w_); 
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void
EigenValue::matrix_matrix_multiply(
  const double (&A)[3][3], const double (&B)[3][3], double (&C)[3][3])
{
  //C = A*B
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double sum = 0;
      for (int k = 0; k < 3; ++k) {
        sum = sum + A[i][k] * B[k][j];
      }
      C[i][j] = sum;
    }
  }
}

void
EigenValue::check_matrix_vector(
  const double (&D)[3][3], const double (&Q)[3][3])
{
  // A = Q*D*QT
  double QT[3][3];
  double B[3][3];
  double A[3][3];

  // compute QT
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      QT[j][i] = Q[i][j];	
    }
  }

  output_matrix(QT, "Summary for QT");

  //mat-vec, B = Q*D
  matrix_matrix_multiply(Q,D,B);

  // mat-vec, A = (Q*D)*QT = B*QT
  matrix_matrix_multiply(B,QT,A);

  output_matrix(A, "Check matrix-matrix multiply based on A = Q*D*QT");
}

void 
EigenValue::sort(const double (&D)[3][3], int (&row)[3]) {

  // fill in data to sort
  double data[3] = {D[0][0], D[1][1], D[2][2]};
  row[0] = 0;
  row[1] = 1;
  row[2] = 2;

  int j = 0;
  double tmp = 0;
  for(int i=0; i<3; ++i){
    j = i;
    for(int k = i; k<3; ++k){
      if(data[j]<data[k]){
        j = k;
      }
    }
    tmp = data[i];
    int tmpI = row[i];
    // deal with data
    data[i] = data[j];
    data[j] = tmp;
    // now row mapping
    row[i] = row[j];
    row[j] = tmpI;
  }
}

void 
EigenValue::diagonalize(
  const double (&A)[3][3], double (&Q)[3][3], double (&D)[3][3])
{
  /*
    obtained from: 
    http://stackoverflow.com/questions/4372224/
    fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition

    A must be a symmetric matrix.
    returns Q and D such that 
    Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
  */

  const int maxsteps=24;
  int k0, k1, k2;
  double o[3], m[3];
  double q [4] = {0.0,0.0,0.0,1.0};
  double jr[4];
  double sqw, sqx, sqy, sqz;
  double tmp1, tmp2, mq;
  double AQ[3][3];
  double thet, sgn, t, c;
  for(int i=0;i < maxsteps;++i) {
    // quat to matrix
    sqx      = q[0]*q[0];
    sqy      = q[1]*q[1];
    sqz      = q[2]*q[2];
    sqw      = q[3]*q[3];
    Q[0][0]  = ( sqx - sqy - sqz + sqw);
    Q[1][1]  = (-sqx + sqy - sqz + sqw);
    Q[2][2]  = (-sqx - sqy + sqz + sqw);
    tmp1     = q[0]*q[1];
    tmp2     = q[2]*q[3];
    Q[1][0]  = 2.0 * (tmp1 + tmp2);
    Q[0][1]  = 2.0 * (tmp1 - tmp2);
    tmp1     = q[0]*q[2];
    tmp2     = q[1]*q[3];
    Q[2][0]  = 2.0 * (tmp1 - tmp2);
    Q[0][2]  = 2.0 * (tmp1 + tmp2);
    tmp1     = q[1]*q[2];
    tmp2     = q[0]*q[3];
    Q[2][1]  = 2.0 * (tmp1 + tmp2);
    Q[1][2]  = 2.0 * (tmp1 - tmp2);
    
    // AQ = A * Q
    AQ[0][0] = Q[0][0]*A[0][0]+Q[1][0]*A[0][1]+Q[2][0]*A[0][2];
    AQ[0][1] = Q[0][1]*A[0][0]+Q[1][1]*A[0][1]+Q[2][1]*A[0][2];
    AQ[0][2] = Q[0][2]*A[0][0]+Q[1][2]*A[0][1]+Q[2][2]*A[0][2];
    AQ[1][0] = Q[0][0]*A[0][1]+Q[1][0]*A[1][1]+Q[2][0]*A[1][2];
    AQ[1][1] = Q[0][1]*A[0][1]+Q[1][1]*A[1][1]+Q[2][1]*A[1][2];
    AQ[1][2] = Q[0][2]*A[0][1]+Q[1][2]*A[1][1]+Q[2][2]*A[1][2];
    AQ[2][0] = Q[0][0]*A[0][2]+Q[1][0]*A[1][2]+Q[2][0]*A[2][2];
    AQ[2][1] = Q[0][1]*A[0][2]+Q[1][1]*A[1][2]+Q[2][1]*A[2][2];
    AQ[2][2] = Q[0][2]*A[0][2]+Q[1][2]*A[1][2]+Q[2][2]*A[2][2];
    // D = Qt * AQ
    D[0][0] = AQ[0][0]*Q[0][0]+AQ[1][0]*Q[1][0]+AQ[2][0]*Q[2][0]; 
    D[0][1] = AQ[0][0]*Q[0][1]+AQ[1][0]*Q[1][1]+AQ[2][0]*Q[2][1]; 
    D[0][2] = AQ[0][0]*Q[0][2]+AQ[1][0]*Q[1][2]+AQ[2][0]*Q[2][2]; 
    D[1][0] = AQ[0][1]*Q[0][0]+AQ[1][1]*Q[1][0]+AQ[2][1]*Q[2][0]; 
    D[1][1] = AQ[0][1]*Q[0][1]+AQ[1][1]*Q[1][1]+AQ[2][1]*Q[2][1]; 
    D[1][2] = AQ[0][1]*Q[0][2]+AQ[1][1]*Q[1][2]+AQ[2][1]*Q[2][2]; 
    D[2][0] = AQ[0][2]*Q[0][0]+AQ[1][2]*Q[1][0]+AQ[2][2]*Q[2][0]; 
    D[2][1] = AQ[0][2]*Q[0][1]+AQ[1][2]*Q[1][1]+AQ[2][2]*Q[2][1]; 
    D[2][2] = AQ[0][2]*Q[0][2]+AQ[1][2]*Q[1][2]+AQ[2][2]*Q[2][2];
    o[0]    = D[1][2];
    o[1]    = D[0][2];
    o[2]    = D[0][1];
    m[0]    = std::abs(o[0]);
    m[1]    = std::abs(o[1]);
    m[2]    = std::abs(o[2]);
    
    k0      = (m[0] > m[1] && m[0] > m[2])?0: (m[1] > m[2])? 1 : 2; // index of largest element of offdiag
    k1      = (k0+1)%3;
    k2      = (k0+2)%3;
    if (o[k0]==0.0) {
      break;  // diagonal already
    }
    thet    = (D[k2][k2]-D[k1][k1])/(2.0*o[k0]);
    sgn     = (thet > 0.0)?1.0:-1.0;
    thet   *= sgn; // make it positive
    t       = sgn /(thet +((thet < 1.E6)? std::sqrt(thet*thet+1.0):thet)) ; // sign(T)/(|T|+sqrt(T^2+1))
    c       = 1.0/std::sqrt(t*t+1.0); //  c= 1/(t^2+1) , t=s/c 
    if(c==1.0) {
      break;  // no room for improvement - reached machine precision.
    }
    jr[0 ]  = jr[1] = jr[2] = jr[3] = 0.0;
    jr[k0]  = sgn*std::sqrt((1.0-c)/2.0);  // using 1/2 angle identity sin(a/2) = std::sqrt((1-cos(a))/2)  
    jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
    jr[3 ]  = std::sqrt(1.0f - jr[k0] * jr[k0]);
    if(jr[3]==1.0) {
      break; // reached limits of floating point precision
    }
    q[0]    = (q[3]*jr[0] + q[0]*jr[3] + q[1]*jr[2] - q[2]*jr[1]);
    q[1]    = (q[3]*jr[1] - q[0]*jr[2] + q[1]*jr[3] + q[2]*jr[0]);
    q[2]    = (q[3]*jr[2] + q[0]*jr[1] - q[1]*jr[0] + q[2]*jr[3]);
    q[3]    = (q[3]*jr[3] - q[0]*jr[0] - q[1]*jr[1] - q[2]*jr[2]);
    mq      = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    q[0]   /= mq;
    q[1]   /= mq;
    q[2]   /= mq;
    q[3]   /= mq;
  }
}

int main()
{
  EigenValue *ev = new EigenValue();

  double A[3][3] = {};
  double Q[3][3] = {};
  double D[3][3] = {};
  int row[3] = {};
   
  const double sqrtThreeOverTwo = std::sqrt(3.0)/2.0;

  const int numProbes = 100;
  std::vector<std::vector<double> >Data;
  
  std::ifstream inFile;
  inFile.open("lineOfSite.dat");

  std::ofstream outFile;
  outFile.open("barycentric.dat");
  outFile << "# cx xB yB uP vP wP ";

  std::string line;
  std::getline(inFile, line);
  std::cout << line << std::endl;

  // Time cx cy cz ux uy uz rs0 rs1 rs2 rs3 rs4 rs5  
  double Time, cx, cy, cz, ux, uy, uz, rs0, rs1, rs2, rs3, rs4, rs5;
  while ( !inFile.eof() ) {
    inFile >> Time >> cx >> cy >> cz 
              >> ux >> uy >> uz 
              >> rs0 >> rs1 >> rs2 >> rs3 >> rs4 >> rs5;
    std::vector<double> tmpData;
    tmpData.push_back(Time);  // 0
    tmpData.push_back(cx);    // 1
    tmpData.push_back(cy);    // 2
    tmpData.push_back(cz);    // 3
    tmpData.push_back(ux);    // 4
    tmpData.push_back(uy);    // 5
    tmpData.push_back(uz);    // 6
    tmpData.push_back(rs0);   // 7
    tmpData.push_back(rs1);   // 8
    tmpData.push_back(rs2);   // 9
    tmpData.push_back(rs3);   // 10
    tmpData.push_back(rs4);   // 11
    tmpData.push_back(rs5);   // 12
    Data.push_back(tmpData);
  }

  const size_t totalSize = Data.size();
  const size_t startIndex = totalSize - numProbes;
  for ( size_t k = startIndex; k < totalSize; ++k ) {

    std::vector<double> tmpData = Data[k];

    const double  Cx = tmpData[1];
    const double  Cy = tmpData[1];
    const double  Cz = tmpData[1];
    const double  ra0 = tmpData[7];
    const double  ra1 = tmpData[8];
    const double  ra2 = tmpData[9];
    const double  ra3 = tmpData[10];
    const double  ra4 = tmpData[11];
    const double  ra5 = tmpData[12];

    A[0][0] = ra0;
    A[0][1] = ra1;
    A[0][2] = ra2;
   
    A[1][0] = ra1;
    A[1][1] = ra3;
    A[1][2] = ra4;
   
    A[2][0] = ra2;
    A[2][1] = ra4;
    A[2][2] = ra5;

    // tke is the 1/2*(Rxx + Ryy + Rzz)
    const double tke = 0.5*(A[0][0] + A[1][1] + A[2][2]);
    
    // scale
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        const double fac = ( i == j ) ? 1.0/3.0 : 0.0;
        A[i][j] = A[i][j]/(2.0*tke) - fac;
      }
    }    

    // eigenvector
    double Q[3][3];
    
    // diagonal, or eigenvalue
    double D[3][3];
    
    // factorize it
    ev->diagonalize(A, Q, D);
    
    // sort eigenvalues
    ev->sort(D, row);
    
    int r0 = row[0];
    int r1 = row[1];
    int r2 = row[2];

    const double lambdaOne = D[r0][r0];
    const double lambdaTwo = D[r1][r1];
    const double lambdaThree = D[r2][r2];

    const double Cone = lambdaOne - lambdaTwo;
    const double Ctwo = 2.0*(lambdaTwo - lambdaThree);
    const double Cthree = 3.0*lambdaThree + 1.0;
    
    const double xB = Cone + Cthree*0.5;
    const double yB = Cthree*sqrtThreeOverTwo;

    std::cout << Cx << " "  
              << ra0 << " " << ra1 << " " << ra2 << " " << ra3 << " " << ra4 << " " << ra5 << " "
              << lambdaOne << " " << lambdaTwo << " " << lambdaThree << " " 
              << lambdaOne + lambdaTwo + lambdaThree << " " 
              << xB << " " << yB << std::endl;
    
    outFile << Cx << " " << xB << " " << yB << " " << std::sqrt(ra0) << " " << std::sqrt(ra3) << std::sqrt(ra5) << std::endl;
  }
  
  inFile.close();
  outFile.close();

  // delete it
  delete ev;

  return(0);
}

