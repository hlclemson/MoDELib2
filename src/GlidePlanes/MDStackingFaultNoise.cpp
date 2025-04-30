/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MDStackingFaultNoise_cpp
#define model_MDStackingFaultNoise_cpp

#include <MDStackingFaultNoise.h>
//#include <thread>

namespace model
{

MDStackingFaultNoise::MDStackingFaultNoise(const PolycrystallineMaterialBase& mat,
                                           const std::string &tag,
                                           const std::string &correlationFile_in,
                                           const int &seed,
                                           const GridSizeType &gridSize,
                                           const GridSpacingType &gridSpacing
                                           ) :
  /* init */ GlidePlaneNoiseBase<1>("MDStackingFaultNoise"+tag,seed,gridSize,gridSpacing)
  /* init */,correlationFile(correlationFile_in)
{
  // allocate correlation array with zero padding in real space
  REAL_SCALAR *Rr = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);

  // allocate correlation array in fourier space
  //COMPLEX *Rk = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
  Rk = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);

  //fftw_plan plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);
  // fft plans
  //fftw_plan plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);
  fftw_plan plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);


  // read the dimension of the original correlation
  const auto originalDimensions(readVTKfileDimension(correlationFile.c_str()));
  const int originalNX = originalDimensions(0);
  const int originalNY = originalDimensions(1);
  if(originalDimensions(2)!=1)
  {
      throw std::runtime_error("vtk correlationFiles 'DIMENSIONS' should have 3rd component == 1.");
  }

  //std::cout << "original NX = " << originalNX << std::endl;
  //std::cout << "origianl NY = " << originalNY << std::endl;

  REAL_SCALAR *Rr_original = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*originalNX*originalNY); // correlation in real space

  // populate Rr_original with the correlation data
  StackingFaultCorrelationReader(correlationFile, Rr_original);

  std::cout << "read sf correlation file." << std::endl;

  // initialize with zeros
  for (int i = 0; i < this->NY; ++i)
  {
      for (int j = 0; j < this->NX; ++j)
      {
          Rr[i*this->NX + j] = 0;
      }
  }

  int start_y = (NY - originalNX) / 2;
  int start_x = (NX - originalNY) / 2;

  // 0-pading from centere
  for (int i = 0; i < originalNY; ++i)
  {
    for (int j = 0; j < originalNX; ++j)
    {
      Rr[(start_y + i) * this->NX + (start_x + j)] = Rr_original[i * originalNX + j];
    }
  }


  std::cout << "padded with zeros" << std::endl;

  // unit conversion (from J^2/m^4 to unitless)
  for (int i = 0; i < NR; ++i)
  {
    //Rr[i] /= (b_SI*mu_SI);
    const double unitconvert = mat.b_SI*mat.mu_SI;
    Rr[i] /= (unitconvert*unitconvert);
  }

  std::cout << "convert unit" << std::endl;

  //plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr_original, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);
  //plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);
  //fftw_plan plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);

  // Execute FFTW plans to populate Rk_xz and Rk_yz
  fftw_execute(plan_R_r2c);

  std::cout << "fftw" << std::endl;

  // Normalize the FFT output
  for (int i = 0; i < this->NX; ++i)
  {
      for (int j = 0; j < (this->NY/2 + 1); ++j)
      {
          Rk[i * (this->NY/2 + 1) + j] /= (this->NX * this->NY);
      }
  }

  std::cout << "normalize fftw output" << std::endl;

  std::cout << "Rk[0] = " << Rk[0] << std::endl;

  //for (int i=0; i<NK; ++i) 
  //{
  //  std::cout << "Rk[i] = " << Rk[i] << std::endl;
  //}

  // correct the standard deviation (this should be done after the noise is sampled and ifft back to real space)
  //const double NRO = static_cast<double>(originalNX)*static_cast<double>(originalNY);
  //for (int i=0; i<NR; ++i) 
  //{
  //    fr[i] = std::sqrt((NRO+(static_cast<double>(NR)-NRO))/NRO) * fr[i];
  //}

  // Destroy FFTW plans
  fftw_destroy_plan(plan_R_r2c);

  // Free allocated memory
  fftw_free(Rr_original);

  // Free allocated memory
  fftw_free(Rr);
}


std::array<MDStackingFaultNoise::COMPLEX,1> MDStackingFaultNoise::kCorrelations(const Eigen::Matrix<double, 3, 1> &kv, const Eigen::Matrix<int, 3, 1> &index) const
{
    std::array<MDStackingFaultNoise::COMPLEX,1> temp;
    std::cout << "temp " << std::endl;
    int idx=(this->NY/2+1)*NZ*index(0) + index(1)*NZ + index(2);
    std::cout << "idx = " << idx << std::endl;
    std::cout << "trying to grep Rk[idx]" << std::endl;
    std::cout << "Rk[idx]" << Rk[0] << std::endl;
    temp[0] = Rk[idx];
    std::cout << "correlation " << std::endl;
    return temp;
}

//int MDStackingFaultNoise::testVec() const
//Eigen::Matrix<double,2,2> MDStackingFaultNoise::nonOrthogonalBasis()
Eigen::Matrix<double,2,2> MDStackingFaultNoise::invTransitionMatrix() 
{
    Eigen::Matrix<double,2,2> transitionMatrix = nonOrthogonalBasisReader(correlationFile);
    // Check for singularity before inversion
    double det = transitionMatrix.determinant();
    // hardcoded stability condition, change this later
    if (std::abs(det) < 1e-10) {
        throw std::runtime_error("non-orthogonal basis matrix is singular.");
    }
    // return the inverse of transition matrix
    //return nonOrthoBasisMatrix.transpose().inverse();
    return transitionMatrix.inverse();

    //const Eigen::Matrix<double,2,2> basis = nonOrthogonalBasisReader(correlationFile);
    //return nonOrthogonalBasisReader(correlationFile);
}

Eigen::Matrix<double,2,2> MDStackingFaultNoise::nonOrthogonalBasisReader(const std::string& fileName_vtk) const
{
    typedef Eigen::Matrix<double,2,1> VectorDimD;
    std::deque<VectorDimD> basis1;
    std::deque<VectorDimD> basis2;
    std::ifstream vtkFile(fileName_vtk); //access vtk file
    // error check
    if (!vtkFile.is_open()) {
        throw std::runtime_error("Error opening stacking fault VTK correlation file!");
    }
    // begin parsing structured_grid vtk file for lines
    std::string line;
    while (std::getline(vtkFile, line)) 
    {
        if(line.find("VECTORS")!=std::string::npos) 
        {
            std::getline(vtkFile, line);
            std::stringstream ss(line);
            double x, y, z;
            if(basis1.empty())
            {
                ss >> x >> y >> z;
                VectorDimD vec = (VectorDimD() << x, y).finished();
                vec.normalize(); // Normalize the vector
                basis1.push_back(vec);
            }
            else
            {
                ss >> x >> y >> z;
                VectorDimD vec = (VectorDimD() << x, y).finished();
                vec.normalize(); // Normalize the vector
                basis2.push_back(vec);
            }

        }
    }

    Eigen::Matrix<double,2,2> nonOrthoBasisMatrix;
    // fill the nonOrthoBasisMatrix with a block that has the size of <2,1> at the beginning of the first column
    nonOrthoBasisMatrix.block<2,1>(0,0) = basis1.front();
    // fill the nonOrthoBasisMatrix with a block that has the size of <2,1> at the beginning of the second column
    nonOrthoBasisMatrix.block<2,1>(0,1) = basis2.front();

    // transpose the matrix to create a transition matrix
    return nonOrthoBasisMatrix.transpose();
}

void MDStackingFaultNoise::StackingFaultCorrelationReader(const std::string &correlationFile, REAL_SCALAR *Rr)
{
  std::cout << "Reading stacking fault correlation" << std::endl;

  std::ifstream vtkFile(correlationFile); // access vtk file
  // error check
  if (!vtkFile.is_open())
  {
    throw std::runtime_error("Error opening stacking fault VTK correlation file!");
  }

  std::string line;
  while (std::getline(vtkFile, line))
  {
    // if the "POINT_DATA" string is read, read the following data
    if (line.find("POINT_DATA") != std::string::npos)
    {
      const size_t numOfHeaders = 2;
      // get the number of points in the file
      const size_t firstSpace(line.find(' '));
      const size_t numOfPoints = std::atoi(line.substr(firstSpace+1).c_str());
      // read the point coordinates
      for (size_t n = 0; n < numOfPoints + numOfHeaders; ++n)
      {
        std::getline(vtkFile, line);
        // ignore the headers that come right after point_data
        if (n < numOfHeaders)
          continue;
        const int ind = n - numOfHeaders;
        // correlationCoeffs.push_back(std::atoi(line.c_str()));
        if (ind >= NR) 
        {
          throw std::runtime_error("Index out of bounds while populating the original correlation array.");
        }

        try
        {
          double value = std::stod(line);
          Rr[ind] = value;
        }
        catch(const std::invalid_argument& e) 
        {
          std::cerr << "Invalid correlation data in line: " << line << std::endl;
        }
        catch(const std::out_of_range& e)
        {
          std::cerr << "Out of range correlation value in line: " << line << std::endl;
        }
        //Rr[ind] = std::atof(line.c_str());
      }
    }
  }
  vtkFile.close(); // Close the file after reading
}

typename MDStackingFaultNoise::GridSizeType MDStackingFaultNoise::readVTKfileDimension(const char *fname)
{
    FILE *InFile=fopen(fname,"r");

    if (InFile == NULL)
    {
        fprintf(stderr, "Can't open stacking fault correlation VTK file %s\n",fname);
        exit(1);
    }
    // return the 5th line of the vtk file
    char line[200];
    for(int i=0;i<4;i++)
    {
        fgets(line, 200, InFile);
    }
    // scan the returned line
    int NXX, NYY, NZZ;
    fscanf(InFile, "%s %d %d %d\n", line, &(NXX), &(NYY), &(NZZ));
    return (GridSizeType()<<NXX,NYY,NZZ).finished();
}

}
#endif
