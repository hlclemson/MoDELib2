/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MDStackingFaultNoise_cpp
#define model_MDStackingFaultNoise_cpp

#include <MDStackingFaultNoise.h>
#include <cmath>
//#include <thread>

namespace model
{

MDStackingFaultNoise::MDStackingFaultNoise(const PolycrystallineMaterialBase& mat,
                                           const std::string& tag,
                                           const std::string& correlationFile_in,
                                           const int& outputNoise,
                                           const std::string& noiseFile,
                                           const int& testNoiseSampling,
                                           const int& seed,
                                           const GridSizeType& gridSize,
                                           const GridSpacingType& gridSpacing
                                           ) :
  /* init */ GlidePlaneNoiseBase<1>("MDStackingFaultNoise"+tag,seed,gridSize,gridSpacing)
  /* init */,correlationFile(correlationFile_in)
{
  // read the dimension of the original correlation
  const auto originalDimensions(readVTKfileDimension(correlationFile.c_str()));
  const int originalNX = originalDimensions(0);
  const int originalNY = originalDimensions(1);
  if(originalDimensions(2)!=1)
  {
      throw std::runtime_error("vtk correlationFiles 'DIMENSIONS' should have 3rd component == 1.");
  }

  REAL_SCALAR *Rr_original = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*originalNX*originalNY); // correlation in real space

  // populate Rr_original with the correlation data
  StackingFaultCorrelationReader(correlationFile, Rr_original, this->NR);

  // unit conversion (from J^2/m^4 to unitless)
  for (int i = 0; i < originalNX*originalNY; ++i)
  {
    const REAL_SCALAR unitconvert = mat.b_SI*mat.mu_SI;
    Rr_original[i] /= (unitconvert*unitconvert);
  }

  // allocate correlation array with zero padding in real space
  REAL_SCALAR *Rr = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
  // explicitly initialize with zeros
  std::memset(Rr, 0, sizeof(REAL_SCALAR) * this->NR);

  // allocate correlation array in fourier space
  Rk = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
  // explicitly initialize the array to zero
  std::memset(Rk, 0, sizeof(COMPLEX) * this->NK);

  // if the assigned grid size is the same as the actual data, just copy the data
  if(this->NR == originalNX*originalNY)
  {
    for (int i = 0; i < this->NY; ++i)
    {
      for (int j = 0; j < this->NX; ++j)
      {
        Rr[i*this->NX + j] = Rr_original[i*this->NX + j];
      }
    }
  }
  // pad with zeros if the assigned grid size is larger than the actual data grid size
  else
  {
    int start_y = (this->NY - originalNY) / 2;
    int start_x = (this->NX - originalNX) / 2;
    // 0-pading from centere
    for (int i = 0; i < originalNY; ++i)
    {
      for (int j = 0; j < originalNX; ++j)
      {
        Rr[(start_y + i) * this->NX + (start_x + j)] = Rr_original[i * originalNX + j];
      }
    }
  }

  // Execute FFTW plans to populate Rk
  fftw_plan plan_R_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr, reinterpret_cast<fftw_complex*>(Rk), FFTW_ESTIMATE);
  fftw_execute(plan_R_r2c);

  // Normalize the FFT output
  for (int i = 0; i < this->NY; ++i)
  {
      for (int j = 0; j < (this->NX/2 + 1); ++j)
      {
          Rk[i * (this->NX/2 + 1) + j] /= static_cast<double>(this->NR);
      }
  }

  // noise output for debugging
  if (outputNoise) { Write_field_slice(mat, seed, noiseFile.c_str()); };
  if (testNoiseSampling) { 
    sampleNoiseRepeatedly(mat, 10);
    sampleNoiseRepeatedly(mat, 100);
    sampleNoiseRepeatedly(mat, 1000);
  };

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
    int idx=(this->NY/2+1)*NZ*index(0) + index(1)*NZ + index(2);
    temp[0] = Rk[idx];
    return temp;
}

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
    return transitionMatrix.inverse();
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

void MDStackingFaultNoise::StackingFaultCorrelationReader(const std::string &correlationFile, REAL_SCALAR *Rr, int NR)
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
          //float value = std::atof(line.c_str());
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

void MDStackingFaultNoise::sampleNoiseLocalInKspace(const PolycrystallineMaterialBase& mat, const int& localSeed, COMPLEX* localRk)
{
  const int J_MAX = this->NY/2 + 1;
  const int K_MAX = 1;

  std::default_random_engine generator(localSeed);
  std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);
  for(int i=0; i<this->NX; i++)
  {
      for(int j=0; j<J_MAX; j++)
      {
          for(int k=0; k<K_MAX; k++)
          {
              const int ind = J_MAX*K_MAX*i + j*K_MAX + k;

              REAL_SCALAR kx = 2.*M_PI/this->LX*REAL_SCALAR(i);
              if(i>this->NX/2)
              {
                  kx = 2.*M_PI/this->LX*REAL_SCALAR(i-this->NX);
              }

              REAL_SCALAR ky = 2*M_PI/this->LY*REAL_SCALAR(j);
              if(j>this->NY/2)
              {
                  ky = 2.*M_PI/this->LY*REAL_SCALAR(j-this->NY);
              }

              REAL_SCALAR kz = 2.*M_PI/this->LZ*REAL_SCALAR(k);

              // random numbers
              REAL_SCALAR Nk = distribution(generator);
              REAL_SCALAR Mk = distribution(generator);

              // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform 
              const double kCorrFactor((j==0 || j==this->NY/2)? 1.0 : 2.0);
              //kNoisyCorrelations[ind]=std::sqrt(Rk[ind]/kCorrFactor)*(Nk+Mk*COMPLEX(0.0,1.0));
              localRk[ind]=std::sqrt(Rk[ind]/kCorrFactor)*(Nk+Mk*COMPLEX(0.0,1.0));
          }
      }
  }
}

// shift the correlation to the center
void MDStackingFaultNoise::circularShift(REAL_SCALAR* rCorrelation)
{
  REAL_SCALAR* temp = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
  const int shiftY = this->NY/2;
  const int shiftX = this->NX/2;
  for (int y = 0; y < this->NY; ++y)
  {
    int newY = (y + shiftY) % this->NY;
    for (int x = 0; x < this->NX; ++x)
    {
      int newX = (x + shiftX) % this->NX;
      temp[newY*this->NX + newX] = rCorrelation[y*this->NX + x];
    }
  }
  // overwrite
  for (int i = 0; i < this->NR; ++i)
  {
    rCorrelation[i] = temp[i];
  }
  fftw_free(temp);
}

/*
*/
void MDStackingFaultNoise::sampleNoiseRepeatedly(const PolycrystallineMaterialBase& mat, const int& realizationNum)
{
  COMPLEX* kCorrEnsemble = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
  // explicitly initialize the array to zero
  std::memset(kCorrEnsemble, 0, sizeof(COMPLEX) * this->NK);

  // allocate a vector for all sampled noise (for noise distribution check)
  std::vector<double> sampledNoiseAll;

  for (int i=1;i<realizationNum+1;i++)
  {
    COMPLEX* kNoise= (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
    // explicitly initialize the array to zero
    std::memset(kNoise, 0, sizeof(COMPLEX) * this->NK);

    //unsigned int randomSeed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()); // loop too fucking fast
    unsigned int randomSeed = i;
    sampleNoiseLocalInKspace(mat, randomSeed, kNoise);

    // accumulate calculated correlation values
    for (int j=0; j<this->NK; ++j)
    {
      // Wiener-Khinchin theorem
      kCorrEnsemble[j] += kNoise[j]*std::conj(kNoise[j]);
    }

    // sample noise in real space (for noise distribution check)
    REAL_SCALAR* noise = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
    // explicitly initialize the array to zero
    std::memset(noise, 0, sizeof(REAL_SCALAR) * this->NR);
    fftw_plan sampleNoise = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kNoise), noise, FFTW_ESTIMATE);
    // populate noise array
    fftw_execute(sampleNoise);
    // save the sampled noises
    for (int k=0; k<this->NR; ++k)
    {
      sampledNoiseAll.push_back(noise[k]);
    }
    // memory cleanup
    fftw_destroy_plan(sampleNoise);
    fftw_free(noise);
    fftw_free(kNoise);
  }

  // average the ensemble
  for (int i=0; i<this->NK; ++i)
  {
    kCorrEnsemble[i] /= static_cast<double>(realizationNum);
  }

  REAL_SCALAR* rCorrEnsemble = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
  // explicitly initialize the array to zero
  std::memset(rCorrEnsemble, 0, sizeof(REAL_SCALAR) * this->NR);
  fftw_plan sampleEnsembleCorr = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kCorrEnsemble), rCorrEnsemble, FFTW_ESTIMATE);
  fftw_execute(sampleEnsembleCorr);

  // circular shift
  // shift the correlation to the center
  circularShift(rCorrEnsemble);

  // save the sampled noises in InputFiles dir
  std::string noiseDistFname = std::format(
    "{}/noiseDistributionR{}.txt",
    std::filesystem::path(mat.materialFile).parent_path().string(),
    realizationNum
  );
  std::ofstream noiseDistOut(noiseDistFname, std::ios::out); //overwrite
  for (const auto& value : sampledNoiseAll) {
    noiseDistOut << value << std::endl;
  }

  // save the ensemble correlation in InputFiles dir
  std::string ensembleCorrFname = std::format(
    "{}/ensembledCorrelationR{}.txt",
    std::filesystem::path(mat.materialFile).parent_path().string(),
    realizationNum
  );
  std::ofstream ensembleCorrOut(ensembleCorrFname, std::ios::out); //overwrite
  for (int i=0; i<this->NR; ++i)
  {
    ensembleCorrOut << rCorrEnsemble[i] << std::endl;
  }

  // memory cleanup
  fftw_destroy_plan(sampleEnsembleCorr);
  fftw_free(rCorrEnsemble);
  fftw_free(kCorrEnsemble);
}


// write noise patch in DDD unit (unitless)
void MDStackingFaultNoise::Write_field_slice(const PolycrystallineMaterialBase& mat,
                                             const int& seed,
                                             const char *fname)
{
    const double DX = this->gridSpacing(0);
    const double DY = this->gridSpacing(1);
    const double DZ = this->gridSpacing(2);
    FILE *OutFile=fopen(fname,"w");

    COMPLEX* kNoise= (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
    // explicitly initialize with zeros
    std::memset(kNoise, 0, sizeof(COMPLEX) * this->NK);
    // populate the noise in fourier space
    sampleNoiseLocalInKspace(mat, seed, kNoise);
    REAL_SCALAR* noise = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
    // explicitly initialize with zeros
    std::memset(noise, 0, sizeof(REAL_SCALAR) * this->NR);
    fftw_plan sampleNoise = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kNoise), noise, FFTW_ESTIMATE);
    fftw_execute(sampleNoise);

    fprintf(OutFile,"# vtk DataFile Version 2.0\n");
    fprintf(OutFile,"iter %d\n",0);
    fprintf(OutFile,"BINARY\n");
    fprintf(OutFile,"DATASET STRUCTURED_POINTS\n");
    fprintf(OutFile,"ORIGIN \t %f %f %f\n",0.,0.,0.);
    fprintf(OutFile,"SPACING \t %f %f %f\n", DX, DY, DZ);
    fprintf(OutFile,"DIMENSIONS \t %d %d %d\n", this->NX, this->NY, 1);
    fprintf(OutFile,"POINT_DATA \t %d\n",this->NX*this->NY);
    fprintf(OutFile,"SCALARS \t volume_scalars double 1\n");
    fprintf(OutFile,"LOOKUP_TABLE \t default\n");

    for(int i=0;i<this->NX;i++)
    {
        for(int j=0;j<this->NY;j++)
        {
            const int k=0;
            const int ind = this->NY*this->NZ*i + j*this->NZ + k;
            //const double temp=NoiseTraitsBase::ReverseDouble(double(rNoisyCorrelations[ind]));
            const double temp=NoiseTraitsBase::ReverseDouble(double(noise[ind]));
            fwrite(&temp, sizeof(REAL_SCALAR), 1, OutFile);
        }
    }
    fclose(OutFile);

    // memory cleanup
    fftw_destroy_plan(sampleNoise);
    fftw_free(kNoise);
    fftw_free(noise);
}

}
#endif
