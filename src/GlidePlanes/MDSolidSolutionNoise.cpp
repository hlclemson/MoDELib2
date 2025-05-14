/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MDSolidSolutionNoise_cpp
#define model_MDSolidSolutionNoise_cpp

#include <MDSolidSolutionNoise.h>

namespace model
{

    MDSolidSolutionNoise::MDSolidSolutionNoise(const PolycrystallineMaterialBase& mat,
                                               const std::string& tag,
                                               const std::string& correlationFile_L,
                                               const std::string& correlationFile_T,
                                               const int& outputNoise,
                                               const int& testNoiseSampling,
                                               const std::string& noiseFile_L,
                                               const std::string& noiseFile_T,
                                               const int& seed,
                                               const GridSizeType& gridSize,
                                               const GridSpacingType& gridSpacing,
                                               const double& a_Cai_in
                                               ) :
    /*init*/ GlidePlaneNoiseBase<2>("MDSolidSolutionNoise"+tag,seed,gridSize,gridSpacing)
    /*init*/,a_cai(a_Cai_in)
    {
        // allocate
        REAL_SCALAR *Rr_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR); //correlation in real space
        REAL_SCALAR *Rr_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR); //correlation in real space
                
        Rk_xz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK); //correlation in fourier space
        Rk_yz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK); //correlation in fourier space
        
        // Initialize FFTW plans as member variables
        fftw_plan plan_R_xz_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr_xz, reinterpret_cast<fftw_complex*>(Rk_xz), FFTW_ESTIMATE);
        fftw_plan plan_R_yz_r2c = fftw_plan_dft_r2c_2d(this->NY, this->NX, Rr_yz, reinterpret_cast<fftw_complex*>(Rk_yz), FFTW_ESTIMATE);
        
        const auto originalDimensionsL(readVTKfileDimension(correlationFile_L.c_str()));
        const auto originalDimensionsT(readVTKfileDimension(correlationFile_T.c_str()));
        if((originalDimensionsL-originalDimensionsT).matrix().squaredNorm()>0)
        {
            throw std::runtime_error("correlationFile_L and correlationFile_T have different grid sizes.");
        }
        const int originalNX = originalDimensionsL(0);
        const int originalNY = originalDimensionsL(1);
        if(originalDimensionsL(2)!=1)
        {
            throw std::runtime_error("vtk correlationFiles 'DIMENSIONS' should have 3rd component == 1.");
        }
        
        REAL_SCALAR *Rr_xz_original = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*originalNX*originalNY); //correlation in real space
        REAL_SCALAR *Rr_yz_original = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*originalNX*originalNY); //correlation in real space
                
        // populate Rr_xy_original with the correlation data
        SolidSolutionCorrelationReader(correlationFile_L, Rr_xz_original, this->NR);
        SolidSolutionCorrelationReader(correlationFile_T, Rr_yz_original, this->NR);
                
        // Divide the values in Rr_xz_original and Rr_yz_original by mat.mu^2
        for (int i = 0; i < originalNX * originalNY; ++i)
        {
            Rr_xz_original[i] /= (mat.mu_SI*mat.mu_SI); // divide by mu^2
            Rr_yz_original[i] /= (mat.mu_SI*mat.mu_SI); // divide by mu^2
        }
        
        int start_y = (this->NY - originalNY) / 2;
        int start_x = (this->NX - originalNX) / 2;
        
        // // initialize with zeros
        for (int i = 0; i < this->NY; ++i)
        {
            for (int j = 0; j < this->NX; ++j)
            {
                Rr_xz[i*this->NX + j] = 0;
                Rr_yz[i*this->NX + j] = 0;
            }
        }
        
        // 0-pading from centere
        for (int i = 0; i < originalNY; ++i)
        {
            for (int j = 0; j < originalNX; ++j)
            {
                Rr_xz[(start_y + i) * this->NX + (start_x + j)] = Rr_xz_original[i * originalNX + j];
                Rr_yz[(start_y + i) * this->NX + (start_x + j)] = Rr_yz_original[i * originalNX + j];
            }
        }
                
        // Execute FFTW plans to populate Rk_xz and Rk_yz
        fftw_execute(plan_R_xz_r2c);
        fftw_execute(plan_R_yz_r2c);
        
        // Normalize the FFT output
        for (int i = 0; i < this->NX; ++i)
        {
            for (int j = 0; j < (this->NY/2 + 1); ++j)
            {
                Rk_xz[i * (this->NY/2 + 1) + j] /= (this->NX * this->NY);
                Rk_yz[i * (this->NY/2 + 1) + j] /= (this->NX * this->NY);
            }
        }
        
        // noise output for debugging
        if (outputNoise) { Write_field_slice(mat, seed, noiseFile_L.c_str(), noiseFile_T.c_str()); };
        if (testNoiseSampling) { 
          sampleNoiseRepeatedly(mat, 10);
          sampleNoiseRepeatedly(mat, 100);
          sampleNoiseRepeatedly(mat, 1000);
        };

        // Destroy FFTW plans
        fftw_destroy_plan(plan_R_xz_r2c);
        fftw_destroy_plan(plan_R_yz_r2c);
        
        // Free allocated memory
        fftw_free(Rr_xz_original);
        fftw_free(Rr_yz_original);
        
        // Free allocated memory
        fftw_free(Rr_xz);
        fftw_free(Rr_yz);
    }

    std::array<MDSolidSolutionNoise::COMPLEX, 2> MDSolidSolutionNoise::kCorrelations(const Eigen::Matrix<double, 3, 1> &kv, const Eigen::Matrix<int, 3, 1> &index) const
    {
        std::array<MDSolidSolutionNoise::COMPLEX, 2> temp;
        int idx=(this->NY/2+1)*NZ*index(0) + index(1)*NZ + index(2);
        temp[0] = Rk_xz[idx];
        temp[1] = Rk_yz[idx];
        if(a_cai>0.0)
        {
            const double wkc2(this->Wk_Cai_squared(kv(0),kv(1),kv(2), a_cai)); // using the square because this is before the square root
            temp[0]*=wkc2;
            temp[1]*=wkc2;
        }
        return temp;
    }

    void MDSolidSolutionNoise::SolidSolutionCorrelationReader(const std::string& correlationFile, REAL_SCALAR *Rr, int NR)
    {
        std::cout << "Reading solid solution correlation" << std::endl;
        
        std::ifstream vtkFile(correlationFile); //access vtk file
        // error check
        if (!vtkFile.is_open())
        {
            throw std::runtime_error("Error opening solid solution VTK correlation file!");
        }
        
        std::string line;
        while (std::getline(vtkFile, line))
        {
            //if the "POINT_DATA" string is read, read the following data
            if(line.find("POINT_DATA")!=std::string::npos)
            {
                const size_t numOfHeaders = 2;
                // get the number of points in the file
                const size_t firstSpace(line.find(' '));
                const size_t numOfPoints = std::atoi(line.substr(firstSpace+1).c_str());
                // read the point coordinates
                for(size_t n=0; n<numOfPoints+numOfHeaders; ++n)
                {
                    std::getline(vtkFile, line);
                    // ignore the headers that come right after point_data
                    if(n<numOfHeaders)
                        continue;
                    const int ind = n-numOfHeaders;

                    //Rr[ind] = std::atof(line.c_str());
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
                }
            }
        }
        vtkFile.close();
    }

    typename MDSolidSolutionNoise::GridSizeType MDSolidSolutionNoise::readVTKfileDimension(const char *fname)
    {
        FILE *InFile=fopen(fname,"r");
        
        if (InFile == NULL)
        {
            fprintf(stderr, "Can't open solid solution correlation VTK file %s\n",fname);
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

    void MDSolidSolutionNoise::sampleNoiseLocalInKspace(const PolycrystallineMaterialBase& mat, const int& localSeed, COMPLEX* localRk_xz, COMPLEX* localRk_yz)
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
                  REAL_SCALAR Nk_yz = distribution(generator);
                  REAL_SCALAR Mk_yz = distribution(generator);
                  REAL_SCALAR Nk_xz, Mk_xz;
                  if(kx*ky>=0)
                  {
                      Nk_xz = Nk_yz;
                      Mk_xz = Mk_yz;
                  }
                  else
                  {
                      Nk_xz = -Nk_yz;
                      Mk_xz = -Mk_yz;
                  }
                  // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform 
                  const double kCorrFactor((j==0 || j==this->NY/2)? 1.0 : 2.0);
                  //double temp_xz[ind]=std::sqrt(Rk_xz[ind]/kCorrFactor)*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
                  //double temp_yz[ind]=std::sqrt(Rk_yz[ind]/kCorrFactor)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
                  COMPLEX temp_xz = std::sqrt(Rk_xz[ind]/kCorrFactor)*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
                  COMPLEX temp_yz = std::sqrt(Rk_yz[ind]/kCorrFactor)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
                  if(a_cai>0.0)
                  {
                      const double wkc2(this->Wk_Cai_squared(kx, ky, kz, a_cai)); // using the square because this is before the square root
                      temp_xz*=wkc2;
                      temp_yz*=wkc2;
                  }
                  localRk_xz[ind] = temp_xz;
                  localRk_yz[ind] = temp_yz;
              }
          }
      }
    }

    // shift the correlation to the center
    void MDSolidSolutionNoise::circularShift(REAL_SCALAR* rCorrelation)
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
    void MDSolidSolutionNoise::sampleNoiseRepeatedly(const PolycrystallineMaterialBase& mat, const int& realizationNum)
    {
      COMPLEX* kCorrEnsemble_xz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
      COMPLEX* kCorrEnsemble_yz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
      // explicitly initialize the array to zero
      std::memset(kCorrEnsemble_xz, 0, sizeof(COMPLEX) * this->NK);
      std::memset(kCorrEnsemble_yz, 0, sizeof(COMPLEX) * this->NK);

      // allocate a vector for all sampled noise (for noise distribution check)
      std::vector<double> sampledNoiseAll_xz;
      std::vector<double> sampledNoiseAll_yz;

      for (int i=1;i<realizationNum+1;i++)
      {
        COMPLEX* kNoise_xz= (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
        COMPLEX* kNoise_yz= (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
        // explicitly initialize the array to zero
        std::memset(kNoise_xz, 0, sizeof(COMPLEX) * this->NK);
        std::memset(kNoise_yz, 0, sizeof(COMPLEX) * this->NK);

        //unsigned int randomSeed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()); // loop too fucking fast
        unsigned int randomSeed = i;
        sampleNoiseLocalInKspace(mat, randomSeed, kNoise_xz, kNoise_yz);
        //sampleNoiseLocalInKspace(mat, randomSeed, kNoise_yz);

        // accumulate calculated correlation values
        for (int j=0; j<this->NK; ++j)
        {
          // Wiener-Khinchin theorem
          kCorrEnsemble_xz[j] += kNoise_xz[j]*std::conj(kNoise_xz[j]);
          kCorrEnsemble_yz[j] += kNoise_yz[j]*std::conj(kNoise_yz[j]);
        }

        // sample noise in real space (for noise distribution check)
        REAL_SCALAR* noise_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
        REAL_SCALAR* noise_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
        // explicitly initialize the array to zero
        std::memset(noise_xz, 0, sizeof(REAL_SCALAR) * this->NR);
        std::memset(noise_yz, 0, sizeof(REAL_SCALAR) * this->NR);
        fftw_plan sampleNoise_xz = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kNoise_xz), noise_xz, FFTW_ESTIMATE);
        fftw_plan sampleNoise_yz = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kNoise_yz), noise_yz, FFTW_ESTIMATE);
        // populate noise array
        fftw_execute(sampleNoise_xz);
        fftw_execute(sampleNoise_yz);
        // save the sampled noises
        for (int k=0; k<this->NR; ++k)
        {
          sampledNoiseAll_xz.push_back(noise_xz[k]);
          sampledNoiseAll_yz.push_back(noise_yz[k]);
        }
        // memory cleanup
        fftw_destroy_plan(sampleNoise_xz);
        fftw_destroy_plan(sampleNoise_yz);
        fftw_free(noise_xz);
        fftw_free(noise_yz);
        fftw_free(kNoise_xz);
        fftw_free(kNoise_yz);
      }

      // average the ensemble
      for (int i=0; i<this->NK; ++i)
      {
        kCorrEnsemble_xz[i] /= static_cast<double>(realizationNum);
        kCorrEnsemble_yz[i] /= static_cast<double>(realizationNum);
      }

      REAL_SCALAR* rCorrEnsemble_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
      REAL_SCALAR* rCorrEnsemble_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
      // explicitly initialize the array to zero
      std::memset(rCorrEnsemble_xz, 0, sizeof(REAL_SCALAR) * this->NR);
      std::memset(rCorrEnsemble_yz, 0, sizeof(REAL_SCALAR) * this->NR);
      fftw_plan sampleEnsembleCorr_xz = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kCorrEnsemble_xz), rCorrEnsemble_xz, FFTW_ESTIMATE);
      fftw_plan sampleEnsembleCorr_yz = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kCorrEnsemble_yz), rCorrEnsemble_yz, FFTW_ESTIMATE);
      fftw_execute(sampleEnsembleCorr_xz);
      fftw_execute(sampleEnsembleCorr_yz);

      // circular shift
      // shift the correlation to the center
      circularShift(rCorrEnsemble_xz);
      circularShift(rCorrEnsemble_yz);

      // save the sampled noises xz in InputFiles dir
      std::string noiseDistFname_xz = std::format(
        "{}/noiseDistribution_xz_R{}.txt",
        std::filesystem::path(mat.materialFile).parent_path().string(),
        realizationNum
      );
      std::ofstream noiseDistOut_xz(noiseDistFname_xz, std::ios::out); //overwrite
      for (const auto& value : sampledNoiseAll_xz) {
        noiseDistOut_xz << value << std::endl;
      }
      // save the sampled noises yz in InputFiles dir
      std::string noiseDistFname_yz = std::format(
        "{}/noiseDistribution_yz_R{}.txt",
        std::filesystem::path(mat.materialFile).parent_path().string(),
        realizationNum
      );
      std::ofstream noiseDistOut_yz(noiseDistFname_yz, std::ios::out); //overwrite
      for (const auto& value : sampledNoiseAll_yz) {
        noiseDistOut_yz << value << std::endl;
      }

      // save the ensemble correlation in InputFiles dir
      std::string ensembleCorrFname_xz = std::format(
        "{}/ensembledCorrelation_xz_R{}.txt",
        std::filesystem::path(mat.materialFile).parent_path().string(),
        realizationNum
      );
      std::ofstream ensembleCorrOut_xz(ensembleCorrFname_xz, std::ios::out); //overwrite
      for (int i=0; i<this->NR; ++i)
      {
        ensembleCorrOut_xz << rCorrEnsemble_xz[i] << std::endl;
      }
      // save the ensemble correlation in InputFiles dir
      std::string ensembleCorrFname_yz = std::format(
        "{}/ensembledCorrelation_yz_R{}.txt",
        std::filesystem::path(mat.materialFile).parent_path().string(),
        realizationNum
      );
      std::ofstream ensembleCorrOut_yz(ensembleCorrFname_yz, std::ios::out); //overwrite
      for (int i=0; i<this->NR; ++i)
      {
        ensembleCorrOut_yz << rCorrEnsemble_yz[i] << std::endl;
      }

      // memory cleanup
      fftw_destroy_plan(sampleEnsembleCorr_xz);
      fftw_destroy_plan(sampleEnsembleCorr_yz);
      fftw_free(rCorrEnsemble_xz);
      fftw_free(rCorrEnsemble_yz);
      fftw_free(kCorrEnsemble_xz);
      fftw_free(kCorrEnsemble_yz);
    }


    // write noise patch in DDD unit (unitless)
    void MDSolidSolutionNoise::Write_field_slice(const PolycrystallineMaterialBase& mat,
                                                const int& seed,
                                                const char *fname_L,
                                                const char *fname_T)
    {
        const double DX = this->gridSpacing(0);
        const double DY = this->gridSpacing(1);
        const double DZ = this->gridSpacing(2);

        COMPLEX* kNoise_xz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
        COMPLEX* kNoise_yz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*this->NK);
        // explicitly initialize with zeros
        std::memset(kNoise_xz, 0, sizeof(COMPLEX) * this->NK);
        std::memset(kNoise_yz, 0, sizeof(COMPLEX) * this->NK);
        // populate the noise in fourier space
        sampleNoiseLocalInKspace(mat, seed, kNoise_xz, kNoise_yz);
        REAL_SCALAR* noise_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
        REAL_SCALAR* noise_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*this->NR);
        // explicitly initialize with zeros
        std::memset(noise_xz, 0, sizeof(REAL_SCALAR) * this->NR);
        std::memset(noise_yz, 0, sizeof(REAL_SCALAR) * this->NR);
        fftw_plan sampleNoise_xz = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kNoise_xz), noise_xz, FFTW_ESTIMATE);
        fftw_plan sampleNoise_yz = fftw_plan_dft_c2r_2d(this->NY, this->NX, reinterpret_cast<fftw_complex*>(kNoise_yz), noise_yz, FFTW_ESTIMATE);
        fftw_execute(sampleNoise_xz);
        fftw_execute(sampleNoise_yz);

        // print xz
        FILE *OutFile_L=fopen(fname_L,"w");
        fprintf(OutFile_L,"# vtk DataFile Version 2.0\n");
        fprintf(OutFile_L,"iter %d\n",0);
        fprintf(OutFile_L,"BINARY\n");
        fprintf(OutFile_L,"DATASET STRUCTURED_POINTS\n");
        fprintf(OutFile_L,"ORIGIN \t %f %f %f\n",0.,0.,0.);
        fprintf(OutFile_L,"SPACING \t %f %f %f\n", DX, DY, DZ);
        fprintf(OutFile_L,"DIMENSIONS \t %d %d %d\n", this->NX, this->NY, 1);
        fprintf(OutFile_L,"POINT_DATA \t %d\n",this->NX*this->NY);
        fprintf(OutFile_L,"SCALARS \t volume_scalars double 1\n");
        fprintf(OutFile_L,"LOOKUP_TABLE \t default\n");

        for(int i=0;i<this->NX;i++)
        {
            for(int j=0;j<this->NY;j++)
            {
                const int k=0;
                const int ind = this->NY*this->NZ*i + j*this->NZ + k;
                const double temp=NoiseTraitsBase::ReverseDouble(double(noise_xz[ind]));
                fwrite(&temp, sizeof(REAL_SCALAR), 1, OutFile_L);
            }
        }
        fclose(OutFile_L);

        // print yz
        FILE *OutFile_T=fopen(fname_T,"w");
        fprintf(OutFile_T,"# vtk DataFile Version 2.0\n");
        fprintf(OutFile_T,"iter %d\n",0);
        fprintf(OutFile_T,"BINARY\n");
        fprintf(OutFile_T,"DATASET STRUCTURED_POINTS\n");
        fprintf(OutFile_T,"ORIGIN \t %f %f %f\n",0.,0.,0.);
        fprintf(OutFile_T,"SPACING \t %f %f %f\n", DX, DY, DZ);
        fprintf(OutFile_T,"DIMENSIONS \t %d %d %d\n", this->NX, this->NY, 1);
        fprintf(OutFile_T,"POINT_DATA \t %d\n",this->NX*this->NY);
        fprintf(OutFile_T,"SCALARS \t volume_scalars double 1\n");
        fprintf(OutFile_T,"LOOKUP_TABLE \t default\n");

        for(int i=0;i<this->NX;i++)
        {
            for(int j=0;j<this->NY;j++)
            {
                const int k=0;
                const int ind = this->NY*this->NZ*i + j*this->NZ + k;
                const double temp=NoiseTraitsBase::ReverseDouble(double(noise_yz[ind]));
                fwrite(&temp, sizeof(REAL_SCALAR), 1, OutFile_T);
            }
        }
        fclose(OutFile_T);

        // memory cleanup
        fftw_destroy_plan(sampleNoise_xz);
        fftw_destroy_plan(sampleNoise_yz);
        fftw_free(kNoise_xz);
        fftw_free(kNoise_yz);
        fftw_free(noise_xz);
        fftw_free(noise_yz);
    }

}
#endif
