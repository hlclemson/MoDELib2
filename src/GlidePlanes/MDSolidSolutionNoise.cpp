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
        if (outputNoise) { Write_field_slice(Rk_xz, Rk_yz, seed, this->NX, this->NY, this->NZ, this->NK, this->NR, this->LX, this->LY, this->LZ, gridSize, gridSpacing, noiseFile_L.c_str(), noiseFile_T.c_str()); };
        
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

    void MDSolidSolutionNoise::Write_field_slice(MDSolidSolutionNoise::COMPLEX *Rk_xz, 
                                                 MDSolidSolutionNoise::COMPLEX *Rk_yz,  
                                                 const int& seed,
                                                 const int& nx,
                                                 const int& ny,
                                                 const int& nz,
                                                 const int& nk,
                                                 const int& nr,
                                                 const int& lx,
                                                 const int& ly,
                                                 const int& lz,
                                                 MDSolidSolutionNoise::GridSizeType gridSize,
                                                 MDSolidSolutionNoise::GridSpacingType gridSpacing,
                                                 const char *fname_xz,
                                                 const char *fname_yz)
    {
        const float dx = gridSpacing(0);
        const float dy = gridSpacing(1);
        const float dz = gridSpacing(2);
        FILE *OutFile_xz=fopen(fname_xz,"w");
        FILE *OutFile_yz=fopen(fname_yz,"w");

        MDSolidSolutionNoise::COMPLEX* kNoisyCorrelations_xz = (MDSolidSolutionNoise::COMPLEX*) fftw_malloc(sizeof(MDSolidSolutionNoise::COMPLEX)*nk);
        MDSolidSolutionNoise::COMPLEX* kNoisyCorrelations_yz = (MDSolidSolutionNoise::COMPLEX*) fftw_malloc(sizeof(MDSolidSolutionNoise::COMPLEX)*nk);
        MDSolidSolutionNoise::REAL_SCALAR* rNoisyCorrelations_xz = (MDSolidSolutionNoise::REAL_SCALAR*) fftw_malloc(sizeof(MDSolidSolutionNoise::REAL_SCALAR)*nr);
        MDSolidSolutionNoise::REAL_SCALAR* rNoisyCorrelations_yz = (MDSolidSolutionNoise::REAL_SCALAR*) fftw_malloc(sizeof(MDSolidSolutionNoise::REAL_SCALAR)*nr);
        
        const int J_MAX = ny/2 + 1;
        const int K_MAX = 1;

        std::default_random_engine generator(seed);
        std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);
        for(int i=0; i<nx; i++)
        {
            for(int j=0; j<J_MAX; j++)
            {
                for(int k=0; k<K_MAX; k++)
                {
                    const int ind = J_MAX*K_MAX*i + j*K_MAX + k;

                    REAL_SCALAR kx = 2.*M_PI/lx*REAL_SCALAR(i);
                    if(i>nx/2)
                    {
                        kx = 2.*M_PI/lx*REAL_SCALAR(i-nx);
                    }
                    
                    REAL_SCALAR ky = 2*M_PI/ly*REAL_SCALAR(j);
                    if(j>ny/2)
                    {
                        ky = 2.*M_PI/ly*REAL_SCALAR(j-ny);
                    }
                    
                    REAL_SCALAR kz = 2.*M_PI/lz*REAL_SCALAR(k);
                    
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
                    
                    const double kCorrFactor((j==0 || j==ny/2)? 1.0 : 2.0); // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
                    kNoisyCorrelations_xz[ind]=sqrt(Rk_xz[ind]/kCorrFactor)*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
                    kNoisyCorrelations_yz[ind]=sqrt(Rk_yz[ind]/kCorrFactor)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
                }
            }
        }


        fftw_plan nPlan_xz = fftw_plan_dft_c2r_2d(ny, nx, reinterpret_cast<fftw_complex*>(kNoisyCorrelations_xz), rNoisyCorrelations_xz, FFTW_ESTIMATE);
        fftw_plan nPlan_yz = fftw_plan_dft_c2r_2d(ny, nx, reinterpret_cast<fftw_complex*>(kNoisyCorrelations_yz), rNoisyCorrelations_yz, FFTW_ESTIMATE);
        fftw_execute(nPlan_xz);
        fftw_execute(nPlan_yz);

        fprintf(OutFile_xz,"# vtk DataFile Version 2.0\n");
        fprintf(OutFile_xz,"iter %d\n",0);
        fprintf(OutFile_xz,"BINARY\n");
        fprintf(OutFile_xz,"DATASET STRUCTURED_POINTS\n");
        fprintf(OutFile_xz,"ORIGIN \t %f %f %f\n",0.,0.,0.);
        fprintf(OutFile_xz,"SPACING \t %f %f %f\n", dx, dy, dz);
        fprintf(OutFile_xz,"DIMENSIONS \t %d %d %d\n", nx, ny, 1);
        fprintf(OutFile_xz,"POINT_DATA \t %d\n",nx*ny);
        fprintf(OutFile_xz,"SCALARS \t volume_scalars double 1\n");
        fprintf(OutFile_xz,"LOOKUP_TABLE \t default\n");

        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                const int k=0;
                const int ind = ny*nz*i + j*nz + k;
                //const double temp=NoiseTraitsBase::ReverseDouble(double(F[ind]));
                const double temp=NoiseTraitsBase::ReverseDouble(double(rNoisyCorrelations_xz[ind]));
                fwrite(&temp, sizeof(double), 1, OutFile_xz);
            }
        }
        fclose(OutFile_xz);

        fprintf(OutFile_yz,"# vtk DataFile Version 2.0\n");
        fprintf(OutFile_yz,"iter %d\n",0);
        fprintf(OutFile_yz,"BINARY\n");
        fprintf(OutFile_yz,"DATASET STRUCTURED_POINTS\n");
        fprintf(OutFile_yz,"ORIGIN \t %f %f %f\n",0.,0.,0.);
        fprintf(OutFile_yz,"SPACING \t %f %f %f\n", dx, dy, dz);
        fprintf(OutFile_yz,"DIMENSIONS \t %d %d %d\n", nx, ny, 1);
        fprintf(OutFile_yz,"POINT_DATA \t %d\n",nx*ny);
        fprintf(OutFile_yz,"SCALARS \t volume_scalars double 1\n");
        fprintf(OutFile_yz,"LOOKUP_TABLE \t default\n");

        for(int i=0;i<nx;i++)
        {
            for(int j=0;j<ny;j++)
            {
                const int k=0;
                const int ind = ny*nz*i + j*nz + k;
                //const double temp=NoiseTraitsBase::ReverseDouble(double(F[ind]));
                const double temp=NoiseTraitsBase::ReverseDouble(double(rNoisyCorrelations_yz[ind]));
                fwrite(&temp, sizeof(double), 1, OutFile_yz);
            }
        }
        fclose(OutFile_yz);
    }
}
#endif
