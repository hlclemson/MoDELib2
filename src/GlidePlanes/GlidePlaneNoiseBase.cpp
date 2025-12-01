/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneNoiseBase_cpp
#define model_GlidePlaneNoiseBase_cpp

#include <chrono>
#include <GlidePlaneNoiseBase.h>
#include <TerminalColors.h>
#include <cstdint>

namespace model
{

// Cai's doubly-convoluted spreading function in Fourier space
template <int N>
typename GlidePlaneNoiseBase<N>::REAL_SCALAR GlidePlaneNoiseBase<N>::Wk_Cai(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz, REAL_SCALAR a)
{
  REAL_SCALAR k = sqrt(kx*kx + ky*ky + kz*kz);
  if(k>0.0)
  {
    return a*k*sqrt(0.5*boost::math::cyl_bessel_k(2,a*k));
  }
  else
  {
    return 1.;
  }
}

// Square of Cai's doubly-convoluted spreading function in Fourier space
template <int N>
typename GlidePlaneNoiseBase<N>::REAL_SCALAR GlidePlaneNoiseBase<N>::Wk_Cai_squared(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz, REAL_SCALAR a)
{
  const REAL_SCALAR k2(kx*kx + ky*ky + kz*kz);
  const REAL_SCALAR k(sqrt(k2));
  if(k2>0.0)
  {
    return a*a*k2*0.5*boost::math::cyl_bessel_k(2,a*k);
  }
  else
  {
    return 1.;
  }
}

// Cai spreading function
template <int N>
typename GlidePlaneNoiseBase<N>::REAL_SCALAR GlidePlaneNoiseBase<N>::W_Cai(REAL_SCALAR r2, REAL_SCALAR a)
{
  return 15.*a*a*a*a/(8.*M_PI*pow(r2+a*a,7./2.));
}

template <int N>
typename GlidePlaneNoiseBase<N>::REAL_SCALAR GlidePlaneNoiseBase<N>::W_t_Cai(REAL_SCALAR r2, REAL_SCALAR a)
{
  return 0.3425*W_Cai(r2,0.9038*a) + 0.6575*W_Cai(r2,0.5451*a);
}

template <int N>
GlidePlaneNoiseBase<N>::GlidePlaneNoiseBase(const std::string& tag_in,
                                            const int& seed_in,
                                            const NoiseTraitsBase::GridSizeType& gridSize_in,
                                            const NoiseTraitsBase::GridSpacingType& gridSpacing_in,
                                            const Eigen::Matrix<double,2,2>& latticeBasis,
                                           const double& effsroAve_in):
  /*init*/ UniformPeriodicGrid<2>(gridSize_in.template segment<2>(0),gridSpacing_in.template segment<2>(0))
  /*init*/,tag(tag_in)
  /*init*/,seed(seed_in)
  /*init*/,gridSize(gridSize_in)
  /*init*/,gridSpacing(gridSpacing_in)
  /*init*/,gridLength(gridSize.template cast<double>()*gridSpacing)
  /*init*/,invTransposeLatticeBasis(latticeBasis.transpose().inverse())
  /*init*/,NX(gridSize(0))
  /*init*/,NY(gridSize(1))
  /*init*/,NZ(gridSize(2))
  /*init*/,NK(NZ>1 ? NX*NY*(NZ/2+1) : NX*(NY/2+1))
  /*init*/,NR(NX*NY*NZ)
  /*init*/,LX(gridLength(0))
  /*init*/,LY(gridLength(1))
  /*init*/,LZ(gridLength(2))
  /* int */,effsroAve(effsroAve_in)
{
}

template <int N>
const typename GlidePlaneNoiseBase<N>::NoiseContainerType& GlidePlaneNoiseBase<N>::noiseVector() const
{
  return *this;
}

template <int N>
typename GlidePlaneNoiseBase<N>::NoiseContainerType& GlidePlaneNoiseBase<N>::noiseVector()
{
  return *this;
}

template <int N>
int GlidePlaneNoiseBase<N>::storageIndex(const int& i,const int& j) const
{/*!\param[in] localPos the  position vector on the grid
      * \returns The grid index periodically wrapped within the gridSize bounds
      */
  return gridSize(1)*i+j;
}

#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_

template <int N>
std::vector<typename GlidePlaneNoiseBase<N>::COMPLEX*> GlidePlaneNoiseBase<N>::originalKCorrelations() const
{
  std::vector<COMPLEX*> okCorr;
  for(int n=0;n<N;++n)
  {
    okCorr.push_back((COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK));
    // explicitly initialize to zeros, avoid garbarge propagation
    std::fill(okCorr[n], okCorr[n]+NK, COMPLEX{0.0, 0.0});
  }

  const int J_MAX = (NZ > 1) ? NY : (NY/2 + 1);
  const int K_MAX = (NZ > 1) ? (NZ/2 + 1) : 1 ;
  for(int i=0; i<NX; i++)
  {
    for(int j=0; j<J_MAX; j++)
    {
      for(int k=0; k<K_MAX; k++)
      {
        //                const int ind = NY*(NZ/2+1)*i + j*(NZ/2+1) + k;
        const int ind = J_MAX*K_MAX*i + j*K_MAX + k;

        REAL_SCALAR kx = 2.*M_PI/LX*REAL_SCALAR(i);
        if(i>NX/2)
        {
          kx = 2.*M_PI/LX*REAL_SCALAR(i-NX);
        }

        REAL_SCALAR ky = 2*M_PI/LY*REAL_SCALAR(j);
        if(j>NY/2)
        {
          ky = 2.*M_PI/LY*REAL_SCALAR(j-NY);
        }

        REAL_SCALAR kz = 2.*M_PI/LZ*REAL_SCALAR(k);

        const Eigen::Matrix<double,3,1> kv((Eigen::Matrix<double,3,1>()<<kx,ky,kz).finished());
        const Eigen::Matrix<int,3,1> kvID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
        //const double kCorrFactor(NZ>1 ? ((k==0 || k==NZ/2)? 1.0 : 2.0) : ((j==0 || j==NY/2)? 1.0 : 2.0)); // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
        const auto kCorr(kCorrelations(kv,kvID));
        for(int n=0;n<N;++n)
        {
          okCorr[n][ind]=sqrt(kCorr[n]);
        }
      }
    }
  }
  return okCorr;
}

template <int N>
std::vector<typename GlidePlaneNoiseBase<N>::COMPLEX*> GlidePlaneNoiseBase<N>::noisyKCorrelations(const std::vector<COMPLEX*>& okCorr,
                                                                                                  std::default_random_engine& generator,
                                                                                                  std::normal_distribution<REAL_SCALAR>& distribution) const
{
  //std::vector<COMPLEX*> nkCorr(okCorr);
  std::vector<COMPLEX*> nkCorr;
  for(int n=0;n<N;++n)
  {
    nkCorr.push_back((COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK));
    // explicitly initialize to zeros, avoid garbarge propagation
    std::fill(nkCorr[n], nkCorr[n]+NK, COMPLEX{0.0, 0.0});
  }

  for(int ind=0;ind<NK;++ind)
  {
    const REAL_SCALAR Nk = distribution(generator);
    const REAL_SCALAR Mk = distribution(generator);
    for(int n=0;n<N;++n)
    {
      nkCorr[n][ind]=okCorr[n][ind]*COMPLEX(Nk,Mk)/std::numbers::sqrt2;
    }
  }

  for(int n=0;n<N;++n)
  {
    nkCorr[n][0]=0;
  }

  return nkCorr;
}

template <int N>
std::vector<typename GlidePlaneNoiseBase<N>::REAL_SCALAR*> GlidePlaneNoiseBase<N>::realNoise(const std::vector<COMPLEX*>& nkCorr) const
{
  std::vector<REAL_SCALAR*> rNoise;
  for(int n=0;n<N;++n)
  {
    rNoise.push_back((REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR));
    // explicitly initialize to zeros, avoid garbarge propagation
    std::fill(rNoise[n], rNoise[n]+NR, REAL_SCALAR{0.0});
  }

  for(int n=0;n<N;++n)
  {
    fftw_plan nPlan = (NZ>1 ? fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(nkCorr[n]), rNoise[n], FFTW_ESTIMATE)
      //: fftw_plan_dft_c2r_2d(NX, NY, reinterpret_cast<fftw_complex*>(nkCorr[n]), rNoise[n], FFTW_ESTIMATE));
      : fftw_plan_dft_c2r_2d(NY, NX, reinterpret_cast<fftw_complex*>(nkCorr[n]), rNoise[n], FFTW_ESTIMATE));
    fftw_execute(nPlan);
  }

  return rNoise;
}
#endif

template <int N>
void GlidePlaneNoiseBase<N>::computeRealNoise()
{
  std::cout<<greenBoldColor<<"Computing "+tag<<defaultColor<<std::flush;
  const auto t0= std::chrono::system_clock::now();

#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_

  std::default_random_engine generator(seed);
  std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);

  std::vector<REAL_SCALAR*> rNoise(realNoise(noisyKCorrelations(originalKCorrelations(),generator,distribution)));

  // Store Real noise
  this->reserve(NX*NY);
  const int k=0;
  for(int i=0;i<NX;i++)
  {
    for(int j=0;j<NY;j++)
    {
      const int ind = NY*NZ*i + j*NZ + k;
      Eigen::Matrix<REAL_SCALAR,1,N> temp;
      for(int n=0;n<N;++n)
      {
        temp(n)=rNoise[n][ind];
      }
      this->push_back(NoiseTraits<N>::fromMatrix(temp));
    }
  }

#else
  this->resize(NX*NY,NoiseTraits<N>::Zero());
#endif
  std::cout<<greenColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
}


template <int N>
std::vector<std::vector<double>> GlidePlaneNoiseBase<N>::averageNoiseCorrelation(const int& reps) const
{
#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_

  std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);

  // store original correlations once
  std::vector<COMPLEX*> okCorr(originalKCorrelations());

  // init storage for average correlation
  std::vector<COMPLEX*> avgkCorr;
  for(int n=0; n<N; ++n) 
  {
    avgkCorr.push_back((COMPLEX*) fftw_malloc(sizeof(COMPLEX) * NK));
    // Zero-initialize explicitly
    std::fill(avgkCorr[n], avgkCorr[n]+NK, COMPLEX{0.0, 0.0});
  }

  for(int k=0;k<reps;++k)
  {
    unsigned int seed = k+1;
    std::default_random_engine generator(seed);
    // nkCorr is a sampled (randomized) noise in k space
    const std::vector<COMPLEX*> nkCorr = noisyKCorrelations(okCorr,generator,distribution);
    for(int n=0;n<N;++n)
    {
      for(int ind=0;ind<NK;++ind)
      {
        // accumulate power spectrum and ensemble average calculated correaltion
        const COMPLEX I = nkCorr[n][ind]*std::conj(nkCorr[n][ind]);
        avgkCorr[n][ind] += I / static_cast<double>(reps);
      }
    }
  }

  for(int n=0;n<N;++n)
  {
    for(int ind=0;ind<NK;++ind)
    {
      // Normalize the output after iFFT
      avgkCorr[n][ind] /= static_cast<double>(NR);
    }
  }

  // do IDFT from correlation in k-space to real space
  //std::vector<REAL_SCALAR*> avgrCorr(realNoise(avgkCorr));
  std::vector<REAL_SCALAR*> avgrCorr;
  for(int n=0;n<N;++n)
  {
    avgrCorr.push_back((REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR));
    // explicitly initialize to zeros, avoid garbarge propagation
    std::fill(avgrCorr[n], avgrCorr[n]+NR, REAL_SCALAR{0.0});
  }

  for(int n=0;n<N;++n)
  {
    fftw_plan fplan = (NZ>1 ? fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(avgkCorr[n]), avgrCorr[n], FFTW_ESTIMATE)
      //: fftw_plan_dft_c2r_2d(NX, NY, reinterpret_cast<fftw_complex*>(nkCorr[n]), rNoise[n], FFTW_ESTIMATE));
      : fftw_plan_dft_c2r_2d(NY, NX, reinterpret_cast<fftw_complex*>(avgkCorr[n]), avgrCorr[n], FFTW_ESTIMATE));
    fftw_execute(fplan);
  }

  // initialize a flattened array for python API
  std::vector<std::vector<double>> avgrCorrFlat;

  const int k=0;
  for(int n=0;n<N;++n)
  {
    std::vector<double> temp;
    for(int i=0;i<NX;i++)
    {
      for(int j=0;j<NY;j++)
      {
        const int ind = NY*NZ*i + j*NZ + k;
        temp.push_back(avgrCorr[n][ind]);
      }
    }
    avgrCorrFlat.push_back(temp);
  }

  return avgrCorrFlat;

#endif
}

template <int N>
std::vector<std::vector<double>> GlidePlaneNoiseBase<N>::sampleAverageNoise(const int& reps) const
{
#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_

  std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);

  // store original correlations once
  std::vector<COMPLEX*> okCorr(originalKCorrelations());

  // init storage for sampled noises
  std::vector<std::vector<double>> rNoises;

  for(int k=0;k<reps;++k)
  {
    unsigned int seed = k+1;
    std::default_random_engine generator(seed);
    for(int n=0;n<N;++n)
    {
      // nkCorr is a sampled (randomized) noise in k space
      const std::vector<COMPLEX*> nkNoise = noisyKCorrelations(okCorr,generator,distribution);
      // sampled noise to real space
      const std::vector<REAL_SCALAR*> nrNoise(realNoise(nkNoise));

      const int k=0;
      std::vector<double> temp;
      for(int i=0;i<NX;i++)
      {
        for(int j=0;j<NY;j++)
        {
          const int ind = NY*NZ*i + j*NZ + k;
          temp.push_back(nrNoise[n][ind]);
        }
      }
      rNoises.push_back(temp);
    }
  }

  return rNoises;

#endif
}

template <int N>
//void GlidePlaneNoiseBase<N>::computeRealNoiseStatistics(const PolycrystallineMaterialBase& mat) const
void GlidePlaneNoiseBase<N>::computeRealNoiseStatistics() const
{
  // Compute Statistics
  NoiseType ave(NoiseTraits<N>::Zero());
  for(const auto& valArr: noiseVector())
  {
    ave+=valArr;
  }
  ave/=noiseVector().size();

  NoiseType var(NoiseTraits<N>::Zero());
  for(const auto& valArr: noiseVector())
  {
    var+=NoiseTraits<N>::squared(valArr-ave);
  }
  var/=noiseVector().size();

  std::cout<<"gridSize= "<<gridSize.transpose()<<std::endl;
  std::cout<<"gridSpacing= "<<gridSpacing.transpose()<<std::endl;

  // these are incorrect for MD stacking fault noises (J/m2)
  //std::cout<<"noiseAverage="<<ave*mat.mu_SI<<" [Pa]"<<std::endl;
  //std::cout<<"noiseVariance="<<var*std::pow(mat.mu_SI,2)<<" [Pa^2]"<<std::endl;
  std::cout<<"noiseAverage="<<ave<<std::endl;
  std::cout<<"noiseVariance="<<var<<std::endl;
}

template <int N>
void GlidePlaneNoiseBase<N>::write_field_slice(const PolycrystallineMaterialBase& mat) const
{
  const std::string matPath = std::filesystem::path(mat.materialFile).parent_path().string();
  const double DX = gridSpacing(0);
  const double DY = gridSpacing(1);
  const double DZ = gridSpacing(2);

  std::vector<std::vector<double>> noiseArrays(N, std::vector<double>(NR, 0.0));
  //const int k=0;
  for(int i=0;i<NX;i++)
  {
    for(int j=0;j<NY;j++)
    {
      // input array is 2D, so index calculation
      // must be changed
      //const int ind = NY*NZ*i + j*NZ + k;
      const int ind = NY*i + j;
      const NoiseType& row = noiseVector()[ind];
      for(int n=0;n<N;++n)
      {
        if constexpr (N==1)
        {
          noiseArrays[n][ind] = row;
        }
        else
        {
          noiseArrays[n][ind] = row(n);
        }
      }
    }
  }

  for(int n=0;n<N;++n)
  {
    std::string fname="noise_patch_"+std::to_string(n)+".vtk";
    std::string fnamePath = matPath+"/"+fname;
    //FILE *OutFile=fopen(fname.c_str(),"w");
    FILE *OutFile=fopen(fnamePath.c_str(),"w");

    fprintf(OutFile,"# vtk DataFile Version 2.0\n");
    fprintf(OutFile,"iter %d\n",0);
    //fprintf(OutFile,"ASCII\n");
    fprintf(OutFile,"BINARY\n");
    fprintf(OutFile,"DATASET STRUCTURED_POINTS\n");
    fprintf(OutFile,"ORIGIN \t %f %f %f\n",0.,0.,0.);
    fprintf(OutFile,"SPACING \t %f %f %f\n", DX, DY, DZ);
    fprintf(OutFile,"DIMENSIONS \t %d %d %d\n", NX, NY, 1);
    fprintf(OutFile,"POINT_DATA \t %d\n",NX*NY);
    fprintf(OutFile,"SCALARS \t volume_scalars double 1\n");
    fprintf(OutFile,"LOOKUP_TABLE \t default\n");

    for(int i=0;i<NX;i++)
    {
      for(int j=0;j<NY;j++)
      {
        const int ind = NY*i + j;
        //const double temp=noiseArrays[n][ind]; // for ascii
        const double temp=NoiseTraitsBase::ReverseDouble(noiseArrays[n][ind]);
        fwrite(&temp, sizeof(double), 1, OutFile);
        //fprintf(OutFile, "%.18e ", temp); // for ascii
      }
    }
    fclose(OutFile);
  }
}

template struct GlidePlaneNoiseBase<1>;
template struct GlidePlaneNoiseBase<2>;

}
#endif

