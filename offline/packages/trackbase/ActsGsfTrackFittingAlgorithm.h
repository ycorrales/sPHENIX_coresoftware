#include "ActsTrackFittingAlgorithm.h"

#include <Acts/Definitions/TrackParametrization.hpp>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/TrackFitting/BetheHeitlerApprox.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#pragma GCC diagnostic push // needed for local Act compilation
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <Acts/TrackFitting/GaussianSumFitter.hpp>
#pragma GCC diagnostic pop
#include <Acts/TrackFitting/GsfMixtureReduction.hpp>
#include <Acts/Utilities/Helpers.hpp>

enum class MixtureReductionAlgorithm
{
  weightCut,
  KLDistance
};

namespace
{
  using BetheHeitlerApprox = Acts::AtlasBetheHeitlerApprox<6, 5>;
  using MultiStepper = Acts::MultiEigenStepperLoop<>;
  using Propagator = Acts::Propagator<MultiStepper, Acts::Navigator>;
  using DirectPropagator = Acts::Propagator<MultiStepper, Acts::DirectNavigator>;

  using Fitter =
      Acts::GaussianSumFitter<Propagator,
                              BetheHeitlerApprox,
                              Acts::VectorMultiTrajectory>;
  using DirectFitter =
      Acts::GaussianSumFitter<DirectPropagator,
                              BetheHeitlerApprox,
                              Acts::VectorMultiTrajectory>;
  using TrackContainer =
      Acts::TrackContainer<Acts::VectorTrackContainer,
                           Acts::VectorMultiTrajectory, std::shared_ptr>;



  struct GsfFitterFunctionImpl
    : public ActsTrackFittingAlgorithm::TrackFitterFunction
  {
    Fitter fitter;

    Acts::GainMatrixUpdater updater;

    std::size_t maxComponents = 0;
    double weightCutoff = 0;
    bool abortOnError = false;
    bool disableAllMaterialHandling = false;
    MixtureReductionAlgorithm reductionAlg =
        MixtureReductionAlgorithm::KLDistance;
    Acts::ComponentMergeMethod mergeMethod =
        Acts::ComponentMergeMethod::eMaxWeight;

    ActsSourceLink::SurfaceAccessor m_slSurfaceAccessor;

    GsfFitterFunctionImpl(Fitter&& f,
                          const Acts::TrackingGeometry& trkGeo)
      : fitter(std::move(f))
      , m_slSurfaceAccessor{trkGeo}
    {
    }

    template <typename calibrator_t>
    auto makeGsfOptions(
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
        const calibrator_t& calibrator)
        const
    {
      Acts::GsfExtensions<Acts::VectorMultiTrajectory> extensions;
      // cppcheck-suppress constStatement
      extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&updater);

      Acts::GsfOptions<Acts::VectorMultiTrajectory> gsfOptions{
          options.geoContext,
          options.magFieldContext,
          options.calibrationContext,
          extensions,
          options.propOptions,
          &(*options.referenceSurface),
          maxComponents,
          weightCutoff,
          abortOnError,
          disableAllMaterialHandling};
      gsfOptions.componentMergeMethod = mergeMethod;
      gsfOptions.extensions.calibrator.connect<&calibrator_t::calibrate>(
          &calibrator);
      gsfOptions.extensions.surfaceAccessor.connect<&ActsSourceLink::SurfaceAccessor::operator()>(&m_slSurfaceAccessor);
      switch (reductionAlg)
      {
      case MixtureReductionAlgorithm::weightCut:
      {
        gsfOptions.extensions.mixtureReducer
            .connect<&Acts::reduceMixtureLargestWeights>();
      }
      break;
      case MixtureReductionAlgorithm::KLDistance:
      {
        gsfOptions.extensions.mixtureReducer
            .connect<&Acts::reduceMixtureWithKLDistance>();
      }
      break;
      }

      return gsfOptions;
    }

    ActsTrackFittingAlgorithm::TrackFitterResult operator()(
        const std::vector<Acts::SourceLink>& sourceLinks,
        const ActsTrackFittingAlgorithm::TrackParameters& initialParameters,
        const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
        const CalibratorAdapter& calibrator,
        ActsTrackFittingAlgorithm::TrackContainer& tracks) const override
    {
      const auto gsfOptions = makeGsfOptions(options, calibrator);
      using namespace Acts::GsfConstants;
      if (not tracks.hasColumn(Acts::hashString(kFinalMultiComponentStateColumn)))
      {
        std::string key(kFinalMultiComponentStateColumn);
        tracks.template addColumn<FinalMultiComponentState>(key);
      }

      return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                        gsfOptions, tracks);
    }
  };

}  // namespace

// Have a separate class befriend the main class to ensure that GSF specific
// track fitting headers only stay here to avoid library clashes
class ActsGsfTrackFittingAlgorithm
{
 public:
  friend class ActsTrackFittingAlgorithm;

  std::shared_ptr<ActsTrackFittingAlgorithm::TrackFitterFunction>
  makeGsfFitterFunction(
      const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      BetheHeitlerApprox betheHeitlerApprox, std::size_t maxComponents,
      double weightCutoff,
      MixtureReductionAlgorithm finalReductionMethod, bool abortOnError,
      bool disableAllMaterialHandling, const Acts::Logger& logger = *Acts::getDefaultLogger("GSF", Acts::Logging::FATAL));
};
