#pragma once

#include "nmopt/application/chapter6.hpp"
#include "nmopt/application/dealii/centroid_split_simplex_mesh.hpp"
#include "nmopt/application/dealii/scalar_function.hpp"
#include "nmopt/application/runner.hpp"
#include "nmopt/compiler/v1/dealii_compiler.hpp"
#include "nmopt/experiment/reduced_envelope.hpp"

#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/tensor_function_parser.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::application::chapter6::dealii
{
  using Backend = nmopt::dealii_backend::SerialBackend;
  using Envelope =
    nmopt::experiment::ReducedSearchExperimentEnvelopeT<Backend>;
  using Evidence =
    nmopt::application::benchmark::BenchmarkExecutionEvidenceT<Envelope>;

  template <int dim>
  class B2DesiredStateFunction final : public ::dealii::Function<dim>
  {
  public:
    static_assert(dim >= 2,
                  "The Chapter 6 B2 target requires at least two coordinates");

    explicit B2DesiredStateFunction(
      const ScalarFunctionDefinition &target_definition)
      : target_definition_(target_definition)
      , target_(::nmopt::application::dealii_support::make_scalar_function<dim>(
                 target_definition_, "B2 target"))
    {}

    double
    value(const ::dealii::Point<dim> &point,
          const unsigned int          component = 0) const override
    {
      return target_->value(point, component);
    }

    const ScalarFunctionDefinition &
    target_definition() const
    {
      return target_definition_;
    }

  private:
    ScalarFunctionDefinition                 target_definition_;
    std::unique_ptr<::dealii::Function<dim>> target_;
  };

  namespace detail
  {
    template <int dim>
    std::unique_ptr<::dealii::TensorFunctionParser<1, dim>>
    make_rank_one_vector_function(
      const RankOneVectorFunctionDefinition &definition,
      const std::string_view                  description)
    {
      validate_rank_one_vector_function_definition(definition, description);
      std::ostringstream variable_names;
      for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
        {
          if (coordinate != 0)
            variable_names << ',';
          variable_names << 'x' << coordinate;
        }

      auto function =
        std::make_unique<::dealii::TensorFunctionParser<1, dim>>();
      try
        {
          function->initialize(
            variable_names.str(),
            definition.expression,
            {{"e", ::dealii::numbers::E}, {"pi", ::dealii::numbers::PI}});
          ::dealii::Point<dim> validation_point;
          for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
            validation_point[coordinate] = 0.371 + 0.113 * coordinate;
          (void)function->value(validation_point);
        }
      catch (const std::exception &exception)
        {
          throw std::invalid_argument(std::string(description) +
                                      " expression is invalid: " +
                                      exception.what());
        }
      return function;
    }

    template <typename Backend>
    contract::PrimalBlockT<Backend>
    make_b2_uniform_control(const contract::LayoutPtr &layout,
                            const double              value)
    {
      if (!std::isfinite(value))
        throw std::invalid_argument(
          "B2 initial independent control value must be finite");
      std::vector<typename Backend::Vector> blocks;
      blocks.reserve(layout->n_blocks());
      for (std::size_t block = 0; block < layout->n_blocks(); ++block)
        {
          auto values = Backend::zeros(layout->dimension(block));
          for (std::size_t index = 0; index < layout->dimension(block); ++index)
            Backend::set_value(values, index, value);
          blocks.push_back(std::move(values));
        }
      return contract::PrimalBlockT<Backend>(layout, std::move(blocks));
    }
  } // namespace detail

  template <int dim>
  struct B2RuntimeDataT
  {
    B2RuntimeDataT(
      const ::dealii::Function<dim> &             forcing_function,
      const ::dealii::Function<dim> &             desired_state_function,
      const ::dealii::Function<dim> &             fixed_temperature_function,
      const ::dealii::TensorFunction<1, dim> &transport_function)
      : forcing(forcing_function)
      , desired_state(desired_state_function)
      , fixed_temperature(fixed_temperature_function)
      , conservative_transport(transport_function)
    {}

    const ::dealii::Function<dim> &             forcing;
    const ::dealii::Function<dim> &             desired_state;
    const ::dealii::Function<dim> &             fixed_temperature;
    const ::dealii::TensorFunction<1, dim> &conservative_transport;
  };

  template <int dim>
  class B2ManufacturedDataT final
  {
  public:
    explicit B2ManufacturedDataT(
      ScalarFunctionDefinition observation_definition =
        b2_manufactured_wings_observation(),
      ScalarFunctionDefinition fixed_dirichlet_data =
        b2_manufactured_fixed_temperature(),
      ScalarFunctionDefinition forcing_definition =
        b2_manufactured_zero_forcing(),
      ScalarFunctionDefinition target_definition =
        b2_manufactured_constant_target(),
      RankOneVectorFunctionDefinition conservative_transport =
        b2_manufactured_graetz_transport())
      : observation_definition_(std::move(observation_definition))
      , fixed_dirichlet_data_(std::move(fixed_dirichlet_data))
      , fixed_temperature_(
          ::nmopt::application::dealii_support::make_scalar_function<dim>(
            fixed_dirichlet_data_, "B2 fixed Dirichlet data"))
      , forcing_definition_(std::move(forcing_definition))
      , forcing_(::nmopt::application::dealii_support::make_scalar_function<dim>(
          forcing_definition_, "B2 forcing"))
      , desired_state_(target_definition)
      , conservative_transport_definition_(std::move(conservative_transport))
      , conservative_transport_(detail::make_rank_one_vector_function<dim>(
          conservative_transport_definition_, "B2 conservative transport"))
      {}

    const ScalarFunctionDefinition &
    observation_definition() const
    {
      return observation_definition_;
    }

    B2RuntimeDataT<dim>
    runtime_data() const
    {
      return {*forcing_,
              desired_state_,
              *fixed_temperature_,
              *conservative_transport_};
    }

    const ScalarFunctionDefinition &
    fixed_dirichlet_data() const
    {
      return fixed_dirichlet_data_;
    }

    const ScalarFunctionDefinition &
    forcing_definition() const
    {
      return forcing_definition_;
    }

    const ScalarFunctionDefinition &
    target_definition() const
    {
      return desired_state_.target_definition();
    }

    const RankOneVectorFunctionDefinition &
    conservative_transport() const
    {
      return conservative_transport_definition_;
    }

  private:
    ScalarFunctionDefinition                 observation_definition_;
    ScalarFunctionDefinition                 fixed_dirichlet_data_;
    std::unique_ptr<::dealii::Function<dim>> fixed_temperature_;
    ScalarFunctionDefinition                 forcing_definition_;
    std::unique_ptr<::dealii::Function<dim>> forcing_;
    B2DesiredStateFunction<dim>             desired_state_;
    RankOneVectorFunctionDefinition          conservative_transport_definition_;
    std::unique_ptr<::dealii::TensorFunctionParser<1, dim>>
      conservative_transport_;
  };

  template <int dim>
  B2RuntimeDataT<dim>
  make_b2_manufactured_runtime_data(
    const B2Scenario &scenario,
    const B2ManufacturedDataT<dim> &data)
  {
    validate_b2(scenario);
    if (scenario.compile.mesh.dimension != dim)
      throw std::invalid_argument(
        "B2 runtime data dimension does not match the scenario mesh");
    if (data.observation_definition() !=
        selected_scalar_function_definition(
          scenario.problem.observation_region_catalog,
          "B2 observation catalog"))
      throw std::invalid_argument(
        "B2 manufactured observation definition does not match the scenario");
    if (data.fixed_dirichlet_data() !=
        scenario.problem.fixed_dirichlet_data)
      throw std::invalid_argument(
        "B2 manufactured fixed Dirichlet data does not match the scenario");
    if (data.forcing_definition() != scenario.problem.forcing)
      throw std::invalid_argument(
        "B2 manufactured forcing does not match the scenario");
    if (data.target_definition() !=
        b2_target_definition(scenario.problem.target_catalog))
      throw std::invalid_argument(
        "B2 manufactured target definition does not match the scenario");
    if (data.conservative_transport() !=
        scenario.problem.conservative_transport)
      throw std::invalid_argument(
        "B2 manufactured conservative transport does not match the scenario");
    return data.runtime_data();
  }

  template <int dim>
  compiler::v1::DealiiDataBindings<dim>
  make_b2_data_bindings(const B2ProblemParameters & parameters,
                        const B2RuntimeDataT<dim> &runtime)
  {
    validate_runtime_data(parameters.data);
    if (parameters.data.fixed_dirichlet_data_provenance.empty() ||
        parameters.data.conservative_transport_provenance.empty())
      throw std::invalid_argument(
        "B2 runtime bindings need fixed-data and transport provenance");

    return {runtime.forcing,
            runtime.desired_state,
            std::optional<double>{parameters.data.diffusion},
            parameters.data.reaction,
            parameters.data.regularisation_weight,
            {parameters.data.forcing_provenance,
             parameters.data.desired_state_provenance,
             parameters.data.fixed_dirichlet_data_provenance},
            std::cref(runtime.fixed_temperature),
            std::nullopt,
            std::nullopt,
            compiler::v1::DealiiConservativeTransportDataBindings<dim>{
              runtime.conservative_transport,
              {parameters.data.conservative_transport_provenance}}};
  }

  inline void
  validate_b2_dealii_boundary_options(const B2BoundaryOptions &options)
  {
    validate_b2_boundary_options(options);
    if (options.fixed_id == ::dealii::numbers::invalid_boundary_id ||
        options.control_id == ::dealii::numbers::invalid_boundary_id ||
        options.outflow_id == ::dealii::numbers::invalid_boundary_id)
      throw std::invalid_argument(
        "B2 boundary IDs must be valid deal.II boundary IDs");
  }

  template <int dim>
  std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>>
  make_b2_compilation_session(const B2Scenario &scenario)
  {
    validate_b2(scenario);
    if (dim != 2 || scenario.compile.mesh.dimension != dim)
      throw std::invalid_argument(
        "B2 compilation sessions require a two-dimensional matching mesh");
    if (scenario.problem.recipe.observed_material_id == 0)
      throw std::invalid_argument(
        "B2 needs a nonzero observed material id");
    validate_b2_dealii_boundary_options(scenario.problem.boundary);

    const auto &observation_definition = selected_scalar_function_definition(
      scenario.problem.observation_region_catalog, "B2 observation catalog");
    std::unique_ptr<::dealii::Function<dim>> observation_indicator;
    try
      {
        observation_indicator =
          ::nmopt::application::dealii_support::make_scalar_function<dim>(
            observation_definition, "B2 observation region");
      }
    catch (const std::exception &exception)
      {
        throw std::invalid_argument(
          std::string("B2 observation region expression is invalid: ") +
          exception.what());
      }

    auto mesh = std::make_unique<::dealii::Triangulation<dim>>();
    ::dealii::Point<dim> lower;
    ::dealii::Point<dim> upper;
    for (unsigned int coordinate = 0; coordinate < dim; ++coordinate)
      {
        lower[coordinate] = scenario.compile.mesh.lower[coordinate];
        upper[coordinate] = scenario.compile.mesh.upper[coordinate];
      }
    switch (scenario.compile.mesh.generation)
      {
        case MeshGeneration::framework_native:
          ::dealii::GridGenerator::hyper_rectangle(*mesh, lower, upper);
          mesh->refine_global(scenario.compile.mesh.refinement);
          break;
        case MeshGeneration::structured_simplex:
          if constexpr (dim == 2)
            ::dealii::GridGenerator::
              subdivided_hyper_rectangle_with_simplices(
                *mesh,
                scenario.compile.mesh.axis_subdivisions,
                lower,
                upper,
                false);
          else
            throw std::invalid_argument(
              "B2 simplex meshes are implemented only in two dimensions");
          break;
        case MeshGeneration::centroid_split_simplex:
          if constexpr (dim == 2)
            {
              ::dealii::Triangulation<2> base_mesh;
              ::dealii::GridGenerator::
                subdivided_hyper_rectangle_with_simplices(
                  base_mesh,
                  scenario.compile.mesh.axis_subdivisions,
                  lower,
                  upper,
                  false);
              mesh = detail::make_centroid_split_simplex_mesh(
                base_mesh,
                scenario.compile.mesh.centroid_splits,
                scenario.compile.mesh.selection_seed);
            }
          else
            throw std::invalid_argument(
              "B2 simplex meshes are implemented only in two dimensions");
          break;
      }

    const auto is_close = [](const double first, const double second) {
      return std::abs(first - second) < 1.0e-12;
    };
    for (auto cell = mesh->begin_active(); cell != mesh->end(); ++cell)
      {
        const auto center = cell->center();
        const double indicator = observation_indicator->value(center);
        if (!std::isfinite(indicator))
          throw std::invalid_argument(
            "B2 observation region indicator returned a non-finite value");
        const bool observed = indicator > 0.0;
        cell->set_material_id(
          observed ? scenario.problem.recipe.observed_material_id : 0);
        for (unsigned int face = 0; face < cell->n_faces(); ++face)
          if (cell->face(face)->at_boundary())
            {
              const auto face_center = cell->face(face)->center();
              const bool left_boundary = is_close(face_center[0], lower[0]);
              const bool right_boundary = is_close(face_center[0], upper[0]);
              const bool horizontal_boundary =
                is_close(face_center[1], lower[1]) ||
                is_close(face_center[1], upper[1]);
              if (!left_boundary && !right_boundary && !horizontal_boundary)
                throw std::invalid_argument(
                  "B2 mesh contains an unclassified exterior boundary face");
              const auto boundary_id =
                left_boundary
                  ? scenario.problem.boundary.fixed_id
                  : right_boundary
                    ? scenario.problem.boundary.outflow_id
                    : face_center[0] <
                        scenario.problem.boundary.upstream_transition
                      ? scenario.problem.boundary.fixed_id
                      : scenario.problem.boundary.control_id;
              cell->face(face)->set_boundary_id(boundary_id);
            }
      }
    return std::make_shared<compiler::v1::DealiiCompilationSession<dim>>(
      std::move(mesh), scenario.compile.mesh.mesh_provenance);
  }

  inline compiler::v1::DealiiDiscretisationPolicy
  make_b2_discretisation_policy(const CompileOptions &options)
  {
    validate_common_compile_options(options);
    if (options.execution != ExecutionSelection::assembled)
      throw std::invalid_argument(
        "B2 runtime execution is limited to assembled operators");

    compiler::v1::DealiiDiscretisationPolicy policy;
    policy.state_degree = options.state_degree;
    policy.execution =
      compiler::v1::DealiiDiscretisationPolicy::Execution::assembled;
    if (!options.volume_observation)
      throw std::invalid_argument(
        "B2 runtime execution needs volume-observation options");
    const auto target_realisation = [&]() {
      switch (options.volume_observation->target_realisation)
        {
          case VolumeObservationTargetRealisation::analytic_quadrature:
            return compiler::v1::VolumeObservationTargetRealisation::
              analytic_quadrature;
          case VolumeObservationTargetRealisation::state_fe_interpolation:
            return compiler::v1::VolumeObservationTargetRealisation::
              state_fe_interpolation;
        }
      throw std::invalid_argument(
        "B2 runtime execution received an unknown observation target realisation");
    }();
    policy.volume_observation =
      compiler::v1::VolumeObservationDiscretisationPolicy{
        options.volume_observation->quadrature_order,
        target_realisation};
    policy.state_solve = {options.state_solve.maximum_iterations,
                          options.state_solve.relative_tolerance,
                          options.state_solve.absolute_tolerance};
    policy.adjoint_solve = {options.adjoint_solve.maximum_iterations,
                            options.adjoint_solve.relative_tolerance,
                            options.adjoint_solve.absolute_tolerance};
    policy.control_metric_solve = {
      options.control_metric_solve.maximum_iterations,
      options.control_metric_solve.relative_tolerance,
      options.control_metric_solve.absolute_tolerance};
    return policy;
  }

  inline compiler::v1::CompilationProduct
  make_b2_compilation_product(const CompileOptions &options)
  {
    if (options.product != ProductSelection::reduced_dto)
      throw std::invalid_argument(
        "B2 runtime execution is limited to the reduced DTO product");
    return compiler::v1::CompilationProduct::reduced_dto;
  }

  inline std::string
  b2_number(const double value)
  {
    std::ostringstream output;
    output << std::setprecision(17) << value;
    return output.str();
  }

  inline std::string
  b2_coordinates(const std::vector<double> &coordinates)
  {
    std::ostringstream output;
    output << std::setprecision(17);
    for (std::size_t index = 0; index < coordinates.size(); ++index)
      {
        if (index != 0)
          output << ',';
        output << coordinates[index];
      }
    return output.str();
  }

  inline double
  b2_relative_reduction(const double initial, const double final)
  {
    return (initial - final) /
           std::max(std::abs(initial), std::numeric_limits<double>::epsilon());
  }

  inline std::string
  b2_history(const std::vector<double> &history)
  {
    std::ostringstream output;
    output << std::setprecision(17);
    for (std::size_t index = 0; index < history.size(); ++index)
      {
        if (index != 0)
          output << ',';
        output << history[index];
      }
    return output.str();
  }

  template <typename Backend>
  double
  b2_block_values_l2_norm(const contract::BlockValuesT<Backend> &values)
  {
    double squared_norm = 0.0;
    for (std::size_t block = 0; block < values.n_blocks(); ++block)
      squared_norm += Backend::dot(values.block(block), values.block(block));
    return std::sqrt(std::max(0.0, squared_norm));
  }

  template <typename Backend>
  std::pair<double, double>
  b2_block_extrema(const contract::BlockValuesT<Backend> &values,
                   const std::size_t                      block = 0)
  {
    contract::require(block < values.n_blocks(),
                      "B2 extrema requested an absent block");
    contract::require(values.block(block).size() > 0,
                      "B2 extrema need a nonempty coefficient block");
    double minimum = values.block(block)[0];
    double maximum = minimum;
    for (std::size_t index = 1; index < values.block(block).size(); ++index)
      {
        minimum = std::min(minimum, values.block(block)[index]);
        maximum = std::max(maximum, values.block(block)[index]);
      }
    contract::require(std::isfinite(minimum) && std::isfinite(maximum),
                      "B2 extrema contain a non-finite coefficient");
    return {minimum, maximum};
  }

  template <typename Backend>
  struct B2DerivativeEvidenceT
  {
    double residual_jvp_error;
    double residual_vjp_error;
    double reduced_gradient_finite_difference_error;
    double reduced_taylor_step;
    double reduced_taylor_error;
    double reduced_taylor_order;
    double initial_objective;
    double initial_gradient_norm;
  };

  template <typename Backend>
  contract::PrimalBlockT<Backend>
  b2_shifted_point(const contract::PrimalBlockT<Backend> &point,
                   const contract::PrimalBlockT<Backend> &direction,
                   const double                            step)
  {
    contract::require(point.layout()->compatible_with(*direction.layout()),
                      "B2 evidence uses incompatible full-point directions");
    contract::PrimalBlockT<Backend> shifted = point;
    for (std::size_t block = 0; block < shifted.n_blocks(); ++block)
      shifted.add_scaled_block(block, step, direction.block(block));
    return shifted;
  }

  template <typename Backend>
  contract::PrimalBlockT<Backend>
  b2_control_shifted(const contract::PrimalBlockT<Backend> &control,
                     const contract::PrimalBlockT<Backend> &direction,
                     const double                            step)
  {
    contract::require(control.layout()->compatible_with(*direction.layout()),
                      "B2 evidence uses incompatible control directions");
    contract::PrimalBlockT<Backend> shifted = control;
    shifted.add_scaled_block(0, step, direction.block(0));
    return shifted;
  }

  template <typename Backend>
  B2DerivativeEvidenceT<Backend>
  make_b2_derivative_evidence(
    const compiler::v1::CompiledProblemT<Backend> &problem,
    const contract::ReducedDTOT<Backend> &         reduced,
    const contract::StateControlPartitionT<Backend> &partition,
    const contract::PrimalBlockT<Backend> &        initial_control)
  {
    using Primal = contract::PrimalBlockT<Backend>;
    using Covector = contract::CovectorBlockT<Backend>;
    using Vector = typename Backend::Vector;

    const auto &model = problem.executable_model();
    const Primal zero_point = Primal::zeros(model.variable_layout());
    Vector state_direction_values(model.variable_layout()->dimension(0));
    Vector control_direction_values(model.variable_layout()->dimension(1));
    for (std::size_t index = 0;
         index < state_direction_values.size();
         ++index)
      state_direction_values[index] =
        (index % 2 == 0 ? 0.02 : -0.015) * static_cast<double>(index + 1);
    for (std::size_t index = 0;
         index < control_direction_values.size();
         ++index)
      control_direction_values[index] =
        (index % 2 == 0 ? -0.03 : 0.025) * static_cast<double>(index + 1);
    const Primal full_direction(model.variable_layout(),
                                {std::move(state_direction_values),
                                 std::move(control_direction_values)});

    Vector test_seed_values(model.test_layout()->dimension(0));
    for (std::size_t index = 0; index < test_seed_values.size(); ++index)
      test_seed_values[index] =
        (index % 3 == 0 ? 0.01 : -0.02) * static_cast<double>(index + 1);
    const Primal test_seed(model.test_layout(), {std::move(test_seed_values)});

    constexpr double residual_step = 1.0e-6;
    const Covector residual_zero = model.residual(zero_point);
    const Covector residual_plus = model.residual(
      b2_shifted_point(zero_point, full_direction, residual_step));
    const Covector residual_jvp =
      model.residual_jvp(zero_point, full_direction);
    Covector residual_finite_difference = residual_plus;
    residual_finite_difference.add_scaled_block(
      0, -1.0, residual_zero.block(0));
    residual_finite_difference.scale_block(0, 1.0 / residual_step);
    residual_finite_difference.add_scaled_block(
      0, -1.0, residual_jvp.block(0));
    const double residual_jvp_error =
      b2_block_values_l2_norm(residual_finite_difference);

    const Covector residual_transpose =
      model.residual_vjp(zero_point, test_seed);
    const double residual_vjp_error = std::abs(
      contract::pair(residual_jvp, test_seed) -
      contract::pair(residual_transpose, full_direction));

    Vector reduced_direction_values(partition.control_layout()->dimension(0));
    for (std::size_t index = 0; index < reduced_direction_values.size(); ++index)
      reduced_direction_values[index] =
        (index % 2 == 0 ? 0.03 : -0.02) * static_cast<double>(index + 1);
    const Primal reduced_direction(partition.control_layout(),
                                   {std::move(reduced_direction_values)});

    constexpr double taylor_step = 1.0e-4;
    const auto base_evaluation = reduced.evaluate(initial_control);
    const auto plus_evaluation = reduced.evaluate(b2_control_shifted(
      initial_control, reduced_direction, taylor_step));
    const auto minus_evaluation = reduced.evaluate(b2_control_shifted(
      initial_control, reduced_direction, -taylor_step));
    const auto half_plus_evaluation = reduced.evaluate(b2_control_shifted(
      initial_control, reduced_direction, 0.5 * taylor_step));
    const double directional_derivative = contract::pair(
      base_evaluation.reduced_derivative, reduced_direction);
    const double central_derivative =
      (plus_evaluation.objective_value - minus_evaluation.objective_value) /
      (2.0 * taylor_step);
    const double reduced_gradient_finite_difference_error =
      std::abs(central_derivative - directional_derivative);
    const double taylor_error = std::abs(
      plus_evaluation.objective_value - base_evaluation.objective_value -
      taylor_step * directional_derivative);
    const double half_taylor_error = std::abs(
      half_plus_evaluation.objective_value - base_evaluation.objective_value -
      0.5 * taylor_step * directional_derivative);
    const double taylor_order =
      taylor_error > std::numeric_limits<double>::epsilon() &&
          half_taylor_error > std::numeric_limits<double>::epsilon()
        ? std::log(taylor_error / half_taylor_error) / std::log(2.0)
        : std::numeric_limits<double>::infinity();

    contract::require(std::isfinite(residual_jvp_error) &&
                        std::isfinite(residual_vjp_error) &&
                        std::isfinite(reduced_gradient_finite_difference_error) &&
                        std::isfinite(taylor_error) &&
                        std::isfinite(taylor_order),
                      "B2 derivative evidence produced a non-finite value");
    contract::require(residual_jvp_error <= 1.0e-7,
                      "B2 residual JVP finite-difference evidence failed");
    contract::require(residual_vjp_error <= 1.0e-7,
                      "B2 residual VJP transpose evidence failed");
    contract::require(reduced_gradient_finite_difference_error <= 1.0e-7,
                      "B2 reduced gradient finite-difference evidence failed");
    contract::require(taylor_order >= 1.5,
                      "B2 reduced Taylor evidence failed");

    return {residual_jvp_error,
            residual_vjp_error,
            reduced_gradient_finite_difference_error,
            taylor_step,
            taylor_error,
            taylor_order,
            base_evaluation.objective_value,
            b2_block_values_l2_norm(base_evaluation.reduced_derivative)};
  }

  template <int dim>
  class B2ReducedExecutionAdapterT final
  {
  public:
    B2ReducedExecutionAdapterT(
      const B2RuntimeDataT<dim> &runtime,
      std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>> session,
      experiment::RunEnvironmentRecord environment,
      std::filesystem::path          native_output_directory = {})
      : runtime_(&runtime)
      , session_(std::move(session))
      , environment_(std::move(environment))
      , native_output_directory_(std::move(native_output_directory))
    {
      if (runtime_ == nullptr || !session_)
        throw std::invalid_argument(
          "B2 execution needs runtime data and an owned compilation session");
    }

    Evidence
    operator()(const semantic::v1::ProblemSpec &specification,
               const B2Scenario &                 scenario) const
    {
      validate_b2(scenario);
      if (specification.id != "scalar_convection_neumann_subdomain_control")
        throw std::invalid_argument(
          "B2 execution received a non-Graetz scalar ProblemSpec");
      if (scenario.solver.method != ReducedMethod::bfgs)
        throw std::invalid_argument("B2 execution requires full BFGS");

      const auto policy = make_b2_discretisation_policy(scenario.compile);
      const auto product = make_b2_compilation_product(scenario.compile);
      const auto bindings = make_b2_data_bindings(scenario.problem, *runtime_);
      compiler::v1::DealiiCompiler compiler;
      auto compilation = compiler.compile(specification,
                                          session_,
                                          bindings,
                                          policy,
                                          std::nullopt,
                                          std::nullopt,
                                          product);
      if (!compilation.succeeded() || !compilation.problem)
        {
          std::ostringstream message;
          message << "B2 compilation failed for the selected recipe and bindings";
          for (const auto &diagnostic : compilation.diagnostics.diagnostics())
            message << " [" << diagnostic.component_id << ":"
                    << diagnostic.capability << ":" << diagnostic.remedy << "]";
          throw std::runtime_error(message.str());
        }

      const auto &volume_observation = *scenario.compile.volume_observation;
      const std::string volume_observation_target =
        volume_observation_target_realisation_name(
          volume_observation.target_realisation);
      const auto &target_profile = scenario.problem.target_profile;
      const auto &target_definition =
        b2_target_definition(scenario.problem.target_catalog);
      const auto &forcing_definition = scenario.problem.forcing;
      const auto &observation_definition = selected_scalar_function_definition(
        scenario.problem.observation_region_catalog, "B2 observation catalog");
      const auto &manifest = compilation.problem->manifest();
      contract::require(
        manifest.observation_realisation.find(
          "target=" + volume_observation_target) != std::string::npos &&
          manifest.observation_realisation.find(
            "(" + std::to_string(volume_observation.quadrature_order) + ")") !=
            std::string::npos,
        "B2 compilation manifest does not match the selected volume observation");

      using Model = compiler::v1::detail::NeumannBoundaryControlModel<dim>;
      const auto *model =
        dynamic_cast<const Model *>(&compilation.problem->executable_model());
      if (model == nullptr)
        throw std::runtime_error(
          "B2 execution needs the Neumann boundary model");

      const auto reduced = compilation.problem->make_reduced_dto();
      const contract::StateControlPartitionT<Backend> partition(
        compilation.problem->executable_model(), 0, 1);
      const auto initial_control =
        detail::make_b2_uniform_control<Backend>(
          partition.control_layout(), scenario.solver.initial_control_value);
      const auto initial_evaluation = reduced.evaluate(initial_control);
      const auto derivative_evidence = make_b2_derivative_evidence(
        *compilation.problem, reduced, partition, initial_control);
      const auto report = [&]() {
        switch (scenario.solver.globalization)
          {
            case ReducedGlobalization::armijo:
              return solvers::ReducedFullBfgsSolverT<Backend>(
                       reduced,
                       compilation.problem->metric(),
                       scenario.solver.parameters,
                       solvers::FullBfgsDirectionPolicyT<Backend>(
                         scenario.solver.full_bfgs))
                .solve(initial_control);
            case ReducedGlobalization::fixed_step:
              return solvers::ReducedFixedStepFullBfgsSolverT<Backend>(
                       reduced,
                       compilation.problem->metric(),
                       scenario.solver.parameters,
                       solvers::FullBfgsDirectionPolicyT<Backend>(
                         scenario.solver.full_bfgs))
                .solve(initial_control);
          }
        throw std::invalid_argument(
          "B2 execution received an unknown reduced globalization");
      }();
      const std::string expected_solver_policy =
        scenario.solver.globalization == ReducedGlobalization::armijo ?
          "armijo" :
          "fixed_step";
      contract::require(report.policy_name == expected_solver_policy &&
                          report.policy_parameters.policy_name ==
                            expected_solver_policy,
                        "B2 solver report does not match the selected globalization");
      if (scenario.experiment.retain_fields &&
          !native_output_directory_.empty())
        {
          model->write_native_output(native_output_directory_,
                                     report.final_evaluation.state,
                                     report.control,
                                     report.final_evaluation.adjoint,
                                     &initial_evaluation.state,
                                     &runtime_->forcing,
                                     &runtime_->desired_state);
        }
      contract::require(!report.objective_history.empty() &&
                          !report.gradient_norm_history.empty(),
                        "B2 solver did not retain objective and gradient histories");
      const double initial_objective = report.objective_history.front();
      const double final_objective = report.objective_history.back();
      const double initial_gradient_norm =
        report.gradient_norm_history.front();
      const double final_gradient_norm = report.gradient_norm_history.back();

      const auto initial_objective_components =
        model->objective_components(initial_evaluation.full_point);
      const auto final_objective_components =
        model->objective_components(report.final_evaluation.full_point);
      contract::require(
        std::abs(initial_objective_components.state_tracking +
                   initial_objective_components.control_regularisation -
                 initial_objective) <=
          1.0e-12 * std::max(1.0, std::abs(initial_objective)) &&
          std::abs(final_objective_components.state_tracking +
                     final_objective_components.control_regularisation -
                   final_objective) <=
            1.0e-12 * std::max(1.0, std::abs(final_objective)),
        "B2 objective components do not sum to the reported objective");

      std::vector<double> coefficient_derivative_norm_history{
        b2_block_values_l2_norm(initial_evaluation.reduced_derivative)};
      coefficient_derivative_norm_history.reserve(
        report.iteration_records.size() + 1);
      for (const auto &record : report.iteration_records)
        coefficient_derivative_norm_history.push_back(
          b2_block_values_l2_norm(
            record.common.accepted_evaluation.reduced_derivative));
      contract::require(
        coefficient_derivative_norm_history.size() ==
          report.objective_history.size(),
        "B2 coefficient-derivative and objective histories disagree");

      const auto uncontrolled_state_extrema =
        b2_block_extrema(initial_evaluation.state);
      const auto optimized_state_extrema =
        b2_block_extrema(report.final_evaluation.state);
      const auto control_extrema = b2_block_extrema(report.control);
      const auto adjoint_extrema =
        b2_block_extrema(report.final_evaluation.adjoint);

      const auto &triangulation = session_->triangulation();
      std::size_t boundary_face_count = 0;
      std::size_t fixed_boundary_face_count = 0;
      std::size_t control_boundary_face_count = 0;
      std::size_t outflow_boundary_face_count = 0;
      double      observation_measure = 0.0;
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->material_id() ==
              scenario.problem.recipe.observed_material_id)
            observation_measure += cell->measure();
          for (unsigned int face = 0; face < cell->n_faces(); ++face)
            if (cell->face(face)->at_boundary())
              {
                ++boundary_face_count;
                if (cell->face(face)->boundary_id() ==
                    scenario.problem.boundary.fixed_id)
                  ++fixed_boundary_face_count;
                else if (cell->face(face)->boundary_id() ==
                         scenario.problem.boundary.control_id)
                  ++control_boundary_face_count;
                else if (cell->face(face)->boundary_id() ==
                         scenario.problem.boundary.outflow_id)
                  ++outflow_boundary_face_count;
                else
                  throw std::runtime_error(
                    "B2 evidence found an unclassified boundary face");
              }
        }
      contract::require(
        boundary_face_count == fixed_boundary_face_count +
                                 control_boundary_face_count +
                                 outflow_boundary_face_count &&
          observation_measure > 0.0,
        "B2 structural evidence is incomplete");
      const auto solver_policy =
        experiment::make_reduced_search_policy_snapshot(report);

      std::vector<std::string> selected_fields{
        "objective_history",
        "gradient_norm_history",
        "step_length_history",
        "objective_change_history",
        "solve_counts"};
      if (scenario.experiment.retain_fields)
        {
          selected_fields.insert(selected_fields.end(),
                                 {"state",
                                  "state_uncontrolled",
                                  "control",
                                  "adjoint",
                                  "negative_adjoint",
                                  "target",
                                  "forcing",
                                  "observation_region"});
        }

      std::vector<benchmark::ArtifactField> fields{
        {"b2.graetz_case",
         b2_case_name(scenario.problem.observation_region, target_profile)},
        {"b2.control_discretisation",
         chapter5::neumann_control_discretisation_name(
           scenario.problem.recipe.control_discretisation)},
        {"b2.volume_observation_quadrature_order",
         std::to_string(volume_observation.quadrature_order)},
        {"b2.volume_observation_target_realisation",
         volume_observation_target},
        {"b2.observation_region",
         scenario.problem.observation_region},
        {"b2.observation_definition", observation_definition.id},
        {"b2.observation_kind", scalar_function_kind_name(
                                   observation_definition.kind)},
        {"b2.observation_value", b2_number(observation_definition.value)},
        {"b2.observation_expression", observation_definition.expression},
        {"b2.observation_provenance", observation_definition.provenance},
        {"b2.observation_realisation", "cell-center-indicator"},
        {"b2.observed_material_id",
         std::to_string(scenario.problem.recipe.observed_material_id)},
        {"b2.fixed_boundary_id",
         std::to_string(scenario.problem.boundary.fixed_id)},
        {"b2.control_boundary_id",
         std::to_string(scenario.problem.boundary.control_id)},
        {"b2.outflow_boundary_id",
         std::to_string(scenario.problem.boundary.outflow_id)},
        {"b2.upstream_transition",
         b2_number(scenario.problem.boundary.upstream_transition)},
        {"b2.target_profile", target_profile},
        {"b2.target_definition", target_definition.id},
        {"b2.target_kind", scalar_function_kind_name(target_definition.kind)},
        {"b2.target_value", b2_number(target_definition.value)},
        {"b2.target_expression", target_definition.expression},
        {"b2.target_provenance", target_definition.provenance},
        {"b2.forcing_definition", forcing_definition.id},
        {"b2.forcing_kind", scalar_function_kind_name(forcing_definition.kind)},
        {"b2.forcing_value", b2_number(forcing_definition.value)},
        {"b2.forcing_expression", forcing_definition.expression},
        {"b2.forcing_provenance", forcing_definition.provenance},
        {"b2.fixed_dirichlet_data",
         scenario.problem.fixed_dirichlet_data.id},
        {"b2.fixed_dirichlet_data_kind",
         scalar_function_kind_name(scenario.problem.fixed_dirichlet_data.kind)},
        {"b2.fixed_dirichlet_data_value",
         b2_number(scenario.problem.fixed_dirichlet_data.value)},
        {"b2.fixed_dirichlet_data_expression",
         scenario.problem.fixed_dirichlet_data.expression},
        {"b2.fixed_dirichlet_data_provenance",
         scenario.problem.fixed_dirichlet_data.provenance},
        {"b2.conservative_transport",
         scenario.problem.conservative_transport.id},
        {"b2.conservative_transport_expression",
         scenario.problem.conservative_transport.expression},
        {"b2.conservative_transport_provenance",
         scenario.problem.conservative_transport.provenance},
        {"b2.regularisation_weight",
         b2_number(scenario.problem.data.regularisation_weight)},
        {"benchmark.mesh_vertices",
         std::to_string(triangulation.n_vertices())},
        {"benchmark.mesh_active_cells",
         std::to_string(triangulation.n_active_cells())},
        {"benchmark.mesh_lower", b2_coordinates(scenario.compile.mesh.lower)},
        {"benchmark.mesh_upper", b2_coordinates(scenario.compile.mesh.upper)},
        {"benchmark.boundary_face_count",
         std::to_string(boundary_face_count)},
        {"benchmark.fixed_boundary_face_count",
         std::to_string(fixed_boundary_face_count)},
        {"benchmark.control_boundary_face_count",
         std::to_string(control_boundary_face_count)},
        {"benchmark.outflow_boundary_face_count",
         std::to_string(outflow_boundary_face_count)},
        {"benchmark.observation_measure", b2_number(observation_measure)},
        {"benchmark.state_dimension",
         std::to_string(report.final_evaluation.state.block(0).size())},
        {"benchmark.state_physical_dimension",
         std::to_string(model->physical_state_dimension())},
        {"benchmark.state_independent_dimension",
         std::to_string(model->independent_state_dimension())},
        {"benchmark.control_dimension",
         std::to_string(report.control.block(0).size())},
        {"benchmark.control_physical_dimension",
         std::to_string(model->physical_control_dimension())},
        {"benchmark.control_independent_dimension",
         std::to_string(model->independent_control_dimension())},
        {"benchmark.adjoint_dimension",
         std::to_string(report.final_evaluation.adjoint.block(0).size())},
        {"benchmark.adjoint_physical_dimension",
         std::to_string(model->physical_state_dimension())},
        {"benchmark.adjoint_independent_dimension",
         std::to_string(model->independent_state_dimension())},
        {"b2.derivative_evidence", "finite_difference_and_taylor"},
        {"b2.residual_jvp_error",
         b2_number(derivative_evidence.residual_jvp_error)},
        {"b2.residual_vjp_error",
         b2_number(derivative_evidence.residual_vjp_error)},
        {"b2.reduced_gradient_finite_difference_error",
         b2_number(derivative_evidence.reduced_gradient_finite_difference_error)},
        {"b2.reduced_taylor_step",
         b2_number(derivative_evidence.reduced_taylor_step)},
        {"b2.reduced_taylor_error",
         b2_number(derivative_evidence.reduced_taylor_error)},
        {"b2.reduced_taylor_order",
         b2_number(derivative_evidence.reduced_taylor_order)},
        {"b2.derivative_evidence_passed", "true"},
        {"b2.initial_objective", b2_number(initial_objective)},
        {"b2.initial_tracking_objective",
         b2_number(initial_objective_components.state_tracking)},
        {"b2.initial_control_regularisation_objective",
         b2_number(initial_objective_components.control_regularisation)},
        {"b2.final_objective", b2_number(final_objective)},
        {"b2.final_tracking_objective",
         b2_number(final_objective_components.state_tracking)},
        {"b2.final_control_regularisation_objective",
         b2_number(final_objective_components.control_regularisation)},
        {"b2.relative_objective_reduction",
         b2_number(b2_relative_reduction(initial_objective, final_objective))},
        {"b2.initial_gradient_norm", b2_number(initial_gradient_norm)},
        {"b2.initial_metric_gradient_norm",
         b2_number(initial_gradient_norm)},
        {"b2.final_gradient_norm", b2_number(final_gradient_norm)},
        {"b2.final_metric_gradient_norm", b2_number(final_gradient_norm)},
        {"b2.relative_gradient_reduction",
         b2_number(b2_relative_reduction(initial_gradient_norm,
                                         final_gradient_norm))},
        {"b2.initial_coefficient_derivative_norm",
         b2_number(coefficient_derivative_norm_history.front())},
        {"b2.final_coefficient_derivative_norm",
         b2_number(coefficient_derivative_norm_history.back())},
        {"b2.relative_coefficient_derivative_reduction",
         b2_number(b2_relative_reduction(
           coefficient_derivative_norm_history.front(),
           coefficient_derivative_norm_history.back()))},
        {"b2.state_l2_norm",
         b2_number(b2_block_values_l2_norm(report.final_evaluation.state))},
        {"b2.control_l2_norm", b2_number(b2_block_values_l2_norm(report.control))},
        {"b2.uncontrolled_state_min",
         b2_number(uncontrolled_state_extrema.first)},
        {"b2.uncontrolled_state_max",
         b2_number(uncontrolled_state_extrema.second)},
        {"b2.optimized_state_min", b2_number(optimized_state_extrema.first)},
        {"b2.optimized_state_max", b2_number(optimized_state_extrema.second)},
        {"b2.control_min", b2_number(control_extrema.first)},
        {"b2.control_max", b2_number(control_extrema.second)},
        {"b2.adjoint_min", b2_number(adjoint_extrema.first)},
        {"b2.adjoint_max", b2_number(adjoint_extrema.second)},
        {"solver.method", "bfgs"},
        {"solver.globalization",
         reduced_globalization_name(scenario.solver.globalization)},
        {"solver.initial_control_value",
         b2_number(scenario.solver.initial_control_value)},
        {"solver.full_bfgs_curvature_tolerance",
         b2_number(scenario.solver.full_bfgs.curvature_tolerance)},
        {"solver.policy", report.policy_name},
        {"solver.stopping_reason",
         solvers::reduced_stopping_reason_name(report.stopping_reason)},
        {"solver.accepted_iterations",
         std::to_string(report.accepted_iterations)},
        {"solver.line_search_trial_count",
         std::to_string(report.line_search_trial_count)},
        {"solver.state_solve_count", std::to_string(report.state_solve_count)},
        {"solver.adjoint_solve_count",
         std::to_string(report.adjoint_solve_count)},
        {"solver.metric_solve_count",
         std::to_string(report.metric_solve_count)},
        {"solver.hessian_action_count",
         std::to_string(report.hessian_action_count)},
        {"solver.direction_reset_count",
         std::to_string(report.direction_reset_count)},
        {"solver.objective_history", b2_history(report.objective_history)},
        {"solver.gradient_norm_history",
         b2_history(report.gradient_norm_history)},
        {"solver.metric_gradient_norm_history",
         b2_history(report.gradient_norm_history)},
        {"solver.coefficient_derivative_norm_history",
         b2_history(coefficient_derivative_norm_history)},
        {"solver.step_length_history", b2_history(report.step_length_history)},
        {"solver.step_norm_history", b2_history(report.step_norm_history)},
        {"solver.objective_change_history",
         b2_history(report.objective_change_history)}};

      Envelope envelope{compilation.problem->manifest(),
                        solver_policy,
              report,
              environment_};
      return {std::move(envelope),
              std::move(compilation.diagnostics),
              {},
              std::move(selected_fields),
              std::move(fields)};
    }

  private:
    const B2RuntimeDataT<dim> *runtime_;
    std::shared_ptr<compiler::v1::DealiiCompilationSession<dim>> session_;
    experiment::RunEnvironmentRecord environment_;
    std::filesystem::path native_output_directory_;
  };
} // namespace nmopt::application::chapter6::dealii
