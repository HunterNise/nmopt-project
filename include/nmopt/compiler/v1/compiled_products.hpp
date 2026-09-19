#pragma once

#include "nmopt/compiler/v1/compilation_manifest.hpp"
#include "nmopt/compiler/v1/compiled_application_view.hpp"
#include "nmopt/contract/reduced_dto.hpp"
#include "nmopt/contract/reduced_hessian.hpp"
#include "nmopt/contract/quadratic_kkt.hpp"
#include "nmopt/contract/pdas.hpp"
#include "nmopt/contract/supplied_otd.hpp"
#include "nmopt/semantic/v1/validation.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1
{

  // One immutable compiled source of truth for coefficientwise cellwise box
  // data. Projection and complementarity products retain this token so a
  // benchmark cannot silently pair independently reconstructed bounds.
  template <typename Backend>
  class CompiledCellwiseBoxDataT final
  {
  public:
    using Metric = contract::MetricT<Backend>;
    using Primal = contract::PrimalBlockT<Backend>;

    CompiledCellwiseBoxDataT(
      contract::LayoutPtr           layout,
      Primal                        lower,
      Primal                        upper,
      std::shared_ptr<const Metric> metric,
      std::string                   semantic_id,
      std::string                   data_provenance)
      : layout_(std::move(layout))
      , lower_(std::move(lower))
      , upper_(std::move(upper))
      , metric_(std::move(metric))
      , semantic_id_(std::move(semantic_id))
      , data_provenance_(std::move(data_provenance))
    {
      contract::require(static_cast<bool>(layout_),
                        "Compiled cellwise box data needs a layout");
      contract::require(layout_->n_blocks() == 1,
                        "Compiled cellwise box data supports one control block");
      contract::require(static_cast<bool>(metric_),
                        "Compiled cellwise box data needs a metric owner");
      contract::require(layout_->compatible_with(*metric_->layout()),
                        "Compiled cellwise box data metric has an incompatible layout");
      contract::require(static_cast<bool>(lower_.layout()) &&
                          static_cast<bool>(upper_.layout()) &&
                          lower_.layout()->compatible_with(*layout_) &&
                          upper_.layout()->compatible_with(*layout_),
                        "Compiled cellwise box data bounds have an incompatible layout");
      contract::require(!semantic_id_.empty() && !data_provenance_.empty(),
                        "Compiled cellwise box data needs provenance");
      for (std::size_t index = 0; index < layout_->dimension(0); ++index)
        {
          const double lower_value = Backend::value(lower_.block(0), index);
          const double upper_value = Backend::value(upper_.block(0), index);
          contract::require(std::isfinite(lower_value) &&
                              std::isfinite(upper_value),
                            "Compiled cellwise box data bounds must be finite");
          contract::require(lower_value <= upper_value,
                            "Compiled cellwise box lower bound exceeds upper bound");
        }
      layout_signature_ = make_layout_signature(*layout_);
      bounds_digest_ = make_bounds_digest(lower_, upper_);
      token_id_ = semantic_id_ + "|" + metric_->id() + "|" +
                  layout_signature_ + "|" + bounds_digest_ + "|" +
                  data_provenance_;
      token_ = std::make_shared<const Token>();
    }

    const contract::LayoutPtr &
    layout() const
    {
      return layout_;
    }

    const Primal &
    lower() const
    {
      return lower_;
    }

    const Primal &
    upper() const
    {
      return upper_;
    }

    const std::shared_ptr<const Metric> &
    metric_owner() const
    {
      return metric_;
    }

    const std::string &
    semantic_id() const
    {
      return semantic_id_;
    }

    const std::string &
    data_provenance() const
    {
      return data_provenance_;
    }

    const std::string &
    layout_signature() const
    {
      return layout_signature_;
    }

    const std::string &
    bounds_digest() const
    {
      return bounds_digest_;
    }

    const std::string &
    token_id() const
    {
      return token_id_;
    }

    const std::shared_ptr<const void> &
    token() const
    {
      return token_;
    }

    bool
    compatible_with(const CompiledCellwiseBoxDataT &other) const
    {
      return token_id_ == other.token_id_ &&
             layout_->compatible_with(*other.layout_) &&
             metric_->realisation_witness().matches(
               other.metric_->realisation_witness());
    }

  private:
    struct Token
    {};

    static std::string
    make_layout_signature(const contract::BlockLayout &layout)
    {
      std::ostringstream result;
      result << layout.label() << ":";
      for (std::size_t block = 0; block < layout.n_blocks(); ++block)
        result << layout.space(block).value << "[" << layout.dimension(block)
               << "];";
      return result.str();
    }

    static void
    update_digest(std::uint64_t &state, const std::uint64_t value)
    {
      for (unsigned int byte = 0; byte < sizeof(value); ++byte)
        {
          state ^= (value >> (byte * 8U)) & 0xffU;
          state *= 1099511628211ULL;
        }
    }

    static std::string
    make_bounds_digest(const Primal &lower, const Primal &upper)
    {
      std::uint64_t state = 1469598103934665603ULL;
      for (const Primal *values : {&lower, &upper})
        for (std::size_t block = 0; block < values->n_blocks(); ++block)
          for (std::size_t index = 0;
               index < values->layout()->dimension(block);
               ++index)
            {
              std::uint64_t bits = 0;
              const double value = Backend::value(values->block(block), index);
              std::memcpy(&bits, &value, sizeof(value));
              update_digest(state, bits);
            }
      std::ostringstream result;
      result << std::hex << std::setw(16) << std::setfill('0') << state;
      return result.str();
    }

    contract::LayoutPtr           layout_;
    Primal                        lower_;
    Primal                        upper_;
    std::shared_ptr<const Metric> metric_;
    std::string                   semantic_id_;
    std::string                   data_provenance_;
    std::string                   layout_signature_;
    std::string                   bounds_digest_;
    std::string                   token_id_;
    std::shared_ptr<const void>   token_;
  };

  template <typename Backend>
  class CompiledQuadraticKKTProblemT final
  {
  public:
    using Product = contract::EqualityConstrainedQuadraticKKTProductT<Backend>;

    CompiledQuadraticKKTProblemT(std::shared_ptr<const Product> product,
                                 CompilationManifest manifest,
                                 std::shared_ptr<const void> lifetime_owner = {})
      : lifetime_owner_(std::move(lifetime_owner))
      , product_(std::move(product))
      , manifest_(std::move(manifest))
    {
      contract::require(static_cast<bool>(product_),
                        "A compiled KKT problem needs a product");
      contract::require(
        manifest_.resolved_decision.kkt_record.present,
        "A compiled KKT problem needs a resolved KKT manifest record");
    }

    const Product &
    product() const
    {
      return *product_;
    }

    const CompilationManifest &
    manifest() const
    {
      return manifest_;
    }

  private:
    std::shared_ptr<const void> lifetime_owner_;
    std::shared_ptr<const Product> product_;
    CompilationManifest            manifest_;
  };

  template <typename Backend>
  class CompiledPDASProblemT final
  {
  public:
    using Product = contract::EqualityConstrainedQuadraticKKTProductT<Backend>;
    using Complementarity = contract::BoxComplementarityT<Backend>;
    using Metric = contract::MetricT<Backend>;
    using Solver = contract::PDASSolverT<Backend>;
    using BoxData = CompiledCellwiseBoxDataT<Backend>;
    using Constraint = contract::ConstraintT<Backend>;

    CompiledPDASProblemT(
      std::shared_ptr<const Product>          product,
      std::shared_ptr<const Complementarity> complementarity,
      std::shared_ptr<const Metric>           metric,
      const std::size_t                       control_block,
      contract::PDASPolicy                   pdas_policy,
      contract::QuadraticKKTSolverPolicy     kkt_solver_policy,
      CompilationManifest                    manifest,
      std::shared_ptr<const void>             lifetime_owner = {},
      std::shared_ptr<const BoxData>          box_data = {},
      std::shared_ptr<const Constraint>       projection_constraint = {})
      : lifetime_owner_(std::move(lifetime_owner))
      , metric_(std::move(metric))
      , box_data_(std::move(box_data))
      , product_(std::move(product))
      , complementarity_(std::move(complementarity))
      , projection_constraint_(std::move(projection_constraint))
      , control_block_(control_block)
      , pdas_policy_(std::move(pdas_policy))
      , kkt_solver_policy_(std::move(kkt_solver_policy))
      , manifest_(std::move(manifest))
    {
      contract::require(static_cast<bool>(product_),
                        "A compiled PDAS problem needs a KKT product");
      contract::require(static_cast<bool>(complementarity_),
                        "A compiled PDAS problem needs box complementarity");
      contract::require(static_cast<bool>(metric_),
                        "A compiled PDAS problem needs its metric owner");
      contract::require(static_cast<bool>(box_data_),
                        "A compiled PDAS problem needs shared box data");
      contract::require(static_cast<bool>(projection_constraint_),
                        "A compiled PDAS problem needs its projection constraint");
      contract::require(
        manifest_.resolved_decision.pdas_record.present,
        "A compiled PDAS problem needs a resolved PDAS manifest record");
      contract::require(control_block_ < product_->layout().primal->n_blocks(),
                        "A compiled PDAS control block is outside the KKT layout");
      contract::require(product_->layout().primal->dimension(control_block_) ==
                          complementarity_->layout()->dimension(0),
                        "A compiled PDAS complementarity layout does not match its control block");
      contract::require(metric_->layout()->compatible_with(
                          *complementarity_->layout()),
                        "A compiled PDAS metric does not match its box layout");
      contract::require(box_data_->layout()->compatible_with(
                          *complementarity_->layout()),
                        "A compiled PDAS shared box does not match its complementarity");
      contract::require(box_data_->metric_owner()->realisation_witness().matches(
                          metric_->realisation_witness()),
                        "A compiled PDAS shared box does not match its metric");
      contract::require(projection_constraint_->layout()->compatible_with(
                          *box_data_->layout()),
                        "A compiled PDAS projection constraint does not match its shared box");
      contract::require(projection_constraint_->box_data_token() ==
                          box_data_->token(),
                        "A compiled PDAS projection constraint does not use its shared box token");
      contract::require(complementarity_->box_data_token() == box_data_->token(),
                        "A compiled PDAS complementarity does not use its shared box token");
      if (lifetime_owner_)
        product_ = std::make_shared<const Product>(
          product_->with_lifetime_owner(lifetime_owner_));
    }

    const Product &
    product() const
    {
      return *product_;
    }

    const Complementarity &
    complementarity() const
    {
      return *complementarity_;
    }

    const Metric &
    metric() const
    {
      return *metric_;
    }

    const BoxData &
    box_data() const
    {
      return *box_data_;
    }

    const Constraint *
    constraint() const
    {
      return projection_constraint_.get();
    }

    const contract::PDASPolicy &
    pdas_policy() const
    {
      return pdas_policy_;
    }

    const contract::QuadraticKKTSolverPolicy &
    kkt_solver_policy() const
    {
      return kkt_solver_policy_;
    }

    Solver
    make_solver(typename Solver::SolveAction solve_action) const
    {
      return Solver(*product_,
                    *complementarity_,
                    control_block_,
                    std::move(solve_action));
    }

    const CompilationManifest &
    manifest() const
    {
      return manifest_;
    }

  private:
    // The metric is declared before the complementarity so callback-held
    // metric references are destroyed after the complementarity object.
    std::shared_ptr<const void>             lifetime_owner_;
    std::shared_ptr<const Metric>            metric_;
    std::shared_ptr<const BoxData>           box_data_;
    std::shared_ptr<const Product>           product_;
    std::shared_ptr<const Complementarity>  complementarity_;
    std::shared_ptr<const Constraint>        projection_constraint_;
    std::size_t                              control_block_;
    contract::PDASPolicy                     pdas_policy_;
    contract::QuadraticKKTSolverPolicy       kkt_solver_policy_;
    CompilationManifest                      manifest_;
  };

  // The compiler product contains the backend-neutral executable ports plus
  // an optional compiler-path compiled application view. A concrete lowerer
  // stays private to its compiler and supplies the model, metric, optional
  // constraint, formulation services, and explicitly typed application seam
  // below; the compiled view is not part of ExecutableModelT.
  template <typename Backend>
  class CompiledProblemT final
  {
  public:
    using Model = contract::ExecutableModelT<Backend>;
    using Metric = contract::MetricT<Backend>;
    using Constraint = contract::ConstraintT<Backend>;
    using ReducedHessian = contract::ReducedHessianT<Backend>;
    using BoxData = CompiledCellwiseBoxDataT<Backend>;
    using CompiledApplicationView = CompiledApplicationViewT<Backend>;

    CompiledProblemT(std::shared_ptr<const Model>             executable,
                     std::shared_ptr<const Metric>            metric,
                     std::shared_ptr<const Constraint>        constraint,
                     contract::StateAdjointSolversT<Backend>   solvers,
                     CompilationManifest                       manifest,
                     std::shared_ptr<const void>               lifetime_owner = {},
                     std::shared_ptr<const ReducedHessian>     reduced_hessian = {},
                     std::shared_ptr<const BoxData>            box_data = {},
                     std::shared_ptr<const CompiledApplicationView>
                       compiled_application_view = {})
      : executable_(std::move(executable))
      , compiled_application_view_(std::move(compiled_application_view))
      , metric_(std::move(metric))
      , constraint_(std::move(constraint))
      , solvers_(std::move(solvers))
      , manifest_(std::move(manifest))
      , lifetime_owner_(std::move(lifetime_owner))
      , reduced_hessian_(std::move(reduced_hessian))
      , box_data_(std::move(box_data))
    {
      contract::require(static_cast<bool>(executable_),
                        "A compiled problem needs an executable model");
      contract::require(static_cast<bool>(metric_),
                        "A compiled problem needs a search metric");
      contract::require(static_cast<bool>(solvers_.solve_state),
                        "A compiled problem needs a state solve service");
      contract::require(static_cast<bool>(solvers_.solve_adjoint),
                        "A compiled problem needs an adjoint solve service");
      contract::require(!manifest_.resolved_decision.semantic_problem_id.empty(),
                        "A compiled problem manifest needs a semantic identifier");
      const contract::StateControlPartitionT<Backend> partition(*executable_, 0, 1);
      contract::require(metric_->layout()->compatible_with(
                          *partition.control_layout()),
                        "A compiled problem metric does not match its control block");
      if (constraint_)
        {
          contract::require(constraint_->layout()->compatible_with(*metric_->layout()),
                            "A compiled problem constraint does not match its metric");
          contract::require(constraint_->supports_projection_in(*metric_),
                            "A compiled problem constraint cannot project in its metric");
        }
      if (box_data_)
        {
          contract::require(static_cast<bool>(constraint_),
                            "A compiled problem shared box needs a projection constraint");
          contract::require(box_data_->layout()->compatible_with(
                              *metric_->layout()),
                            "A compiled problem shared box does not match its metric");
          contract::require(box_data_->metric_owner()->realisation_witness().matches(
                              metric_->realisation_witness()),
                            "A compiled problem shared box does not match its metric owner");
          contract::require(constraint_->layout()->compatible_with(
                              *box_data_->layout()),
                            "A compiled problem projection constraint does not match its shared box");
          contract::require(constraint_->box_data_token() == box_data_->token(),
                            "A compiled problem projection constraint does not use its shared box token");
        }
      if (reduced_hessian_)
        contract::require(
          reduced_hessian_->layout()->compatible_with(*metric_->layout()),
          "A compiled problem Hessian does not match its metric");
    }

    const Model &
    executable_model() const
    {
      return *executable_;
    }

    const CompiledApplicationView *
    compiled_application_view() const noexcept
    {
      return compiled_application_view_ ? compiled_application_view_.get() : nullptr;
    }

    const Metric &
    metric() const
    {
      return *metric_;
    }

    const Constraint *
    constraint() const
    {
      return constraint_ ? constraint_.get() : nullptr;
    }

    const BoxData *
    box_data() const
    {
      return box_data_ ? box_data_.get() : nullptr;
    }

    const ReducedHessian *
    reduced_hessian() const
    {
      return reduced_hessian_ ? reduced_hessian_.get() : nullptr;
    }

    const contract::StateAdjointSolversT<Backend> &
    state_adjoint_solvers() const
    {
      return solvers_;
    }

    contract::ReducedDTOT<Backend>
    make_reduced_dto() const
    {
      return contract::ReducedDTOT<Backend>(
        executable_,
        contract::StateControlPartitionT<Backend>(*executable_, 0, 1),
        solvers_,
        lifetime_owner_);
    }

    const CompilationManifest &
    manifest() const
    {
      return manifest_;
    }

  private:
    std::shared_ptr<const Model>           executable_;
    std::shared_ptr<const CompiledApplicationView> compiled_application_view_;
    std::shared_ptr<const Metric>          metric_;
    std::shared_ptr<const Constraint>      constraint_;
    contract::StateAdjointSolversT<Backend> solvers_;
    CompilationManifest                     manifest_;
    std::shared_ptr<const void>             lifetime_owner_;
    std::shared_ptr<const ReducedHessian>   reduced_hessian_;
    std::shared_ptr<const BoxData>          box_data_;
  };

  // A supplied-OTD compilation is a distinct product.  It intentionally does
  // not expose the reduced DTO service bundle used by CompiledProblemT.
  template <typename Backend>
  class CompiledSuppliedOTDProblemT final
  {
  public:
    using System = contract::SuppliedOTDSystemT<Backend>;

    CompiledSuppliedOTDProblemT(std::shared_ptr<const System> system,
                                CompilationManifest           manifest,
                                std::shared_ptr<const void>    lifetime_owner = {})
      : lifetime_owner_(std::move(lifetime_owner))
      , system_(std::move(system))
      , manifest_(std::move(manifest))
    {
      contract::require(static_cast<bool>(system_),
                        "A supplied-OTD problem needs an executable system");
      contract::require(
        manifest_.resolved_decision.formulation_record.kind ==
            semantic::v1::FormulationKind::all_at_once &&
          manifest_.resolved_decision.formulation_record.provenance ==
            semantic::v1::FormulationProvenance::supplied_otd,
        "A supplied-OTD problem needs an all-at-once supplied-OTD manifest");
      contract::require(manifest_.resolved_decision.supplied_otd_record.present,
                        "A supplied-OTD problem needs typed block provenance");
      contract::require(
        manifest_.resolved_decision.supplied_otd_record.declaration.has_value(),
        "A supplied-OTD problem needs its typed formulation declaration");
    }

    const System &
    system() const
    {
      return *system_;
    }

    const CompilationManifest &
    manifest() const
    {
      return manifest_;
    }

  private:
    // The owner must be declared first so reverse member destruction tears
    // down callback-held systems before the session and its mesh.
    std::shared_ptr<const void>    lifetime_owner_;
    std::shared_ptr<const System> system_;
    CompilationManifest           manifest_;
  };

  template <typename Backend>
  struct CompilationResultT
  {
    semantic::v1::ValidationReport             diagnostics;
    std::shared_ptr<const CompiledProblemT<Backend>> problem;
    std::shared_ptr<const CompiledQuadraticKKTProblemT<Backend>>
      kkt_problem;
    std::shared_ptr<const CompiledPDASProblemT<Backend>> pdas_problem;
    std::shared_ptr<const CompiledSuppliedOTDProblemT<Backend>>
      supplied_otd_problem;

    bool
    succeeded() const
    {
      return diagnostics.valid() &&
             (static_cast<bool>(problem) ||
              static_cast<bool>(kkt_problem) ||
              static_cast<bool>(pdas_problem) ||
              static_cast<bool>(supplied_otd_problem));
    }
  };
} // namespace nmopt::compiler::v1
