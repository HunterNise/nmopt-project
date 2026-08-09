#pragma once

#include "nmopt/semantic/v1/types.hpp"

namespace nmopt::compiler::v1
{
  // This small built-in registry is intentionally explicit. A semantic kind
  // is executable only when the v1 compiler has registered the corresponding
  // lowerer capability; it never falls back to a nearby term.
  class DealiiCapabilityRegistryV1 final
  {
  public:
    bool
    has_residual_term_lowerer(const semantic::v1::ResidualTermKind kind) const
    {
      switch (kind)
        {
          case semantic::v1::ResidualTermKind::unspecified:
            return false;
          case semantic::v1::ResidualTermKind::diffusion_reaction:
          case semantic::v1::ResidualTermKind::tensor_diffusion:
          case semantic::v1::ResidualTermKind::conservative_transport:
          case semantic::v1::ResidualTermKind::advective_transport:
          case semantic::v1::ResidualTermKind::reaction:
          case semantic::v1::ResidualTermKind::parameter_diffusion_reaction:
          case semantic::v1::ResidualTermKind::volume_source:
          case semantic::v1::ResidualTermKind::volume_control:
          case semantic::v1::ResidualTermKind::neumann_control:
          case semantic::v1::ResidualTermKind::robin_bilinear:
          case semantic::v1::ResidualTermKind::robin_source:
            return true;
        }
      return false;
    }

    bool
    has_observation_lowerer(const semantic::v1::ObservationKind kind) const
    {
      switch (kind)
        {
          case semantic::v1::ObservationKind::unspecified:
            return false;
          case semantic::v1::ObservationKind::volume_restriction:
          case semantic::v1::ObservationKind::boundary_trace:
          case semantic::v1::ObservationKind::boundary_restriction:
            return true;
        }
      return false;
    }

    bool
    has_loss_lowerer(const semantic::v1::LossKind kind) const
    {
      switch (kind)
        {
          case semantic::v1::LossKind::unspecified:
            return false;
          case semantic::v1::LossKind::quadratic_tracking:
          case semantic::v1::LossKind::quadratic_control_regularisation:
          case semantic::v1::LossKind::quadratic_h1_control_regularisation:
          case semantic::v1::LossKind::quadratic_parameter_regularisation:
            return true;
        }
      return false;
    }

    bool
    has_metric_lowerer(const semantic::v1::MetricKind kind) const
    {
      return kind == semantic::v1::MetricKind::l2 ||
             kind == semantic::v1::MetricKind::h1;
    }

    bool
    has_constraint_lowerer(const semantic::v1::ConstraintKind kind) const
    {
      return kind == semantic::v1::ConstraintKind::cellwise_box ||
             kind == semantic::v1::ConstraintKind::facewise_box;
    }

    bool
    has_transformation_lowerer(
      const semantic::v1::TransformationKind kind) const
    {
      return kind == semantic::v1::TransformationKind::fixed_dirichlet_reconstruction ||
             kind == semantic::v1::TransformationKind::dirichlet_control_lifting;
    }
  };
} // namespace nmopt::compiler::v1
