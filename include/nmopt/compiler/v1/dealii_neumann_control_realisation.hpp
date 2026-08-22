#pragma once

#include "nmopt/compiler/v1/dealii_reference_cell.hpp"
#include "nmopt/contract/executable_model.hpp"
#include "nmopt/dealii/facewise_box_constraint.hpp"
#include "nmopt/dealii/mass_metric.hpp"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <locale>
#include <memory>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nmopt::compiler::v1::detail
{
  template <int dim>
  class FacewiseNeumannControlRealisation final
  {
  public:
    using Vector = dealii::Vector<double>;

    FacewiseNeumannControlRealisation(
      const dealii::DoFHandler<dim> &              state_dof_handler,
      const dealii::FiniteElement<dim> &            state_fe,
      const std::vector<bool> &                     constrained_state_dofs,
      std::set<dealii::types::boundary_id>          control_boundary_ids,
      const unsigned int                            quadrature_order,
      const double                                  coupling_scale)
      : control_boundary_ids_(std::move(control_boundary_ids))
    {
      contract::require(!control_boundary_ids_.empty(),
                        "Facewise Neumann control needs a boundary selection");
      contract::require(quadrature_order > 0,
                        "Facewise Neumann control needs a positive quadrature order");
      contract::require(std::isfinite(coupling_scale) && coupling_scale > 0.0,
                        "Facewise Neumann control needs a positive finite coupling scale");

      enumerate_faces(state_dof_handler);
      contract::require(!faces_.empty(),
                        "The selected Neumann boundary has no active boundary faces");
      layout_ = std::make_shared<const contract::BlockLayout>(
        "control",
        std::vector<contract::SpaceId>{{"control"}},
        std::vector<std::size_t>{faces_.size()});
      initialise_matrices(state_dof_handler,
                          state_fe,
                          constrained_state_dofs);
      assemble(state_dof_handler,
               state_fe,
               constrained_state_dofs,
               quadrature_order,
               coupling_scale);
    }

    std::size_t
    dimension() const
    {
      return faces_.size();
    }

    const contract::LayoutPtr &
    layout() const
    {
      return layout_;
    }

    const std::vector<dealii::Point<dim>> &
    coordinates() const
    {
      return coordinates_;
    }

    bool
    is_control_face(
      const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
      const unsigned int face) const
    {
      return cell->face(face)->at_boundary() &&
             control_boundary_ids_.count(cell->face(face)->boundary_id()) != 0;
    }

    Vector
    coupling_action(const Vector &control) const
    {
      require_control(control);
      Vector value(control_coupling_.m());
      control_coupling_.vmult(value, control);
      return value;
    }

    Vector
    coupling_transpose_action(const Vector &state_covector) const
    {
      contract::require(state_covector.size() == control_coupling_.m(),
                        "Facewise Neumann coupling transpose received the wrong state dimension");
      Vector value(control_coupling_.n());
      control_coupling_.Tvmult(value, state_covector);
      return value;
    }

    double
    regularisation_objective(const Vector &control,
                             const double  regularisation_weight) const
    {
      require_control(control);
      Vector mass_times_control(dimension());
      control_mass_->vmult(mass_times_control, control);
      return 0.5 * regularisation_weight * (control * mass_times_control);
    }

    Vector
    regularisation_derivative(const Vector &control,
                              const double  regularisation_weight) const
    {
      require_control(control);
      Vector value(dimension());
      control_mass_->vmult(value, control);
      value *= regularisation_weight;
      return value;
    }

    dealii_backend::MassMetric
    l2_metric(
      dealii_backend::MassMetricSolveParameters solve_parameters = {}) const
    {
      return dealii_backend::MassMetric("l2_facewise",
                                        layout_,
                                        control_mass_,
                                        solve_parameters);
    }

    dealii_backend::FacewiseBoxConstraint
    l2_box_constraint(
      Vector                              lower,
      Vector                              upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::FacewiseBoxConstraint(layout_,
                                                    std::move(lower),
                                                    std::move(upper),
                                                    projection_metric);
    }

    dealii_backend::FacewiseBoxConstraint
    l2_box_constraint(
      const double                        lower,
      const double                        upper,
      const dealii_backend::MassMetric & projection_metric) const
    {
      return dealii_backend::FacewiseBoxConstraint(layout_,
                                                    lower,
                                                    upper,
                                                    projection_metric);
    }

    void
    write_native_output(const std::filesystem::path &path,
                        const Vector &                control) const
    {
      static_assert(dim == 2,
                    "Facewise Neumann control output currently supports two dimensions");
      require_control(control);

      std::size_t point_count = 0;
      for (const auto &face : faces_)
        point_count += face.vertices.size();

      std::ofstream output(path);
      if (!output)
        throw std::runtime_error("could not open Neumann boundary-control output");
      output.imbue(std::locale::classic());
      output << "<?xml version=\"1.0\"?>\n"
             << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
                "byte_order=\"LittleEndian\">\n"
             << "<UnstructuredGrid>\n"
             << "<Piece NumberOfPoints=\"" << point_count
             << "\" NumberOfCells=\"" << faces_.size() << "\">\n"
             << "<PointData>\n</PointData>\n"
             << "<CellData Scalars=\"control\">\n"
             << "<DataArray type=\"Float64\" Name=\"control\" "
                "format=\"ascii\">\n";
      output.precision(17);
      for (typename Vector::size_type index = 0; index < control.size(); ++index)
        output << control[index] << ' ';
      output << "\n</DataArray>\n</CellData>\n"
             << "<Points>\n"
             << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                "format=\"ascii\">\n";
      for (const auto &face : faces_)
        for (const auto &point : face.vertices)
          output << point[0] << ' ' << point[1] << " 0 ";
      output << "\n</DataArray>\n</Points>\n"
             << "<Cells>\n"
             << "<DataArray type=\"Int32\" Name=\"connectivity\" "
                "format=\"ascii\">\n";
      std::size_t point_index = 0;
      for (const auto &face : faces_)
        for (std::size_t vertex = 0; vertex < face.vertices.size(); ++vertex)
          output << point_index++ << ' ';
      output << "\n</DataArray>\n"
             << "<DataArray type=\"Int32\" Name=\"offsets\" "
                "format=\"ascii\">\n";
      point_index = 0;
      for (const auto &face : faces_)
        {
          point_index += face.vertices.size();
          output << point_index << ' ';
        }
      output << "\n</DataArray>\n"
             << "<DataArray type=\"UInt8\" Name=\"types\" "
                "format=\"ascii\">\n";
      for (std::size_t face = 0; face < faces_.size(); ++face)
        output << "3 ";
      output << "\n</DataArray>\n</Cells>\n"
             << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
      if (!output)
        throw std::runtime_error("could not write Neumann boundary-control output");
    }

  private:
    struct FaceRecord
    {
      typename dealii::DoFHandler<dim>::active_cell_iterator cell;
      unsigned int                                           face = 0;
      std::vector<dealii::Point<dim>>                         vertices;
    };

    void
    enumerate_faces(const dealii::DoFHandler<dim> &state_dof_handler)
    {
      for (auto cell = state_dof_handler.begin_active();
           cell != state_dof_handler.end();
           ++cell)
        for (unsigned int face = 0; face < cell->n_faces(); ++face)
          if (is_control_face(cell, face))
            {
              FaceRecord record{cell, face, {}};
              for (unsigned int vertex = 0;
                   vertex < cell->face(face)->n_vertices();
                   ++vertex)
                record.vertices.push_back(cell->face(face)->vertex(vertex));
              faces_.push_back(std::move(record));
              coordinates_.push_back(cell->face(face)->center());
            }
    }

    void
    initialise_matrices(const dealii::DoFHandler<dim> &state_dof_handler,
                        const dealii::FiniteElement<dim> &state_fe,
                        const std::vector<bool> &constrained_state_dofs)
    {
      dealii::DynamicSparsityPattern coupling_dsp(state_dof_handler.n_dofs(),
                                                   dimension());
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe.dofs_per_cell);
      for (std::size_t control_index = 0;
           control_index < faces_.size();
           ++control_index)
        {
          faces_[control_index].cell->get_dof_indices(state_indices);
          for (const auto state_index : state_indices)
            if (!constrained_state_dofs.at(state_index))
              coupling_dsp.add(state_index, control_index);
        }
      control_sparsity_.copy_from(coupling_dsp);
      control_coupling_.reinit(control_sparsity_);

      dealii::DynamicSparsityPattern mass_dsp(dimension(), dimension());
      for (std::size_t index = 0; index < dimension(); ++index)
        mass_dsp.add(index, index);
      control_mass_sparsity_.copy_from(mass_dsp);
      control_mass_ = std::make_shared<dealii::SparseMatrix<double>>();
      control_mass_->reinit(control_mass_sparsity_);
    }

    void
    assemble(const dealii::DoFHandler<dim> &state_dof_handler,
             const dealii::FiniteElement<dim> &state_fe,
             const std::vector<bool> &constrained_state_dofs,
             const unsigned int quadrature_order,
             const double coupling_scale)
    {
      const auto face_quadrature = make_gauss_face_quadrature(
        state_dof_handler.get_triangulation(),
        quadrature_order,
        "The facewise Neumann control realization");
      dealii::FEFaceValues<dim> face_values(state_fe,
                                            *face_quadrature,
                                            dealii::update_values |
                                              dealii::update_JxW_values);
      std::vector<dealii::types::global_dof_index> state_indices(
        state_fe.dofs_per_cell);
      for (std::size_t control_index = 0;
           control_index < faces_.size();
           ++control_index)
        {
          const auto &record = faces_[control_index];
          face_values.reinit(record.cell, record.face);
          record.cell->get_dof_indices(state_indices);
          double control_measure = 0.0;
          for (unsigned int q = 0; q < face_quadrature->size(); ++q)
            {
              const double weight = face_values.JxW(q);
              control_measure += weight;
              for (unsigned int i = 0; i < state_fe.dofs_per_cell; ++i)
                {
                  const auto global_i = state_indices[i];
                  if (!constrained_state_dofs.at(global_i))
                    control_coupling_.add(global_i,
                                          control_index,
                                          coupling_scale *
                                            face_values.shape_value(i, q) *
                                            weight);
                }
            }
          contract::require(control_measure > 0.0,
                            "A Neumann control face has zero measure");
          control_mass_->add(control_index, control_index, control_measure);
        }
    }

    void
    require_control(const Vector &control) const
    {
      contract::require(control.size() == dimension(),
                        "Facewise Neumann control has an incompatible dimension");
    }

    const std::set<dealii::types::boundary_id> control_boundary_ids_;
    std::vector<FaceRecord>                    faces_;
    std::vector<dealii::Point<dim>>            coordinates_;
    dealii::SparsityPattern                    control_sparsity_;
    dealii::SparsityPattern                    control_mass_sparsity_;
    dealii::SparseMatrix<double>               control_coupling_;
    std::shared_ptr<dealii::SparseMatrix<double>> control_mass_;
    contract::LayoutPtr                           layout_;
  };
} // namespace nmopt::compiler::v1::detail
