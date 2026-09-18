#include "dealii_backend_support.hpp"

namespace
{
  void
  run_hminus1_metric_contract_test()
  {
    const auto layout = std::make_shared<const contract::BlockLayout>(
      "hminus1_metric_test",
      std::vector<contract::SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    dealii::DynamicSparsityPattern dsp(2, 2);
    dsp.add(0, 0);
    dsp.add(0, 1);
    dsp.add(1, 0);
    dsp.add(1, 1);
    dealii::SparsityPattern sparsity;
    sparsity.copy_from(dsp);

    auto mass = std::make_shared<dealii::SparseMatrix<double>>();
    mass->reinit(sparsity);
    mass->set(0, 0, 2.0);
    mass->set(1, 1, 3.0);
    auto laplace = std::make_shared<dealii::SparseMatrix<double>>();
    laplace->reinit(sparsity);
    laplace->set(0, 0, 4.0);
    laplace->set(0, 1, 1.0);
    laplace->set(1, 0, 1.0);
    laplace->set(1, 1, 2.0);

    dealii_backend::MetricSolveParameters solve_parameters;
    solve_parameters.maximum_iterations = 100;
    solve_parameters.relative_tolerance = 1e-13;
    solve_parameters.absolute_tolerance = 1e-15;
    const dealii_backend::Hminus1Metric metric(
      "hminus1_continuous", layout, mass, laplace, solve_parameters);

    dealii::Vector<double> primal_values(2);
    primal_values[0] = 1.0;
    primal_values[1] = -2.0;
    const Primal primal(layout, {std::move(primal_values)});
    const Covector applied = metric.apply(primal);
    require_close(applied.block(0)[0],
                  20.0 / 7.0,
                  1e-12,
                  "H-1 metric M K^-1 M first component");
    require_close(applied.block(0)[1],
                  -78.0 / 7.0,
                  1e-12,
                  "H-1 metric M K^-1 M second component");

    const Primal recovered = metric.inverse_apply(applied);
    dealii::Vector<double> recovery_error = recovered.block(0);
    recovery_error.add(-1.0, primal.block(0));
    require_close(recovery_error.l2_norm(),
                  0.0,
                  1e-11,
                  "H-1 metric apply/inverse pairing");

    dealii::Vector<double> second_values(2);
    second_values[0] = 0.25;
    second_values[1] = 0.5;
    const Primal second(layout, {std::move(second_values)});
    require_close(contract::pair(metric.apply(primal), second),
                  contract::pair(metric.apply(second), primal),
                  1e-12,
                  "H-1 metric symmetry");
    contract::require(metric.id() == "hminus1_continuous" &&
                        metric.solve_parameters().maximum_iterations == 100,
                      "H-1 metric omitted its identity or solve policy");
  }

  void
  run_quadratic_form_contract_test()
  {
    dealii::DynamicSparsityPattern dsp(2, 2);
    dsp.add(0, 0);
    dsp.add(0, 1);
    dsp.add(1, 0);
    dsp.add(1, 1);
    dealii::SparsityPattern sparsity;
    sparsity.copy_from(dsp);

    dealii::SparseMatrix<double> matrix(sparsity);
    matrix.set(0, 0, 2.0);
    matrix.set(0, 1, 1.0);
    matrix.set(1, 0, 1.0);
    matrix.set(1, 1, 3.0);

    dealii::Vector<double> linear(2);
    linear[0] = 1.0;
    linear[1] = -2.0;
    dealii::Vector<double> coordinates(2);
    coordinates[0] = 2.0;
    coordinates[1] = -1.0;

    const dealii_backend::QuadraticForm form(matrix, linear, 3.0);
    require_close(form.value(coordinates), 1.0, 1e-14,
                  "Quadratic-form value");
    const auto gradient = form.gradient(coordinates);
    require_close(gradient[0], 2.0, 1e-14,
                  "Quadratic-form gradient first component");
    require_close(gradient[1], 1.0, 1e-14,
                  "Quadratic-form gradient second component");
    const auto action = form.hessian_action(coordinates);
    require_close(action[0], 3.0, 1e-14,
                  "Quadratic-form Hessian action first component");
    require_close(action[1], -1.0, 1e-14,
                  "Quadratic-form Hessian action second component");

    const dealii_backend::QuadraticForm homogeneous(matrix);
    require_close(homogeneous.value(coordinates), 3.5, 1e-14,
                  "Homogeneous quadratic-form value");

    dealii::Vector<double> bad_linear(1);
    test_support::require_contract_error(
      [&]() {
        const dealii_backend::QuadraticForm invalid(matrix, bad_linear, 0.0);
        (void)invalid;
      },
      "Quadratic-form affine load has an incompatible dimension",
      "Quadratic-form affine-load validation");
    dealii::Vector<double> bad_coordinates(1);
    test_support::require_contract_error(
      [&]() { (void)form.value(bad_coordinates); },
      "Quadratic-form action received an incompatible dimension",
      "Quadratic-form coordinate validation");
  }

  void
  run_projection_compatibility_contract_test()
  {
    const auto layout = std::make_shared<const contract::BlockLayout>(
      "projection_compatibility",
      std::vector<contract::SpaceId>{{"control"}},
      std::vector<std::size_t>{2});
    dealii::DynamicSparsityPattern dynamic_pattern(2, 2);
    for (std::size_t row = 0; row < 2; ++row)
      for (std::size_t column = 0; column < 2; ++column)
        dynamic_pattern.add(row, column);
    dealii::SparsityPattern pattern;
    pattern.copy_from(dynamic_pattern);

    auto diagonal_matrix =
      std::make_shared<dealii::SparseMatrix<double>>(pattern);
    diagonal_matrix->set(0, 0, 2.0);
    diagonal_matrix->set(1, 1, 2.0);
    const dealii_backend::MassMetric cellwise_metric(
      "l2_cellwise", layout, diagonal_matrix);
    const dealii_backend::CellwiseBoxConstraint cellwise_box(
      layout, 0.0, 1.0, cellwise_metric);
    contract::require(cellwise_box.supports_projection_in(cellwise_metric),
                      "Cellwise box rejected its coupled diagonal metric");

    auto non_diagonal_matrix =
      std::make_shared<dealii::SparseMatrix<double>>(pattern);
    non_diagonal_matrix->set(0, 0, 2.0);
    non_diagonal_matrix->set(0, 1, 1.0);
    non_diagonal_matrix->set(1, 0, 1.0);
    non_diagonal_matrix->set(1, 1, 2.0);
    const dealii_backend::MassMetric spoofed_cellwise_metric(
      "l2_cellwise", layout, non_diagonal_matrix);
    contract::require(
      !cellwise_box.supports_projection_in(spoofed_cellwise_metric),
      "A non-diagonal deal.II metric obtained cellwise clipping by reusing the l2_cellwise display identifier");
    test_support::require_contract_error(
      [&layout, &spoofed_cellwise_metric]() {
        (void)dealii_backend::CellwiseBoxConstraint(
          layout, 0.0, 1.0, spoofed_cellwise_metric);
      },
      "Cellwise box projection needs a positive diagonal metric realization",
      "non-diagonal cellwise projection coupling");

    const dealii_backend::MassMetric facewise_metric(
      "l2_facewise", layout, diagonal_matrix);
    const dealii_backend::FacewiseBoxConstraint facewise_box(
      layout, 0.0, 1.0, facewise_metric);
    contract::require(facewise_box.supports_projection_in(facewise_metric),
                      "Facewise box rejected its coupled diagonal metric");
    const dealii_backend::MassMetric spoofed_facewise_metric(
      "l2_facewise", layout, non_diagonal_matrix);
    contract::require(
      !facewise_box.supports_projection_in(spoofed_facewise_metric),
      "A non-diagonal deal.II metric obtained facewise clipping by reusing the l2_facewise display identifier");
    const dealii_backend::MassMetric h1_metric(
      "h1_continuous", layout, non_diagonal_matrix);
    contract::require(!cellwise_box.supports_projection_in(h1_metric),
                      "Cellwise clipping accepted an H1 metric realization");
  }

  void
  run_backend_size_conversion_contract_test()
  {
    const std::size_t maximum =
      dealii_backend::SerialBackend::maximum_native_size();
    contract::require(
      static_cast<std::size_t>(
        dealii_backend::SerialBackend::checked_native_size(maximum)) == maximum,
      "Serial deal.II size conversion rejected its native maximum");

    if constexpr (dealii_backend::SerialBackend::native_size_is_narrower)
      test_support::require_contract_error(
        [maximum]() {
          (void)dealii_backend::SerialBackend::checked_native_size(maximum + 1);
        },
        "Serial deal.II vector size exceeds its native range",
        "oversized serial deal.II vector dimension");
  }

  template <int dim>
  void
  run_continuous_neumann_control_metric_contract_test()
  {
    for (unsigned int family = 0; family < 2; ++family)
      {
        const bool simplex = family == 1;
        dealii::Triangulation<dim> triangulation;
        if (simplex)
          dealii::GridGenerator::subdivided_hyper_cube_with_simplices(
            triangulation,
            2);
        else
          dealii::GridGenerator::subdivided_hyper_cube(triangulation, 2);
        for (auto cell = triangulation.begin_active();
             cell != triangulation.end();
             ++cell)
          for (unsigned int face = 0; face < cell->n_faces(); ++face)
            if (cell->face(face)->at_boundary())
              {
                const auto center = cell->face(face)->center();
                cell->face(face)->set_boundary_id(
                  center[0] > 1.0 - 1e-12 || center[1] > 1.0 - 1e-12 ?
                    1 :
                    0);
              }

        const compiler::v1::detail::ContinuousNeumannControlRealisation<dim>
          realisation(triangulation, {1}, 3);
        const auto &coordinates = realisation.coordinates();
        contract::require(realisation.dimension() == 5 &&
                            realisation.physical_dimension() == 5 &&
                            realisation.independent_dimension() == 5 &&
                            coordinates.size() == 5,
                          "Continuous Neumann trace has the wrong dimension");
        const auto coordinate_count = [&coordinates](
                                        const dealii::Point<dim> &target) {
          return std::count_if(
            coordinates.begin(),
            coordinates.end(),
            [&target](const dealii::Point<dim> &coordinate) {
              return coordinate.distance(target) < 1e-14;
            });
        };
        contract::require(
          coordinate_count(dealii::Point<dim>(1.0, 1.0)) == 1,
          "Continuous Neumann trace did not deduplicate its connected corner");
        contract::require(
          coordinate_count(dealii::Point<dim>(0.0, 1.0)) == 1 &&
            coordinate_count(dealii::Point<dim>(1.0, 0.0)) == 1,
          "Continuous Neumann trace did not include selected-boundary endpoints");

        const dealii_backend::MassMetric metric(
          "l2_neumann_trace",
          realisation.layout(),
          realisation.control_mass_matrix());
        contract::require(
          metric.id() == "l2_neumann_trace" &&
            metric.layout()->dimension(0) == realisation.dimension(),
          "Continuous Neumann trace metric lost its identity or layout");
        dealii::Vector<double> unit_values(realisation.dimension());
        unit_values = 1.0;
        const Primal unit(realisation.layout(), {unit_values});
        require_close(contract::pair(metric.apply(unit), unit),
                      2.0,
                      1e-13,
                      "Continuous Neumann trace mass changed boundary measure");

        std::size_t corner_index = coordinates.size();
        for (std::size_t index = 0; index < coordinates.size(); ++index)
          if (coordinates[index].distance(dealii::Point<dim>(1.0, 1.0)) <
              1e-14)
            corner_index = index;
        contract::require(corner_index < coordinates.size(),
                          "Continuous Neumann trace omitted its connected corner");
        dealii::Vector<double> basis_values(realisation.dimension());
        basis_values[corner_index] = 1.0;
        const Primal basis(realisation.layout(), {basis_values});
        const Covector basis_action = metric.apply(basis);
        bool has_off_diagonal_mass = false;
        for (std::size_t index = 0; index < realisation.dimension(); ++index)
          if (index != corner_index &&
              std::abs(basis_action.block(0)[index]) > 1e-14)
            has_off_diagonal_mass = true;
        contract::require(
          has_off_diagonal_mass &&
            !metric.supports_coefficientwise_box_projection(),
          "Continuous Neumann trace mass is not the consistent nodal matrix");

        dealii::Vector<double> first_values(realisation.dimension());
        dealii::Vector<double> second_values(realisation.dimension());
        for (std::size_t index = 0; index < realisation.dimension(); ++index)
          {
            first_values[index] =
              (index % 2 == 0 ? 0.2 : -0.15) * static_cast<double>(index + 1);
            second_values[index] =
              (index % 2 == 0 ? -0.1 : 0.25) * static_cast<double>(index + 1);
          }
        const Primal first(realisation.layout(), {first_values});
        const Primal second(realisation.layout(), {second_values});
        const Covector first_action = metric.apply(first);
        const Covector second_action = metric.apply(second);
        require_close(contract::pair(first_action, second),
                      contract::pair(second_action, first),
                      1e-13,
                      "Continuous Neumann trace mass is not symmetric");
        contract::require(contract::pair(first_action, first) > 1e-12,
                          "Continuous Neumann trace mass is not positive definite");
        const Primal recovered = metric.inverse_apply(first_action);
        dealii::Vector<double> recovery_error = recovered.block(0);
        recovery_error.add(-1.0, first.block(0));
        require_close(recovery_error.l2_norm(),
                      0.0,
                      1e-11,
                      "Continuous Neumann trace metric apply/inverse mismatch");
      }
  }

  template <int dim>
  void
  run_continuous_neumann_control_output_contract_test()
  {
    for (unsigned int family = 0; family < 2; ++family)
      {
        const bool simplex = family == 1;
        dealii::Triangulation<dim> triangulation;
        if (simplex)
          dealii::GridGenerator::subdivided_hyper_cube_with_simplices(
            triangulation,
            2);
        else
          dealii::GridGenerator::subdivided_hyper_cube(triangulation, 2);
        for (auto cell = triangulation.begin_active();
             cell != triangulation.end();
             ++cell)
          for (unsigned int face = 0; face < cell->n_faces(); ++face)
            if (cell->face(face)->at_boundary())
              {
                const auto center = cell->face(face)->center();
                cell->face(face)->set_boundary_id(
                  center[0] > 1.0 - 1e-12 || center[1] > 1.0 - 1e-12 ?
                    1 :
                    0);
              }

        const compiler::v1::detail::ContinuousNeumannControlRealisation<dim>
          realisation(triangulation, {1}, 3);
        dealii::Vector<double> output_values(realisation.dimension());
        for (std::size_t index = 0; index < output_values.size(); ++index)
          output_values[index] =
            (index % 2 == 0 ? 0.125 : -0.25) *
            static_cast<double>(index + 1);
        const auto output_path =
          std::filesystem::temp_directory_path() /
          (simplex ? "nmopt-continuous-neumann-simplex.vtu" :
                     "nmopt-continuous-neumann-hypercube.vtu");
        std::filesystem::remove(output_path);
        realisation.write_native_output(output_path, output_values);
        std::ifstream output(output_path);
        contract::require(
          static_cast<bool>(output),
          "Continuous Neumann trace output file was not created");
        const std::string document((std::istreambuf_iterator<char>(output)),
                                   std::istreambuf_iterator<char>());
        output.close();

        const auto data_array_values = [&document](const std::string &marker) {
          const auto marker_position = document.find(marker);
          contract::require(
            marker_position != std::string::npos,
            "Continuous Neumann trace output omitted a data array");
          const auto values_begin = document.find('\n', marker_position);
          const auto values_end = document.find("</DataArray>", values_begin);
          contract::require(
            values_begin != std::string::npos &&
              values_end != std::string::npos,
            "Continuous Neumann trace output contains a malformed data array");
          std::istringstream values_stream(
            document.substr(values_begin + 1,
                            values_end - values_begin - 1));
          std::vector<double> values;
          for (double value = 0.0; values_stream >> value;)
            values.push_back(value);
          return values;
        };
        const auto written_control =
          data_array_values("Name=\"control\"");
        const auto written_points =
          data_array_values("NumberOfComponents=\"3\"");
        const auto written_connectivity =
          data_array_values("Name=\"connectivity\"");
        const auto written_offsets =
          data_array_values("Name=\"offsets\"");
        const auto written_types = data_array_values("Name=\"types\"");
        contract::require(
          document.find(
            "<Piece NumberOfPoints=\"5\" NumberOfCells=\"4\">") !=
              std::string::npos &&
            document.find("<PointData Scalars=\"control\">") !=
              std::string::npos &&
            document.find("<CellData>\n</CellData>") !=
              std::string::npos &&
            written_control.size() == output_values.size() &&
            written_points.size() == 3 * output_values.size() &&
            written_connectivity.size() == 8 &&
            written_offsets.size() == 4 && written_types.size() == 4,
          "Continuous Neumann trace output has the wrong topology or data association");
        for (std::size_t index = 0; index < output_values.size(); ++index)
          require_close(written_control[index],
                        output_values[index],
                        1e-15,
                        "Continuous Neumann trace output changed a nodal value");
        std::vector<bool> used_points(output_values.size(), false);
        for (const double value : written_connectivity)
          {
            contract::require(
              value >= 0.0 && value == std::floor(value) &&
                value < static_cast<double>(used_points.size()),
              "Continuous Neumann trace output has invalid connectivity");
            used_points[static_cast<std::size_t>(value)] = true;
          }
        contract::require(
          static_cast<std::size_t>(
            std::count(used_points.begin(), used_points.end(), true)) ==
            used_points.size(),
          "Continuous Neumann trace output does not connect every nodal point");
        for (std::size_t cell = 0; cell < written_offsets.size(); ++cell)
          contract::require(
            written_offsets[cell] == 2.0 * static_cast<double>(cell + 1) &&
              written_types[cell] == 3.0,
            "Continuous Neumann trace output does not contain VTK line cells");
        std::filesystem::remove(output_path);
      }
  }

  void
  run_serial_spd_reporting_contract_test()
  {
    dealii::DynamicSparsityPattern dynamic_pattern(2, 2);
    for (std::size_t row = 0; row < 2; ++row)
      for (std::size_t column = 0; column < 2; ++column)
        dynamic_pattern.add(row, column);
    dealii::SparsityPattern sparsity;
    sparsity.copy_from(dynamic_pattern);
    dealii::SparseMatrix<double> matrix(sparsity);
    matrix.set(0, 0, 4.0);
    matrix.set(0, 1, 1.0);
    matrix.set(1, 0, 1.0);
    matrix.set(1, 1, 3.0);
    dealii::Vector<double> right_hand_side(2);
    right_hand_side[0] = 1.0;
    right_hand_side[1] = 2.0;
    dealii::Vector<double> approximate_solution(2);
    const auto failed_report = dealii_backend::solve_serial_spd(
      matrix,
      approximate_solution,
      right_hand_side,
      dealii_backend::SPDLinearSolvePolicy{1, 1e-15, 1e-15});
    contract::require(!failed_report.converged() &&
                        failed_report.maximum_iterations == 1,
                      "serial SPD service did not report deliberate nonconvergence");

    dealii::Vector<double> zero_right_hand_side(2);
    dealii::Vector<double> zero_solution(2);
    const auto zero_report = dealii_backend::solve_serial_spd(
      matrix, zero_solution, zero_right_hand_side, {});
    contract::require(zero_report.converged() &&
                        zero_solution.l2_norm() == 0.0,
                      "serial SPD service did not handle a zero right-hand side");
  }

} // namespace

int
main(const int argc, char **argv)
{
  try
    {
      const std::vector<test_support::Scenario> scenarios{
        {"hminus1_metric",
         "nmopt.dealii.hminus1_metric",
         {"dealii", "contract", "metric"},
         30,
         run_hminus1_metric_contract_test},
        {"quadratic_form",
         "nmopt.dealii.quadratic_form",
         {"dealii", "contract", "objective"},
         30,
         run_quadratic_form_contract_test},
        {"continuous_neumann_control_metric",
         "nmopt.dealii.continuous_neumann_control_metric",
         {"dealii", "contract", "metric"},
         30,
         []() { run_continuous_neumann_control_metric_contract_test<2>(); }},
        {"continuous_neumann_control_output",
         "nmopt.dealii.continuous_neumann_control_output",
         {"dealii", "contract", "output"},
         30,
         []() { run_continuous_neumann_control_output_contract_test<2>(); }},
        {"projection_compatibility",
         "nmopt.dealii.projection_compatibility",
         {"dealii", "contract", "constraint"},
         60,
         run_projection_compatibility_contract_test},
        {"serial_spd_reporting",
         "nmopt.dealii.serial_spd_reporting",
         {"dealii", "contract", "solver"},
         30,
         run_serial_spd_reporting_contract_test},
        {"backend_size_conversion",
         "nmopt.dealii.backend_size_conversion",
         {"dealii", "contract", "adapter"},
         30,
         run_backend_size_conversion_contract_test}};
      const auto result = test_support::run_requested_scenarios(
        argc, argv, scenarios, std::cout);
      if (!result.listed)
        std::cout << "deal.II diffusion DTO contract scenario passed: "
                  << result.executed << '\n';
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "deal.II diffusion DTO contract test failed: "
                << exception.what() << '\n';
      return 1;
    }
}
