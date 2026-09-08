#include "tutorial_application.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int
main()
{
  try
    {
      TutorialApplication application;
      TutorialApplication::Vector control(application.control_dimension());
      control = 0.0;

      const auto state = application.solve_state(control);
      if (!(application.residual(state, control).l2_norm() < 1e-10))
        throw std::runtime_error("tutorial binding state residual is too large");

      const auto derivative = application.objective_derivative(state, control);
      if (derivative.state.size() != application.state_dimension() ||
          derivative.control.size() != application.control_dimension())
        throw std::runtime_error(
          "tutorial binding objective derivative has the wrong layout");

      const auto metric = application.control_metric_matrix();
      if (!metric || metric->m() != application.control_dimension() ||
          metric->n() != application.control_dimension())
        throw std::runtime_error("tutorial binding metric has the wrong size");

      const auto output_directory =
        std::filesystem::path("runs/external-dealii/step-4/binding");
      application.write_native_output(output_directory, state);
      if (!std::filesystem::exists(output_directory / "fields-volume.vtu"))
        throw std::runtime_error("tutorial binding native output is missing");

      std::cout << "external tutorial binding smoke passed (state="
                << application.state_dimension()
                << ", control=" << application.control_dimension() << ")\n";
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "external tutorial binding smoke failed: "
                << exception.what() << '\n';
      return 1;
    }
}
