#include "external_poisson_fixture.hpp"
#include "../support/scoped_temporary_directory.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int
main()
{
  try
    {
      external_poisson_application::PoissonControlApplication application;
      external_poisson_application::PoissonControlApplication::Vector control(
        application.control_dimension());
      const auto state = application.solve_state(control);

      if (!(application.residual_norm(state, control) < 1e-10))
        throw std::runtime_error(
          "standalone external Poisson state residual is too large");

      const nmopt::test_support::ScopedTemporaryDirectory temporary_directory(
        "external-poisson-standalone");
      const auto &output_directory = temporary_directory.path();
      application.write_native_output(output_directory, state);
      if (!std::filesystem::exists(output_directory / "fields-volume.vtu"))
        throw std::runtime_error(
          "standalone external Poisson output is missing");
      std::cout << "standalone external Poisson application passed\n";
      return 0;
    }
  catch (const std::exception &exception)
    {
      std::cerr << "standalone external Poisson application failed: "
                << exception.what() << '\n';
      return 1;
    }
}
