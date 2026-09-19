#pragma once

// Compatibility aggregate for the first v1 semantic slice. New users may
// include only types.hpp, validation.hpp, or problem_library.hpp when they do
// not need the complete public surface.
#include "nmopt/semantic/v1/problem_library.hpp"
#include "nmopt/semantic/v1/resolved_problem.hpp"
#include "nmopt/semantic/v1/types.hpp"
#include "nmopt/semantic/v1/validation.hpp"
