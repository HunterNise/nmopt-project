#pragma once

#include <deal.II/grid/tria.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace nmopt::application::chapter6::dealii::detail
{
  inline std::uint64_t
  centroid_split_selection_score(const unsigned int triangle_index,
                                 const unsigned int selection_seed)
  {
    std::uint64_t value =
      static_cast<std::uint64_t>(triangle_index) ^
      (static_cast<std::uint64_t>(selection_seed) << 32);
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  }

  inline std::unique_ptr<::dealii::Triangulation<2>>
  make_centroid_split_simplex_mesh(
    const ::dealii::Triangulation<2> &base_mesh,
    const unsigned int                centroid_splits,
    const unsigned int                selection_seed)
  {
    if (!base_mesh.all_reference_cells_are_simplex())
      throw std::invalid_argument(
        "centroid splitting requires a triangular base mesh");
    if (centroid_splits == 0 ||
        centroid_splits > base_mesh.n_active_cells())
      throw std::invalid_argument(
        "centroid splitting needs a valid positive split count");

    struct BaseTriangle
    {
      std::array<unsigned int, 3> vertices;
      ::dealii::types::material_id material_id;
    };

    std::vector<BaseTriangle> base_triangles;
    base_triangles.reserve(base_mesh.n_active_cells());
    for (const auto &cell : base_mesh.active_cell_iterators())
      {
        BaseTriangle triangle{{{cell->vertex_index(0),
                                cell->vertex_index(1),
                                cell->vertex_index(2)}},
                              cell->material_id()};
        const auto &a = base_mesh.get_vertices()[triangle.vertices[0]];
        const auto &b = base_mesh.get_vertices()[triangle.vertices[1]];
        const auto &c = base_mesh.get_vertices()[triangle.vertices[2]];
        const double signed_area =
          (b[0] - a[0]) * (c[1] - a[1]) -
          (b[1] - a[1]) * (c[0] - a[0]);
        if (signed_area == 0.0)
          throw std::invalid_argument(
            "centroid splitting requires nondegenerate base triangles");
        if (signed_area < 0.0)
          std::swap(triangle.vertices[1], triangle.vertices[2]);
        base_triangles.push_back(triangle);
      }

    std::vector<std::pair<std::uint64_t, unsigned int>> ranked_triangles;
    ranked_triangles.reserve(base_triangles.size());
    for (unsigned int index = 0; index < base_triangles.size(); ++index)
      ranked_triangles.emplace_back(
        centroid_split_selection_score(index, selection_seed), index);
    std::sort(ranked_triangles.begin(), ranked_triangles.end());

    std::vector<bool> split(base_triangles.size(), false);
    for (unsigned int rank = 0; rank < centroid_splits; ++rank)
      split[ranked_triangles[rank].second] = true;

    auto vertices = base_mesh.get_vertices();
    std::vector<::dealii::CellData<2>> cells;
    cells.reserve(base_triangles.size() + 2 * centroid_splits);
    const auto add_triangle = [&cells](const unsigned int a,
                                       const unsigned int b,
                                       const unsigned int c,
                                       const ::dealii::types::material_id material_id) {
        ::dealii::CellData<2> cell(3);
        cell.vertices = {a, b, c};
        cell.material_id = material_id;
        cells.push_back(std::move(cell));
      };

    for (unsigned int index = 0; index < base_triangles.size(); ++index)
      {
        const auto &triangle = base_triangles[index];
        const auto  a = triangle.vertices[0];
        const auto  b = triangle.vertices[1];
        const auto  c = triangle.vertices[2];
        if (!split[index])
          {
            add_triangle(a, b, c, triangle.material_id);
            continue;
          }

        const auto centroid =
          (vertices[a] + vertices[b] + vertices[c]) / 3.0;
        const auto centroid_index = static_cast<unsigned int>(vertices.size());
        vertices.push_back(centroid);
        add_triangle(a, b, centroid_index, triangle.material_id);
        add_triangle(b, c, centroid_index, triangle.material_id);
        add_triangle(c, a, centroid_index, triangle.material_id);
      }

    auto mesh = std::make_unique<::dealii::Triangulation<2>>();
    mesh->create_triangulation(vertices, cells, ::dealii::SubCellData{});
    return mesh;
  }
} // namespace nmopt::application::chapter6::dealii::detail
