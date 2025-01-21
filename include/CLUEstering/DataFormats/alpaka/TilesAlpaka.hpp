
#pragma once

#include <alpaka/core/Common.hpp>
#include <alpaka/alpaka.hpp>
#include <cstddef>
#include <cstdint>
#include <stdint.h>

#include "../../AlpakaCore/alpakaWorkDiv.hpp"
#include "../../AlpakaCore/alpakaConfig.hpp"
#include "../../AlpakaCore/alpakaMemory.hpp"
#include "AlpakaVecArray.hpp"

using clue::VecArray;

constexpr uint32_t max_tile_depth{1 << 10};
constexpr uint32_t max_n_tiles{1 << 15};

namespace ALPAKA_ACCELERATOR_NAMESPACE_CLUE {

  template <uint8_t Ndim>
  class CoordinateExtremes {
  private:
    float m_data[2 * Ndim];

  public:
    CoordinateExtremes() = default;

    ALPAKA_FN_HOST_ACC const float* data() const { return m_data; }
    ALPAKA_FN_HOST_ACC float* data() { return m_data; }

    ALPAKA_FN_HOST_ACC float min(int i) const { return m_data[2 * i]; }
    ALPAKA_FN_HOST_ACC float& min(int i) { return m_data[2 * i]; }
    ALPAKA_FN_HOST_ACC float max(int i) const { return m_data[2 * i + 1]; }
    ALPAKA_FN_HOST_ACC float& max(int i) { return m_data[2 * i + 1]; }
    ALPAKA_FN_HOST_ACC float range(int i) const { return m_data[2 * i + 1] - m_data[2 * i]; }
    ALPAKA_FN_HOST_ACC float& range(int i) { return m_data[2 * i + 1] - m_data[2 * i]; }
  };

  template <typename TAcc, uint8_t Ndim>
  class TilesAlpaka {
  public:
    TilesAlpaka() {};

    ALPAKA_FN_HOST_ACC inline void setAcc(const TAcc& Acc) { acc = Acc; };

    ALPAKA_FN_HOST_ACC inline constexpr const float* minMax() const {
      return min_max.data();
    }
    ALPAKA_FN_HOST_ACC inline constexpr float* minMax() { return min_max.data(); }

    ALPAKA_FN_HOST_ACC inline constexpr const float* tileSize() const {
      return tile_size;
    }
    ALPAKA_FN_HOST_ACC inline constexpr float* tileSize() { return tile_size; }

    ALPAKA_FN_HOST_ACC void resizeTiles(std::size_t nTiles, int nPerDim) {
      this->n_tiles = nTiles;
      this->n_tiles_per_dim = nPerDim;

      this->m_tiles.resize(nTiles);
    }

    ALPAKA_FN_HOST_ACC inline constexpr int getBin(float coord_,
                                                   int dim_) const {
      int coord_Bin;
      if (wrapped[dim_]) {
        coord_Bin =
            static_cast<int>((normalizeCoordinate(coord_, dim_) - min_max.min(dim_)) / tile_size[dim_]);
      } else {
        coord_Bin =
            static_cast<int>((coord_ - min_max.min(dim_)) / tile_size[dim_]);
      }


      // Address the cases of underflow and overflow
      coord_Bin = alpaka::math::min(acc, coord_Bin, n_tiles_per_dim - 1);
      coord_Bin = alpaka::math::max(acc, coord_Bin, 0);

      return coord_Bin;
    }

    ALPAKA_FN_HOST_ACC inline constexpr int getGlobalBin(const float* coords) const {
      int globalBin = 0;
      for (int dim = 0; dim != Ndim - 1; ++dim) {
        globalBin += alpaka::math::pow(acc, n_tiles_per_dim, Ndim - dim - 1) *
                     getBin(acc, coords[dim], dim);
      }
      globalBin += getBin(acc, coords[Ndim - 1], Ndim - 1);
      return globalBin;
    }

    ALPAKA_FN_HOST_ACC inline constexpr int getGlobalBinByBin(
        const VecArray<uint32_t, Ndim>& Bins) const {
      uint32_t globalBin = 0;
      for (int dim = 0; dim != Ndim - 1; ++dim) {
        auto bin_i = wrapped[dim] ? (Bins[dim] % n_tiles_per_dim) : Bins[dim];
        globalBin += alpaka::math::pow(acc, n_tiles_per_dim, Ndim - dim - 1) * bin_i;
      }
      globalBin += wrapped[Ndim] ? (Bins[Ndim - 1] % n_tiles_per_dim) : Bins[Ndim - 1];
      return globalBin;
    }

    ALPAKA_FN_ACC inline constexpr void fill(const float* coords,
                                             int i) {
      m_tiles[getGlobalBin(acc, coords)].push_back(acc, i);
    }

    ALPAKA_FN_ACC inline void searchBox(
        const VecArray<VecArray<float, 2>, Ndim>& sb_extremes,
        VecArray<VecArray<uint32_t, 2>, Ndim>* search_box) {
      for (int dim{}; dim != Ndim; ++dim) {
        VecArray<uint32_t, 2> dim_sb;
        auto infBin = getBin(acc, sb_extremes[dim][0], dim);
        auto supBin = getBin(acc, sb_extremes[dim][1], dim);
        if (wrapped[dim] and infBin > supBin)
          supBin += n_tiles_per_dim;
        dim_sb.push_back_unsafe(infBin);
        dim_sb.push_back_unsafe(supBin);

        search_box->push_back_unsafe(dim_sb);
      }
    }

    ALPAKA_FN_HOST_ACC inline constexpr float distance(const float* coord_i, const float* coord_j) {
      float dist_sq = 0.f;
      for (int dim = 0; dim != Ndim; ++dim) {
        if (wrapped[dim])
          dist_sq += normalizeCoordinate(coord_i - coord_j, dim) * normalizeCoordinate(coord_i - coord_j, dim);
        else
          dist_sq += (coord_i - coord_j) * (coord_i - coord_j);
      }
      return dist_sq;
    }

    ALPAKA_FN_HOST_ACC inline constexpr auto size() { return n_tiles; }

    ALPAKA_FN_HOST_ACC inline constexpr int nPerDim() const { return n_tiles_per_dim; }

    ALPAKA_FN_HOST_ACC inline constexpr void clear() {
      for (int i{}; i < n_tiles; ++i) {
        m_tiles[i].reset();
      }
    }

    ALPAKA_FN_HOST_ACC inline constexpr void clear(uint32_t i) { m_tiles[i].reset(); }

    ALPAKA_FN_HOST_ACC inline constexpr VecArray<uint32_t, max_tile_depth>& operator[](
        int globalBinId) {
      return m_tiles[globalBinId];
    }

  private:

    float normalizeCoordinate(float coord, int dim){
      const float range = min_max.range(dim);
      float remainder = coord - static_cast<int>(coord / range) * range;
      if (remainder >= min_max.max(dim))
        remainder -= range;
      else if (remainder < min_max.min(dim))
        remainder += range;
      return remainder;
    }

    TAcc acc;
    std::size_t n_tiles;
    int n_tiles_per_dim;
    int wrapped[Ndim];
    CoordinateExtremes<Ndim> min_max;
    float tile_size[Ndim];
    VecArray<VecArray<uint32_t, max_tile_depth>, max_n_tiles> m_tiles;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE_CLUE
