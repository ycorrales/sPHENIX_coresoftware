#ifndef MVTX_NOISEMAP_H
#define MVTX_NOISEMAP_H

/*!
 * @file MvtxNoiseMap.h
 * @author Y. Corrales <ycmorales@bnl.gov>
 * @date Apr 2024
 * @brief NoiseMap class for the MVTX
 */

#include <TNamed.h>
#include <sys/types.h>
#include "Rtypes.h"

#include <cassert>
#include <climits>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

class MvtxNoiseMap : public TNamed
{
 public:
  using Map = std::map<int, int>;

  //! Constructor, initializing values for position, charge and readout frame
  //! ctor
  MvtxNoiseMap() = default;
  explicit MvtxNoiseMap(std::vector<Map>& noise)
    : TNamed("mvtx_noiseMap", "mvtx_noiseMap")
  {
    mNoisyPixels.swap(noise);
  }
  explicit MvtxNoiseMap(int nchips)
    : TNamed("mvtx_noiseMap", "mvtx_noiseMap")
  {
    mNoisyPixels.assign(nchips, Map());
  }

  //! cp/mv ctor
  MvtxNoiseMap(MvtxNoiseMap&) = default;
  MvtxNoiseMap(MvtxNoiseMap&&) = default;

  //! cp/mv assignment
  MvtxNoiseMap& operator=(const MvtxNoiseMap&) = default;
  MvtxNoiseMap& operator=(MvtxNoiseMap&&) = default;

  //! dtor
  ~MvtxNoiseMap() override = default;

  //! Get the noise level for this pixels
  float getNoiseLevel(int chip, int row, int col) const
  {
    assert(chip < (int) mNoisyPixels.size());
    const auto keyIt = mNoisyPixels[chip].find(static_cast<int>(getKey(row, col)));
    if (keyIt != mNoisyPixels[chip].end())
    {
      return static_cast<float>(keyIt->second);
    }
    return 0;
  }

  void increaseNoiseCount(int lay, int stv, int chipId, int row, int col)
  {
    auto chip = getChipId(lay, stv, chipId);
    increaseNoiseCount(chip, row, col);
  }

  void increaseNoiseCount(int chip, int row, int col)
  {
    assert(chip < (int) mNoisyPixels.size());
    mNoisyPixels[chip][static_cast<int>(getKey(row, col))]++;
  }

  void increaseNoiseCount(int chip, const std::vector<int>& rowcolKey)
  {
    assert(chip < (int) mNoisyPixels.size());
    auto& ch = mNoisyPixels[chip];
    for (const auto k : rowcolKey)
    {
      ch[k]++;
    }
  }

  int dumpAboveThreshold(int t = 3) const
  {
    int n = 0;
    auto chipID = mNoisyPixels.size();
    while (chipID--)
    {
      const auto& map = mNoisyPixels[chipID];
      for (const auto& pair : map)
      {
        if (pair.second <= t)
        {
          continue;
        }
        n++;
        auto key = pair.first;
        auto row = key2Row(key);
        auto col = key2Col(key);
        std::cout << "chip, row, col, noise: " << chipID << ' ' << row << ' ' << col << ' ' << pair.second << '\n';
      }
    }
    return n;
  }
  int dumpAboveProbThreshold(float p = 1e-7) const
  {
    return dumpAboveThreshold(std::ceil(p * static_cast<float>(mNumOfStrobes)));
  }

  void applyProbThreshold(float t, long int n, float relErr = 0.2F, int minChipID = 0, int maxChipID = 431)
  {
    // Remove from the maps all pixels with the firing probability below the threshold
    // Apply the cut only for chips between minChipID and maxChipID (included)
    if (n < 1)
    {
      std::cerr << "Cannot apply threshold with " << n << " ROFs scanned" << std::endl;
      return;
    }
    mProbThreshold = t;
    mNumOfStrobes = n;
    float minFiredForErr = 0.F;
    if (relErr > 0)
    {
      minFiredForErr = relErr * relErr - 1.F / static_cast<float>(n);
      if (minFiredForErr <= 0.F)
      {
        std::cerr << "Noise threshold " << t << " with relative error " << relErr << " is not reachable with " << n << " ROFs processed, mask all permanently fired pixels" << std::endl;
        minFiredForErr = static_cast<float>(n);
      }
      else
      {
        minFiredForErr = 1.F / minFiredForErr;
      }
    }
    int minFired = std::ceil(std::max(t * static_cast<float>(mNumOfStrobes), minFiredForErr));  // min number of fired pixels exceeding requested threshold
    auto req = getMinROFs(t, relErr);
    if (n < req)
    {
      mProbThreshold = float(minFired) / static_cast<float>(n);
      std::cerr << "Requested relative error " << relErr << " with prob.threshold " << t << " needs > " << req << " ROFs, " << n << " provided: pixels with noise >" << mProbThreshold << " will be masked" << std::endl;
    }

    int currChipID = 0;
    for (auto& map : mNoisyPixels)
    {
      if (currChipID < minChipID || currChipID > maxChipID)
      {  // chipID range
        currChipID++;
        continue;
      }
      for (auto it = map.begin(); it != map.end();)
      {
        if (it->second < minFired)
        {
          it = map.erase(it);
        }
        else
        {
          ++it;
        }
      }
      currChipID++;
    }
  }
  float getProbThreshold() const { return mProbThreshold; }
  long int getNumOfStrobes() const { return mNumOfStrobes; }

  bool isNoisy(int chip, int row, int col) const
  {
    assert(chip < (int) mNoisyPixels.size());
    return (mNoisyPixels[chip].find(static_cast<int>(getKey(row, col))) != mNoisyPixels[chip].end());
  }

  bool isNoisyOrFullyMasked(int chip, int row, int col) const
  {
    assert(chip < (int) mNoisyPixels.size());
    return isNoisy(chip, row, col) || isFullChipMasked(chip);
  }

  bool isNoisy(int chip) const
  {
    assert(chip < (int) mNoisyPixels.size());
    return !mNoisyPixels[chip].empty();
  }

  // Methods required by the calibration framework
  void print();
  //  void fill(const gsl::span<const CompClusterExt> data);
  //  void merge(const MvtxNoiseMap* prev) {}

  const Map* getChipMap(int chip) const { return chip < (int) mNoisyPixels.size() ? &mNoisyPixels[chip] : nullptr; }
  Map& getChip(int chip) { return mNoisyPixels[chip]; }
  const Map& getChip(int chip) const { return mNoisyPixels[chip]; }

  void maskFullChip(int chip, bool cleanNoisyPixels = false)
  {
    if (cleanNoisyPixels)
    {
      resetChip(chip);
    }
    increaseNoiseCount(chip, -1, -1);
  }

  bool isFullChipMasked(int chip) const
  {
    return isNoisy(chip, -1, -1);
  }

  void resetChip(int chip)
  {
    assert(chip < (int) mNoisyPixels.size());
    mNoisyPixels[chip].clear();
  }

  static long getMinROFs(float t, float relErr)
  {
    // calculate min number of ROFs needed to reach threshold t with relative error relErr
    relErr = relErr >= 0.F ? relErr : 0.1F;
    t = t >= 0.F ? t : 1e-6F;
    return std::ceil((1. + 1. / t) / (relErr * relErr));
  }

  size_t size() const { return mNoisyPixels.size(); }
  void setNumOfStrobes(long n) { mNumOfStrobes = n; }
  void addStrobes(long n) { mNumOfStrobes += n; }
  static uint32_t getKey(int row, int col) { return (static_cast<uint32_t>(row) << SHIFT) + static_cast<uint32_t>(col); }
  static uint32_t key2Row(uint32_t key) { return key >> SHIFT; }
  static uint32_t key2Col(uint32_t key) { return key & MASK; }
  static uint32_t getChipId(uint8_t layId, uint8_t stvId, uint8_t chipId)
  {
    return chipZeroPerLayer[layId] + (stvId * 9) + chipId;
  }

  void read(const std::string& /*filename*/);
  void write(const std::string& /*filename*/) const;

 private:
  static constexpr uint8_t chipZeroPerLayer[3] = {0, 12 * 9, 16 * 9};
  static constexpr uint32_t SHIFT = 10, MASK = (0x1U << SHIFT) - 1U;
  std::vector<Map> mNoisyPixels;  ///< Internal noise map representation
  long int mNumOfStrobes = 0;     ///< Accumulated number of ALPIDE strobes
  float mProbThreshold = 0;       ///< Probability threshold for noisy pixels

  ClassDefOverride(MvtxNoiseMap, 1);
};

#endif /* MVTX_NOISEMAP_H */
