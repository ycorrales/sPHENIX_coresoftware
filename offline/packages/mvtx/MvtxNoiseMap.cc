/*!
 * @file MvtxNoiseMap.h
 * @author Y. Corrales
 * @date Apr 2024
 * @brief  Definition NoiseMap class
 */

#include "MvtxNoiseMap.h"

#include <cdbobjects/CDBTTree.h>
#include <fstream>
#include <string>

ClassImp(MvtxNoiseMap);

namespace
{
  const std::string m_threshold("noiseThres");
  const std::string m_NumStrobes("numStrobes");
  const std::string m_total_entries_key("totalHotPixels");
  const std::string m_chip_key = "chip";
  const std::string m_px_key = "key";
  const std::string m_count_key = "count";
}  // namespace

//+-----------------------------------------------------------------------+
void MvtxNoiseMap::read(const std::string& filename)
{
  std::cout << "MvtxNoiseMap::read - filename: " << filename << std::endl;
  // clear existing data
  for (auto& map : mNoisyPixels)
  {
    map.clear();
  }
  // make sure file exists before loading, otherwise crashes
  if (!std::ifstream(filename.c_str()).good())
  {
    std::cout << "MvtxNoiseMap::read -"
              << " filename: " << filename << " does not exist."
              << " No calibration loaded" << std::endl;
    return;
  }
  // use generic CDBTree to load
  CDBTTree cdbttree(filename);
  cdbttree.LoadCalibrations();

  mNumOfStrobes = static_cast<int>(cdbttree.GetSingleUInt64Value(m_NumStrobes));
  mProbThreshold = cdbttree.GetSingleFloatValue(m_threshold);

  // read total number of hot channels
  const int m_total_entries = cdbttree.GetSingleIntValue(m_total_entries_key);
  for (int i = 0; i < m_total_entries; ++i)
  {
    // read channel id
    const int px_chip = cdbttree.GetIntValue(i, m_chip_key);
    const int px_key = cdbttree.GetIntValue(i, m_px_key);
    const int px_count = cdbttree.GetIntValue(i, m_count_key);
    if (std::isnan(px_chip) || std::isnan(px_key) || std::isnan(px_count))
    {
      continue;
    }
    mNoisyPixels[px_chip][px_key] = px_count;
  }
  std::cout << "MvtxNoiseMap::read - total entries: " << mNoisyPixels.size() << std::endl;
}

//+-----------------------------------------------------------------------+
void MvtxNoiseMap::write(const std::string& filename) const
{
  std::cout << "MvtxNoiseMap::write - filename: " << filename << std::endl;
  if (mNoisyPixels.empty())
  {
    return;
  }
  // use generic CDBTree to load
  CDBTTree cdbttree(filename);

  cdbttree.SetSingleUInt64Value(m_NumStrobes, mNumOfStrobes);
  cdbttree.SetSingleFloatValue(m_threshold, mProbThreshold);

  int chip_index = -1;
  int index = 0;
  for (const auto& map : mNoisyPixels)
  {
    ++chip_index;
    if (map.empty())
    {
      continue;
    }
    for (const auto& [px_key, counts] : map)
    {
      cdbttree.SetIntValue(index, m_chip_key, chip_index);
      cdbttree.SetIntValue(index, m_px_key, static_cast<int>(px_key));
      cdbttree.SetIntValue(index, m_count_key, static_cast<int>(counts));
      ++index;
    }
  }
  cdbttree.SetSingleIntValue(m_total_entries_key, index);

  // commit and write
  cdbttree.Commit();
  cdbttree.CommitSingle();
  cdbttree.WriteCDBTTree();
}

//+-----------------------------------------------------------------------+
// void NoiseMap::print()
//{
//   int nc = 0, np = 0, nm = 0;
//   for (const auto& map : mNoisyPixels) {
//     if (!map.empty()) {
//       nc++;
//     }
//     np += map.size();
//     if (map.find(getKey(-1, -1)) != map.end()) {
//       nm++;
//       nc--;
//       np--;
//     }
//   }
//   LOG(info) << "Number of fully maske chips " << nm;
//   LOG(info) << "Number of noisy chips: " << nc;
//   LOG(info) << "Number of noisy pixels: " << np;
//   LOG(info) << "Number of of strobes: " << mNumOfStrobes;
//   LOG(info) << "Probability threshold: " << mProbThreshold;
// }
//
// void NoiseMap::fill(const gsl::span<const CompClusterExt> data)
//{
//   for (const auto& c : data) {
//     if (c.getPatternID() != o2::itsmft::CompCluster::InvalidPatternID) {
//       // For the noise calibration, we use "pass1" clusters...
//       continue;
//     }
//
//     auto id = c.getSensorID();
//     auto row = c.getRow();
//     auto col = c.getCol();
//
//     // A simplified 1-pixel calibration
//     increaseNoiseCount(id, row, col);
//   }
// }
