#ifndef MVTX_MVTXSAPRDFUNPACKER_H
#define MVTX_MVTXSAPRDFUNPACKER_H

/*!
 * @file MvtxSaPrdfUnpacker.h
 * @author Y. Corrales <ycorrale@cern.ch>
 * @date Mar-2025
 * @brief Prdf unpacker for Mvtx calibration
 */

#include "fun4all/SubsysReco.h"

// #include <algorithm>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>
// #include <vector>

class Packet;
class Event;
class mvtx_pool;
class MvtxNoiseMap;

namespace mvtx
{
  using _hit = struct _hit
  {
    uint16_t feeid = 0xFF;
    uint8_t chip_id = 0xF;
    uint16_t row_pos = 0xFFFF;
    uint16_t col_pos = 0xFFFF;
  };
};  // namespace mvtx

class MvtxSaPrdfUnpacker : public SubsysReco
{
 public:
  explicit MvtxSaPrdfUnpacker(const std::string &name = "MvtxSaPrdfUnpacker");

  MvtxSaPrdfUnpacker(const MvtxSaPrdfUnpacker &) = default;
  MvtxSaPrdfUnpacker &operator=(const MvtxSaPrdfUnpacker &) = default;
  MvtxSaPrdfUnpacker(MvtxSaPrdfUnpacker &&) = default;
  MvtxSaPrdfUnpacker &operator=(MvtxSaPrdfUnpacker &&) = default;

  ~MvtxSaPrdfUnpacker() override;

  //! global initialization
  int Init(PHCompositeNode * /*dummy*/) override;

  //! run initialization
  int InitRun(PHCompositeNode * /*dummy*/) override;

  //! event processing
  int process_event(PHCompositeNode * /*dummy*/) override;

  //! end of processing
  int End(PHCompositeNode * /*dummy*/) override;

  void SetOutFileName(const std::string &name)
  {
    m_outputfile = name;
  }

  std::string GetOutFileName() const
  {
    return m_outputfile;
  }

 protected:
 private:
  int GetNodes(PHCompositeNode *);

  Packet **plist{nullptr};
  Event *m_event{nullptr};

  std::map<int, mvtx_pool *> poolmap;

  std::string m_outputfile;
  std::map<uint64_t, size_t> m_feeidMap;
  std::map<uint64_t, std::vector<mvtx::_hit>> m_ChipMap;
};

#endif  // MVTX_MVTXSAPRDFUNPACKER_H
