#ifndef FUN4ALLRAW_SINGLEMVTXPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMVTXPOOLINPUT_H

#include "SingleStreamingInput.h"

#include <algorithm>
#include <map>
#include <vector>

class MvtxRawHit;
class Packet;
class mvtx_pool;

class SingleMvtxPoolInput : public SingleStreamingInput
{
 public:
  explicit SingleMvtxPoolInput(const std::string &name);
  ~SingleMvtxPoolInput() override;
  void FillPool(const uint64_t minBCO) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  unsigned int GetBcoRange() const { return m_BcoRange; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

  std::set<int> &getFeeIdSet(const uint64_t &bco) { return m_BeamClockFEE[bco]; };
  std::set<uint64_t>& getGtmL1BcoSet() { return m_gtmL1BcoSetRef; }
  const std::map<int, std::set<uint64_t>>& getFeeGTML1BCOMap() const { return m_FeeGTML1BCOMap; }
  void clearGtmL1BcoSet() { m_gtmL1BcoSetRef.clear(); }
  void clearFeeGTML1BCOMap() { 
    for(auto& [key, set] : m_FeeGTML1BCOMap)
    {
      set.clear();
    }
    m_FeeGTML1BCOMap.clear(); 
  }
 protected:

 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};

  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MvtxRawHit *>> m_MvtxRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::map<int, uint64_t> m_FeeStrobeMap;
  std::set<uint64_t> m_BclkStack;
  std::set<uint64_t> gtmL1BcoSet;  // GTM L1 BCO
  std::set<uint64_t> m_gtmL1BcoSetRef;
  std::map<int, std::set<uint64_t>> m_FeeGTML1BCOMap;
  std::map<int, mvtx_pool *> poolmap;
};

#endif
