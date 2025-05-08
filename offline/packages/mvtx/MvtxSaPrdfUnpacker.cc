/*!
 * @file MvtxSaPrdfUnpacker.cc
 * @author Y. Corrales <ycorrale@cern.ch>
 * @date Mar-2025
 * @brief Prdf unpacker for Mvtx calibration
 */

#include "MvtxSaPrdfUnpacker.h"

// #include "MvtxNoiseMap.h"
#include "fun4all/Fun4AllReturnCodes.h"
#include "fun4all/PHTFileServer.h"
#include "fun4allraw/MvtxRawDefs.h"
#include "fun4allraw/mvtx_pool.h"
#include "phool/PHCompositeNode.h"
#include "phool/getClass.h"

#include "Event/Event.h"
#include "Event/EventTypes.h"

#include <iostream>

//========================================================================+
MvtxSaPrdfUnpacker::MvtxSaPrdfUnpacker(const std::string &name)
  : SubsysReco(name)
  , plist(new Packet *[2])
{
  // m_rawHitContainerName = "MVTXRAWHIT";
}

//========================================================================+
MvtxSaPrdfUnpacker::~MvtxSaPrdfUnpacker()
{
  delete[] plist;
  for (auto &iter : poolmap)
  {
    if (Verbosity() > 2)
    {
      std::cout << "deleting mvtx pool for id " << iter.first << std::endl;
    }
    delete (iter.second);
  }
}

//========================================================================+
int MvtxSaPrdfUnpacker::Init(PHCompositeNode * /*topnode*/)
{
  // Create the output file and trees
  PHTFileServer::open(m_outputfile, "RECREATE");

  // m_noiseMap = new MvtxNoiseMap(48 * 9);

  return Fun4AllReturnCodes::EVENT_OK;
}

//========================================================================+
int MvtxSaPrdfUnpacker::InitRun(PHCompositeNode * /*topnode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//========================================================================+
int MvtxSaPrdfUnpacker::GetNodes(PHCompositeNode *topNode)
{
  m_event = findNode::getClass<Event>(topNode, "PRDF");
  if (!m_event)
  {
    std::cout << PHWHERE << " Event not found" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//========================================================================+
int MvtxSaPrdfUnpacker::process_event(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  if (m_event->getEvtType() != DATAEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
    // RunNumber(m_event->getRunNumber());
  }
  int EventSequence = m_event->getEvtSequence();
  int npackets = m_event->getPacketList(plist, 2);

  if (npackets > 2)
  {
    exit(1);
  }
  for (int i = 0; i < npackets; i++)
  {
    // Ignoring packet not from MVTX detector
    if (Verbosity() > 1)
    {
      plist[i]->identify();
    }

    if (poolmap.find(plist[i]->getIdentifier()) == poolmap.end())
    {
      if (Verbosity() > 1)
      {
        std::cout << "starting new mvtx pool for packet " << plist[i]->getIdentifier() << std::endl;
      }
      poolmap[plist[i]->getIdentifier()] = new mvtx_pool();
    }
    poolmap[plist[i]->getIdentifier()]->addPacket(plist[i]);
    delete plist[i];
  }

  for (auto &iter : poolmap)
  {
    mvtx_pool *pool = iter.second;
    int num_feeId = pool->get_feeidSet_size();
    if (Verbosity() > 1)
    {
      std::cout << "Number of feeid in RCDAQ events: " << num_feeId << " for packet "
                << iter.first << std::endl;
    }
    if (num_feeId > 0)
    {
      for (int i_fee{0}; i_fee < num_feeId; ++i_fee)
      {
        auto feeId = pool->get_feeid(i_fee);
        auto link = MvtxRawDefs::decode_feeid(feeId);

        auto num_strobes = pool->get_strbSet_size(feeId);
        auto num_L1Trgs = pool->get_trgSet_size(feeId);
        for (int iL1 = 0; iL1 < num_L1Trgs; ++iL1)
        {
          // auto l1Trg_bco = pool->get_L1_IR_BCO(feeId, iL1);
        }
        for (int i_strb{0}; i_strb < num_strobes; ++i_strb)
        {
          auto strb_detField = pool->get_TRG_DET_FIELD(feeId, i_strb);
          uint64_t strb_bco = pool->get_TRG_IR_BCO(feeId, i_strb);
          if (strb_detField != 0)
          {
            if (Verbosity() > 3)
            {
              std::cout << "Skipping strobe from trigger ramp for feeid: " << feeId << std::endl;
            }
            continue;
          }
          ++m_feeidMap[strb_bco];
          // for (const auto &chipId : {0, 1, 2})
          // {
          //   auto chip = MvtxNoiseMap::getChipId(link.layer, link.stave, MvtxRawDefs::gbtChipId_to_staveChipId[link.gbtid][chipId]);
          //   m_noiseMap->addStrobePerChipId(chip, 1);
          // }
          // auto strb_bc = pool->get_TRG_IR_BC(feeId, i_strb);
          auto num_hits = pool->get_TRG_NR_HITS(feeId, i_strb);

          if (Verbosity() > 4)
          {
            std::cout << "evtno: " << EventSequence << ", Fee: " << feeId;
            std::cout << " Layer: " << link.layer << " Stave: " << link.stave;
            std::cout << " GBT: " << link.gbtid << ", bco: 0x" << std::hex << strb_bco << std::dec;
            std::cout << ", n_hits: " << num_hits << std::endl;
          }
          auto hits = pool->get_hits(feeId, i_strb);
          for (const auto &hit : hits)
          {
            mvtx::_hit new_hit;
            new_hit.feeid = feeId;
            new_hit.chip_id = hit->chip_id;
            new_hit.row_pos = hit->row_pos;
            new_hit.col_pos = hit->col_pos;

            m_ChipMap[strb_bco].push_back(new_hit);
          }  //!<! hits
        }  //!<! strobes
      }  //!<! feeId
    }  //!<! if feeids != 0
  }  //!<! poolmap
  // if (m_feeidMap.size() > 1000)
  // {
  //   for (const auto &[bco, numFeeids] : m_feeidMap)
  //   {
  //     if (numFeeids == 24)
  //     {
  //       for (const auto &chip : m_ChipMap[bco])
  //       {
  //         auto link = MvtxRawDefs::decode_feeid(chip.feeid);
  //         auto layer = link.layer;
  //         auto stave = link.stave;
  //         auto chip_id = MvtxRawDefs::gbtChipId_to_staveChipId[link.gbtid][chip.chip_id];
  //         auto col = chip.col_pos;
  //         auto row = chip.row_pos;
  //       }
  //     }
  //     while (m_feeidMap.begin()->first < bco)
  //     {
  //       m_ChipMap[bco].clear();
  //       m_ChipMap.erase(m_feeidMap.begin()->first);
  //       m_feeidMap.erase(m_feeidMap.begin()->first);
  //     }
  //   }
  // }
  return Fun4AllReturnCodes::EVENT_OK;
}

//========================================================================+
int MvtxSaPrdfUnpacker::End(PHCompositeNode * /*topnode*/)
{
  // while (!m_feeidMap.empty())
  // {
  //   auto bco = m_feeidMap.begin()->first;
  //   auto numFeeids = m_feeidMap.begin()->second;
  //   if (numFeeids == 24)
  //   {
  //     for (const auto &chip : m_ChipMap[bco])
  //     {
  //       auto link = MvtxRawDefs::decode_feeid(chip.feeid);
  //       auto layer = link.layer;
  //       auto stave = link.stave;
  //       auto chip_id = MvtxRawDefs::gbtChipId_to_staveChipId[link.gbtid][chip.chip_id];
  //       auto col = chip.col_pos;
  //       auto row = chip.row_pos;
  //     }
  //   }
  //   m_ChipMap[bco].clear();
  //   m_ChipMap.erase(m_feeidMap.begin()->first);
  //   m_feeidMap.erase(m_feeidMap.begin()->first);
  // }
  PHTFileServer::cd(m_outputfile);
  PHTFileServer::write(m_outputfile);
  PHTFileServer::close();
  return Fun4AllReturnCodes::EVENT_OK;
}
