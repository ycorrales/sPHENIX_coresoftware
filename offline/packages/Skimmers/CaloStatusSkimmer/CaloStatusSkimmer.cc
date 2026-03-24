#include "CaloStatusSkimmer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <qautils/QAHistManagerDef.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <cassert>
#include <cstdint>
#include <iostream>

#include <TH1F.h>

//____________________________________________________________________________..
CaloStatusSkimmer::CaloStatusSkimmer(const std::string &name)
    : SubsysReco(name)
{
  std::cout << "CaloStatusSkimmer::CaloStatusSkimmer(const std::string &name) ""Calling ctor" << std::endl;
}


//____________________________________________________________________________..
int CaloStatusSkimmer::Init([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "CaloStatusSkimmer::Init(PHCompositeNode *topNode) This is Init..." << std::endl;

  if (b_produce_QA_histograms)
  {
    auto* hm = QAHistManagerDef::getHistoManager();
    assert(hm);

    h_EMC_nTowers_notinstr = new TH1F("h_EMC_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in EMCal; nNotInstrTowers; Counts", 193, -0.5, 192.5);
    h_EMC_nTowers_notinstr->SetDirectory(nullptr);
    h_HCal_nTowers_notinstr = new TH1F("h_HCal_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in HCal; nNotInstrTowers; Counts", 193, -0.5, 192.5);
    h_HCal_nTowers_notinstr->SetDirectory(nullptr);
    h_sEPD_nTowers_notinstr = new TH1F("h_sEPD_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in sEPD; nNotInstrTowers; Counts", 17, -0.5, 16.5);
    h_sEPD_nTowers_notinstr->SetDirectory(nullptr);
    h_ZDC_nTowers_notinstr = new TH1F("h_ZDC_nTowers_notinstr", "Number of not-instrumented(empty/missing pckt) towers in ZDC; nNotInstrTowers; Counts", 5, -0.5, 4.5);
    h_ZDC_nTowers_notinstr->SetDirectory(nullptr);

    h_EMC_nEvents = new TH1F("h_EMC_nEvents", "Number of events", 2, 0.5, 2.5);
    h_EMC_nEvents->GetXaxis()->SetBinLabel(1, "Total events processed");
    h_EMC_nEvents->GetXaxis()->SetBinLabel(2, "Total events skimmed");
    h_EMC_nEvents->SetDirectory(nullptr); 

    h_HCal_nEvents = new TH1F("h_HCal_nEvents", "Number of events", 2, 0.5, 2.5);
    h_HCal_nEvents->GetXaxis()->SetBinLabel(1, "Total events processed");
    h_HCal_nEvents->GetXaxis()->SetBinLabel(2, "Total events skimmed");
    h_HCal_nEvents->SetDirectory(nullptr);

    h_sEPD_nEvents = new TH1F("h_sEPD_nEvents", "Number of events", 2, 0.5, 2.5);
    h_sEPD_nEvents->GetXaxis()->SetBinLabel(1, "Total events processed");
    h_sEPD_nEvents->GetXaxis()->SetBinLabel(2, "Total events skimmed");
    h_sEPD_nEvents->SetDirectory(nullptr);

    h_ZDC_nEvents = new TH1F("h_ZDC_nEvents", "Number of events", 2, 0.5, 2.5);
    h_ZDC_nEvents->GetXaxis()->SetBinLabel(1, "Total events processed");
    h_ZDC_nEvents->GetXaxis()->SetBinLabel(2, "Total events skimmed");
    h_ZDC_nEvents->SetDirectory(nullptr);

    hm->registerHisto(h_EMC_nTowers_notinstr);
    hm->registerHisto(h_HCal_nTowers_notinstr);
    hm->registerHisto(h_sEPD_nTowers_notinstr);
    hm->registerHisto(h_ZDC_nTowers_notinstr);

    hm->registerHisto(h_EMC_nEvents);
    hm->registerHisto(h_HCal_nEvents);
    hm->registerHisto(h_sEPD_nEvents);
    hm->registerHisto(h_ZDC_nEvents);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloStatusSkimmer::process_event(PHCompositeNode *topNode)
{
  n_eventcounter++;
  if (m_EMC_skim_threshold > 0)
  {
    TowerInfoContainer *towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
    if (!towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_CEMC" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const uint32_t ntowers = towers->size();
    uint16_t notinstr_count = 0;
    for (uint32_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_count;
      }
    }
    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in EMCal = " << ntowers << ", not-instrumented(empty/missing pckt) towers in EMCal = " << notinstr_count << std::endl;
    }

    if (b_produce_QA_histograms)
    {
      h_EMC_nTowers_notinstr->Fill(notinstr_count);
    }

    if (notinstr_count >= m_EMC_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (m_HCal_skim_threshold > 0)
  {
    TowerInfoContainer *hcalin_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
    TowerInfoContainer *hcalout_towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALOUT");
    if (!hcalin_towers || !hcalout_towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_HCALIN or TOWERS_HCALOUT" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    const uint32_t ntowers_hcalin = hcalin_towers->size();
    uint16_t notinstr_count_hcalin = 0;
    for (uint32_t ch = 0; ch < ntowers_hcalin; ++ch)
    {
      TowerInfo *tower_in = hcalin_towers->get_tower_at_channel(ch);
      if (tower_in->get_isNotInstr())
      {
        ++notinstr_count_hcalin;
      }
    }

    const uint32_t ntowers_hcalout = hcalout_towers->size();
    uint16_t notinstr_count_hcalout = 0;
    for (uint32_t ch = 0; ch < ntowers_hcalout; ++ch)
    {
      TowerInfo *tower_out = hcalout_towers->get_tower_at_channel(ch);
      if (tower_out->get_isNotInstr())
      {
        ++notinstr_count_hcalout;
      }
    }

    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in HCalIn = " << ntowers_hcalin << ", not-instrumented(empty/missing pckt) towers in HCalIn = " << notinstr_count_hcalin << ", ntowers in HCalOut = " << ntowers_hcalout << ", not-instrumented(empty/missing pckt) towers in HCalOut = " << notinstr_count_hcalout << std::endl;
    }

    if (b_produce_QA_histograms)
    {
      h_HCal_nTowers_notinstr->Fill(notinstr_count_hcalin);
      h_HCal_nTowers_notinstr->Fill(notinstr_count_hcalout);
    }

    if (notinstr_count_hcalin >= m_HCal_skim_threshold ||
        notinstr_count_hcalout >= m_HCal_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (m_sEPD_skim_threshold > 0)
  {
    TowerInfoContainer *sepd_towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_SEPD");
    if (!sepd_towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_SEPD" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const uint32_t ntowers = sepd_towers->size();
    uint16_t notinstr_count = 0;
    for (uint32_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = sepd_towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_count;
      }
    }

    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in sEPD = " << ntowers << ", not-instrumented(empty/missing pckt) towers in sEPD = " << notinstr_count << std::endl;
    }

    if(b_produce_QA_histograms)
    {
      h_sEPD_nTowers_notinstr->Fill(notinstr_count);
    }

    if (notinstr_count >= m_sEPD_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (m_ZDC_skim_threshold > 0)
  {
    TowerInfoContainer *zdc_towers =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
    if (!zdc_towers)
    {
      n_notowernodecounter++;
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << "CaloStatusSkimmer::process_event: missing TOWERS_ZDC" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    const uint32_t ntowers = zdc_towers->size();
    uint16_t notinstr_count = 0;
    for (uint32_t ch = 0; ch < ntowers; ++ch)
    {
      TowerInfo *tower = zdc_towers->get_tower_at_channel(ch);
      if (tower->get_isNotInstr())
      {
        ++notinstr_count;
      }
    }

    if (Verbosity() > 9)
    {
      std::cout << "CaloStatusSkimmer::process_event: event " << n_eventcounter << ", ntowers in ZDC = " << ntowers << ", not-instrumented(empty/missing pckt) towers in ZDC = " << notinstr_count << std::endl;
    }

    if(b_produce_QA_histograms)
    {
      h_ZDC_nTowers_notinstr->Fill(notinstr_count);
    }

    if (notinstr_count >= m_ZDC_skim_threshold)
    {
      n_skimcounter++;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloStatusSkimmer::End([[maybe_unused]] PHCompositeNode *topNode)
{
  std::cout << "CaloStatusSkimmer::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  std::cout << "CaloStatusSkimmer::End Total events processed: " << n_eventcounter << std::endl;
  std::cout << "CaloStatusSkimmer::End Total events skimmed: " << n_skimcounter << std::endl;
  std::cout << "CaloStatusSkimmer::End Total events with missing tower nodes: " << n_notowernodecounter << std::endl;

  if (b_produce_QA_histograms)
  {
    h_EMC_nEvents->SetBinContent(1, n_eventcounter);
    h_EMC_nEvents->SetBinContent(2, n_skimcounter);

    h_HCal_nEvents->SetBinContent(1, n_eventcounter);
    h_HCal_nEvents->SetBinContent(2, n_skimcounter);

    h_sEPD_nEvents->SetBinContent(1, n_eventcounter);
    h_sEPD_nEvents->SetBinContent(2, n_skimcounter);

    h_ZDC_nEvents->SetBinContent(1, n_eventcounter);
    h_ZDC_nEvents->SetBinContent(2, n_skimcounter);

  }

  return Fun4AllReturnCodes::EVENT_OK;
}
