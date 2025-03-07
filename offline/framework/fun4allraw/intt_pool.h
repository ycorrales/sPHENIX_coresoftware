#ifndef FUN4ALLRAW_INTT_POOL_H
#define FUN4ALLRAW_INTT_POOL_H

#include <Event/packet.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>


class  intt_pool  {


public:
  intt_pool( const unsigned int required_depth=30000, const unsigned int low_mark =1000);
  // delete copy ctor and assignment operator (cppcheck)
  explicit intt_pool(const intt_pool &) = delete;
  intt_pool &operator=(const intt_pool &) = delete;

  virtual ~intt_pool();

  virtual int addPacket( Packet *p);

  virtual void drain() { _low_mark = 0;};

  virtual unsigned int  rawValue(const int fee, const int index);


  virtual int iValue(const int hit, const int field);
  virtual int iValue(const int hit, const char *what);
  virtual int iValue(const int, const int, const char *what);
  virtual long long lValue(const int hit, const int field);
  virtual long long lValue(const int hit, const char *what);
  virtual long long lValue(const int, const int, const char * what);


  virtual unsigned int min_depth() const; // the lowest vector length
  virtual bool depth_ok() const;
  virtual int next();

  virtual void  dump ( OSTREAM& os = std::cout);
  virtual int getIdentifier() const {return _myPacketid;};

  

  void Verbosity(const int i) {verbosity = i;}
  void Name(const std::string &n) {name = n;}
  const std::string &Name() const {return name;}


protected:
  int intt_decode ();

  int intt_decode_hitlist (std::vector<unsigned int> & /*hitlist*/ , const int /*fee*/);

  unsigned long long calcBCO(unsigned int *hitlist) const;

  static const int MAX_FEECOUNT {14};
  
  int verbosity {0};
  int _is_decoded {0};
  
  unsigned int _required_depth;
  unsigned int _low_mark;
  int _myPacketid {-1};  // we are not locked in yet

  int currentpos {0};
  int writeindex {0};

  struct intt_hit
  {
    uint64_t bco;
    uint16_t fee;
    uint16_t channel_id;
    uint16_t chip_id;
    uint16_t adc;
    uint16_t FPHX_BCO;
    uint16_t full_FPHX;
    uint16_t full_ROC;
    uint16_t amplitude;
    uint16_t full_fphx;
    uint32_t event_counter;
    uint32_t word;
  };

    
  //std::vector<unsigned int> packetData;
  unsigned int * packetData;
  int _allocated_size{0};

  std::vector<unsigned int> fee_data[MAX_FEECOUNT];
  std::vector<intt_hit *> intt_hits;

  std::array<unsigned int,MAX_FEECOUNT> last_index{};
  std::map<unsigned int, uint64_t> last_bco;
  std::string name;

  std::set<unsigned int> FEE_List;
  std::set<unsigned long long> BCO_List;
  // a list of FEEs with data belonging to a given BCO
  std::map<unsigned long long, std::set<unsigned int>> FEEs_by_BCO;
  // a list of BCOs with data belonging to a given FEE
  std::map<unsigned int, std::set<unsigned long long>> BCOs_by_FEE;

};

 
#endif 
