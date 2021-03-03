#if defined(__NVCC__) && !defined(UPCXX_SERIALIZED_FIELDS)
#    define UPCXX_SERIALIZED_FIELDS(...)
#endif
#include "helper.hpp"

struct ReadSeq {
  std::string read_id;
  std::string seq;
  std::string quals;
  UPCXX_SERIALIZED_FIELDS(read_id, seq, quals);
};

struct CtgWithReads {
  int64_t cid;
  std::string seq;
  double depth;
  unsigned max_reads;
  std::vector<ReadSeq> reads_left;
  std::vector<ReadSeq> reads_right;
};
