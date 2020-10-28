#include "main.hpp"

bool _verbose = false;

StageTimers stage_timers = {
    .merge_reads = new IntermittentTimer(__FILENAME__ + string(":") + "Merge reads", "Merging reads"),
    .cache_reads = new IntermittentTimer(__FILENAME__ + string(":") + "Load reads into cache", "Loading reads into cache"),
    .load_ctgs = new IntermittentTimer(__FILENAME__ + string(":") + "Load contigs", "Loading contigs"),
    .analyze_kmers = new IntermittentTimer(__FILENAME__ + string(":") + "Analyze kmers", "Analyzing kmers"),
    .dbjg_traversal = new IntermittentTimer(__FILENAME__ + string(":") + "Traverse deBruijn graph", "Traversing deBruijn graph"),
    .alignments = new IntermittentTimer(__FILENAME__ + string(":") + "Alignments", "Aligning reads to contigs"),
    .kernel_alns = new IntermittentTimer(__FILENAME__ + string(":") + "Kernel alignments", ""),
    .localassm = new IntermittentTimer(__FILENAME__ + string(":") + "Local assembly", "Locally extending ends of contigs"),
    .cgraph = new IntermittentTimer(__FILENAME__ + string(":") + "Traverse contig graph", "Traversing contig graph"),
    .dump_ctgs = new IntermittentTimer(__FILENAME__ + string(":") + "Dump contigs", "Dumping contigs"),
    .compute_kmer_depths = new IntermittentTimer(__FILENAME__ + string(":") + "Compute kmer depths", "Computing kmer depths")};

#if defined(ENABLE_GASNET_STATS)

// We may be compiling with debug-mode GASNet with optimization.
// GASNet has checks to prevent users from blindly doing this,
// because it's a bad idea to run that way in production.
// However in this case we know what we are doing...
#undef NDEBUG
#undef __OPTIMIZE__
#include <gasnetex.h>
#include <gasnet_tools.h>
string _gasnet_stats_stage = "";
void mhm2_trace_set_mask(const char *newmask) { GASNETT_TRACE_SETMASK(newmask); }

#endif
