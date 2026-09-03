#ifndef REGIONAL_LD_H
#define REGIONAL_LD_H

#include "bmediator.h"
#include "plink_ld.h"

namespace bmediator {

void compute_multisignal_regional_evidence(ProteinData& protein,
                                           const PlinkData& plink,
                                           const Options& opts);

} // namespace bmediator

#endif // REGIONAL_LD_H
