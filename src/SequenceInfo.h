#ifndef SEQUENCEINFO_H
#define SEQUENCEINFO_H
#include <string>
#include <tuple>
#include <sstream>
#include <vector>
#include <algorithm>
#include "fmt/core.h"
#include "myUtility.h"


using trim_end = std::tuple<int, float, float>; // std::tuple<end_target_len, align_percent, align_identity>
class SequenceInfo {
public:
    /*LSK114; NBD114*/
    SequenceInfo(
        const std::string& name,
        const std::string& top5end_query,
        const trim_end& top5end,
        const std::string& top3end_query,
        const trim_end& top3end
    );


    /*RAD114; RBK114; ULK114*/
    SequenceInfo(
        const std::string& name,
        const std::string& top5end_query,
        const trim_end& top5end
    );

    /*PCS114; PCB114*/
    SequenceInfo(
        const std::string& name,

        const std::string& top5end_query,
        const trim_end& top5end,

        const std::string& top3end_query,
        const trim_end& top3end,

        const std::string& bot5end_query,
        const trim_end& bot5end,

        const std::string& bot3end_query,
        const trim_end& bot3end
    );
    /*Custom Sequence*/
    SequenceInfo(
        const std::string& forward,
        const std::string& reversed);

    SequenceInfo(const SequenceInfo&) = default;
    SequenceInfo(SequenceInfo&&) = delete;
    SequenceInfo& operator=(const SequenceInfo&) = delete;
    SequenceInfo& operator=(SequenceInfo&&) = delete;

public:
    const std::string m_name;
    // const std::string barcode;
    /*
    query is the adapter or barcode or primers.
    target is the read from sequencer.
    we align query to the front or rear N bases of target reads and ignore the middle of target/
    top is one strand of dsDNA and bottom is another reverse complement strand.

    each sequence kit has its own unique parameters, such as kit_name, adapter, primers that is fixed by Kit
    And searching length of target read, align_percent, align_identity, etc. that could be specified by user in command line

    * align_percent = align_length/query_length
    * align_identity = match_bases/align_length
    * if real align_percent > this kit align_percent and real align_identity > this kit align_identity, we think we found
      the right target, and the bases before the alignment stop idx of target for 5'end and the bases after the
      alignment start idx of target for 3'end will be trimmed;
    */
    const std::string m_top5end_query;
    trim_end m_top5end;
    // const int m_top_5end_target_len;
    // const float m_top_5end_align_percent;
    // const float m_top_5end_align_identity;

    const std::string m_top3end_query{};
    trim_end m_top3end;

    // the bottom parameter may be empty.
    // Example for reads from Rapid Kit. We ignore bottom strand for this situation
    // Example for reads from LSK Kit. In ths case, the top strand and bottom strand adapters are same
    // This situation is most suitable for PCS114 and PCB114
    const std::string m_bot5end_query{};
    trim_end m_bot5end;

    const std::string m_bot3end_query{};
    trim_end m_bot3end;

public:
    std::string seq_info();
    SequenceInfo& update_sequence_info(int top5end_len, float top5end_percent,
                                       float top5end_identity, int top3end_len, float top3end_percent,
                                       float top3end_identity, int bot5end_len, float bot5end_percent,
                                       float bot5end_identity, int bot3end_len, float bot3end_percent,
                                       float bot3end_identity);
};


#endif //SEQUENCEINFO_H
