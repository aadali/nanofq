#include <sstream>
#include <vector>
#include <algorithm>
#include "SequenceInfo.h"

#include <fmt/core.h>

#include "myUtility.h"

// LSK114; NBD114
SequenceInfo::SequenceInfo(
    const std::string& name,
    const std::string& top5end_query,
    const trim_end& top5end,
    const std::string& top3end_query,
    const trim_end& top3end)
    : m_name(name),
      m_top5end_query(top5end_query),
      m_top5end(top5end),
      m_top3end_query(top3end_query),
      m_top3end(top3end) {}

// RAD114; RBK114; ULK114
SequenceInfo::SequenceInfo(
    const std::string& name,
    const std::string& top5end_query,
    const trim_end& top5end)
    : m_name(name),
      m_top5end_query(top5end_query),
      m_top5end(top5end) {}

// PCS114, PCB114
SequenceInfo::SequenceInfo(
    const std::string& name,
    const std::string& top5end_query,
    const trim_end& top5end,
    const std::string& top3end_query,
    const trim_end& top3end,
    const std::string& bot5end_query,
    const trim_end& bot5end,
    const std::string& bot3end_query,
    const trim_end& bot3end)
    : m_name(name),
      m_top5end_query(top5end_query),
      m_top5end(top5end),
      m_top3end_query(top3end_query),
      m_top3end(top3end),
      m_bot5end_query(bot5end_query),
      m_bot5end(bot5end),
      m_bot3end_query(bot3end_query),
      m_bot3end(bot3end) {}

// Custom Sequence
SequenceInfo::SequenceInfo(const std::string& forward,
                           const std::string& reversed)
    : m_name("Custom"),
      m_top5end_query(forward),
      m_top5end(180, 0.8, 0.8),
      m_top3end_query(reversed),
      m_top3end(180, 0.8, 0.8),
      m_bot5end_query(myUtility::rev_com(reversed)),
      m_bot5end(180, 0.8, 0.8),
      m_bot3end_query(myUtility::rev_com(forward)),
      m_bot3end(180, 0.8, 0.8) {}

std::string SequenceInfo::seq_info() {
    std::stringstream info;
    std::vector<std::string> no_barcode_kits{"SQK-LSK114", "SQK-RAD114", "SQK-ULK114", "SQK-PCS114"};
    if (std::ranges::find(no_barcode_kits, m_name) == no_barcode_kits.end()) {
        info << fmt::format("#Kit: {}\n", m_name);
    } else {
        info << fmt::format("#Kit-Barcode: {}\n", m_name);
    }
    info << "#[Description]: (QuerySeq, QuerySeqLen, SearchLenFromReadEnd, QueryAlignPercent, QueryAlignIdentity)\n";
    if (!m_top5end_query.empty()) {
        info << fmt::format("#[Expect sequence found in read 5end]: ({}, {}, {}, {}, {})\n",
                            m_top5end_query, m_top5end_query.size(), get<0>(m_top5end), get<1>(m_top5end),
                            get<2>(m_top5end));
    }

    if (!m_top3end_query.empty()) {
        info << fmt::format("#[Expect sequence found in read 3end]: ({}, {}, {}, {}, {})\n",
                            m_top3end_query, m_top3end_query.size(), std::get<0>(m_top3end), std::get<1>(m_top3end),
                            std::get<2>(m_top3end));
    }
    if (!m_bot5end_query.empty()) {
        info << fmt::format(
            "#[Expect sequence found in read 5end if it's reversed complemented]: ({}, {}, {}, {}, {})\n",
            m_bot5end_query, m_bot5end_query.size(), std::get<0>(m_bot5end), std::get<1>(m_bot5end),
            std::get<2>(m_bot5end));
    }
    if (!m_bot3end_query.empty()) {
        info << fmt::format(
            "#[Expect sequence found in read 3end if it's reversed complemented]: ({}, {}, {}, {}, {})\n",
            m_bot3end_query, m_bot3end_query.size(), std::get<0>(m_bot3end), std::get<1>(m_bot3end),
            std::get<2>(m_bot3end));
    }
    return info.str();
}

SequenceInfo& SequenceInfo::update_sequence_info(int top5end_len, float top5end_percent, float top5end_identity,
                                                 int top3end_len, float top3end_percent, float top3end_identity,
                                                 int bot5end_len, float bot5end_percent, float bot5end_identity,
                                                 int bot3end_len, float bot3end_percent, float bot3end_identity) {
    /*
* parameter: *len, *percent, *identity should be checked at the stage of parsing argument from CLI
* all of them couldn't be negative
* *percent, *identity must >0 and <1.0
* if parameter not be set,
default *len = std::numeric_limits<int>::max()
default *percent = DEFAULT_FLOAT // 3.14
default *identity = DEFAULT_FLOAT // 3.14
*/
    if (top5end_len != DEFAULT_INT) {
        get<0>(m_top5end) = top5end_len;
    }
    if (static_cast<int>(top5end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(m_top5end) = top5end_percent;
    }
    if (static_cast<int>(top5end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(m_top5end) = top5end_identity;
    }

    if (top3end_len != DEFAULT_INT) {
        get<0>(m_top3end) = top3end_len;
    }
    if (static_cast<int>(top3end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(m_top3end) = top3end_percent;
    }
    if (static_cast<int>(top3end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(m_top3end) = top3end_identity;
    }


    if (bot5end_len != DEFAULT_INT) {
        get<0>(m_bot5end) = bot5end_len;
    }
    if (static_cast<int>(bot5end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(m_bot5end) = bot5end_percent;
    }
    if (static_cast<int>(bot5end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(m_bot5end) = bot5end_identity;
    }

    if (bot3end_len != DEFAULT_INT) {
        get<0>(m_bot3end) = bot3end_len;
    }
    if (static_cast<int>(bot3end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(m_bot3end) = bot3end_percent;
    }
    if (static_cast<int>(bot3end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(m_bot3end) = bot3end_identity;
    }
    return *this;
}
