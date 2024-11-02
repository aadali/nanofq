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

std::string SequenceInfo::seq_info() {
    std::stringstream info;
    std::vector<std::string> no_barcode_kits{"SQK-LSK114", "SQK-RAD114", "SQK-ULK114", "SQK-PCS114"};
    if (std::ranges::find(no_barcode_kits, m_name) == no_barcode_kits.end()) {
        info << fmt::format("#Kit: {}\n", m_name);
    } else {
        info << fmt::format("#Kit-Barcode: {}\n", m_name);
    }
    info << "#[Description]: (QuerySequence)(QueryLength,SearchLengthFromReadEnd,QueryAlignPercent,QueryAlignIdentity)\n";
    if (!m_top5end_query.empty()) {
        info << fmt::format("#[Expect sequence found in read 5end]: ({})({},{},{},{})\n",
        m_top5end_query, m_top5end_query.size(),get<0>(m_top5end), get<1>(m_top5end), get<2>(m_top5end));
        // info << fmt::format("#Search length(target5) in read 5end: {}\n", get<0>(m_top5end));
        // info << fmt::format("#Align percent of query5 for query5 align to target5: {}\n", get<1>(m_top5end));
        // info << fmt::format("#Align identity of query5: (match bases / query5 size): {}\n", get<2>(m_top5end));
    }

    if (!m_top3end_query.empty()) {
        info << fmt::format("#[Expect sequence found in read 3end]: ({})({},{},{},{})\n",
        m_top3end_query, m_top3end_query.size(), std::get<0>(m_top3end), std::get<1>(m_top3end), std::get<2>(m_top3end));
        // info << fmt::format("#Search length(target3) in read 3end: {}\n", get<0>(m_top3end));
        // info << fmt::format("#Align percent of query3 for query align to target3: {}\n", get<1>(m_top3end));
        // info << fmt::format("#Align identity of query3: (match bases / query3 size): {}\n", get<2>(m_top3end));
    }
    if (!m_bot5end_query.empty()) {
        info << fmt::format("#[Expect sequence found in read 5end if it's reversed complemented]: ({})({},{},{},{})\n",
                            m_bot5end_query, m_bot5end_query.size(),std::get<0>(m_bot5end), std::get<1>(m_bot5end), std::get<2>(m_bot5end));
        // info << fmt::format("#Search length(target5) in this reverse complemented read 5end: {}\n", get<0>(m_bot5end));
        // info << fmt::format("#Align percent of rev_com_query5 for rev_com_query5 align to target5: {}\n",
        //                     get<1>(m_bot5end));
        // info << fmt::format("#Align identity of rev_com_query5: (match bases / rev_com_query5 size): {}\n",
        //                     get<2>(m_bot5end));
    }
    if (!m_bot3end_query.empty()) {
        info << fmt::format("#[Expect sequence found in read 3end if it's reversed complemented]: ({})({},{},{},{})\n",
                            m_bot3end_query, m_bot3end_query.size(),std::get<0>(m_bot3end), std::get<1>(m_bot3end), std::get<2>(m_bot3end));
        // info << fmt::format("#Search length(target3) in this reverse complemented read 3end: {}\n", get<0>(m_bot3end));
        // info << fmt::format("#Align percent of rev_com_query3 for rev_com_query3 align to target3: {}\n",
        //                     get<1>(m_bot3end));
        // info << fmt::format("#Align identity of rev_com_query3: (match bases / rev_com_query3 size): {}\n",
        //                     get<2>(m_bot3end));
    }
    return info.str();
}
