#include "SequenceInfo.h"

// LSK114; NBD114
SequenceInfo::SequenceInfo(
    const std::string& name,
    const std::string& top_5end_query,
    const trim_end& top_5end,
    const std::string& top_3end_query,
    const trim_end& top_3end)
    : m_name(name),
      m_top_5end_query(top_5end_query),
      m_top_5end(top_5end),
      m_top_3end_query(top_3end_query),
      m_top_3end(top_3end) {}

// RAD114; RBK114; ULK114
SequenceInfo::SequenceInfo(
    const std::string& name,
    const std::string& top_5end_query,
    const trim_end& top_5end)
    : m_name(name),
      m_top_5end_query(top_5end_query),
      m_top_5end(top_5end) {}

// PCS114, PCB114
SequenceInfo::SequenceInfo(
    const std::string& name,
    const std::string& top_5end_query,
    const trim_end& top_5end,
    const std::string& top_3end_query,
    const trim_end& top_3end,
    const std::string& bottom_5end_query,
    const trim_end& bottom_5end,
    const std::string& bottom_3end_query,
    const trim_end& bottom_3end)
    : m_name(name),
      m_top_5end_query(top_5end_query),
      m_top_5end(top_5end),
      m_top_3end_query(top_3end_query),
      m_top_3end(top_3end),
      m_bottom_5end_query(bottom_5end_query),
      m_bottom_5end(bottom_5end),
      m_bottom_3end_query(bottom_3end_query),
      m_bottom_3end(bottom_3end) {}
