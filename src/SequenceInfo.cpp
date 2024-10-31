#include "SequenceInfo.h"

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
