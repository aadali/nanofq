use std::fmt::{Display, Formatter};
use std::io::{BufRead, BufReader, Read};

pub type EachStats = (Box<String>, usize, f64, Option<f64>);

pub struct FqRecord {
    id: Box<String>,
    desc: Option<String>,
    seq: String,
    quality: String,
}

impl FqRecord {
    pub fn new() -> Self {
        Self {
            id: Box::new(String::new()),
            desc: None,
            seq: String::new(),
            quality: String::new(),
        }
    }

    fn is_empty(&self) -> bool {
        self.id.len() == 0 && self.desc == None && self.seq.len() == 0 && self.quality.len() == 0
    }
    fn clear(&mut self) {
        // self.id = Rc::new(String::new());
        self.id.clear();
        self.desc = None;
        self.seq.clear();
        self.quality.clear();
    }

    #[inline]
    fn id(&self) -> &str {
        &self.id
    }

    #[inline]
    fn len(&self) -> usize {
        self.seq.len()
    }

    #[inline]
    fn seq(&self) -> &str {
        &self.seq
    }

    #[inline]
    fn quality(&self) -> &str {
        &self.quality
    }

    #[inline]
    fn calculate_read_quality(&self) -> f64 {
        (self
            .quality
            .as_bytes()
            .iter()
            .map(|x| 10.0f64.powf((x - 33) as f64 / -10.0))
            .sum::<f64>()
            / self.len() as f64)
            .log10()
            * -10.0
    }

    #[inline]
    fn gc_count(&self) -> f64 {
        let gc_number = self
            .seq
            .as_bytes()
            .iter()
            .filter_map(|x| {
                if x == &b'G' || x == &b'C' {
                    Some(1)
                } else {
                    None
                }
            })
            .count();
        gc_number as f64 / self.len() as f64
    }

    pub fn stats(self, stats_gc: bool) -> EachStats {
        let len = self.len();
        let read_quality = self.calculate_read_quality();
        let gc = if stats_gc {
            Some(self.gc_count())
        } else {
            None
        };
        (self.id, len, read_quality, gc)
    }

    fn is_passed(
        &self,
        min_len: usize,
        max_len: usize,
        min_read_qvalue: f64,
        filter_gc: bool,
        min_gc: f64,
        max_gc: f64,
    ) -> bool {
        let gc_passed = if filter_gc {
            let gc = self.gc_count();
            gc > min_gc && gc < max_gc
        } else {
            true
        };
        self.len() >= min_len
            && self.len() <= max_len
            && self.calculate_read_quality() > min_read_qvalue
            && gc_passed
    }
}

pub struct FqReader<T> {
    reader: BufReader<T>,
    line_buffer: String,
    line_number: u32,
    find_last_idx: fn(char) -> bool,
}

impl<T> FqReader<T>
where
    T: Read,
{
    pub fn new(in_read: T) -> Self {
        FqReader {
            reader: std::io::BufReader::new(in_read),
            line_buffer: String::new(),
            line_number: 0u32,
            find_last_idx: |c| c != '\n' && c != '\r',
        }
    }

    pub fn get_record(&mut self, record: &mut FqRecord) -> Result<bool, String> {
        if !self.line_buffer.is_empty() {
            self.line_buffer.clear()
        }
        if record.is_empty() {
            record.clear()
        }
        let header_line_size = self
            .reader
            .read_line(&mut self.line_buffer)
            .expect(&format!(
                "Error when read header of record at line number: {}",
                self.line_number
            ));
        if header_line_size == 0 {
            return Ok(false);
        }
        if !self.line_buffer.starts_with('@') {
            return Err(format!(
                "Expect first char is @ for header line at line number: {}",
                self.line_number
            ));
        }
        let mut fields = self
            .line_buffer
            .trim_end()
            .splitn(2, |x| char::is_ascii_whitespace(&x));
        record.id.push_str(&fields.next().unwrap()[1..]);
        record.desc = fields.next().map(|x| x.to_string());
        if record.id.len() < 1 {
            return Err(format!(
                "Expect fastq record id at line number: {}",
                self.line_number
            ));
        }
        self.line_number += 1;
        self.line_buffer.clear();

        self.reader
            .read_line(&mut self.line_buffer)
            .expect(&format!(
                "Error when read seq line at line number: {}",
                self.line_number
            ));
        record.seq.push_str(self.line_buffer.trim_end());
        self.line_number += 1;

        self.reader
            .read_line(&mut self.line_buffer)
            .expect(&format!(
                "Error when read line number: {}",
                self.line_number
            ));
        self.line_buffer.clear();

        self.reader
            .read_line(&mut self.line_buffer)
            .expect(&format!(
                "Error when read quality line at line number: {}",
                self.line_number
            ));
        record.quality.push_str(self.line_buffer.trim_end());

        if record.quality.len() != record.seq.len() {
            return Err(format!(
                "Not equal length between seq and quality for {}",
                record.id
            ));
        }

        if record.seq.len() == 0 {
            return Err(format!("Length is zero for {}", record.id));
        }
        Ok(true)
    }

    pub fn stats(&mut self, stats_gc: bool) -> Vec<EachStats> {
        let mut stats_vec = vec![];
        for record_res in self {
            if let Ok(fq_record) = record_res {
                stats_vec.push(fq_record.stats(stats_gc))
            } else {
                panic!("{:?}", record_res.err().unwrap())
            }
        }
        stats_vec
    }
}

impl<T> Iterator for FqReader<T>
where
    T: Read,
{
    type Item = Result<FqRecord, String>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut record = FqRecord::new();
        match self.get_record(&mut record) {
            Ok(read_ok) => {
                if read_ok {
                    return Some(Ok(record));
                }
                None
            }
            Err(error) => Some(Err(error)),
        }
    }
}

impl Display for FqRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let header = match &self.desc {
            None => &self.id,
            Some(desc) => &format!("{} {}", &self.id, desc),
        };
        write!(f, "@{}\n{}\n+\n{}", header, self.seq, self.quality)
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn add_one() {
        assert_eq!(45, 45)
    }
}
