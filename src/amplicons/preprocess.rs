use crate::fastq2::{FastqRecord, read_fastq};
use crate::primer_barcode::{Barcode, PO, Primer, get_myers_from_primers};
use ahash::{HashMap, HashSet, RandomState};
use bio::alignment::distance::levenshtein;
use bio::alphabets::dna::revcomp;
use bio::pattern_matching::myers::Myers;
use log::info;
use std::cmp::Reverse;
use std::io::BufWriter;

pub struct ReadsCollector<'a> {
    fastq_file: &'a str,
    barcode: &'a Barcode,
    left_range: usize,  // 150
    right_range: usize, // 180
    max_distance: u8,
}

impl<'a> ReadsCollector<'a> {
    pub fn new(
        fastq_file: &'a str,
        barcode: &'a Barcode,
        left_range: usize,
        right_range: usize,
        max_distance: u8,
    ) -> Self {
        ReadsCollector {
            fastq_file,
            barcode,
            left_range,
            right_range,
            max_distance,
        }
    }

    pub fn collect_fastqs(&self) -> HashMap<usize, FastqRecord> {
        let raw_all_reads = read_fastq(self.fastq_file, true);
        let mut all_reads =
            HashMap::with_capacity_and_hasher(raw_all_reads.len(), RandomState::new());
        let mut front_bar_myers = self.barcode.front_myers();
        let mut rear_bar_myers = self.barcode.rear_myers();
        let mut length0_number = 0;
        let mut total_reads = 0usize;
        for (idx, mut read) in raw_all_reads.into_iter().enumerate() {
            total_reads += 1;
            if read.split_off_front_barcode_end(
                &mut front_bar_myers,
                self.left_range,
                self.max_distance,
            ) {
                if read.truncate_at_rear_barcode_start(
                    &mut rear_bar_myers,
                    self.right_range,
                    self.max_distance,
                ) {
                    all_reads.insert(idx, read);
                } else {
                    length0_number += 1;
                }
            } else {
                length0_number += 1;
            }
        }
        info!("{} reads found in {}", total_reads, self.fastq_file);
        info!("{} barcoded trimmed records", all_reads.len());
        info!("{length0_number} dropped cause zero length after trimmed");
        all_reads
    }
}

pub struct ReadsClassifier {
    pub all_reads: HashMap<usize, FastqRecord>,
    // known_primers: &'a HashMap<String, Primer>,
    lead_length: usize,
    left_range: usize,  // 100
    right_range: usize, // 100
    max_distance: u8,
    pub analysis_name: String,
}

impl ReadsClassifier {
    pub fn new(
        all_reads: HashMap<usize, FastqRecord>,
        lead_length: usize,
        left_range: usize,
        right_range: usize,
        max_distance: u8,
        analysis_name: String,
    ) -> Self {
        ReadsClassifier {
            all_reads,
            lead_length,
            left_range,
            right_range,
            max_distance,
            analysis_name,
        }
    }

    pub fn save_fastq(&self, output_file: &str) {
        let mut file = BufWriter::new(
            std::fs::File::create(output_file)
                .expect(&format!("Failed to create file: {output_file}")),
        );
        for (_, read) in &self.all_reads {
            read.write(&mut file).unwrap();
        }
    }

    pub fn detect_one_primer(
        &self,
        guess_reads_number: usize,
        min_rev_match_reads_number: usize,
        rev_primer_found_ratio: f64,
        primer_idx: usize,
        merge_similar_lead: bool,
        min_lead_supported: usize,
    ) -> Option<Primer> {
        let lead_freq = self.lead_stats(merge_similar_lead, min_lead_supported);
        info!(
            "First lead seq is: {} with {} reads",
            str::from_utf8(&lead_freq[0].0.as_slice()).unwrap(),
            &lead_freq[0].1.len()
        );
        if lead_freq.len() > 1 {
            info!(
                "Second lead seq is: {} with {} reads",
                str::from_utf8(&lead_freq[1].0.as_slice()).unwrap(),
                &lead_freq[1].1.len()
            );
        }
        for (fwd_primer_seq, reads_idx) in lead_freq.iter() {
            for (rev_primer_seq, _) in lead_freq.iter() {
                if fwd_primer_seq == rev_primer_seq {
                    continue;
                }
                let mut rev_rc_primer = Myers::<u64>::new(revcomp(rev_primer_seq));
                let mut match_reads_number = 0;
                for (idx, read_idx) in reads_idx.iter().enumerate() {
                    /* TODO optimize search rev primer depending on the read length. Start search from reads with mean lengths */
                    let read = self.all_reads.get(read_idx).unwrap();
                    if read.find_rev_primer(&mut rev_rc_primer, self.right_range, self.max_distance)
                    {
                        match_reads_number += 1;
                        if idx == guess_reads_number - 1 {
                            if (match_reads_number as f64) / (guess_reads_number as f64)
                                > rev_primer_found_ratio
                            {
                                return Some(Primer::new(
                                    &format!("{}_primer{primer_idx}", self.analysis_name),
                                    fwd_primer_seq,
                                    rev_primer_seq,
                                ));
                            }
                        }
                    }
                }
                if match_reads_number < min_rev_match_reads_number {
                    return None;
                } else {
                    if (match_reads_number as f64) / (guess_reads_number as f64)
                        > rev_primer_found_ratio
                    {
                        return Some(Primer::new(
                            &format!("{}_primer{primer_idx}", self.analysis_name),
                            fwd_primer_seq,
                            rev_primer_seq,
                        ));
                    }
                }
            }
        }
        None
    }

    fn lead_stats(
        &self,
        merge_similar: bool,
        min_lead_supported: usize,
    ) -> Vec<(Vec<u8>, Vec<usize>)> {
        let mut lead_freq = HashMap::with_hasher(RandomState::new());
        for (idx, read) in &self.all_reads {
            if (read.len() as usize) < self.lead_length {
                continue;
            }
            let lead_seq = &read.seq[..self.lead_length];
            match lead_freq.get_mut(lead_seq) {
                None => {
                    lead_freq.insert(lead_seq.to_vec(), vec![*idx]);
                }
                Some(reads_idx) => reads_idx.push(*idx),
            }
        }
        if !merge_similar {
            let mut lead_freq: Vec<_> = lead_freq.into_iter().collect();
            lead_freq.sort_by_key(|(_, reads_idx)| Reverse(reads_idx.len()));
            lead_freq
        } else {
            let mut lead_freq_vec = lead_freq
                .iter()
                .map(|(lead_seq, reads_idx)| (lead_seq.to_vec(), reads_idx.len()))
                .collect::<Vec<_>>();
            lead_freq_vec.sort_by_key(|x| Reverse(x.1));
            let mut seen = HashSet::with_hasher(RandomState::new());
            for (outer_lead_seq, _) in &lead_freq_vec {
                if seen.contains(outer_lead_seq) {
                    continue;
                }
                let mut similar_reads_index: HashSet<usize> =
                    HashSet::with_hasher(RandomState::new());
                for (inner_lead_seq, _) in &lead_freq_vec {
                    if outer_lead_seq == inner_lead_seq || seen.contains(inner_lead_seq) {
                        continue;
                    }
                    let dis = levenshtein(outer_lead_seq, inner_lead_seq);
                    if dis < (self.max_distance as u32) {
                        similar_reads_index
                            .extend(lead_freq.get(inner_lead_seq.as_slice()).unwrap());
                        seen.insert(inner_lead_seq);
                        lead_freq.remove(inner_lead_seq.as_slice());
                    }
                }
                lead_freq
                    .get_mut(outer_lead_seq)
                    .unwrap()
                    .extend(similar_reads_index);
            }
            lead_freq.retain(|_, reads_idx| reads_idx.len() > min_lead_supported);
            let mut lead_freq: Vec<_> = lead_freq.into_iter().collect();
            lead_freq.sort_by_key(|(_, reads_idx)| Reverse(reads_idx.len()));
            lead_freq
        }
    }

    pub fn classify_reads_with_known_primers(
        &mut self,
        known_primers: &HashMap<String, Primer>,
    ) -> HashMap<String, Vec<usize>> {
        let mut primer_name2reads: HashMap<String, Vec<usize>> =
            HashMap::with_hasher(RandomState::new());
        let primer_seq2name = Primer::primer_seq2_primer_name(known_primers);
        let mut primer_name2myers = get_myers_from_primers(known_primers);
        let lead_length = primer_seq2name
            .keys()
            .min_by_key(|x| x.len())
            .unwrap()
            .len();
        for (idx, read) in self.all_reads.iter_mut() {
            if (read.len() as usize) < lead_length {
                continue;
            }
            let read_lead_seq = &read.seq[..lead_length];
            match primer_seq2name.get(read_lead_seq) {
                Some((po, primer_name)) => {
                    if *po == PO::F {
                        let rev_rc_primer_pat =
                            &mut primer_name2myers.get_mut(*primer_name).unwrap()[1];
                        if read.truncate_at_rev_primer_start(
                            rev_rc_primer_pat,
                            self.right_range,
                            self.max_distance,
                        ) {
                            match primer_name2reads.get_mut(*primer_name) {
                                None => {
                                    primer_name2reads.insert(primer_name.to_string(), vec![*idx]);
                                }
                                Some(good_reads_idx) => good_reads_idx.push(*idx),
                            }
                            // primer_name2reads
                            //     .entry(*primer_name)
                            //     .and_modify(|x| x.push(*idx))
                            //     .or_insert(vec![]);
                        }
                    } else {
                        let rev_rc_primer_pat =
                            &mut primer_name2myers.get_mut(*primer_name).unwrap()[3];
                        if read.truncate_at_rev_primer_start(
                            rev_rc_primer_pat,
                            self.right_range,
                            self.max_distance,
                        ) {
                            read.reversed();
                            match primer_name2reads.get_mut(*primer_name) {
                                None => {
                                    primer_name2reads.insert(primer_name.to_string(), vec![*idx]);
                                }
                                Some(good_reads_idx) => {
                                    good_reads_idx.push(*idx);
                                }
                            }
                            // primer_name2reads
                            //     .entry(*primer_name)
                            //     .and_modify(|x| x.push(*idx))
                            //     .or_insert(vec![]);
                        }
                    }
                }
                None => {
                    'search_each_primer: for (primer_name, _) in known_primers {
                        let [
                            fwd_primer_pat,
                            rev_rc_primer_pat,
                            rev_primer_pat,
                            fwd_rc_primer_pat,
                        ] = primer_name2myers.get_mut(primer_name).unwrap();
                        if read.split_off_fwd_primer(
                            fwd_primer_pat,
                            self.left_range,
                            self.max_distance,
                        ) {
                            if read.truncate_at_rev_primer_start(
                                rev_rc_primer_pat,
                                self.right_range,
                                self.max_distance,
                            ) {
                                match primer_name2reads.get_mut(primer_name) {
                                    None => {
                                        primer_name2reads
                                            .insert(primer_name.to_string(), vec![*idx]);
                                    }
                                    Some(good_reads_idx) => good_reads_idx.push(*idx),
                                }
                                // primer_name2reads
                                //     .entry(primer_name.as_ref())
                                //     .and_modify(|x| x.push(*idx))
                                //     .or_insert(vec![]);
                                break 'search_each_primer;
                            }
                        }

                        if read.split_off_fwd_primer(
                            rev_primer_pat,
                            self.left_range,
                            self.max_distance,
                        ) {
                            if read.truncate_at_rev_primer_start(
                                fwd_rc_primer_pat,
                                self.right_range,
                                self.max_distance,
                            ) {
                                read.reversed();
                                match primer_name2reads.get_mut(primer_name) {
                                    None => {
                                        primer_name2reads
                                            .insert(primer_name.to_string(), vec![*idx]);
                                    }
                                    Some(good_reads_idx) => good_reads_idx.push(*idx),
                                }
                                //
                                // read.reversed();
                                // primer_name2reads
                                //     .entry(primer_name.as_ref())
                                //     .and_modify(|x| x.push(*idx))
                                //     .or_insert(vec![]);
                                break 'search_each_primer;
                            }
                        }
                    }
                }
            }
        }
        primer_name2reads
    }

    pub fn remove_reads_with_idxes(&mut self, reads_idx: &[usize]) -> HashMap<usize, FastqRecord> {
        let mut records = HashMap::with_capacity_and_hasher(reads_idx.len(), RandomState::new());

        for idx in reads_idx {
            records.insert(
                *idx,
                self.all_reads
                    .remove(idx)
                    .expect(&format!("Could found read with idx: {idx}")),
            );
        }
        records
    }

    pub fn remove_reads_with_names(&mut self, read_names: &HashSet<String>) {
        self.all_reads
            .retain(|_, read| !read_names.contains(&read.name))
    }
}

pub struct ReadsWithPrimer {
    pub reads: HashMap<usize, FastqRecord>,
    pub redundant_reads: HashMap<usize, FastqRecord>,
    // primer: Primer,
    reads_downsample: usize,
}

impl ReadsWithPrimer {
    pub fn new(
        reads: HashMap<usize, FastqRecord>,
        // primer: Primer,
        reads_downsample: usize,
    ) -> ReadsWithPrimer {
        ReadsWithPrimer {
            reads,
            redundant_reads: HashMap::with_hasher(RandomState::new()),
            // primer,
            reads_downsample,
        }
    }

    pub fn filter(&mut self, read_q: f64, length_range: f64) -> usize {
        let mean_length = self
            .reads
            .iter()
            .fold(0usize, |acc, x| acc + x.1.len() as usize) as f64
            / (self.reads.len() as f64);

        let mut failed_reads_idxes = vec![];
        for (read_idx, read) in self.reads.iter() {
            if !(read.qual(false) > (read_q as f32)
                && (read.len() as f64) < (mean_length * (1.0 + length_range))
                && (read.len() as f64) > (mean_length * (1.0 - length_range)))
            {
                failed_reads_idxes.push(*read_idx);
            }
        }
        failed_reads_idxes.iter().for_each(|x| {
            self.reads.remove(x).unwrap();
        });

        if self.reads.len() > self.reads_downsample {
            // Calculate mean_length again, this mean_length is more precisely.
            // Because some longer reads joined by multi individual amplicon has been removed at previous length filtered
            let mean_length =
                self.reads
                    .iter()
                    .fold(0usize, |acc, x| acc + x.1.len() as usize) as f64
                    / (self.reads.len() as f64);

            let (mut longer, mut shorter): (Vec<_>, Vec<_>) = self
                .reads
                .iter()
                .partition(|&(_, read)| read.len() as f64 > mean_length);

            longer.sort_by_key(|&(_, read)| read.len());
            shorter.sort_by_key(|&(_, read)| Reverse(read.len()));

            let mut idx = 0usize;
            let mut candidate_reads_idxes = HashSet::with_hasher(RandomState::new());
            loop {
                match shorter.get(idx) {
                    None => {}
                    Some(&(read_idx, _)) => {
                        candidate_reads_idxes.insert(*read_idx);
                        if candidate_reads_idxes.len() == self.reads_downsample {
                            break;
                        }
                    }
                }
                match longer.get(idx) {
                    None => {}
                    Some(&(read_idx, _)) => {
                        candidate_reads_idxes.insert(*read_idx);
                        if candidate_reads_idxes.len() == self.reads_downsample {
                            break;
                        }
                    }
                }
                idx += 1;
            }

            let mut tmp_reads = candidate_reads_idxes
                .iter()
                .map(|x| (*x, self.reads.remove(x).unwrap()))
                .collect::<HashMap<usize, FastqRecord>>();
            std::mem::swap(&mut tmp_reads, &mut self.reads);
            self.redundant_reads = tmp_reads;
        }
        failed_reads_idxes.len()
    }

    pub fn save_fastq(&self, output_file: &str) {
        let mut outf = BufWriter::new(std::fs::File::create(output_file).unwrap());
        for (_, read) in &self.reads {
            read.write(&mut outf).unwrap()
        }
    }

    pub fn save_redundant_fastq(&self, output_file: &str) {
        let mut outf = BufWriter::new(std::fs::File::create(output_file).unwrap());
        for (_, read) in &self.redundant_reads {
            read.write(&mut outf).unwrap()
        }
    }
}
