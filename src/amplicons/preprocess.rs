use crate::fastq::{FastqRecord, read_fastq};
use crate::primer_barcode::{Barcode, PO, Primer, get_myers_from_primers};
use ahash::{HashMap, HashSet, RandomState};
use bio::alignment::pairwise::Aligner;
use bio::alphabets::dna::revcomp;
use bio::pattern_matching::myers::Myers;
use log::info;
use rayon::prelude::*;
use std::cell::RefCell;
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

    pub fn collect_fastqs(&self, thread: usize) -> HashMap<usize, FastqRecord> {
        let raw_all_reads = read_fastq(self.fastq_file, true);
        let mut all_reads =
            HashMap::with_capacity_and_hasher(raw_all_reads.len(), RandomState::new());
        let front_bar_myers = self.barcode.front_myers();
        let rear_bar_myers = self.barcode.rear_myers();
        let mut empty_reads = 0;
        let total_reads = raw_all_reads.len();
        let mut barcode_trimmed_reads = 0;

        thread_local! {
            static FRONT_BAR_MYERS: RefCell<Option<Myers>> = RefCell::default();
            static REAR_BAR_MYERS: RefCell<Option<Myers>> = RefCell::default();
        }

        rayon::ThreadPoolBuilder::new()
            .num_threads(thread)
            .start_handler(move |_| {
                FRONT_BAR_MYERS.with_borrow_mut(|x| *x = Some(front_bar_myers.clone()));
                REAR_BAR_MYERS.with_borrow_mut(|x| *x = Some(rear_bar_myers.clone()));
            })
            .build_global()
            .unwrap();

        raw_all_reads
            .into_par_iter()
            .enumerate()
            .map(|(idx, mut read)| {
                let mut read_is_empty: bool = false;
                let mut is_trimmed: bool = false;
                FRONT_BAR_MYERS.with_borrow_mut(|x| {
                    let (read_no_len, is_split_off) = read.split_off_front_barcode_end(
                        x.as_mut().unwrap(),
                        self.left_range,
                        self.max_distance,
                    );
                    read_is_empty = read_no_len;
                    is_trimmed = is_split_off;
                });
                if read_is_empty {
                    return (true, idx, None);
                }
                REAR_BAR_MYERS.with_borrow_mut(|x| {
                    let (read_no_len, is_truncated) = read.truncate_at_rear_barcode_start(
                        x.as_mut().unwrap(),
                        self.right_range,
                        self.max_distance,
                    );
                    read_is_empty = read_no_len;
                    if !is_trimmed {
                        is_trimmed = is_truncated;
                    }
                });
                if read_is_empty {
                    return (true, idx, None);
                }
                return (is_trimmed, idx, Some(read));
            })
            .collect::<Vec<_>>()
            .into_iter()
            .for_each(|(read_is_trimmed, read_idx, read_opt)| {
                if read_is_trimmed {
                    barcode_trimmed_reads += 1;
                }
                if read_opt.is_some() {
                    all_reads.insert(read_idx, read_opt.unwrap());
                } else {
                    empty_reads += 1;
                }
            });

        {
            info!("{} reads found in {}", total_reads, self.fastq_file);
            info!(
                "{} reads reserved, {} of them undergone barcode trimmed",
                all_reads.len(),
                barcode_trimmed_reads
            );
            info!("{empty_reads} dropped cause zero length after trimmed");
        }
        all_reads
    }
}
pub struct DetectedPrimer {
    fwd_primer: Vec<u8>,
    pub fwd_primer_reads: Vec<usize>,
    rev_primer: Vec<u8>,
    pub rev_primer_reads: Vec<usize>,
    _guess_reads_number: usize,
    _rev_primer_found: usize,
}

impl DetectedPrimer {
    pub fn generate_primer(&self, primer_name: &str) -> Primer {
        Primer::new(
            primer_name,
            self.fwd_primer.as_slice(),
            self.rev_primer.as_slice(),
        )
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

    pub fn save_clean_fastq(&self, output_file: &str) {
        let mut read_idx = 0usize;
        let mut count = 0;
        let mut file = BufWriter::new(
            std::fs::File::create(output_file)
                .expect(&format!("Failed to create and open file: {output_file}")),
        );
        loop {
            match self.all_reads.get(&mut read_idx) {
                None => {}
                Some(read) => {
                    count += 1;
                    read.write(&mut file).unwrap();
                }
            }
            read_idx += 1;
            if count == self.all_reads.len() {
                break;
            }
        }
        {
            info!(
                "{} barcode trimmed records, saved into {}",
                count, output_file
            );
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

    fn check_rev_rc_primer(
        &self,
        records: &[&FastqRecord],
        rev_rc_primer_myers: &mut Myers,
    ) -> usize {
        let mut rev_rc_primer_found = 0usize;
        for record in records.iter() {
            if record.find_rev_primer(rev_rc_primer_myers, self.right_range, self.max_distance) {
                rev_rc_primer_found += 1;
            }
        }
        rev_rc_primer_found
    }

    pub fn detect_one_primer(
        &self,
        guess_reads_number: usize,
        expect_amplicons: usize,
        merge_similar_lead: bool,
        min_lead_supported: usize,
    ) -> DetectedPrimer {
        let mut lead_freq = self.lead_stats(merge_similar_lead, min_lead_supported);
        lead_freq.truncate(if lead_freq.len() > expect_amplicons * 2 * 10 {
            expect_amplicons * 2 * 10
        } else {
            lead_freq.len()
        });

        {
            info!("top {} lead seqs are: ", lead_freq.iter().take(10).len());
            for (idx, (x, y)) in lead_freq.iter().take(10).enumerate() {
                info!(
                    "{}:\t{} >> {}",
                    idx + 1,
                    str::from_utf8(x).unwrap(),
                    y.len()
                )
            }
        }

        let most_fwd_primer = lead_freq[0].0.clone();
        let most_fwd_primer_reads = lead_freq[0].1.len();
        let mut lead_freq_hashmap = lead_freq
            .into_iter()
            .collect::<HashMap<Vec<u8>, Vec<usize>>>();
        let mut rev_rc_primer_myers = lead_freq_hashmap
            .iter()
            .map(|(primer, _)| (primer.as_slice(), Myers::<u64>::new(revcomp(primer))))
            .collect::<HashMap<&[u8], Myers>>();
        let mut detected_rev_primers = vec![];
        let records_used_to_search_rev_primer = lead_freq_hashmap
            .get::<Vec<u8>>(most_fwd_primer.as_ref())
            .unwrap()
            .iter()
            .take(guess_reads_number)
            .map(|x| self.all_reads.get(x).unwrap())
            .collect::<Vec<_>>();

        for (rev_primer, rev_primer_reads) in lead_freq_hashmap.iter() {
            if rev_primer == &most_fwd_primer {
                continue;
            }
            let find_rev_primer_reads = self.check_rev_rc_primer(
                &records_used_to_search_rev_primer,
                rev_rc_primer_myers.get_mut(rev_primer.as_slice()).unwrap(),
            );
            let estimate_rev_primer_matched = (find_rev_primer_reads as f64
                / (records_used_to_search_rev_primer.len() as f64))
                * (rev_primer_reads.len() as f64);

            detected_rev_primers.push((
                rev_primer,
                rev_primer_reads.len(),
                find_rev_primer_reads,
                estimate_rev_primer_matched,
            ));
        }
        detected_rev_primers.sort_by_key(|x| Reverse((x.3 * 1000.0) as usize));

        {
            info!(
                "top {} candidate rev primers for candidate primer: {} with {}",
                detected_rev_primers.iter().take(10).len(),
                str::from_utf8(most_fwd_primer.as_ref()).unwrap(),
                most_fwd_primer_reads
            );
            info!(
                "Index: LeadSeq >> LeadFreq >> SearchReadsNum >> FwdPriFoundRear >> EstFwdPriFoundRear"
            );
            for (idx, x) in detected_rev_primers.iter().take(10).enumerate() {
                info!(
                    "{}:\t{} >> {} >> {} >> {} >> {}",
                    idx + 1,
                    str::from_utf8(x.0).unwrap(),
                    x.1,
                    records_used_to_search_rev_primer.len(),
                    x.2,
                    x.3
                );
            }
        }

        let (rev_primer, _, rev_primer_found, _) = detected_rev_primers.first().unwrap();
        let rev_primer = rev_primer.to_vec();
        let rev_primer_found = *rev_primer_found;
        let fwd_primer_reads = lead_freq_hashmap.remove(&most_fwd_primer).unwrap();
        let rev_primer_reads = lead_freq_hashmap.remove(&rev_primer).unwrap();
        let detected_paired_primer = DetectedPrimer {
            fwd_primer: most_fwd_primer,
            fwd_primer_reads,
            rev_primer,
            rev_primer_reads,
            _guess_reads_number: records_used_to_search_rev_primer.len(),
            _rev_primer_found: rev_primer_found,
        };
        detected_paired_primer
    }

    pub fn lead_stats(
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
        lead_freq.retain(|_, y| y.len() > min_lead_supported);
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
            let mut aligner =
                Aligner::with_capacity(self.lead_length, self.lead_length, -2, -1, |x, y| {
                    if x == y { 1 } else { -1 }
                });
            for (outer_lead_seq, _) in &lead_freq_vec {
                if seen.contains(outer_lead_seq) {
                    continue;
                }
                let mut outer_lead_seq_myers = Myers::<u64>::new(outer_lead_seq);
                let mut similar_reads_index: HashSet<usize> =
                    HashSet::with_hasher(RandomState::new());
                for (inner_lead_seq, _) in &lead_freq_vec {
                    if outer_lead_seq == inner_lead_seq || seen.contains(inner_lead_seq) {
                        continue;
                    }
                    let mut all_matches =
                        outer_lead_seq_myers.find_all(inner_lead_seq, self.max_distance);

                    if all_matches.next().is_some() {
                        similar_reads_index
                            .extend(lead_freq.get(inner_lead_seq.as_slice()).unwrap());
                        seen.insert(inner_lead_seq);
                        lead_freq.remove(inner_lead_seq.as_slice());
                        continue;
                    }

                    let alignment = aligner.local(outer_lead_seq, inner_lead_seq);
                    if alignment.score >= ((self.lead_length - self.max_distance as usize) as i32) {
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
                seen.insert(outer_lead_seq);
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

pub struct ReadsWithPairedPrimers {
    pub reads: HashMap<usize, FastqRecord>,
    pub redundant_reads: HashMap<usize, FastqRecord>,
    // primer: Primer,
    reads_downsample: usize,
}

impl ReadsWithPairedPrimers {
    pub fn new(
        reads: HashMap<usize, FastqRecord>,
        // primer: Primer,
        reads_downsample: usize,
    ) -> ReadsWithPairedPrimers {
        ReadsWithPairedPrimers {
            reads,
            redundant_reads: HashMap::with_hasher(RandomState::new()),
            // primer,
            reads_downsample,
        }
    }

    pub fn filter(
        &mut self,
        read_q: f64,
        length_range: f64,
        failed_fastq_opt: Option<&str>,
    ) -> usize {
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

        match failed_fastq_opt {
            None => {
                failed_reads_idxes.iter().for_each(|x| {
                    self.reads.remove(x).unwrap();
                });
            }
            Some(failed_fastq) => {
                let mut failed_fastq_writer =
                    BufWriter::new(std::fs::File::create(failed_fastq).unwrap());
                failed_reads_idxes.iter().for_each(|x| {
                    let _ = self
                        .reads
                        .remove(x)
                        .unwrap()
                        .write(&mut failed_fastq_writer);
                });
            }
        }

        if self.reads.len() > self.reads_downsample {
            // Calculate mean_length again, this mean_length is more precisely.
            // Because some longer reads joined by multi individual amplicon had been removed at previous length filtered
            let mean_length =
                self.reads
                    .iter()
                    .fold(0usize, |acc, x| acc + x.1.len() as usize) as f64
                    / (self.reads.len() as f64);

            let (mut longer, mut shorter): (Vec<_>, Vec<_>) = self
                .reads
                .iter()
                .partition(|&(_, read)| read.len() as f64 > mean_length);

            longer.par_sort_by_key(|&(_, read)| {
                (read.len(), Reverse((read.qual(false) * 1000.00) as u32)) //
            });

            shorter.par_sort_by_key(|&(_, read)| {
                (Reverse(read.len()), (read.qual(false) * 1000.00) as u32)
            });

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
