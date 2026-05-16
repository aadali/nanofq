use ahash::{HashMap, RandomState};
use bio::pattern_matching::myers::Myers;
use crate::fastq2::FastqRecord;
use crate::primer_barcode::{Primer, PO};

pub struct ReadsClassifier<'a> {
    all_reads: HashMap<usize, FastqRecord>,
    known_primers: &'a HashMap<String, Primer>,
    left_range: usize,  // 100
    right_range: usize, // 100
    max_distance: u8,
}

impl<'a> ReadsClassifier<'a> {
    fn new(
        all_reads: HashMap<usize, FastqRecord>,
        known_primers: &'a HashMap<String, Primer>,
        left_range: usize,
        right_range: usize,
        max_distance: u8,
    ) -> Self {
        ReadsClassifier {
            all_reads,
            known_primers,
            left_range,
            right_range,
            max_distance,
        }
    }

    fn get_myers(&self) -> HashMap<String, [Myers; 4]> {
        self.known_primers
            .iter()
            .map(|(primer_name, primer)| {
                (
                    primer_name.clone(),
                    [
                        primer.fwd_myers(),
                        primer.rev_rc_myers(),
                        primer.rev_myers(),
                        primer.fwd_rc_myers(),
                    ],
                )
            })
            .collect::<HashMap<_, _>>()
    }

    fn classify_reads(&mut self) -> HashMap<&str, Vec<usize>> {
        let mut primer_name2reads: HashMap<&str, Vec<usize>> =
            HashMap::with_hasher(RandomState::new());
        let primer_seq2name = Primer::primer_seq2_primer_name(self.known_primers);
        let mut primer_name2myers = self.get_myers();
        let lead_length = primer_seq2name
            .keys()
            .min_by_key(|x| x.len())
            .unwrap()
            .len();
        for (idx, read) in self.all_reads.iter_mut() {
            if (read.len() as usize) < lead_length {
                continue;
            }
            let lead_seq = &read.seq[..lead_length];
            match primer_seq2name.get(lead_seq) {
                Some((po, primer_name)) => {
                    if *po == PO::F {
                        let rev_rc_primer_pat =
                            &mut primer_name2myers.get_mut(*primer_name).unwrap()[1];
                        if read.truncate_at_rev_primer(
                            rev_rc_primer_pat,
                            self.right_range,
                            self.max_distance,
                        ) {
                            primer_name2reads
                                .entry(*primer_name)
                                .and_modify(|x| x.push(*idx))
                                .or_insert(vec![]);
                        }
                    } else {
                        let rev_rc_primet_pat =
                            &mut primer_name2myers.get_mut(*primer_name).unwrap()[3];
                        if read.truncate_at_rev_primer(
                            rev_rc_primet_pat,
                            self.right_range,
                            self.max_distance,
                        ) {
                            read.reversed();
                            primer_name2reads
                                .entry(*primer_name)
                                .and_modify(|x| x.push(*idx))
                                .or_insert(vec![]);
                        }
                    }
                }
                None => {
                    'search_each_primer: for (primer_name, primer) in self.known_primers {
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
                            if read.truncate_at_rev_primer(
                                rev_rc_primer_pat,
                                self.right_range,
                                self.max_distance,
                            ) {
                                primer_name2reads
                                    .entry(primer_name.as_ref())
                                    .and_modify(|x| x.push(*idx))
                                    .or_insert(vec![]);
                                break 'search_each_primer;
                            }
                        }

                        if read.split_off_fwd_primer(
                            rev_primer_pat,
                            self.left_range,
                            self.max_distance,
                        ) {
                            if read.truncate_at_rev_primer(
                                fwd_rc_primer_pat,
                                self.right_range,
                                self.max_distance,
                            ) {
                                read.reversed();
                                primer_name2reads
                                    .entry(primer_name.as_ref())
                                    .and_modify(|x| x.push(*idx))
                                    .or_insert(vec![]);
                                break 'search_each_primer;
                            }
                        }
                    }
                }
            }
        }
        primer_name2reads
    }
}
