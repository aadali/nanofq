use ahash::{HashMap, RandomState};
use crate::fastq2::{read_fastq, FastqRecord};
use crate::primer_barcode::Barcode;

pub struct ReadsCollector<'a> {
    fastq_file: &'a str,
    barcode: &'a Barcode,
    left_range: usize,  // 150
    right_range: usize, // 180
    max_distance: u8,
}

impl<'a> ReadsCollector<'a> {
    fn new(
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

    fn collect_fastqs(&self) -> HashMap<usize, FastqRecord> {
        let raw_all_reads = read_fastq(self.fastq_file, true);
        let mut all_reads =
            HashMap::with_capacity_and_hasher(raw_all_reads.len(), RandomState::new());
        let mut front_bar_myers = self.barcode.front_myers();
        let mut rear_bar_myers = self.barcode.rear_myers();
        for (idx, mut read) in raw_all_reads.into_iter().enumerate() {
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
                }
            }
        }
        all_reads
    }
}
