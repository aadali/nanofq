use ansi_term::Color;
use std::cmp::Reverse;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::process::Stdio;
use std::sync::OnceLock;

static DEGE_BASES: OnceLock<HashMap<u8, HashSet<u8>>> = OnceLock::new();
static BASES: OnceLock<HashMap<u8, u8>> = OnceLock::new();
static Q2P_TABLE: OnceLock<[f64; 128]> = OnceLock::new();

pub const DORADO_TRIM_LEADING_BASE_NUMBER: usize = 60;

pub fn quit_with_error(msg: &str) -> ! {
    eprintln!();
    eprintln!("{}", Color::Red.paint(msg));
    std::process::exit(1)
}

pub fn get_dege_bases() -> &'static HashMap<u8, HashSet<u8>> {
    DEGE_BASES.get_or_init(|| {
        HashMap::from([
            (b'R', HashSet::from([b'A', b'G'])),
            (b'Y', HashSet::from([b'C', b'T'])),
            (b'M', HashSet::from([b'C', b'A'])),
            (b'K', HashSet::from([b'G', b'T'])),
            (b'S', HashSet::from([b'C', b'G'])),
            (b'W', HashSet::from([b'A', b'T'])),
            (b'H', HashSet::from([b'A', b'T', b'C'])),
            (b'B', HashSet::from([b'G', b'T', b'C'])),
            (b'V', HashSet::from([b'G', b'A', b'C'])),
            (b'D', HashSet::from([b'G', b'A', b'T'])),
            (b'N', HashSet::from([b'G', b'A', b'T', b'C'])),
        ])
    })
}

pub fn get_bases() -> &'static HashMap<u8, u8> {
    BASES.get_or_init(|| {
        HashMap::from([
            (b'A', b'T'),
            (b'T', b'A'),
            (b'G', b'C'),
            (b'C', b'G'),
            (b'R', b'Y'),
            (b'Y', b'R'),
            (b'M', b'K'),
            (b'K', b'M'),
            (b'S', b'S'),
            (b'W', b'W'),
            (b'H', b'D'),
            (b'B', b'V'),
            (b'V', b'B'),
            (b'D', b'H'),
            (b'a', b'T'),
            (b't', b'A'),
            (b'g', b'C'),
            (b'c', b'G'),
            (b'N', b'N'),
            (b'n', b'N'),
        ])
    })
}

pub fn get_q2p_table() -> &'static [f64; 128] {
    Q2P_TABLE.get_or_init(|| {
        let mut arr = [f64::NAN; 128];
        for q in 33..127usize {
            arr[q] = 10.0f64.powf((q - 33) as f64 / -10.0)
        }
        arr
    })
}

pub fn rev_com(seq: &str) -> String {
    seq.as_bytes()
        .iter()
        .map(|x| *get_bases().get(x).unwrap() as char)
        .rev()
        .collect::<String>()
}

// ref_base from primer or reference can be degenerate base
pub static IS_MATCHED: fn(&u8, &u8) -> bool = |ref_base, read_base| {
    ref_base == read_base
        || get_dege_bases()
            .get(ref_base)
            .map_or(false, |x| x.contains(read_base))
};

pub const SEP_LINE: &str =
    "----------------------------------------------------------------------\n";

pub const ERR_PROB_TABLE: [f64; 127] = [
    f64::NAN,           // 0  00  000 NUL空字符 (Null)
    f64::NAN,           // 1  01  001 SOH标题开始 (Start of Heading)
    f64::NAN,           // 2  02  002 STX正文开始 (Start of Text)
    f64::NAN,           // 3  03  003 ETX正文结束 (End of Text)
    f64::NAN,           // 4  04  004 EOT传输结束 (End of Transmission)
    f64::NAN,           // 5  05  005 ENQ询问 (Enquiry)
    f64::NAN,           // 6  06  006 ACK确认 (Acknowledge)
    f64::NAN,           // 7  07  007 BEL响铃 (Bell)
    f64::NAN,           // 8  08  010 BS退格 (Backspace)
    f64::NAN,           // 9  09  011 HT水平制表符 (Horizontal Tab)
    f64::NAN,           // 10 0A  012 LF换行符 (Line Feed)
    f64::NAN,           // 11 0B  013 VT垂直制表符 (Vertical Tab)
    f64::NAN,           // 12 0C  014 FF换页符 (Form Feed)
    f64::NAN,           // 13 0D  015 CR回车符 (Carriage Return)
    f64::NAN,           // 14 0E  016 SO移出 (Shift Out)
    f64::NAN,           // 15 0F  017 SI移入 (Shift In)
    f64::NAN,           // 16 10  020 DLE数据链路转义 (Data Link Escape)
    f64::NAN,           // 17 11  021 DC1设备控制 1 (Device Control 1)
    f64::NAN,           // 18 12  022 DC2设备控制 2 (Device Control 2)
    f64::NAN,           // 19 13  023 DC3设备控制 3 (Device Control 3)
    f64::NAN,           // 20 14  024 DC4设备控制 4 (Device Control 4)
    f64::NAN,           // 21 15  025 NAK否定确认 (Negative Acknowledge)
    f64::NAN,           // 22 16  026 SYN同步 (Synchronous Idle)
    f64::NAN,           // 23 17  027 ETB传输块结束 (End of Transmission Block)
    f64::NAN,           // 24 18  030 CAN取消 (Cancel)
    f64::NAN,           // 25 19  031 EM介质结束 (End of Medium)
    f64::NAN,           // 26 1A  032 SUB替换 (Substitute)
    f64::NAN,           // 27 1B  033 ESC转义 (Escape)
    f64::NAN,           // 28 1C  034 FS文件分隔符 (File Separator)
    f64::NAN,           // 29 1D  035 GS组分隔符 (Group Separator)
    f64::NAN,           // 30 1E  036 RS记录分隔符 (Record Separator)
    f64::NAN,           // 31 1F  037 US单元分隔符 (Unit Separator)
    f64::NAN,           // 32 20  040 空格空格 (Space)
    1.0000000000000000, // 33   21  041 !   0
    0.7943282347242815, // 34   22  042 "   1
    0.6309573444801932, // 35   23  043 #   2
    0.5011872336272722, // 36   24  044 $   3
    0.3981071705534972, // 37   25  045 %   4
    0.3162277660168379, // 38   26  046 &   5
    0.2511886431509580, // 39   27  047 '   6
    0.1995262314968880, // 40   28  050 (   7
    0.1584893192461113, // 41   29  051 )   8
    0.1258925411794167, // 42   2A  052 *   9
    0.1000000000000000, // 43   2B  053 +   10
    0.0794328234724281, // 44   2C  054 ,   11
    0.0630957344480193, // 45   2D  055 -   12
    0.0501187233627272, // 46   2E  056 .   13
    0.0398107170553497, // 47   2F  057 /   14
    0.0316227766016838, // 48   30  060 0   15
    0.0251188643150958, // 49   31  061 1   16
    0.0199526231496888, // 50   32  062 2   17
    0.0158489319246111, // 51   33  063 3   18
    0.0125892541179417, // 52   34  064 4   19
    0.0100000000000000, // 53   35  065 5   20
    0.0079432823472428, // 54   36  066 6   21
    0.0063095734448019, // 55   37  067 7   22
    0.0050118723362727, // 56   38  070 8   23
    0.0039810717055350, // 57   39  071 9   24
    0.0031622776601684, // 58   3A  072 :   25
    0.0025118864315096, // 59   3B  073 ;   26
    0.0019952623149689, // 60   3C  074 <   27
    0.0015848931924611, // 61   3D  075 =   28
    0.0012589254117942, // 62   3E  076 >   29
    0.0010000000000000, // 63   3F  077 ?   30
    0.0007943282347243, // 64   40  100 @   31
    0.0006309573444802, // 65   41  101 A   32
    0.0005011872336273, // 66   42  102 B   33
    0.0003981071705535, // 67   43  103 C   34
    0.0003162277660168, // 68   44  104 D   35
    0.0002511886431510, // 69   45  105 E   36
    0.0001995262314969, // 70   46  106 F   37
    0.0001584893192461, // 71   47  107 G   38
    0.0001258925411794, // 72   48  110 H   39
    0.0001000000000000, // 73   49  111 I   40
    0.0000794328234724, // 74   4A  112 J   41
    0.0000630957344480, // 75   4B  113 K   42
    0.0000501187233627, // 76   4C  114 L   43
    0.0000398107170553, // 77   4D  115 M   44
    0.0000316227766017, // 78   4E  116 N   45
    0.0000251188643151, // 79   4F  117 O   46
    0.0000199526231497, // 80   50  120 P   47
    0.0000158489319246, // 81   51  121 Q   48
    0.0000125892541179, // 82   52  122 R   49
    0.0000100000000000, // 83   53  123 S   50
    /*
    In practice, Dorado reports a maximum Q-score of 50 which is S (upper case).
    https://software-docs.nanoporetech.com/dorado/latest/basecaller/qscore/#q-string
    */
    0.0000079432823472, // 84   54  124 T
    0.0000063095734448, // 85   55  125 U
    0.0000050118723363, // 86   56  126 V
    0.0000039810717055, // 87   57  127 W
    0.0000031622776602, // 88   58  130 X
    0.0000025118864315, // 89   59  131 Y
    0.0000019952623150, // 90   5A  132 Z
    0.0000015848931925, // 91   5B  133 [
    0.0000012589254118, // 92   5C  134 \
    0.0000010000000000, // 93   5D  135 ]
    0.0000007943282347, // 94   5E  136 ^
    0.0000006309573445, // 95   5F  137 _
    0.0000005011872336, // 96   60  140 `
    0.0000003981071706, // 97   61  141 a
    0.0000003162277660, // 98   62  142 b
    0.0000002511886432, // 99   63  143 c
    0.0000001995262315, // 100  64  144 d
    0.0000001584893192, // 101  65  145 e
    0.0000001258925412, // 102  66  146 f
    0.0000001000000000, // 103  67  147 g
    0.0000000794328235, // 104  68  150 h
    0.0000000630957344, // 105  69  151 i
    0.0000000501187234, // 106  6A  152 j
    0.0000000398107171, // 107  6B  153 k
    0.0000000316227766, // 108  6C  154 l
    0.0000000251188643, // 109  6D  155 m
    0.0000000199526231, // 110  6E  156 n
    0.0000000158489319, // 111  6F  157 o
    0.0000000125892541, // 112  70  160 p
    0.0000000100000000, // 113  71  161 q
    0.0000000079432823, // 114  72  162 r
    0.0000000063095734, // 115  73  163 s
    0.0000000050118723, // 116  74  164 t
    0.0000000039810717, // 117  75  165 u
    0.0000000031622777, // 118  76  166 v
    0.0000000025118864, // 119  77  167 w
    0.0000000019952623, // 120  78  170 x
    0.0000000015848932, // 121  79  171 y
    0.0000000012589254, // 122  7A  172 z
    0.0000000010000000, // 123  7B  173 {
    0.0000000007943282, // 124  7C  174 |
    0.0000000006309573, // 125  7D  175 }
    0.0000000005011872, // 126  7E  176 ~
];

pub fn calculate_read_q<T: AsRef<[u8]>>(quality: T, use_dorado_quality: bool) -> f32 {
    let quality = quality.as_ref();
    debug_assert!(quality.as_ref().len() > 0);
    if use_dorado_quality {
        if quality.len() > DORADO_TRIM_LEADING_BASE_NUMBER {
            calculate_q(&quality[60..])
        } else {
            calculate_q(quality)
        }
    } else {
        calculate_q(quality)
    }
}

fn calculate_q<T: AsRef<[u8]>>(quality: T) -> f32 {
    let quality = quality.as_ref();
    let x = (quality
        .iter()
        .map(|x| ERR_PROB_TABLE[*x as usize])
        .sum::<f64>()
        / (quality.len() as f64))
        .log10()
        * -10.0;
    x as f32
}

pub fn gc<T: AsRef<[u8]>>(sequence: T) -> f32 {
    let sequence = sequence.as_ref();
    debug_assert!(sequence.len() > 0);
    let gc = (sequence
        .iter()
        .map(|&x| {
            if x == b'G' || x == b'C' || x == b'g' || x == b'c' {
                1
            } else {
                0
            }
        })
        .sum::<u32>() as f64)
        / (sequence.len() as f64);
    gc as f32
}

pub fn collect_fastq_dir(path: &str) -> Vec<PathBuf> {
    let all_fqs = Path::new(path)
        .read_dir()
        .expect(&format!("Failed to read directory: {path}"))
        .filter_map(|x| match x {
            Ok(fs) => {
                let fs_path = fs.path();
                let fs_path_str = fs_path.to_str().unwrap();
                return if fs_path_str.ends_with(".fastq")
                    || fs_path_str.ends_with(".fq")
                    || fs_path_str.ends_with(".fastq.gz")
                    || fs_path_str.ends_with(".fq.gz")
                {
                    Some(fs_path)
                } else {
                    None
                };
            }
            Err(_) => None,
        })
        .collect::<Vec<PathBuf>>();
    all_fqs
}

pub type MatchPosition = (usize, usize, u8); // (start, exclude end, distance)

pub fn find_most_right_front(
    all_matches: Vec<MatchPosition>,
    max_dist: u8,
) -> Option<MatchPosition> {
    if all_matches.is_empty() {
        None
    } else {
        let most_right = all_matches.iter().max_by_key(|x| x.1).unwrap().1;
        all_matches
            .into_iter()
            .filter(|x| (most_right - x.1) <= (max_dist as usize))
            .min_by_key(|x| x.2) // use min distance found
    }
}

pub fn find_most_left_rear(all_matches: Vec<MatchPosition>, max_dist: u8) -> Option<MatchPosition> {
    if all_matches.is_empty() {
        None
    } else {
        let most_left = all_matches.iter().min_by_key(|x| x.0).unwrap().0;
        all_matches
            .into_iter()
            .filter(|x| (x.0 - most_left) <= (max_dist as usize))
            .min_by_key(|x| x.2) // use min distance found
    }
}

pub fn complement(n: u8) -> u8 {
    match n {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b'U' | b't' | b'u' => b'A',
        other => quit_with_error(&format!("Illegal base char found {}", other as char)),
    }
}

/// used to check minimap2 and samtools or other external tools
pub fn check_program<'a>(program: &'static str, program_path: Option<&'a str>) -> &'a str {
    let program_path = program_path.unwrap_or(program);
    let mm_path = Path::new(program_path);
    let cmd = std::process::Command::new(mm_path)
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .status()
        .expect(&format!("Check {program} failed"));
    if cmd.success() {
        return program_path;
    }
    quit_with_error(&format!("check {program} path error"))
}

pub fn run_minimap2_and_index(
    work_dir: &str,
    fastq_file: &str,
    draft: &str,
    prefix: &str,
    minimap2: Option<&str>,
    samtools: Option<&str>,
) {
    let minimap2 = check_program("minimap2", minimap2);
    let samtools = check_program("samtools", samtools);
    let work_path = std::path::Path::new(work_dir);
    if !work_path.exists() {
        std::fs::create_dir_all(work_path).expect(&format!("Failed to create dir: {work_dir}"));
    }
    let mm2_child = std::process::Command::new(minimap2)
        .current_dir(work_dir)
        .args(["-a", "-x", "map-ont", "-t", "4", draft, fastq_file])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("minimap2 failed to start");

    let view_child = std::process::Command::new(samtools)
        .current_dir(work_dir)
        .stdin(mm2_child.stdout.expect("Failed to get minimap2 output"))
        .args([
            "view",
            "-",
            "-b",
            "-S",
            "-o",
            &format!("{prefix}.raw.bam"),
            "--threads",
            "2",
        ])
        .status();
    if !view_child
        .expect(&format!(
            "Failed to complete samtools view for {prefix}.raw.bam"
        ))
        .success()
    {
        quit_with_error("Samtools view failed")
    }

    let sort_child = std::process::Command::new(samtools)
        .current_dir(work_dir)
        .args([
            "sort",
            "--threads",
            "2",
            "-o",
            &format!("{prefix}.sorted.bam"),
            &format!("{prefix}.raw.bam"),
        ])
        .status();
    if !sort_child
        .expect(&format!(
            "Failed to complete samtools sort for {prefix}.raw.bam"
        ))
        .success()
    {
        quit_with_error("Samtools sort failed")
    }

    let index_child = std::process::Command::new(samtools)
        .current_dir(work_dir)
        .args(["index", &format!("{prefix}.sorted.bam")])
        .status();
    if !index_child
        .expect(&format!("Failed to index {prefix}.sorted.bam"))
        .success()
    {
        quit_with_error("Samtools index failed")
    }
}

pub fn run_abpoa(fastq_file: &str, work_dir: &str, amplicon_name: &str, abpoa: Option<&str>) {
    let abpoa = check_program("abpoa", abpoa);
    check_and_create_dir(work_dir);
    let abpoa_child = std::process::Command::new(abpoa)
        .current_dir(work_dir)
        .arg(fastq_file)
        .output()
        .expect(&format!("Failed to run abpoa {fastq_file}"));
    if !abpoa_child.status.success() {
        quit_with_error("Run abpoa failed")
    }
    let consensus = String::from_utf8(abpoa_child.stdout).unwrap();
    let abpoa_log = String::from_utf8(abpoa_child.stderr).unwrap();
    let sequence = consensus.lines().last().unwrap();
    let consensus_output = format!("{work_dir}/{amplicon_name}.draft_consensus.fasta");
    let consensus_log = format!("{work_dir}/{amplicon_name}.log");
    
    std::fs::write(
        &consensus_output,
        format!(">{amplicon_name}_{}\n{}", sequence.len(), sequence),
    )
    .expect(&format!(
        "Failed to write draft consensus into {consensus_output}"
    ));

    std::fs::write(&consensus_log, abpoa_log)
        .expect(&format!("Failed to write abpoa log into {consensus_log}"));
}

pub fn check_and_create_dir(target_dir: &str) {
    let target_dir_path = Path::new(target_dir);
    if !target_dir_path.exists() {
        std::fs::create_dir_all(target_dir_path)
            .expect(&format!("Failed to create dir: {target_dir}"))
    } else {
        if !target_dir_path.is_dir() {
            quit_with_error(&format!("{target_dir} already exists and it's not a directory"))
        }
    }
}

#[cfg(test)]
mod utils_test {
    use super::*;
    #[test]
    #[ignore]
    fn test_dege_base() {
        assert!(IS_MATCHED(&b'V', &b'A'));
        assert!(IS_MATCHED(&b'A', &b'A'));
        assert!(!IS_MATCHED(&b'C', &b'A'));
        // assert!(IS_MATCHED(&b'G', &b'V'));
        assert!(IS_MATCHED(&b'B', &b'C'));
        assert!(IS_MATCHED(&b'B', &b'T'));
        assert!(IS_MATCHED(&b'B', &b'G'));
        assert!(IS_MATCHED(&b'W', &b'T'));
    }

    #[test]
    fn t() {
        // check_minimap2(Some("/opt/homebrew/bin/minimap2"));
        // check_program("minimap2", Some("/Users/aadali/mm2"));
        // check_program("samtools", None);
    }

    #[test]
    fn minimap2() {
        run_minimap2_and_index(
            "/Users/aadali/projects/RustProjects/nanoamp/test_data",
            "/Users/aadali/projects/RustProjects/nanoamp/test_data/py-barcode04-1600-head1000.fastq",
            "/Users/aadali/projects/RustProjects/nanoamp/test_data/py-barcode04-1600-true-consensus.fasta",
            "test001",
            Some("minimap2"),
            Some("samtools"),
        )
    }

    #[test]
    fn abpoa() {
        run_abpoa(
            "/Users/aadali/projects/RustProjects/nanofq/test_data/py-barcode04-1600-filtered.fastq",
            "/Users/aadali/projects/RustProjects/nanofq/test_data",
            "py-barcode04-1600",
            // Some("abpoa")
            // None
            Some("/Users/aadali/biotools/abPOA-v1.5.6_arm64-macos/abpoa"),
        );
    }
}
