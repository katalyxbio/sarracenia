#[derive(Debug, Clone)]
struct ModificationGroup {
    header: String,       // e.g. "C+m?" or "C+m."
    canonical_base: char, // 'C'
    strand: char,         // '+'
    skips: Vec<usize>,    // skips list
}

fn parse_group(group_str: &str) -> Option<ModificationGroup> {
    if group_str.is_empty() {
        return None;
    }
    let parts: Vec<&str> = group_str.split(',').collect();
    if parts.is_empty() {
        return None;
    }
    let header = parts[0].to_string();
    if header.len() < 2 {
        return None;
    }
    let mut chars = header.chars();
    let canonical_base = chars.next()?;
    let strand = chars.next()?;

    let mut skips = Vec::new();
    for &skip_str in parts.iter().skip(1) {
        if let Ok(skip) = skip_str.trim().parse::<usize>() {
            skips.push(skip);
        }
    }
    Some(ModificationGroup {
        header,
        canonical_base,
        strand,
        skips,
    })
}

fn parse_mm_tag(mm_str: &str) -> Vec<ModificationGroup> {
    mm_str
        .split(';')
        .filter_map(|s| {
            let s = s.trim();
            if s.is_empty() { None } else { parse_group(s) }
        })
        .collect()
}

fn get_canonical_occurrences(seq: &[u8], canonical_base: char, strand: char) -> Vec<usize> {
    let target_base = if strand == '-' {
        match canonical_base {
            'A' | 'a' => b'T',
            'C' | 'c' => b'G',
            'G' | 'g' => b'C',
            'T' | 't' => b'A',
            'U' | 'u' => b'A',
            _ => canonical_base as u8,
        }
    } else {
        canonical_base as u8
    };

    let target_upper = target_base.to_ascii_uppercase();
    let mut indices = Vec::new();

    if strand == '-' {
        // Count in reverse order: from L-1 down to 0
        for i in (0..seq.len()).rev() {
            if seq[i].to_ascii_uppercase() == target_upper {
                indices.push(i);
            }
        }
    } else {
        // Count in forward order: from 0 to L-1
        for i in 0..seq.len() {
            if seq[i].to_ascii_uppercase() == target_upper {
                indices.push(i);
            }
        }
    }
    indices
}

fn complement_base(base: char) -> char {
    match base {
        'A' => 'T',
        'a' => 't',
        'C' => 'G',
        'c' => 'g',
        'G' => 'C',
        'g' => 'c',
        'T' => 'A',
        't' => 'a',
        'U' => 'A',
        'u' => 'a',
        _ => base,
    }
}

fn toggle_strand(strand: char) -> char {
    match strand {
        '+' => '-',
        '-' => '+',
        _ => strand,
    }
}

fn rebuild_header(orig_header: &str, base: char, strand: char) -> String {
    let mut chars: Vec<char> = orig_header.chars().collect();
    if chars.len() >= 2 {
        chars[0] = base;
        chars[1] = strand;
    }
    chars.into_iter().collect()
}

/// Recalculates MM and ML tags for a trimmed and optionally flipped/reverse-complemented read.
pub fn recalculate_base_mods(
    seq: &[u8],
    trimmed_seq: &[u8],
    start: usize,
    end: usize,
    flip: bool,
    mm_str: &str,
    ml_probs: &[u8],
) -> (String, Vec<u8>) {
    let original_groups = parse_mm_tag(mm_str);
    let mut ml_offset = 0;

    let mut new_mm_groups = Vec::new();
    let mut new_ml_probs = Vec::new();

    for group in original_groups {
        let n_sites = group.skips.len();
        if ml_offset + n_sites > ml_probs.len() {
            break;
        }
        let group_probs = &ml_probs[ml_offset..ml_offset + n_sites];
        ml_offset += n_sites;

        // 1. Get original sequence indices of the modified bases
        let orig_occurrences = get_canonical_occurrences(seq, group.canonical_base, group.strand);
        let mut modified_sites = Vec::new();
        let mut current_occurrence_idx = 0;
        for (idx, &skip) in group.skips.iter().enumerate() {
            current_occurrence_idx += skip;
            if current_occurrence_idx < orig_occurrences.len() {
                let seq_pos = orig_occurrences[current_occurrence_idx];
                let prob = group_probs[idx];
                modified_sites.push((seq_pos, prob));
            }
            current_occurrence_idx += 1; // move past modified site
        }

        // 2. Filter and transform positions for trimming and flipping
        let mut kept_sites = Vec::new();
        for (seq_pos, prob) in modified_sites {
            if seq_pos >= start && seq_pos < end {
                let mut new_pos = seq_pos - start;
                let mut new_base = group.canonical_base;
                let mut new_strand = group.strand;

                if flip {
                    let trimmed_len = end - start;
                    new_pos = trimmed_len - 1 - new_pos;
                    new_base = complement_base(new_base);
                    new_strand = toggle_strand(new_strand);
                }

                kept_sites.push((new_pos, prob, new_base, new_strand));
            }
        }

        if kept_sites.is_empty() {
            continue;
        }

        // Since all sites in a group have the same base and strand:
        let group_base = kept_sites[0].2;
        let group_strand = kept_sites[0].3;

        // 3. Sort according to new counting order
        if group_strand == '-' {
            kept_sites.sort_by(|a, b| b.0.cmp(&a.0)); // Descending order
        } else {
            kept_sites.sort_by(|a, b| a.0.cmp(&b.0)); // Ascending order
        }

        // 4. Calculate new skips in the trimmed sequence
        let new_occurrences = get_canonical_occurrences(trimmed_seq, group_base, group_strand);
        let mut new_skips = Vec::new();
        let mut prev_occurrence_idx = None;

        for (new_pos, prob, _, _) in &kept_sites {
            if let Some(occ_idx) = new_occurrences.iter().position(|&p| p == *new_pos) {
                let skip = match prev_occurrence_idx {
                    None => occ_idx,
                    Some(prev) => occ_idx - prev - 1,
                };
                new_skips.push(skip);
                new_ml_probs.push(*prob);
                prev_occurrence_idx = Some(occ_idx);
            }
        }

        if !new_skips.is_empty() {
            let new_header = rebuild_header(&group.header, group_base, group_strand);
            let skips_str = new_skips
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>()
                .join(",");
            new_mm_groups.push(format!("{},{}", new_header, skips_str));
        }
    }

    let new_mm_str = if new_mm_groups.is_empty() {
        "".to_string()
    } else {
        new_mm_groups.join(";") + ";"
    };

    (new_mm_str, new_ml_probs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_recalculate_base_mods_forward() {
        let seq = b"ACGCGCG";
        let trimmed_seq = b"GCGC";
        let start = 2;
        let end = 6;
        let flip = false;
        let mm_str = "C+m?,0,1;";
        let ml_probs = vec![255, 128];

        let (new_mm, new_ml) = recalculate_base_mods(seq, trimmed_seq, start, end, flip, mm_str, &ml_probs);
        assert_eq!(new_mm, "C+m?,1;");
        assert_eq!(new_ml, vec![128]);
    }

    #[test]
    fn test_recalculate_base_mods_reverse() {
        let seq = b"CGCGCGT";
        let trimmed_seq = b"CGCG";
        let start = 2;
        let end = 6;
        let flip = false;
        let mm_str = "C-m,0,1;";
        let ml_probs = vec![255, 128];

        let (new_mm, new_ml) = recalculate_base_mods(seq, trimmed_seq, start, end, flip, mm_str, &ml_probs);
        assert_eq!(new_mm, "C-m,0;");
        assert_eq!(new_ml, vec![255]);
    }
}
