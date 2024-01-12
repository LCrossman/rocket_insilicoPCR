use std::collections::HashMap;
use regex::Regex;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::io::{BufRead, BufReader, Read, Cursor};
use std::io;
use std::str;
use chacha20poly1305::ChaCha20Poly1305; 
use chacha20poly1305::aead::NewAead;
use chacha20poly1305::aead::stream::DecryptorBE32;
use sha256::try_digest;
use anyhow::anyhow;


//author="Lisa Crossman",version,about="InSilicoPCR from DNA sequences, adapted from perl script by Egan Ozer based on a php script by Joseba Bikandi"

#[derive(Debug, Clone)]
struct InSilicoPCRProduct {
    seqid: String,
    pcrid: String,
    seqpos: usize,
    length: usize,
    sequence: String,
}

struct ChunksReader {
    chunks: Vec<Vec<u8>>,
    current_chunk: usize,
}

impl Read for ChunksReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
       if self.current_chunk < self.chunks.len() {
          let remaining = &self.chunks[self.current_chunk];
	  let to_copy = std::cmp::min(buf.len(), remaining.len());
	  buf[..to_copy].copy_from_slice(&remaining[..to_copy]);
	  self.current_chunk +=1;
	  Ok(to_copy)
       } else {
          Ok(0)
       }
    }
}

pub fn in_silico_pcr(filename: &std::path::PathBuf, opti_m: bool, opti_c: bool, opti_l: usize, pfile: &std::path::PathBuf, k2: &[u8; 32], n2: &[u8; 7], k3: &[u8; 32], n3: &[u8; 7], wf: &std::path::PathBuf, ) -> Result<String, anyhow::Error> {
      let opt_s = filename;
    //let opt_a = String::new();
    //let opt_b = String::new();
    //if cli.a.is_some() {
    //   let opt_a = &cli.a.unwrap().to_string();
    //   }
    //if cli.b.is_some() {
    //   let opt_b = &cli.b.unwrap().to_string();
    //   }
      let opt_p = pfile;
      let opt_l = opti_l; //.parse::<usize>().unwrap();
      let opt_m = opti_m;
      let opt_c = opti_c;
      println!("opt_c {} opt_m {}",  opt_c, opt_m);
      let mut w = File::create(wf).unwrap();
      let primey = decrypt_large_file(opt_p, &k2, &n2);
      let prim_flat: Vec<u8> = primey.unwrap().as_slice().into_iter().flatten().map(|&x| x).collect();
      let prim_cursor = Cursor::new(&prim_flat[..]);
      let mut prims: Vec<(String,String,String)> = Vec::new();
      read_primers(BufReader::new(prim_cursor), &mut prims);
      let mut sequences: HashMap<String, String> = HashMap::new();
      let seqchunks_reader =  decrypt_large_file(filename, &k3, &n3);
      let seqchunks_flattened: Vec<u8> = seqchunks_reader.unwrap().as_slice().into_iter().flatten().map(|&y| y).collect();
      let seq_cursor = Cursor::new(&seqchunks_flattened[..]);
      let buf_reader = BufReader::new(seq_cursor);
      read_sequences(buf_reader, &mut sequences);
      let iupac: [(&str, &str); 12] = [
        ("R", "[AG]"),
        ("Y", "[CT]"),
        ("S", "[GC]"),
        ("W", "[AT]"),
        ("K", "[GT]"),
        ("M", "[AC]"),
        ("B", "[CGT]"),
        ("D", "[AGT]"),
        ("H", "[ACT]"),
        ("V", "[ACG]"),
        ("N", "."),
        ("U", "T"),
      ];
      writeln!(&mut w, "AmpId\tSequenceId\tPositionInSequence\tLength").expect("unable to write to file");
      for slice in &prims {
        let (primer1, primer2, pref) = slice.clone();
        let mut pattern1 = primer1.to_uppercase().replace(|c: char| !c.is_ascii_alphabetic(), "");
         let mut pattern2 = primer2.to_uppercase().replace(|c: char| !c.is_ascii_alphabetic(), "");

        if opt_m {
            pattern1 = include_n(&primer1);
            pattern2 = include_n(&primer2);
        }

        let mut _indel1 = String::new();
        let mut _indel2 = String::new();
//        if opt_i {
  //          indel1 = indel(&primer1);
    //        indel2 = indel(&primer2);
      //  }

        for (from, to) in &iupac {
            pattern1 = pattern1.replace(from, to);
            pattern2 = pattern2.replace(from, to);
  //          if opt_i {
    //            indel1 = indel1.replace(from, to);
      //          indel2 = indel2.replace(from, to);
        //    }
        }

        let start_pattern = format!("{}", pattern1);
	let endpattern = format!("{}", pattern2);
	let end_pattern = endpattern.chars().rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'T' => 'A',
                'G' => 'C',
                '[' => ']',
                ']' => '[',
                _ => c,
            })
            .collect::<String>();

        let full_start = start_pattern.clone();
        let full_end = end_pattern.clone();
	//Implement opt_i later if requested
    //    if opt_i {
      //      let start_indel = format!("{}|{}", indel1, indel2);
	    //println!("here is start_indel");
        //    let end_indel = start_indel.chars().rev()
          //      .map(|c| match c {
            //        'A' => 'T',
              //      'C' => 'G',
                //    'T' => 'A',
                  //  'G' => 'C',
   //                 '[' => ']',
     //               ']' => '[',
       //             _ => c,
         //       })
           //     .collect::<String>();
         //   full_start = format!("{}|{}", full_start, start_indel);
          //  full_end = format!("{}|{}", full_end, end_indel);
        //}

        let (for_pattern, rev_pattern, pref) = (pattern1.clone(), pattern2.clone(), pref.clone());
	let mut results_hash: Vec<InSilicoPCRProduct> = Vec::new();
	if !opt_c {
	   results_hash = amplify(&full_start, &full_end, opt_l, &sequences, &for_pattern, &iupac, &opt_s, results_hash);
	   }
	
        if results_hash.len() > 0 {
            let mut count = 0;
            for isp in results_hash.iter() {
                    count += 1;
                    let amp = isp.sequence.chars().collect::<String>();
                    let amp_id = if isp.pcrid.is_empty() {
                        format!("amp_{}", count)
                    } else {
                        format!("{}_amp_{}", isp.pcrid, count)
                    };
		    writeln!(&mut w, 
                        "{}\t{}\t{}\t{}",
                        amp_id,
                        isp.seqid,
                        isp.seqpos,
                        isp.length
                    ).expect("unable to write results to file");
                    writeln!(&mut w, "").expect("unable to write to file");
                    let amp_copy = amp.clone();
		    writeln!(&mut w, ">{}\n{}", amp_id, amp_copy).expect("unable to write results to file");
                }
		return Ok(try_digest(&wf).unwrap());
        } else {
            if opt_c {
	        //use function to get results just using a single primer to the end of the sequence if less than option l for length, good for reads
                let c_array = c_amplify(&for_pattern, &rev_pattern, &pref, &iupac, &sequences); //, opt_i);
                if c_array.len() > 0 {
                    let mut sorted_c_array = c_array.clone();
                    sorted_c_array.sort_by(|a, b| a.0.cmp(&b.0));
                    for slice in sorted_c_array {
                        let (amp_type, id, pos, mut seq) = slice;
            //            if opt_i {
			    //println!("doing opt i");
              //              let check = split_indel_pattern2(&seq, &start_pattern);
                //            if check[0].len() == 1 {
                  //              seq = seq.chars().skip(1).collect();
                    //        } else if check[check.len() - 1].len() == 1 {
                      //          seq = seq.chars().take(seq.len() - 1).collect();
                  //          }
                    //        let check2 = split_indel_pattern2(&seq, &end_pattern);
                      //      if check2[0].len() == 1 {
                        //        seq = seq.chars().skip(1).collect();
                    //        } else if check2[check2.len() - 1].len() == 1 {
                      //          seq = seq.chars().take(seq.len() - 1).collect();
                        //    }
                       // }
                        let seqleng = seq.len();
                        if seqleng > opt_l {
                            seq = seq.chars().take(opt_l).collect();
                        }
			writeln!(&mut w, "{}\t{}\t{}\t{}", amp_type, id, pos, opt_l).expect("unable to write results to file");
			writeln!(&mut w, ">{}\n{}", amp_type, seq).expect("unable to write seq results to file");
                    }
                } else {
                    if !pref.is_empty() {
                        println!("{}\t", pref);
                    }
		    writeln!(&mut w, "No amplification").expect("unable to write to file");
                }
            } else {
                if !pref.is_empty() {
                    println!("{}\t", pref);
                }
		writeln!(&mut w, "No amplification").expect("unable to write to file");
            }
	    }
	    }
	    Ok(try_digest(&wf).unwrap())
}


fn read_sequences<R: Read>(reader: BufReader<R>, sequences: &mut HashMap<String, String>) {
    //function to read the query sequence from file
    let mut id = String::new();
    let mut seq = String::new();

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        if line.starts_with('>') {
            if !id.is_empty() {
                sequences.insert(id.clone(), seq.clone());
                seq.clear();
            }
            id = line[1..].split_whitespace().next().unwrap_or("").to_string();
        } else {
	    println!("adding sequences {:?}", &line);
            seq += &line.to_uppercase();
        }
    }

    if !id.is_empty() {
        sequences.insert(id, seq);
    }
}

pub fn decrypt_large_file<'a>(
    source_file_path: &'a PathBuf, k: &'a [u8; 32], n: &'a [u8; 7],) -> Result<Vec<Vec<u8>>, anyhow::Error> {
    let aead = ChaCha20Poly1305::new(k.as_ref().into());
    let mut stream_decryptor = DecryptorBE32::from_aead(aead, n.as_ref().into());
    const BUFFER_LEN: usize = 500+16;
    let mut buffer = [0u8; BUFFER_LEN];
    let mut chunks_reader: Vec<Vec<u8>> = Vec::new();
    let mut source_file = std::fs::File::open(source_file_path).expect("source file not open");
    loop {
       let read_count = source_file.read(&mut buffer).unwrap();
       if read_count == BUFFER_LEN {
          let ciphertext = stream_decryptor
	     .decrypt_next(buffer.as_slice())
	     .map_err(|err| anyhow!("decrypting file error {}", err))?;
	  chunks_reader.push(ciphertext);
       } else {
          let ciphertext = stream_decryptor
	     .decrypt_last(&buffer[..read_count])
	     .map_err(|err| anyhow!("decrypting file error {}", err))?;
	  chunks_reader.push(ciphertext);
	  break;
      }
    }
    Ok(chunks_reader)
}

fn read_primers<R: Read>(reader: BufReader<R>, prims: &mut Vec<(String, String, String)>) {
    //function to read primers from file if option p specified
    for line in reader.lines() {
        let fline = line.expect("Failed to read line");
        let lines: Vec<&str> = fline.split("\n").collect();
        for fline in lines {
            let fields: Vec<&str> = fline.split('\t').collect();
            if fields.len() >= 3 {
                let for_prim = fields[0].to_string();
                let rev_prim = fields[1].to_string();
                let ipref = fields[2].to_string();
                prims.push((for_prim, rev_prim, ipref));
            }
        }
    }
}

fn include_n(primer: &str) -> String {
    //function to allow one mismatch
    let mut new_pattern = String::new();
    if primer.len() > 2 {
        new_pattern = format!(".{}", &primer[1..]);
        for pos in 1..primer.len() {
            new_pattern += &format!("|{}.{}", &primer[0..pos], &primer[pos + 1..]);
        }
    }
    new_pattern
}

fn indel(primer: &str) -> String {
    //function to allow one indel in the match, currently not implemented
    let mut ins = format!("{}.{}", &primer[0..1], &primer[1..]);
    for pos in 1..(primer.len() - 1) {
        ins += &format!("|{}.{}", &primer[0..pos + 1], &primer[1 + pos..]);
    }

    let mut del = format!("{}{}", &primer[0..1], &primer[2..]);
    for pos in 2..(primer.len() - 1) {
        del += &format!("|{}{}", &primer[0..pos], &primer[1 + pos..]);
    }

    format!("{}|{}", ins, del)
}

fn amplify(
    start_pattern: &str,
    end_pattern: &str,
    maxlength: usize,
    //opt_r: bool,
    sequences: &HashMap<String, String>,
    for_pattern: &str,
    iupac: &[(&str, &str); 12],
    filename: &std::path::Path,
    mut results_hash: Vec<InSilicoPCRProduct>,
) -> Vec<InSilicoPCRProduct> {
    let pattern1_rc: String = start_pattern.chars().rev().map(|c| match c {
        'A' => 'T',
	'C' => 'G',
	'G' => 'C',
	'T' => 'A',
	'[' => ']',
	']' => '[',
	other => other,
	}).collect();
    let pattern2_rc: String = end_pattern.chars().rev().map(|c| match c {
        'A' => 'T',
	'C' => 'G',
	'G' => 'C',
	'T' => 'A',
	'[' => ']',
	']' => '[',
	other => other,
	}).collect();
    let pattern_regex1 = Regex::new(start_pattern).expect("Invalid regex for the forward primer");
    let pattern_regex1_rc = Regex::new(&pattern1_rc).expect("invalid Regex pattern for reverse complement forward primer for_rc");
    let pattern_regex2 = Regex::new(end_pattern).expect("invalid regex pattern for reverse primer");
    let pattern_regex2_rc = Regex::new(&pattern2_rc).expect("invalid regex pattern for reverse complement reverse primer rev_rc");
    
    let regexes = vec![(pattern_regex1.clone(), pattern_regex2.clone()), (pattern_regex1_rc.clone(), pattern_regex2_rc.clone()), (pattern_regex1.clone(), pattern_regex2_rc.clone()),(pattern_regex1_rc.clone(), pattern_regex2.clone()),(pattern_regex2.clone(), pattern_regex1_rc.clone()),(pattern_regex2_rc.clone(), pattern_regex1.clone())];
    let mut r = 0;
    let mut frag_cnt = 0;

    for (regexy_k, regexy_v) in regexes.iter() {
        //loop through the regexes for the primers
        for (id, sequencey) in sequences.iter() {
	    frag_cnt = 0;
            for (_m, fragment) in regexy_k.find_iter(sequencey).enumerate() {
	       let mut position = fragment.start();
	       let mut otherpos = if (sequencey.len()-position) < maxlength { sequencey.len() } else {  maxlength+position };
	       let subfragment_to_maximum = if r > 2 { &sequencey[position..otherpos] } else { &sequencey[position..otherpos] };
	       let mut adjusted_subfragment = subfragment_to_maximum.to_string();
               //implement opt_i later if requested
              // if opt_i {
               //        //use the indel function if indels required 
                //       let check = split_indel_pattern(&adjusted_subfragment, &regexy_v);
                 //      if check[0].len() == 1 {
                  //            fragments2 = adjusted_subfragment.chars().skip(1).collect();
                   //           }  else if check[check.len() - 1].len() == 1 {
                    //          fragments2 = adjusted_subfragment.chars().take(adjusted_subfragment.len() - 1).collect();
                   //           }  
                    //    }
              // else {
	        //    let mut fragments2 = regexy_v.split(&adjusted_subfragment);
                 //   }
	       let mut fragments2 = regexy_v.split(&adjusted_subfragment);
	       let mut capsmatches: Vec<String> = Vec::new();
	       let mut capscnt = 0;
	       for fcaps in regexy_v.captures_iter(&adjusted_subfragment) {
		   let capsgroup = fcaps[0].to_string();
		   println!("capsgroup {:?}", &capsgroup);
	           capsmatches.push(capsgroup.clone());
		   capscnt+=1;
		   }
               if let Some(fragment2_0) = fragments2.next() {
		     let capsy = if capsmatches.len() > 0 { capsmatches.pop().expect("troublepopping").to_string() } else { "".to_string() };
		     let amp = format!("{}{}",fragment2_0.to_string(),capsy); //&sequence[position..position + lenfragment];
                     if let Some(fragment2_1) = fragments2.next() {
                                let lenfragment = amp.len();
				let nam: &str = filename.file_stem().unwrap().to_str().unwrap();
				if amp.len() > 40 {
				    let new_isp = InSilicoPCRProduct { seqid: id.to_string(), pcrid: format!("{}_{}", nam.to_string(),frag_cnt), seqpos: position, length: lenfragment, sequence: amp.to_string() };
				    results_hash.push(new_isp);
				    }
                                }
                     }
                     position += fragment.len() + subfragment_to_maximum.len() + end_pattern.len();
                 }
	    }
	r+=1;
	}
    results_hash
}

fn split_indel_pattern(seq: &str, indel_pattern: &regex::Regex) -> Vec<String> {
    indel_pattern.split(seq).into_iter().map(|s| s.to_string()).collect()

}

fn split_indel_pattern2(seq: &str, indel_pattern: &str) -> Vec<String> {
     seq.split(indel_pattern).into_iter()
        .map(|s| s.to_string())
        .collect()
}

fn c_amplify(
    for_pattern: &str,
    rev_pattern: &str,
    pref: &str,
    iupac: &[(&str, &str); 12],
    sequences: &HashMap<String, String>,
    //opt_i: bool,
) -> Vec<(String, String, usize, String)> {
     let pattern1_rc: String = for_pattern.chars().rev().map(|c| match c {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        '[' => ']',
        ']' => '[',
        other => other,
    }).collect();

    let pattern2_rc: String = rev_pattern.chars().rev().map(|c| match c {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        '[' => ']',
        ']' => '[',
        other => other,
    }).collect();

    let mut results: Vec<(String, String, usize, String)> = Vec::new();
    let mut position = 0;
    let pattern_regex1 = Regex::new(for_pattern).expect("Invalid Regex pattern for forward primer");
    let pattern_regex1_rc = Regex::new(&pattern1_rc).expect("invalid Regex pattern for reverse complement forward primer for_rc");
    let pattern_regex2 = Regex::new(rev_pattern).expect("invalid regex pattern for reverse primer");
    let pattern_regex2_rc = Regex::new(&pattern2_rc).expect("invalid regex pattern for reverse complement reverse primer rev_rc");

    let regexes = vec![pattern_regex1, pattern_regex2, pattern_regex1_rc, pattern_regex2_rc];
    let mut r = 0;
    let mut ptype: HashMap<String, u32> = HashMap::new();
  
    for regexy in regexes {
      //loop through the regexes for the primers
      for (id, sequence) in sequences.iter() {
        position = 0;
	let mut rev_sequence = String::new();
        let mut subfragment: &str = "";
	//find matches for each sequence 
        for (m, fragment) in regexy.find_iter(sequence).enumerate() {
	    
            let outpos = fragment.start();
	    //if reverse complement match  get appropriate sequence, otherwise use the forward sequence
	    let subfragment = if r > 2 { &sequence[0..outpos+fragment.len()] } else { &sequence[outpos..sequence.len()] };
	    let mut adjusted_subfragment = subfragment.to_string();
	    //implement opt_i later if requested
	   // if opt_i {
	       //use the indel function if indels required 
	     //  let check = split_indel_pattern(&adjusted_subfragment, &regexy);
	    //   if check[0].len() == 1 {
	     //     adjusted_subfragment = adjusted_subfragment.chars().skip(1).collect();
//		  }  else if check[check.len() - 1].len() == 1 {
//		       adjusted_subfragment = adjusted_subfragment.chars().take(adjusted_subfragment.len() - 1).collect();
//		       }
//		  }    
		    let keya = vec!["p1", "p2", "p1rc", "p2rc"];
		    if r == 0 {
		        *ptype.entry("p1".to_string()).or_default()+=1;
	            } else if r == 1 {
			*ptype.entry("p2".to_string()).or_default()+=1;
		    } else if r == 2 {
			*ptype.entry("p1rc".to_string()).or_default()+=1;
		    } else if r == 3 {
			*ptype.entry("p2rc".to_string()).or_default()+=1;
		    } else {
		        println!("Unexpected index number for the primer search {:?}", r);
			}
		    position = outpos + fragment.len();
		    let use_variable = if r > 2 { position } else { outpos };
		    let amp_type = format!("{}_{}_{}", pref, keya[r], ptype.get(keya[r]).unwrap());
                    results.push((amp_type, id.to_string(), use_variable, adjusted_subfragment.to_string()));
                    }
       }
       r+=1;
    }
    results
}