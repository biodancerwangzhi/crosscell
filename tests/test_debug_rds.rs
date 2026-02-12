//! Debug test for RDS parsing with detailed output

use std::path::Path;
use std::fs::File;
use std::io::{BufReader, Read};
use flate2::read::GzDecoder;

#[test]
fn test_debug_dgc_bytes() {
    let path = Path::new("tests/data/test_dgc_simple.rds");
    
    if !path.exists() {
        eprintln!("Skipping: {} not found", path.display());
        return;
    }
    
    // Read and decompress
    let file = File::open(path).unwrap();
    let buf_reader = BufReader::new(file);
    let mut gz = GzDecoder::new(buf_reader);
    let mut data = Vec::new();
    gz.read_to_end(&mut data).unwrap();
    
    println!("\n========================================");
    println!("Decompressed RDS file ({} bytes)", data.len());
    println!("========================================\n");
    
    // Print all bytes in hex
    for (i, chunk) in data.chunks(16).enumerate() {
        print!("{:04x}: ", i * 16);
        for b in chunk {
            print!("{:02x} ", b);
        }
        // Pad if needed
        for _ in 0..(16 - chunk.len()) {
            print!("   ");
        }
        // Print ASCII
        print!(" |");
        for b in chunk {
            if *b >= 32 && *b < 127 {
                print!("{}", *b as char);
            } else {
                print!(".");
            }
        }
        println!("|");
    }
    
    // Parse header manually
    println!("\n--- Header Analysis ---");
    println!("Magic: {:02x} {:02x} = '{}{}'", data[0], data[1], data[0] as char, data[1] as char);
    
    let format_version = u32::from_be_bytes([data[2], data[3], data[4], data[5]]);
    println!("Format version: {}", format_version);
    
    let writer_version = &data[6..10];
    println!("Writer version: {:?}", writer_version);
    
    let reader_version = &data[10..14];
    println!("Reader version: {:?}", reader_version);
    
    if format_version >= 3 {
        let encoding_len = u32::from_be_bytes([data[14], data[15], data[16], data[17]]) as usize;
        println!("Encoding length: {}", encoding_len);
        let encoding = String::from_utf8_lossy(&data[18..18+encoding_len]);
        println!("Encoding: {}", encoding);
        
        let obj_start = 18 + encoding_len;
        println!("\n--- Main Object Header (offset {}) ---", obj_start);
        parse_object_at(&data, obj_start, 0);
    }
}

fn parse_object_at(data: &[u8], offset: usize, depth: usize) -> usize {
    let indent = "  ".repeat(depth);
    
    if offset + 4 > data.len() {
        println!("{}ERROR: Not enough bytes at offset {}", indent, offset);
        return offset;
    }
    
    let header = &data[offset..offset+4];
    let sexp_type = header[3];
    let flags = header[2];
    let levels = ((header[0] as u16) << 8) | (header[1] as u16);
    
    let is_object = (flags & 0x01) != 0;
    let has_attributes = (flags & 0x02) != 0;
    let has_tag = (flags & 0x04) != 0;
    
    println!("{}[{}] Header: {:02x} {:02x} {:02x} {:02x} = {} (obj={}, attr={}, tag={})", 
        indent, offset, header[0], header[1], header[2], header[3],
        sexp_type_name(sexp_type), is_object, has_attributes, has_tag);
    
    let mut pos = offset + 4;
    
    match sexp_type {
        25 => { // S4
            if has_attributes {
                println!("{}  Parsing S4 attributes...", indent);
                pos = parse_object_at(data, pos, depth + 1);
            }
        }
        2 => { // LIST (pairlist)
            if has_tag {
                println!("{}  Parsing pairlist tag...", indent);
                pos = parse_object_at(data, pos, depth + 1);
            }
            println!("{}  Parsing pairlist value...", indent);
            pos = parse_object_at(data, pos, depth + 1);
            println!("{}  Parsing pairlist CDR...", indent);
            pos = parse_object_at(data, pos, depth + 1);
        }
        1 => { // SYM
            println!("{}  Parsing symbol name (CHARSXP)...", indent);
            pos = parse_object_at(data, pos, depth + 1);
        }
        9 => { // CHAR
            if pos + 4 > data.len() {
                println!("{}  ERROR: Not enough bytes for string length", indent);
                return pos;
            }
            let strlen = i32::from_be_bytes([data[pos], data[pos+1], data[pos+2], data[pos+3]]);
            pos += 4;
            if strlen == -1 {
                println!("{}  String: NA", indent);
            } else {
                let s = String::from_utf8_lossy(&data[pos..pos+strlen as usize]);
                println!("{}  String: \"{}\" (len={})", indent, s, strlen);
                pos += strlen as usize;
            }
        }
        13 => { // INT
            if pos + 4 > data.len() {
                println!("{}  ERROR: Not enough bytes for int length", indent);
                return pos;
            }
            let len = u32::from_be_bytes([data[pos], data[pos+1], data[pos+2], data[pos+3]]) as usize;
            pos += 4;
            println!("{}  Integer vector length: {}", indent, len);
            pos += len * 4;
        }
        14 => { // REAL
            if pos + 4 > data.len() {
                println!("{}  ERROR: Not enough bytes for real length", indent);
                return pos;
            }
            let len = u32::from_be_bytes([data[pos], data[pos+1], data[pos+2], data[pos+3]]) as usize;
            pos += 4;
            println!("{}  Double vector length: {}", indent, len);
            pos += len * 8;
        }
        19 => { // VEC (list)
            if pos + 4 > data.len() {
                println!("{}  ERROR: Not enough bytes for list length", indent);
                return pos;
            }
            let len = u32::from_be_bytes([data[pos], data[pos+1], data[pos+2], data[pos+3]]) as usize;
            pos += 4;
            println!("{}  List length: {}", indent, len);
            for i in 0..len {
                println!("{}  List element {}:", indent, i);
                pos = parse_object_at(data, pos, depth + 1);
            }
            if has_attributes {
                println!("{}  Parsing list attributes...", indent);
                pos = parse_object_at(data, pos, depth + 1);
            }
        }
        16 => { // STR
            if pos + 4 > data.len() {
                println!("{}  ERROR: Not enough bytes for str length", indent);
                return pos;
            }
            let len = u32::from_be_bytes([data[pos], data[pos+1], data[pos+2], data[pos+3]]) as usize;
            pos += 4;
            println!("{}  String vector length: {}", indent, len);
            for i in 0..len {
                println!("{}  String element {}:", indent, i);
                pos = parse_object_at(data, pos, depth + 1);
            }
            if has_attributes {
                println!("{}  Parsing string vector attributes...", indent);
                pos = parse_object_at(data, pos, depth + 1);
            }
        }
        254 => { // NILVALUE
            println!("{}  NILVALUE", indent);
        }
        255 => { // REF
            let ref_index = ((header[0] as usize) << 16) | ((header[1] as usize) << 8) | (header[2] as usize);
            println!("{}  Reference index: {}", indent, ref_index);
        }
        _ => {
            println!("{}  Unknown type, skipping...", indent);
        }
    }
    
    pos
}

fn sexp_type_name(t: u8) -> &'static str {
    match t {
        0 => "NIL",
        1 => "SYM",
        2 => "LIST",
        3 => "CLO",
        4 => "ENV",
        5 => "PROM",
        6 => "LANG",
        9 => "CHAR",
        10 => "LGL",
        13 => "INT",
        14 => "REAL",
        16 => "STR",
        19 => "VEC",
        25 => "S4",
        254 => "NILVALUE",
        255 => "REF",
        _ => "UNKNOWN",
    }
}
